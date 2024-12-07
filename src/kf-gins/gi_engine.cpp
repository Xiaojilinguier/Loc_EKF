#include "gi_engine.h"
#include "common/earth.h"
#include "common/rotation.h"
#include "insmech.h"
#include "satellite.h"
#include <random>

// 计算飞行器与卫星之间的测距
double GIEngine::computeLaserRange(const Eigen::Vector3d &aircraft_pos, const Eigen::Vector3d &satellite_pos) {
    return (aircraft_pos - satellite_pos).norm(); // 返回飞行器与卫星之间的欧几里得距离
}

// 计算并绘制误差图
void GIEngine::compute_and_print_error(int timeIndex) {
    std::vector<double> real_pos = getRealPos(timeIndex); // 纬经高(deg) 形式
    double x, y, z;
    LLA_to_ecef(real_pos[0], real_pos[1], real_pos[2], x, y, z);
    std::vector<double> real_vel = getRealVel(timeIndex);
    std::vector<double> est_pos  = getEstimatePos();
    std::vector<double> est_vel  = getEstimateVel();
    double error_x               = x - est_pos[0];
    double error_y               = y - est_pos[1];
    double error_z               = z - est_pos[2];
    double error_vx              = real_vel[0] - est_vel[0];
    double error_vy              = real_vel[1] - est_vel[1];
    double error_vz              = real_vel[2] - est_vel[2];
    std::cout << "pos_error(ECEF, m): " << error_x << ", " << error_y << ", " << error_z << std::endl;
    std::cout << "vel_error(NED, m/s): " << error_vx << ", " << error_vy << ", " << error_vz << std::endl;
}

// BLH to ECEF 的Jacobian矩阵
Matrix3d GIEngine::computeBLHtoECEFJacobian(const Vector3d &blh) {
    double phi    = blh[0]; // Latitude
    double lambda = blh[1]; // Longitude
    double h      = blh[2]; // Height

    Eigen::Vector2d rmn = Earth::meridianPrimeVerticalRadius(phi);
    double N            = rmn[1];
    double e2           = 0.00669437999014; // WGS84 eccentricity squared

    Matrix3d jacobian = Matrix3d::Zero();
    // Partial derivatives w.r.t. latitude (phi)
    jacobian(0, 0) = -(N + h) * sin(phi) * cos(lambda);
    jacobian(1, 0) = -(N + h) * sin(phi) * sin(lambda);
    jacobian(2, 0) = (N * (1 - e2) + h) * cos(phi);

    // Partial derivatives w.r.t. longitude (lambda)
    jacobian(0, 1) = -(N + h) * cos(phi) * sin(lambda);
    jacobian(1, 1) = (N + h) * cos(phi) * cos(lambda);
    jacobian(2, 1) = 0;

    // Partial derivatives w.r.t. height (h)
    jacobian(0, 2) = cos(phi) * cos(lambda);
    jacobian(1, 2) = cos(phi) * sin(lambda);
    jacobian(2, 2) = sin(phi);

    return jacobian;
}

GIEngine::GIEngine(GINSOptions &options, std::string tlepath, std::string coverpath, std::string truthpath,
                   double starttime, double endtime) {
    this->starttime_     = starttime;
    this->tlepath_       = tlepath;
    this->coverpath_     = coverpath;
    this->truthpath_     = truthpath;
    this->options_       = options;
    this->dataProcesser_ = new DataProcessing(truthpath, coverpath);
    dataProcesser_->load_position_velocity(starttime, endtime, realPositions_, realVelocities_, truthpath);
    // dataProcesser
    options_.print_options();
    timestamp_        = 0;
    choose_sat_index_ = 0;
    // 设置协方差矩阵，系统噪声阵和系统误差状态矩阵大小
    // resize covariance matrix, system noise matrix, and system error state matrix
    Cov_.resize(RANK, RANK);          // 协方差阵
    Qc_.resize(NOISERANK, NOISERANK); // 系统噪声矩阵
    dx_.resize(RANK, 1);              // 状态误差向量
    Cov_.setZero();
    Qc_.setZero();
    dx_.setZero();

    // 初始化系统噪声阵
    // 方差的计算通过cwiseProduct和asDiagonal生成一个对角阵
    // 六类噪声
    // initialize noise matrix
    auto imunoise                   = options_.imunoise;
    Qc_.block(ARW_ID, ARW_ID, 3, 3) = imunoise.gyr_arw.cwiseProduct(imunoise.gyr_arw).asDiagonal();
    Qc_.block(VRW_ID, VRW_ID, 3, 3) = imunoise.acc_vrw.cwiseProduct(imunoise.acc_vrw).asDiagonal();
    Qc_.block(BGSTD_ID, BGSTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.gyrbias_std.cwiseProduct(imunoise.gyrbias_std).asDiagonal();
    Qc_.block(BASTD_ID, BASTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.accbias_std.cwiseProduct(imunoise.accbias_std).asDiagonal();
    Qc_.block(SGSTD_ID, SGSTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.gyrscale_std.cwiseProduct(imunoise.gyrscale_std).asDiagonal();
    Qc_.block(SASTD_ID, SASTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.accscale_std.cwiseProduct(imunoise.accscale_std).asDiagonal();

    // 设置系统状态(位置、速度、姿态和IMU误差)初值和初始协方差
    // set initial state (position, velocity, attitude and IMU error) and covariance
    initialize(options_.initstate, options_.initstate_std);
}

// 获取卫星ECEF位置
bool GIEngine::getSatellitePosition(satellite *sat, const std::string &now_time, Eigen::Vector3d &satellite_pos) {
    if (!sat->get_position(now_time, satellite_pos)) {
        return false; // 如果获取卫星位置失败，返回 false
    }
    return true;
}

void GIEngine::initialize(const NavState &initstate, const NavState &initstate_std) {

    // 初始化位置、速度、姿态
    // initialize position, velocity and attitude
    pvacur_.pos       = initstate.pos;
    pvacur_.vel       = initstate.vel;
    pvacur_.att.euler = initstate.euler;
    pvacur_.att.cbn   = Rotation::euler2matrix(pvacur_.att.euler);
    pvacur_.att.qbn   = Rotation::euler2quaternion(pvacur_.att.euler);
    // 初始化IMU误差
    // initialize imu error
    imuerror_ = initstate.imuerror;

    // 给上一时刻状态赋同样的初值
    // set the same value to the previous state
    pvapre_ = pvacur_;

    // 初始化协方差
    // initialize covariance
    ImuError imuerror_std            = initstate_std.imuerror;
    Cov_.block(P_ID, P_ID, 3, 3)     = initstate_std.pos.cwiseProduct(initstate_std.pos).asDiagonal();
    Cov_.block(V_ID, V_ID, 3, 3)     = initstate_std.vel.cwiseProduct(initstate_std.vel).asDiagonal();
    Cov_.block(PHI_ID, PHI_ID, 3, 3) = initstate_std.euler.cwiseProduct(initstate_std.euler).asDiagonal();
    Cov_.block(BG_ID, BG_ID, 3, 3)   = imuerror_std.gyrbias.cwiseProduct(imuerror_std.gyrbias).asDiagonal();
    Cov_.block(BA_ID, BA_ID, 3, 3)   = imuerror_std.accbias.cwiseProduct(imuerror_std.accbias).asDiagonal();
    Cov_.block(SG_ID, SG_ID, 3, 3)   = imuerror_std.gyrscale.cwiseProduct(imuerror_std.gyrscale).asDiagonal();
    Cov_.block(SA_ID, SA_ID, 3, 3)   = imuerror_std.accscale.cwiseProduct(imuerror_std.accscale).asDiagonal();
}

void GIEngine::newImuProcess(int truthdatarate) {

    // 当前IMU时间作为系统当前状态时间,
    // set current IMU time as the current state time
    timestamp_ = imucur_.time;

    double miu   = 0;   // 测距误差平均值
    double sigma = 0.1; // 测距误差标准差
    // 随机数生成器
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(miu, sigma);

    // 如果距离上次激光测距的时间超过2秒，则更新卫星位置和激光测距信息
    if (timestamp_ - last_laser_time_ >= 2.0) {
        std::string now_time = dataProcesser_->time_to_mmss(std::to_string(timestamp_));
        std::cout << "now_time: " << now_time << std::endl;

        // 计算并输出定位误差
        int timeindex = (timestamp_ - starttime_) * truthdatarate;
        std::cout << "滤波前误差：" << std::endl;
        compute_and_print_error(timeindex);

        // 获取卫星信息
        std::vector<std::string> satellite_id_list;
        if (!dataProcesser_->find_satellites_by_time(now_time, satellite_id_list)) {
            std::cerr << "No satellites found at time " << timestamp_ << std::endl;
        }

        std::cout << "卫星可覆盖列表：[ ";
        for (int i = 0; i < satellite_id_list.size(); i++) {
            std::cout << satellite_id_list[i];
            if (i != satellite_id_list.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;

        if (satellite_id_list.empty()) {
            std::cout << "No satellite coverage at time " << now_time << std::endl;
            // 只传播导航状态
            // only propagate navigation state
            insPropagation(imupre_, imucur_);
        } else {
            insPropagation(imupre_, imucur_);
            // 随机误差
            double random_error = distribution(generator);

            // 选择卫星并计算激光测距
            std::size_t satellite_index = (choose_sat_index_++) % satellite_id_list.size();
            std::string satellite_id    = satellite_id_list[satellite_index];

            std::cout << "选择卫星：" << satellite_id << std::endl;
            Eigen::Vector3d satellite_pos_ecef;

            // 计算真实卫星位置
            satellite *sat = new satellite(satellite_id, tlepath_);
            if (!getSatellitePosition(sat, now_time, satellite_pos_ecef)) {
                std::cerr << "Failed to get satellite position for " << satellite_id << std::endl;
                return;
            }

            // 计算飞行器与卫星之间的真实测距
            int pos_index = (timestamp_ - starttime_) * truthdatarate;

            Eigen::Vector3d aircraft_real_pos_now_blh = Eigen::Vector3d(
                realPositions_[pos_index][0] * D2R, realPositions_[pos_index][1] * D2R, realPositions_[pos_index][2]);
            Eigen::Vector3d aircraft_real_pos_now_ecef = Earth::blh2ecef(aircraft_real_pos_now_blh);

            double rangeMeasurement = computeLaserRange(aircraft_real_pos_now_ecef, satellite_pos_ecef);

            // 添加高斯误差
            rangeMeasurement += random_error;

            // 添加激光测距数据到 GIEngine 并进行状态更新
            LaserUpdate(rangeMeasurement, satellite_pos_ecef);
            stateFeedback();
            // 计算并输出定位误差
            int timeindex = (timestamp_ - starttime_) * truthdatarate;
            std::cout << "滤波后误差：" << std::endl;
            compute_and_print_error(timeindex);
        }

        // 更新上次激光测距的时间
        last_laser_time_ = timestamp_;
        std::cout << "----------------------------" << std::endl;
    }

    else {
        // 还未开始激光测距，只用IMU传播导航状态
        // only propagate navigation state
        insPropagation(imupre_, imucur_);
    }
    // // 计算并输出定位误差
    // int timeindex = (timestamp_ - starttime_) * truthdatarate;
    // std::cout << "滤波后误差：" << std::endl;
    // compute_and_print_error(timeindex);
    // 检查协方差矩阵对角线元素
    // check diagonal elements of current covariance matrix
    checkCov();

    // 更新上一时刻的状态和IMU数据
    // update system state and imudata at the previous epoch
    pvapre_ = pvacur_;
    imupre_ = imucur_;
}

void GIEngine::imuCompensate(IMU &imu) {

    // 补偿IMU零偏
    // compensate the imu bias
    imu.dtheta -= imuerror_.gyrbias * imu.dt;
    imu.dvel -= imuerror_.accbias * imu.dt;

    // 补偿IMU比例因子
    // compensate the imu scale
    Eigen::Vector3d gyrscale, accscale;
    gyrscale   = Eigen::Vector3d::Ones() + imuerror_.gyrscale;
    accscale   = Eigen::Vector3d::Ones() + imuerror_.accscale;
    imu.dtheta = imu.dtheta.cwiseProduct(gyrscale.cwiseInverse()); // 按元素乘 比例因子的倒数来补偿
    imu.dvel   = imu.dvel.cwiseProduct(accscale.cwiseInverse());
}

void GIEngine::insPropagation(IMU &imupre, IMU &imucur) {

    // 对当前IMU数据(imucur)补偿误差, 上一IMU数据(imupre)已经补偿过了
    // compensate imu error to 'imucur', 'imupre' has been compensated
    imuCompensate(imucur);
    // IMU状态更新(机械编排算法) 位置速度姿态更新
    // update imustate(mechanization)
    INSMech::insMech(pvapre_, pvacur_, imupre, imucur);

    // 系统噪声传播，姿态误差采用phi角误差模型
    // system noise propagate, phi-angle error model for attitude error
    Eigen::MatrixXd Phi, F, Qd, G;

    // 初始化Phi阵(状态转移矩阵)，F阵，Qd阵(传播噪声阵)，G阵(噪声驱动阵)
    // initialize Phi (state transition), F matrix, Qd(propagation noise) and G(noise driven) matrix
    Phi.resizeLike(Cov_);
    F.resizeLike(Cov_);
    Qd.resizeLike(Cov_);
    G.resize(RANK, NOISERANK);
    Phi.setIdentity();
    F.setZero();
    Qd.setZero();
    G.setZero();

    // 使用上一历元状态计算状态转移矩阵
    // compute state transition matrix using the previous state
    Eigen::Vector2d rmrn;
    Eigen::Vector3d wie_n, wen_n;
    double gravity;
    rmrn    = Earth::meridianPrimeVerticalRadius(pvapre_.pos[0]);
    gravity = Earth::gravity(pvapre_.pos);
    wie_n << WGS84_WIE * cos(pvapre_.pos[0]), 0, -WGS84_WIE * sin(pvapre_.pos[0]);
    wen_n << pvapre_.vel[1] / (rmrn[1] + pvapre_.pos[2]), -pvapre_.vel[0] / (rmrn[0] + pvapre_.pos[2]),
        -pvapre_.vel[1] * tan(pvapre_.pos[0]) / (rmrn[1] + pvapre_.pos[2]);

    Eigen::Matrix3d temp;
    Eigen::Vector3d accel, omega;
    double rmh, rnh;

    rmh   = rmrn[0] + pvapre_.pos[2];  // RM+H
    rnh   = rmrn[1] + pvapre_.pos[2];  // RN+H
    accel = imucur.dvel / imucur.dt;   // 根据输入的速度增量计算加速度
    omega = imucur.dtheta / imucur.dt; // 根据输入的角增量计算角速度

    // 位置误差
    // position error
    temp.setZero();
    temp(0, 0)                = -pvapre_.vel[2] / rmh;
    temp(0, 2)                = pvapre_.vel[0] / rmh;
    temp(1, 0)                = pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh;
    temp(1, 1)                = -(pvapre_.vel[2] + pvapre_.vel[0] * tan(pvapre_.pos[0])) / rnh;
    temp(1, 2)                = pvapre_.vel[1] / rnh;
    F.block(P_ID, P_ID, 3, 3) = temp;
    F.block(P_ID, V_ID, 3, 3) = Eigen::Matrix3d::Identity();

    // 速度误差
    // velocity error
    temp.setZero();
    temp(0, 0) = -2 * pvapre_.vel[1] * WGS84_WIE * cos(pvapre_.pos[0]) / rmh -
                 pow(pvapre_.vel[1], 2) / rmh / rnh / pow(cos(pvapre_.pos[0]), 2);
    temp(0, 2) = pvapre_.vel[0] * pvapre_.vel[2] / rmh / rmh - pow(pvapre_.vel[1], 2) * tan(pvapre_.pos[0]) / rnh / rnh;
    temp(1, 0) = 2 * WGS84_WIE * (pvapre_.vel[0] * cos(pvapre_.pos[0]) - pvapre_.vel[2] * sin(pvapre_.pos[0])) / rmh +
                 pvapre_.vel[0] * pvapre_.vel[1] / rmh / rnh / pow(cos(pvapre_.pos[0]), 2);
    temp(1, 2) = (pvapre_.vel[1] * pvapre_.vel[2] + pvapre_.vel[0] * pvapre_.vel[1] * tan(pvapre_.pos[0])) / rnh / rnh;
    temp(2, 0) = 2 * WGS84_WIE * pvapre_.vel[1] * sin(pvapre_.pos[0]) / rmh;
    temp(2, 2) = -pow(pvapre_.vel[1], 2) / rnh / rnh - pow(pvapre_.vel[0], 2) / rmh / rmh +
                 2 * gravity / (sqrt(rmrn[0] * rmrn[1]) + pvapre_.pos[2]);
    F.block(V_ID, P_ID, 3, 3) = temp;
    temp.setZero();
    temp(0, 0)                  = pvapre_.vel[2] / rmh;
    temp(0, 1)                  = -2 * (WGS84_WIE * sin(pvapre_.pos[0]) + pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh);
    temp(0, 2)                  = pvapre_.vel[0] / rmh;
    temp(1, 0)                  = 2 * WGS84_WIE * sin(pvapre_.pos[0]) + pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh;
    temp(1, 1)                  = (pvapre_.vel[2] + pvapre_.vel[0] * tan(pvapre_.pos[0])) / rnh;
    temp(1, 2)                  = 2 * WGS84_WIE * cos(pvapre_.pos[0]) + pvapre_.vel[1] / rnh;
    temp(2, 0)                  = -2 * pvapre_.vel[0] / rmh;
    temp(2, 1)                  = -2 * (WGS84_WIE * cos(pvapre_.pos(0)) + pvapre_.vel[1] / rnh);
    F.block(V_ID, V_ID, 3, 3)   = temp;
    F.block(V_ID, PHI_ID, 3, 3) = Rotation::skewSymmetric(pvapre_.att.cbn * accel);
    F.block(V_ID, BA_ID, 3, 3)  = pvapre_.att.cbn;
    F.block(V_ID, SA_ID, 3, 3)  = pvapre_.att.cbn * (accel.asDiagonal());

    // 姿态误差
    // attitude error
    temp.setZero();
    temp(0, 0) = -WGS84_WIE * sin(pvapre_.pos[0]) / rmh;
    temp(0, 2) = pvapre_.vel[1] / rnh / rnh;
    temp(1, 2) = -pvapre_.vel[0] / rmh / rmh;
    temp(2, 0) = -WGS84_WIE * cos(pvapre_.pos[0]) / rmh - pvapre_.vel[1] / rmh / rnh / pow(cos(pvapre_.pos[0]), 2);
    temp(2, 2) = -pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh / rnh;
    F.block(PHI_ID, P_ID, 3, 3) = temp;
    temp.setZero();
    temp(0, 1)                    = 1 / rnh;
    temp(1, 0)                    = -1 / rmh;
    temp(2, 1)                    = -tan(pvapre_.pos[0]) / rnh;
    F.block(PHI_ID, V_ID, 3, 3)   = temp;
    F.block(PHI_ID, PHI_ID, 3, 3) = -Rotation::skewSymmetric(wie_n + wen_n);
    F.block(PHI_ID, BG_ID, 3, 3)  = -pvapre_.att.cbn;
    F.block(PHI_ID, SG_ID, 3, 3)  = -pvapre_.att.cbn * (omega.asDiagonal());

    // IMU零偏误差和比例因子误差，建模成一阶高斯-马尔科夫过程
    // imu bias error and scale error, modeled as the first-order Gauss-Markov process
    F.block(BG_ID, BG_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(BA_ID, BA_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(SG_ID, SG_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(SA_ID, SA_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();

    // 系统噪声驱动矩阵
    // system noise driven matrix
    G.block(V_ID, VRW_ID, 3, 3)    = pvapre_.att.cbn;
    G.block(PHI_ID, ARW_ID, 3, 3)  = pvapre_.att.cbn;
    G.block(BG_ID, BGSTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(BA_ID, BASTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(SG_ID, SGSTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(SA_ID, SASTD_ID, 3, 3) = Eigen::Matrix3d::Identity();

    // 状态转移矩阵
    // compute the state transition matrix
    Phi.setIdentity();
    Phi = Phi + F * imucur.dt;

    // 计算系统传播噪声
    // compute system propagation noise
    Qd = G * Qc_ * G.transpose() * imucur.dt;
    Qd = (Phi * Qd * Phi.transpose() + Qd) / 2;

    // EKF预测传播系统协方差和系统误差状态
    // do EKF predict to propagate covariance and error state
    EKFPredict(Phi, Qd);
}

int GIEngine::isToUpdate(double imutime1, double imutime2, double updatetime) const {

    if (abs(imutime1 - updatetime) < TIME_ALIGN_ERR) {
        // 更新时间靠近imutime1
        // updatetime is near to imutime1
        return 1;
    } else if (abs(imutime2 - updatetime) <= TIME_ALIGN_ERR) {
        // 更新时间靠近imutime2
        // updatetime is near to imutime2
        return 2;
    } else if (imutime1 < updatetime && updatetime < imutime2) {
        // 更新时间在imutime1和imutime2之间, 但不靠近任何一个
        // updatetime is between imutime1 and imutime2, but not near to either
        return 3;
    } else {
        // 更新时间不在imutimt1和imutime2之间，且不靠近任何一个
        // updatetime is not bewteen imutime1 and imutime2, and not near to either.
        return 0;
    }
}

void GIEngine::EKFPredict(Eigen::MatrixXd &Phi, Eigen::MatrixXd &Qd) {

    assert(Phi.rows() == Cov_.rows());
    assert(Qd.rows() == Cov_.rows());

    // 传播系统协方差和误差状态
    // propagate system covariance and error state
    Cov_ = Phi * Cov_ * Phi.transpose() + Qd;
    dx_  = Phi * dx_;
}

// 新增的添加激光测距数据到系统并进行更新的函数
void GIEngine::LaserUpdate(const double &z, const Eigen::Vector3d &satellite_pos_ecef) {
    // 激光测距数据的测量矩阵 H
    Eigen::MatrixXd H;
    H.resize(1, Cov_.rows());
    H.setZero();

    Eigen::Matrix3d Dr_inv;
    Dr_inv = Earth::DRi(pvacur_.pos);

    Eigen::Vector3d delta_r              = dx_.block(P_ID, 0, 3, 1);
    Eigen::Vector3d ins_pos_feedback_blh = pvacur_.pos - Dr_inv * delta_r; // 单位：rad rad m

    Eigen::Vector3d ins_pos_feedback_ecef = Earth::blh2ecef(ins_pos_feedback_blh);

    //  计算飞行器与卫星之间的测距
    double range_measured  = z;                                                            // 观测值
    double range_predicted = computeLaserRange(ins_pos_feedback_ecef, satellite_pos_ecef); // 计算预测的测距

    // 计算测量残差 dz
    Eigen::VectorXd dz(1);
    dz(0) = range_measured - range_predicted; // 观测值减去预测值

    // 观测矩阵 H 计算（通过位置对测距的影响）
    Eigen::Vector3d diff = ins_pos_feedback_ecef - satellite_pos_ecef;
    double range         = diff.norm();
    Matrix3d J           = computeBLHtoECEFJacobian(ins_pos_feedback_blh); // 计算雅可比矩阵J
    // if (range > 1e-6) {
    //     // 计算H矩阵，反映位置误差对测距的影响,
    //     H.block(0, P_ID, 1, 3) = diff.transpose() / range * J * Dr_inv * (-1);
    // }
    H.block(0, P_ID, 1, 3) = diff.transpose() / range * J * Dr_inv * (-1);
    // 激光测距的测量噪声矩阵 R（假设为常数噪声，具体根据实际情况调整）
    Eigen::MatrixXd R = Eigen::MatrixXd::Identity(1, 1) * 0.1; // 假设测量噪声为0.1米

    // 执行EKF更新
    EKFUpdate(dz, H, R);
}

void GIEngine::EKFUpdate(Eigen::VectorXd &dz, Eigen::MatrixXd &H, Eigen::MatrixXd &R) {

    assert(H.cols() == Cov_.rows());
    assert(dz.rows() == H.rows());
    assert(dz.rows() == R.rows());
    assert(dz.cols() == 1);

    // 计算Kalman增益
    // Compute Kalman Gain
    auto temp         = H * Cov_ * H.transpose() + R;
    Eigen::MatrixXd K = Cov_ * H.transpose() * temp.inverse();

    // 更新系统误差状态和协方差
    // update system error state and covariance
    Eigen::MatrixXd I;
    I.resizeLike(Cov_);
    I.setIdentity();
    I = I - K * H;
    // 如果每次更新后都进行状态反馈，则更新前dx_一直为0，下式可以简化为：dx_ = K * dz;
    // if state feedback is performed after every update, dx_ is always zero before the update
    // the following formula can be simplified as : dx_ = K * dz;
    dx_  = dx_ + K * dz;
    Cov_ = I * Cov_ * I.transpose() + K * R * K.transpose();
    // Cov_ = (Eigen::MatrixXd::Identity(Cov_.rows(), Cov_.rows()) - K * H) * Cov_; // 更新协方差矩阵 (GPT写的)
}

void GIEngine::stateFeedback() { // 误差反馈
    // dx_：位置误差(n系) 速度误差(n系) 等效旋转矢量(相对于n系) 陀螺零偏 加表零偏 陀螺比例因子 加表比例因子

    Eigen::Vector3d vectemp;

    // 位置误差反馈
    // posisiton error feedback
    Eigen::Vector3d delta_r = dx_.block(P_ID, 0, 3, 1);
    Eigen::Matrix3d Dr_inv  = Earth::DRi(pvacur_.pos); // 计算从当前 n系 到 地理坐标系 的变换矩阵
    pvacur_.pos -= Dr_inv * delta_r;                   // 把误差转换到地理系 然后再减

    // 速度误差反馈
    // velocity error feedback
    vectemp = dx_.block(V_ID, 0, 3, 1);
    pvacur_.vel -= vectemp;

    // 姿态误差反馈
    // attitude error feedback
    vectemp                = dx_.block(PHI_ID, 0, 3, 1);           // 等效旋转矢量
    Eigen::Quaterniond qpn = Rotation::rotvec2quaternion(vectemp); // 转换成四元数
    pvacur_.att.qbn        = qpn * pvacur_.att.qbn;
    pvacur_.att.cbn        = Rotation::quaternion2matrix(pvacur_.att.qbn);
    pvacur_.att.euler      = Rotation::matrix2euler(pvacur_.att.cbn);

    // IMU零偏误差反馈
    // IMU bias error feedback
    vectemp = dx_.block(BG_ID, 0, 3, 1);
    imuerror_.gyrbias += vectemp;
    vectemp = dx_.block(BA_ID, 0, 3, 1);
    imuerror_.accbias += vectemp;

    // IMU比例因子误差反馈
    // IMU sacle error feedback
    vectemp = dx_.block(SG_ID, 0, 3, 1);
    imuerror_.gyrscale += vectemp;
    vectemp = dx_.block(SA_ID, 0, 3, 1);
    imuerror_.accscale += vectemp;

    // 误差状态反馈到系统状态后,将误差状态清零
    // set 'dx' to zero after feedback error state to system state
    dx_.setZero();
}

NavState GIEngine::getNavState() {

    NavState state;

    state.pos      = pvacur_.pos;
    state.vel      = pvacur_.vel;
    state.euler    = pvacur_.att.euler;
    state.imuerror = imuerror_;

    return state;
}
