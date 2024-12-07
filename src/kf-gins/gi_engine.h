/*
 * KF-GINS: An EKF-Based GNSS/INS Integrated Navigation System
 *
 * Copyright (C) 2022 i2Nav Group, Wuhan University
 *
 *     Author : Liqiang Wang
 *    Contact : wlq@whu.edu.cn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef GI_ENGINE_H
#define GI_ENGINE_H
#include "DataProcessing.h"
#include "common/types.h"
#include "satellite.h"
#include <Eigen/Dense>
#include <string>
#include <vector>

#include "kf_gins_types.h"

class GIEngine {

public:
    explicit GIEngine(GINSOptions &options, std::string tlepath, std::string coverpath, std::string truthpath,
                      double starttime, double endtime);

    ~GIEngine() = default;

    // 获取卫星位置
    bool getSatellitePosition(satellite *sat, const std::string &now_time, Eigen::Vector3d &satellite_pos);

    /**
     * @brief 添加新的IMU数据，(不)补偿IMU误差
     *        add new imudata, do (not) compensate imu error
     * @param [in] imu        新的IMU原始数据
     *                        new raw imudata
     * @param [in] compensate 是否补偿IMU误差
     *                        if compensate imu error to new imudata
     * */
    void addImuData(const IMU &imu, bool compensate = false) {

        imupre_ = imucur_;
        imucur_ = imu;

        if (compensate) {
            imuCompensate(imucur_);
        }
    }

    /**
     * @brief 添加新的GNSS数据
     *        add new gnssdata
     * @param [in] gnss 新的GNSS数据
     *                  new gnssdata
     * */
    void addGnssData(const GNSS &gnss) {

        gnssdata_ = gnss;
        // 暂不进行数据有效性检查，GNSS数据默认有效
        // do not check the validity of gnssdata, the gnssdata is valid by default
        gnssdata_.isvalid = true;
    }

    /**
     * @brief 处理新的IMU数据
     *        process new imudata
     * */
    void newImuProcess(int imudatarate);

    /**
     * @brief 内插增量形式的IMU数据到指定时刻
     *        interpolate incremental imudata to given timestamp
     * @param [in]     imu1      前一时刻IMU数据
     *                           the previous imudata
     * @param [in,out] imu2      当前时刻IMU数据
     *                           the current imudata
     * @param [in]     timestamp 给定内插到的时刻
     *                           given interpolate timestamp
     * @param [in,out] midimu    输出内插时刻的IMU数据
     *                           output imudata at given timestamp
     * */
    static void imuInterpolate(const IMU &imu1, IMU &imu2, const double timestamp, IMU &midimu) {

        if (imu1.time > timestamp || imu2.time < timestamp) {
            return;
        }

        double lamda = (timestamp - imu1.time) / (imu2.time - imu1.time);

        midimu.time   = timestamp;
        midimu.dtheta = imu2.dtheta * lamda;
        midimu.dvel   = imu2.dvel * lamda;
        midimu.dt     = timestamp - imu1.time;

        imu2.dtheta = imu2.dtheta - midimu.dtheta;
        imu2.dvel   = imu2.dvel - midimu.dvel;
        imu2.dt     = imu2.dt - midimu.dt;
    }

    /**
     * @brief 获取当前时间
     *        get current time
     * */
    double timestamp() const {
        return timestamp_;
    }

    /**
     * @brief 获取当前IMU状态
     *        get current navigation state
     * */
    NavState getNavState();

    /**
     * @brief 获取当前状态协方差
     *        get current state covariance
     * */
    Eigen::MatrixXd getCovariance() {
        return Cov_;
    }

    // 由于加载进内存的真实航迹数据从starttime到endtime 所以需要用index来索引时间对应的下标
    std::vector<double> getRealPos(int index) {
        return realPositions_[index]; // 纬经高(deg)
    }

    std::vector<double> getRealVel(int index) {
        return realVelocities_[index];
    }

    // 获取估计的位置
    std::vector<double> getEstimatePos() {
        double lat = pvacur_.pos[0] * R2D;
        double lon = pvacur_.pos[1] * R2D;
        double alt = pvacur_.pos[2];
        std::vector<double> res;
        double x, y, z;
        LLA_to_ecef(lat, lon, alt, x, y, z);
        res.push_back(x);
        res.push_back(y);
        res.push_back(z);
        return res;
    }

    // 获取估计的北东地速度
    std::vector<double> getEstimateVel() {
        std::vector<double> res;
        res.push_back(pvacur_.vel[0]);
        res.push_back(pvacur_.vel[1]);
        res.push_back(pvacur_.vel[2]);
        return res;
    }

    void LLA_to_ecef(double latitude_deg, double longitude_deg, double altitude_m, double &x, double &y, double &z) {
        const double a  = 6378137.0; // 半长轴（米）
        const double e  = 8.1819190842622e-2;
        const double e2 = e * e; // 偏心率平方
        double lat_rad  = latitude_deg * M_PI / 180.0;
        double lon_rad  = longitude_deg * M_PI / 180.0;

        double N = a / sqrt(1.0 - e2 * sin(lat_rad) * sin(lat_rad));

        x = (N + altitude_m) * cos(lat_rad) * cos(lon_rad);
        y = (N + altitude_m) * cos(lat_rad) * sin(lon_rad);
        z = ((1.0 - e2) * N + altitude_m) * sin(lat_rad);
    }

    // 计算飞行器与卫星之间的测距
    double computeLaserRange(const Eigen::Vector3d &aircraft_pos, const Eigen::Vector3d &satellite_pos);

    void compute_and_print_error(int timeIndex);

    double last_laser_time_;

private:
    /**
     * @brief 初始化系统状态和协方差
     *        initialize state and state covariance
     * @param [in] initstate     初始状态
     *                           initial state
     * @param [in] initstate_std 初始状态标准差
     *                           initial state std
     * */
    void initialize(const NavState &initstate, const NavState &initstate_std);

    /**
     * @brief 当前IMU误差补偿到IMU数据中
     *        componsate imu error to the imudata
     * @param [in,out] imu 需要补偿的IMU数据
     *                     imudata to be compensated
     * */
    void imuCompensate(IMU &imu);

    /**
     * @brief 判断是否需要更新,以及更新哪一时刻系统状态
     *        determine if we should do upate and which navstate to update
     * @param [in] imutime1   上一IMU状态时间
     *                        the last state time
     * @param [in] imutime2   当前IMU状态时间
     *                        the current state time
     * @param [in] updatetime 状态更新的时间
     *                        time to update state
     * @return 0: 不需要更新
     *            donot need update
     *         1: 需要更新上一IMU状态
     *            update the last navstate
     *         2: 需要更新当前IMU状态
     *            update the current navstate
     *         3: 需要将IMU进行内插到状态更新时间
     *            need interpolate imudata to updatetime
     * */
    int isToUpdate(double imutime1, double imutime2, double updatetime) const;

    /**
     * @brief 进行INS状态更新(IMU机械编排算法), 并计算IMU状态转移矩阵和噪声阵
     *        do INS state update(INS mechanization), and compute state transition matrix and noise matrix
     * @param [in,out] imupre 前一时刻IMU数据
     *                        imudata at the previous epoch
     * @param [in,out] imucur 当前时刻IMU数据
     *                        imudata at the current epoch
     * */
    void insPropagation(IMU &imupre, IMU &imucur);

    /**
     * @brief 使用GNSS位置观测更新系统状态
     *        update state using gnss position
     * @param [in,out] gnssdata
     * */
    void gnssUpdate(GNSS &gnssdata);

    /**
     * @brief Kalman 预测,
     *        Kalman Filter Predict process
     * @param [in,out] Phi 状态转移矩阵
     *                     state transition matrix
     * @param [in,out] Qd  传播噪声矩阵
     *                     propagation noise matrix
     * */
    void EKFPredict(Eigen::MatrixXd &Phi, Eigen::MatrixXd &Qd);

    /**
     * @brief Kalman 更新
     *        Kalman Filter Update process
     * @param [in] dz 观测新息
     *                measurement innovation
     * @param [in] H  观测矩阵
     *                measurement matrix
     * @param [in] R  观测噪声阵
     *                measurement noise matrix
     * */
    void EKFUpdate(Eigen::VectorXd &dz, Eigen::MatrixXd &H, Eigen::MatrixXd &R);

    /**
     * @brief 反馈误差状态到当前状态
     *        feedback error state to the current state
     * */
    void stateFeedback();

    /**
     * @brief 检查协方差对角线元素是否都为正
     *        Check if covariance diagonal elements are all positive
     * */
    void checkCov() {

        for (int i = 0; i < RANK; i++) {
            if (Cov_(i, i) < 0) {
                std::cout << "Covariance is negative at " << std::setprecision(10) << timestamp_ << " !" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }

    // 激光测距观测更新
    void LaserUpdate(const double &z, const Eigen::Vector3d &satellite_pos_ned);

    // BLH to ECEF Jacobian
    Matrix3d computeBLHtoECEFJacobian(const Vector3d &blh);

private:
    double starttime_;
    std::string tlepath_;
    std::string coverpath_;
    std::string truthpath_;
    GINSOptions options_;
    DataProcessing *dataProcesser_;
    double timestamp_;
    std::vector<std::vector<double>> realPositions_;  // 飞行器真实位置
    std::vector<std::vector<double>> realVelocities_; // 飞行器真实速度
    int choose_sat_index_;                            // 选星下标
    // 更新时间对齐误差，IMU状态和观测信息误差小于它则认为两者对齐
    // updata time align error
    const double TIME_ALIGN_ERR = 0.001;

    // IMU和GNSS原始数据
    // raw imudata and gnssdata
    IMU imupre_;
    IMU imucur_;
    GNSS gnssdata_;

    // IMU状态（位置、速度、姿态和IMU误差）
    // imu state (position, velocity, attitude and imu error)
    PVA pvacur_;
    PVA pvapre_;
    ImuError imuerror_;

    // Kalman滤波相关
    // ekf variables
    Eigen::MatrixXd Cov_;
    Eigen::MatrixXd Qc_;
    Eigen::MatrixXd dx_;

    const int RANK      = 21;
    const int NOISERANK = 18;

    // 状态ID和噪声ID （位置速度误差为n系 北东地
    // state ID and noise ID
    enum StateID { P_ID = 0, V_ID = 3, PHI_ID = 6, BG_ID = 9, BA_ID = 12, SG_ID = 15, SA_ID = 18 };
    enum NoiseID { VRW_ID = 0, ARW_ID = 3, BGSTD_ID = 6, BASTD_ID = 9, SGSTD_ID = 12, SASTD_ID = 15 };
};

#endif // GI_ENGINE_H
