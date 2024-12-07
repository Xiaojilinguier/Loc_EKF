#include "drawError.h"

bool DrawError::loadData() {
    this->est_data_   = loadNavData(est_file);
    this->truth_data_ = loadNavData(truth_file);
    if (this->est_data_.empty()) {
        std::cout << "导航结果文件读取失败" << std::endl;
        return 0;
    }
    if (this->truth_data_.empty()) {
        std::cout << "真实航迹文件读取失败" << std::endl;
        return 0;
    }
    return !this->est_data_.empty() && !this->truth_data_.empty();
}

void DrawError::computeErrors() {
    for (auto &pair : this->est_data_) {
        const auto &key    = pair.first;
        const NavData &est = pair.second;

        if (this->truth_data_.find(key) != this->truth_data_.end()) {
            const NavData &truth = truth_data_.at(key);

            double est_x, est_y, est_z, truth_x, truth_y, truth_z;
            convertToCartesian(est, est_x, est_y, est_z);
            convertToCartesian(truth, truth_x, truth_y, truth_z);

            pos_errors_x.push_back(fabs(est_x - truth_x));
            pos_errors_y.push_back(fabs(est_y - truth_y));
            pos_errors_z.push_back(fabs(est_z - truth_z));

            vel_errors_x.push_back(fabs(est.vel_n - truth.vel_n));
            vel_errors_y.push_back(fabs(est.vel_e - truth.vel_e));
            vel_errors_z.push_back(fabs(est.vel_d - truth.vel_d));

            att_errors_x.push_back(fabs(est.roll - truth.roll));
            att_errors_y.push_back(fabs(est.pitch - truth.pitch));
            att_errors_z.push_back(fabs(est.yaw - truth.yaw));

            time.push_back(est.seconds);
        }
    }
}

// 绘制误差曲线
void DrawError::plotErrors(int time_step) {
    // 采样间隔为time_step

    // 创建新容器存储采样后的数据
    std::vector<double> sampled_time;
    std::vector<double> sampled_pos_errors_x, sampled_pos_errors_y, sampled_pos_errors_z;
    std::vector<double> sampled_vel_errors_x, sampled_vel_errors_y, sampled_vel_errors_z;
    std::vector<double> sampled_att_errors_x, sampled_att_errors_y, sampled_att_errors_z;
    int j = 0;
    // 遍历原始数据，根据时间间隔取样
    for (size_t i = 0; i < time.size(); ++i) {
        if (i == 0 || time[i] - time[j] >= time_step) {
            j = i;
            sampled_time.push_back(time[i]);
            sampled_pos_errors_x.push_back(pos_errors_x[i]);
            sampled_pos_errors_y.push_back(pos_errors_y[i]);
            sampled_pos_errors_z.push_back(pos_errors_z[i]);
            sampled_vel_errors_x.push_back(vel_errors_x[i]);
            sampled_vel_errors_y.push_back(vel_errors_y[i]);
            sampled_vel_errors_z.push_back(vel_errors_z[i]);
            sampled_att_errors_x.push_back(att_errors_x[i]);
            sampled_att_errors_y.push_back(att_errors_y[i]);
            sampled_att_errors_z.push_back(att_errors_z[i]);
        }
    }

    // 绘制位置误差
    plt::figure();
    plt::named_plot("Position Error X", sampled_time, sampled_pos_errors_x, "r-");
    plt::named_plot("Position Error Y", sampled_time, sampled_pos_errors_y, "g-");
    plt::named_plot("Position Error Z", sampled_time, sampled_pos_errors_z, "b-");
    plt::xlabel("Time (s)");
    plt::ylabel("Position Error (m)");
    plt::legend();
    plt::save("/home/sxh/Loc_EKF/output/position_error_xyz.png");

    // 绘制速度误差
    plt::figure();
    plt::named_plot("Velocity Error X", sampled_time, sampled_vel_errors_x, "r-");
    plt::named_plot("Velocity Error Y", sampled_time, sampled_vel_errors_y, "g-");
    plt::named_plot("Velocity Error Z", sampled_time, sampled_vel_errors_z, "b-");
    plt::xlabel("Time (s)");
    plt::ylabel("Velocity Error (m/s)");
    plt::legend();
    plt::save("/home/sxh/Loc_EKF/output/velocity_error_xyz.png");

    // 绘制姿态误差
    plt::figure();
    plt::named_plot("Attitude Error X (Roll)", sampled_time, sampled_att_errors_x, "r-");
    plt::named_plot("Attitude Error Y (Pitch)", sampled_time, sampled_att_errors_y, "g-");
    plt::named_plot("Attitude Error Z (Yaw)", sampled_time, sampled_att_errors_z, "b-");
    plt::xlabel("Time (s)");
    plt::ylabel("Attitude Error (deg)");
    plt::legend();
    plt::save("/home/sxh/Loc_EKF/output/attitude_error_xyz.png");
}

double roundToThreeDecimalPlaces(double number) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(3) << number;

    // 获取格式化后的字符串
    std::string roundedStr = ss.str();

    // 将字符串转换回 double
    double roundedNum = std::stod(roundedStr);

    // 返回格式化后的字符串和对应的 double 值
    return roundedNum;
}

std::map<double, NavData> DrawError::loadNavData(const std::string &filename) {
    std::ifstream file(filename);
    std::map<double, NavData> nav_data;
    std::string line;
    while (std::getline(file, line)) {
        NavData data           = parseNavLine(line);
        data.seconds           = roundToThreeDecimalPlaces(data.seconds);
        nav_data[data.seconds] = data;
    }
    return nav_data;
}

// 解析导航数据
NavData DrawError::parseNavLine(const std::string &line) {
    std::stringstream ss(line);
    NavData data;
    ss >> data.gnss_week >> data.seconds >> data.latitude >> data.longitude >> data.altitude >> data.vel_n >>
        data.vel_e >> data.vel_d >> data.roll >> data.pitch >> data.yaw;
    return data;
}

// 经纬度转换为笛卡尔坐标
void DrawError::convertToCartesian(const NavData &data, double &x, double &y, double &z) {
    const double a  = 6378137.0;           // WGS-84 长半轴
    const double f  = 1.0 / 298.257223563; // 扁率
    const double e2 = f * (2 - f);         // 偏心率的平方

    double lat_rad = data.latitude * M_PI / 180.0;
    double lon_rad = data.longitude * M_PI / 180.0;

    double N = a / sqrt(1 - e2 * sin(lat_rad) * sin(lat_rad));
    x        = (N + data.altitude) * cos(lat_rad) * cos(lon_rad);
    y        = (N + data.altitude) * cos(lat_rad) * sin(lon_rad);
    z        = (N * (1 - e2) + data.altitude) * sin(lat_rad);
}