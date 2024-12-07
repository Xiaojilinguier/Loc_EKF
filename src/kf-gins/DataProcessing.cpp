#define _CRT_SECURE_NO_WARNINGS
#include "DataProcessing.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

DataProcessing::DataProcessing(const std::string &real_traj_filename, const std::string &satellite_coverage_file)
    : real_traj_filename(real_traj_filename)
    , satellite_coverage_file(satellite_coverage_file) {
}

// 将时间字符串 'MM:SS.00' 转换为总秒数
double DataProcessing::parse_time(const std::string &time_str) {
    double minutes;
    double seconds;
    sscanf(time_str.c_str(), "%lf:%lf", &minutes, &seconds);
    return minutes * 60 + seconds;
}

std::string DataProcessing::time_to_mmss(const std::string &target_time) {
    int minute     = std::__cxx11::stoi(target_time) / 60;
    double seconds = std::__cxx11::stod(target_time) - 60 * minute;
    std::string now_time;
    if (minute < 10) {
        now_time += "0";
    }
    now_time += std::to_string(minute);
    now_time += ":";
    if (seconds < 10) {
        now_time += "0";
    }
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << seconds;
    std::string result = oss.str();
    now_time += result;
    return now_time;
}

// 实现函数
bool DataProcessing::load_position_velocity(const double &start_time, const double &end_time,
                                            std::vector<std::vector<double>> &RealPositions,
                                            std::vector<std::vector<double>> &RealVelocities,
                                            const std::string &truthpath) {

    double duration = end_time - start_time;
    // 打开文件进行读取
    std::ifstream file(truthpath);
    if (!file.is_open()) {
        std::cerr << "Failed to open the real traj file: " << truthpath << std::endl;
        return false;
    }

    // 读取文件中的每一行
    std::string line;
    double current_time = 0.0;
    int line_number     = 0; // 行号，用于错误报告

    while (std::getline(file, line)) {
        line_number++;
        if (line.empty()) {
            // 跳过空行
            continue;
        }

        std::stringstream ss(line);

        // 定义变量存储每一列的数据
        double gnss_week, time_s, lat, lon, height;
        double vel_n, vel_e, vel_d;
        double roll, pitch, yaw;

        // 使用提取运算符读取每一列数据
        if (!(ss >> gnss_week >> time_s >> lat >> lon >> height >> vel_n >> vel_e >> vel_d >> roll >> pitch >> yaw)) {
            std::cerr << "Error parsing line " << line_number << ": " << line << std::endl;
            continue; // 跳过该行，继续读取下一行
        }

        // 判断时间是否在加载时间跨度内
        if (time_s >= start_time && time_s <= start_time + duration) {
            // 存储位置数据
            RealPositions.emplace_back(std::vector<double>{lat, lon, height});

            // 存储速度数据
            RealVelocities.emplace_back(std::vector<double>{vel_n, vel_e, vel_d});

            // 更新当前时间（以供后续检查时间跨度）
            current_time = time_s;
        }
    }

    file.close();
    return true;
}

// // 读取真实位置和速度 并转换到ECEF坐标系下
// bool DataProcessing::load_position_velocity(const double &start_time, const double &end_time,
//                                             std::vector<std::vector<double>> &RealPositions,
//                                             std::vector<std::vector<double>> &RealVelocities,
//                                             const std::string &truthpath) {

//     double duration = end_time - start_time;
//     // 打开文件进行读取
//     std::ifstream file(truthpath); // 请替换为文件的路径
//     if (!file.is_open()) {
//         std::cerr << "Failed to open the real traj file!" << std::endl;
//         return false;
//     }

//     // 读取文件中的每一行
//     std::string line;
//     double current_time = 0.0;
//     while (std::getline(file, line)) {
//         std::stringstream ss(line);
//         std::string token;

//         // 读取每行的数据
//         double gnss_week, time_s, lat, lon, height, vel_n, vel_e, vel_d, roll, pitch, yaw;
//         std::getline(ss, token, '\t');
//         gnss_week = std::stod(token); // GNSS周数
//         std::getline(ss, token, '\t');
//         time_s = std::stod(token); // 时间 (秒)
//         std::getline(ss, token, '\t');
//         lat = std::stod(token); // 纬度
//         std::getline(ss, token, '\t');
//         lon = std::stod(token); // 经度
//         std::getline(ss, token, '\t');
//         height = std::stod(token); // 高度
//         std::getline(ss, token, '\t');
//         vel_n = std::stod(token); // 北向速度 (m/s)
//         std::getline(ss, token, '\t');
//         vel_e = std::stod(token); // 东向速度 (m/s)
//         std::getline(ss, token, '\t');
//         vel_d = std::stod(token); // 垂向速度 (m/s)
//         std::getline(ss, token, '\t');
//         roll = std::stod(token); // 横滚角 (deg)
//         std::getline(ss, token, '\t');
//         pitch = std::stod(token); // 俯仰角 (deg)
//         std::getline(ss, token, '\t');
//         yaw = std::stod(token); // 航向角 (deg)

//         // 判断时间是否在加载时间跨度内
//         if (time_s >= start_time && time_s <= start_time + duration) {
//             // 存储位置数据
//             RealPositions.push_back({lat, lon, height});

//             // 存储速度数据
//             RealVelocities.push_back({vel_n, vel_e, vel_d});

//             // 更新当前时间（以供后续检查时间跨度）
//             current_time = time_s;
//         }
//     }

//     file.close();
//     return true;
// }

// // 读取惯导指示的加速度
// bool DataProcessing::load_accel(const std::string &start_time_str, double duration,
//                                 std::vector<std::vector<double>> &accelerations) {
//     std::ifstream file(ins_csv_path);
//     if (!file.is_open()) {
//         std::cerr << "Unable to open acceleration file." << std::endl;
//         return false;
//     }

//     std::string line;
//     std::getline(file, line); // 跳过表头

//     double start_time_seconds = parse_time(start_time_str);
//     double end_time_seconds   = start_time_seconds + duration;

//     while (std::getline(file, line)) {
//         std::stringstream ss(line);
//         std::string time_str;
//         std::getline(ss, time_str, ',');

//         double total_seconds = parse_time(time_str);
//         if (total_seconds >= start_time_seconds && total_seconds < end_time_seconds) {
//             std::vector<double> accel(3);
//             std::string value;
//             // 跳过不需要的列
//             for (int i = 0; i < 12; ++i) {
//                 std::getline(ss, value, ',');
//             }
//             for (int i = 0; i < 3; ++i) {
//                 std::getline(ss, value, ',');
//                 accel[i] = std::stod(value);
//             }
//             accelerations.push_back(accel);
//         }
//     }

//     return true;
// }

// // 获取校正起始时刻飞行器惯导计算出的初始位置和速度
// bool DataProcessing::get_ini_position_velocity_array(const std::string &time_str,
//                                                      std::vector<double> &position_velocity_ini_array) {
//     std::ifstream file(ins_csv_path);
//     if (!file.is_open()) {
//         std::cerr << "Unable to open initial position and velocity file." << std::endl;
//         return false;
//     }

//     std::string line;
//     std::getline(file, line); // 跳过表头

//     while (std::getline(file, line)) {
//         std::stringstream ss(line);
//         std::string current_time_str;
//         std::getline(ss, current_time_str, ',');

//         if (current_time_str == time_str) {
//             std::string value;
//             // 前三列为位置
//             for (int i = 0; i < 3; ++i) {
//                 std::getline(ss, value, ',');
//                 position_velocity_ini_array.push_back(std::stod(value));
//             }
//             // 跳过不需要的列
//             for (int i = 0; i < 3; ++i) {
//                 std::getline(ss, value, ',');
//             }
//             for (int i = 0; i < 3; ++i) {
//                 std::getline(ss, value, ',');
//                 position_velocity_ini_array.push_back(std::stod(value));
//             }
//             return true;
//         }
//     }

//     std::cerr << "Time not found in initial position and velocity file." << std::endl;
//     return false;
// }

void DataProcessing::build_index() {
    if (!time_to_file_pos.empty())
        return;

    std::ifstream file(satellite_coverage_file);
    if (!file.is_open()) {
        std::cerr << "无法打开卫星覆盖文件。" << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line); // 跳过表头

    while (true) {
        std::streampos pos = file.tellg();
        if (!std::getline(file, line))
            break;

        std::stringstream ss(line);
        std::string time_str;
        std::getline(ss, time_str, ',');

        time_to_file_pos[time_str] = pos;
    }

    file.close();
}

void split_and_trim(const std::string &str, char delimiter, std::vector<std::string> &tokens) {
    size_t start = 0;
    size_t end   = str.find(delimiter);

    while (end != std::string::npos) {
        std::string token = str.substr(start, end - start);
        // 去除首尾空格
        token.erase(0, token.find_first_not_of(' '));
        token.erase(token.find_last_not_of(' ') + 1);
        tokens.push_back(token);

        start = end + 1;
        end   = str.find(delimiter, start);
    }

    // 处理最后一个token
    std::string token = str.substr(start);
    token.erase(0, token.find_first_not_of(' '));
    token.erase(token.find_last_not_of(' ') + 1);
    tokens.push_back(token);
}
bool DataProcessing::find_satellites_by_time(const std::string &now_time, std::vector<std::string> &satellite_list) {
    build_index();

    auto it = time_to_file_pos.find(now_time);
    if (it == time_to_file_pos.end()) {
        return false;
    }

    std::ifstream file(satellite_coverage_file);
    if (!file.is_open()) {
        std::cerr << "无法打开卫星覆盖文件。" << std::endl;
        return false;
    }

    file.seekg(it->second);

    std::string line;
    std::getline(file, line);

    // 解析行内容
    std::stringstream ss(line);
    std::string time_str;
    std::getline(ss, time_str, ',');

    // 验证时间匹配
    if (time_str != now_time) {
        return false;
    }

    std::string satellites_str;
    std::getline(ss, satellites_str);
    satellites_str.erase(std::remove(satellites_str.begin(), satellites_str.end(), '"'), satellites_str.end());
    satellites_str.pop_back(); // 去除最后一个换行符
    split_and_trim(satellites_str, ',', satellite_list);

    file.close();

    return true;
}

// int main() {
//     std::string ins_csv_path = "C:\\Users\\yuanchuang\\Desktop\\simulate_res4e-1.csv";
//     std::string real_traj_filename = "C:\\Users\\yuanchuang\\Desktop\\RealTraj.csv";
//     std::string satellite_coverage_file = "C:\\Users\\yuanchuang\\Desktop\\satellite_cover.txt";
//     DataProcessing dataprocessor(ins_csv_path, real_traj_filename, satellite_coverage_file);
//     std::string target_time = "00:00.00";
//     std::vector<double> position_velocity_ini_array;
//     dataprocessor.get_ini_position_velocity_array(target_time, position_velocity_ini_array);
//
//
//     for (auto ele : position_velocity_ini_array) {
//         std::cout << ele << " ";
//     }
//
//     std::vector<std::vector<double>> accel_vec;
//     dataprocessor.load_accel(target_time, 2.0, accel_vec);
//     for (int i = 0; i < accel_vec.size(); i++) {
//         for (int j = 0; j < accel_vec[0].size(); j++) {
//             std::cout << accel_vec[i][j] << " ";
//         }
//         std::cout << std::endl;
//     }
//
//
//     std::vector<std::vector<double>> RealPositions;
//     std::vector<std::vector<double>> RealVelocities;
//     dataprocessor.load_position_velocity(target_time, 2, RealPositions, RealVelocities);
//     for (int i = 0; i < RealPositions.size(); i++) {
//         for (int j = 0; j < RealPositions[0].size(); j++) {
//             std::cout << RealPositions[i][j] << " ";
//         }
//         std::cout << std::endl;
//     }
//     std::cout << "--------------------------------------------------------------------------" << std::endl;
//
//     for (int i = 0; i < RealVelocities.size(); i++) {
//         for (int j = 0; j < RealVelocities[0].size(); j++) {
//             std::cout << RealVelocities[i][j] << " ";
//         }
//         std::cout << std::endl;
//     }
//
//     return 0;
// }