#pragma once
#ifndef DATAPROCESSING_H
#define DATAPROCESSING_H
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
class DataProcessing {
public:
    // 赋值构造函数（赋值惯导文件、真实航迹、卫星覆盖文件）
    DataProcessing(const std::string &real_traj_filename, const std::string &satellite_coverage_file);

    // 根据开始时间戳与持续时间 读取真实位置和速度
    bool load_position_velocity(const double &start_time, const double &end_time,
                                std::vector<std::vector<double>> &RealPositions,
                                std::vector<std::vector<double>> &RealVelocities, const std::string &truthpath);

    // // 读取惯导指示的加速度
    // bool load_accel(const std::string &start_time_str, double duration,
    //                 std::vector<std::vector<double>> &accelerations);

    // // 获取校正起始时刻飞行器惯导计算出的初始位置和速度
    // bool get_ini_position_velocity_array(const std::string &time_str, std::vector<double>
    // &position_velocity_ini_array);

    // 根据时间从覆盖分析文件中找到能覆盖到的卫星列表
    bool find_satellites_by_time(const std::string &target_time, std::vector<std::string> &satellite_list);

    // 建立时间到文件位置的索引
    void build_index();

    // 将时间字符串 'MM:SS.00' 转换为总秒数
    static double parse_time(const std::string &time_str);

    // 将秒数转变成'MM:SS.00'
    std::string time_to_mmss(const std::string &target_time);

private:
    // std::string ins_csv_path;                                         // 惯导文件
    std::string real_traj_filename;                                   // 真实航迹文件
    std::string satellite_coverage_file;                              // 卫星覆盖文件
    std::unordered_map<std::string, std::streampos> time_to_file_pos; // 时间到文件位置的索引
};
#endif // DATAPROCESSING_H