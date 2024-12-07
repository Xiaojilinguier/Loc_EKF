#pragma once
#ifndef SATELLITE_H
#define SATELLITE_H
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <pybind11/embed.h> // pybind11 嵌入API
#include <pybind11/pybind11.h>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace py = pybind11;

#define M_PI 3.14159265358979323846

class satellite {
public:
    satellite(const std::string &satellite_id, const std::string &tle_file_path);

    // 根据时间和tle数据获取位置（没有对应的库
    bool get_position(const std::string &time_input, Eigen::Vector3d &sat_position);

    // 根据卫星id 查找对应的tle数据
    bool find_satellite_tle(std::string &line1, std::string &line2);

    // 经纬高和地心地固的转换
    void LLA_to_ecef(double latitude_deg, double longitude_deg, double altitude_m, double &x, double &y, double &z);

    // 转换utc时间
    py::object convert_to_skyfield_utc(const std::string &time_str, py::object &ts);

    static py::scoped_interpreter guard;

private:
    // 卫星id
    std::string satellite_id_;

    // tle文件名
    std::string tle_file_path_;
};
#endif