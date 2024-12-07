#include "satellite.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <regex>
// #include <boost/algorithm/string.hpp>
#define M_PI 3.14159265358979323846
// WGS84椭球常数
const double a  = 6378137.0;        // 半长轴（米）
const double e2 = 6.69437999014e-3; // 偏心率平方
py::scoped_interpreter guard{};
satellite::satellite(const std::string &satellite_id, const std::string &tle_file_path)
    : satellite_id_(satellite_id)
    , tle_file_path_(tle_file_path) {
}

// 将输入时间转换为 Skyfield 的 UTC 格式
py::object satellite::convert_to_skyfield_utc(const std::string &time_str, py::object &ts) {
    int year  = 2024;
    int month = 3;
    int day   = 16;
    int hour  = 4; // 4 AM

    // 使用正则表达式分离时间的分钟、秒和毫秒部分
    std::regex time_regex("(\\d+):(\\d+)\\.(\\d+)");
    std::smatch match;
    if (std::regex_search(time_str, match, time_regex)) {
        int minutes       = std::stoi(match[1].str());
        int seconds       = std::stoi(match[2].str());
        double subseconds = std::stod("0." + match[3].str());

        // 调用 Skyfield 的 ts.utc 函数
        return ts.attr("utc")(year, month, day, hour, minutes, seconds + subseconds);
    }

    throw std::runtime_error("Invalid time format");
}

// LLA（经纬度和高度）转换为ECEF坐标
void satellite::LLA_to_ecef(double latitude_deg, double longitude_deg, double altitude_m, double &x, double &y,
                            double &z) {
    const double a  = 6378137.0;        // 半长轴（米）
    const double e2 = 6.69437999014e-3; // 偏心率平方
    double lat_rad  = latitude_deg * M_PI / 180.0;
    double lon_rad  = longitude_deg * M_PI / 180.0;

    double N = a / sqrt(1 - e2 * sin(lat_rad) * sin(lat_rad));

    x = (N + altitude_m) * cos(lat_rad) * cos(lon_rad);
    y = (N + altitude_m) * cos(lat_rad) * sin(lon_rad);
    z = ((1 - e2) * N + altitude_m) * sin(lat_rad);
}

// 根据卫星名称找TLE数据
bool satellite::find_satellite_tle(std::string &line1, std::string &line2) {
    // 从卫星名称中提取数字ID
    std::regex re("\\d+");
    std::smatch match;
    if (!std::regex_search(satellite_id_, match, re)) {
        std::cerr << "Invalid satellite ID." << std::endl;
        return false;
    }

    std::string satellite_number = match.str();

    // 打开TLE数据文件并读取行
    std::ifstream file(tle_file_path_);
    if (!file.is_open()) {
        std::cerr << "Unable to open TLE file." << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {

        std::string sub = line.substr(2, 5);
        if (sub == satellite_id_.substr(4, 5)) {
            line1 = line;
            if (std::getline(file, line2)) {
                return true;
            } else {
                std::cerr << "TLE data incomplete." << std::endl;
                return false;
            }
        }
    }

    std::cerr << "satellite TLE not found." << std::endl;
    return false;
}

// 根据卫星ID和当前时间计算卫星ECEF坐标
bool satellite::get_position(const std::string &time_input, Eigen::Vector3d &sat_position) {
    // py::scoped_interpreter guard{}; // 启动Python解释器

    py::object load  = py::module::import("skyfield.api").attr("load");
    py::object ts    = load.attr("timescale")();
    py::object wgs84 = py::module::import("skyfield.api").attr("wgs84");
    std::string line1;
    std::string line2;
    find_satellite_tle(line1, line2);

    py::object t = convert_to_skyfield_utc(time_input, ts);

    py::object satellite_module = py::module::import("skyfield.api");
    py::object satellite_class  = satellite_module.attr("EarthSatellite");
    py::object satellite        = satellite_class(line1, line2);

    py::object geocentric = satellite.attr("at")(t);
    py::object subpoint   = wgs84.attr("subpoint")(geocentric);

    double latitude  = subpoint.attr("latitude").attr("degrees").cast<double>();
    double longitude = subpoint.attr("longitude").attr("degrees").cast<double>();
    double elevation = subpoint.attr("elevation").attr("m").cast<double>();

    double x, y, z;
    LLA_to_ecef(latitude, longitude, elevation, x, y, z);
    sat_position[0] = x;
    sat_position[1] = y;
    sat_position[2] = z;

    return 1;
}

// int main() {
//     std::string time_input = "15:46.75";
//     std::string satellite_id = "tle-60060";
//     std::string tle_file_path = "C:\\Users\\yuanchuang\\Desktop\\all_tle_two_line.txt";  //
//     std::string line1, line2;
//
//     satellite sat(satellite_id, tle_file_path);
//     sat.find_satellite_tle(line1, line2);
//     std::cout << line1 << std::endl;
//     std::cout << line2 << std::endl;
//     std::vector<double> sat_pos;
//     sat.get_position(time_input, sat_pos);
//     for (auto ele : sat_pos) {
//         std::cout << ele << " ";
//     }
//
//     return 0;
// }