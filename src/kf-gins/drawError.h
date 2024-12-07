#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <matplotlibcpp.h>
#include <sstream>
#include <string>
#include <vector>
namespace plt = matplotlibcpp;

// 结构体存储导航数据
struct NavData {
    double gnss_week;
    double seconds;
    double latitude;
    double longitude;
    double altitude;
    double vel_n, vel_e, vel_d;
    double roll, pitch, yaw;
};

// DrawError 类
class DrawError {
public:
    DrawError(const std::string &estFile, const std::string &truthFile)
        : est_file(estFile)
        , truth_file(truthFile) {
    }

    // 读取文件并解析数据
    bool loadData();

    // 计算位置、速度和姿态的误差
    void computeErrors();

    // 绘制误差曲线
    void plotErrors(int time_step);

private:
    std::string est_file;
    std::string truth_file;
    std::map<double, NavData> est_data_;
    std::map<double, NavData> truth_data_;
    std::vector<double> time;
    std::vector<double> pos_errors_x; // Position error vector
    std::vector<double> pos_errors_y;
    std::vector<double> pos_errors_z;
    std::vector<double> vel_errors_x; // Velocity error vector
    std::vector<double> vel_errors_y;
    std::vector<double> vel_errors_z;
    std::vector<double> att_errors_x; // Attitude error vector
    std::vector<double> att_errors_y;
    std::vector<double> att_errors_z;
    // 加载导航数据
    std::map<double, NavData> loadNavData(const std::string &filename);

    // 解析导航数据
    NavData parseNavLine(const std::string &line);
    // 经纬度转换为笛卡尔坐标
    void convertToCartesian(const NavData &data, double &x, double &y, double &z);
};

// int main() {
//     // 创建 DrawError 对象并加载数据
//     DrawError drawError("/home/sxh/KF-GINS-main/dataset/KF_GINS_Navresult.nav",
//                         "/home/sxh/KF-GINS-main/dataset/truth.nav");

//     if (!drawError.loadData()) {
//         std::cerr << "Error loading data files." << std::endl;
//         return -1;
//     }

//     // 计算误差并绘制误差曲线
//     drawError.computeErrors();
//     drawError.plotErrors();

//     return 0;
// }
