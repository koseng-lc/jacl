#include "mainwindow.h"
#include <QApplication>

int main(int argc, char** argv){
    auto bias(.0);
    auto scale(1.);
    auto dead_zone(.0);
    auto gamma_pos(5.0);
    auto gamma_spd(5.0);
    if(argc == 6){
        bias = std::stod(argv[1]);
        scale = std::stod(argv[2]);
        dead_zone = std::stod(argv[3]);
        gamma_pos = std::stod(argv[4]);
        gamma_spd = std::stod(argv[5]);
        std::cout << "Setup used : " << std::endl;
        std::cout << "Bias : " << bias << std::endl;
        std::cout << "Scale : " << scale << std::endl;
        std::cout << "Dead-zone : " << dead_zone << std::endl;
        std::cout << "Gamma pos : " << gamma_pos << std::endl;
        std::cout << "Gamma spd : " << gamma_spd << std::endl;
    }else{
        std::cerr << "Number of argument is not enough" << std::endl;
        std::cerr << "Example use -> $ control_sim [bias] [scale] [dead_zone] [gamma_pos] [gamma_spd]" << std::endl;
        exit(1);
    }

    QApplication a(argc, argv);
    MainWindow w(bias,scale,dead_zone,gamma_pos,gamma_spd);
    w.show();

    return a.exec();
}
