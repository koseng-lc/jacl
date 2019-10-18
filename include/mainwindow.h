/**
*   @author : koseng (Lintang)
*   @brief : GUI Header
*/

#pragma once

#include <QMainWindow>
#include <QGroupBox>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QGridLayout>
#include <QPushButton>

#include <jacl>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:

    //-- GUI stuff

    Ui::MainWindow *ui;

    QWidget* main_widget_;
    QGridLayout* main_layout_;

    QGroupBox* params_gb_;
    QGridLayout* params_gl_;
    QLabel* params_label_[6];
    QDoubleSpinBox* params_dsb_[6];

    QGroupBox* command_gb_;
    QGridLayout* command_gl_;
    QPushButton* perturb_pb_;
    QPushButton* reset_pb_;

    void setupWidgets();
    void setupActions();

    //-- Control stuff

    enum ParamIndex{
        iBm, // Viscous Friction (N m s / rad)
        iJm, // Rotor Intertia (kg m^2)
        iKi, // Torque Constant (N m / A)
        iLa, // Armature Inductance (H)
        iKb, // Back EMF Constant (V s / rad)
        iRa  // Armature Resistance (Ohm)
    };


    JACL::PhysicalParameter bm;
    JACL::PhysicalParameter jm;
    JACL::PhysicalParameter ki;
    JACL::PhysicalParameter la;
    JACL::PhysicalParameter kb;
    JACL::PhysicalParameter ra;

    using SS = JACL::StateSpace<3, 2, 2,
                                JACL::PhysicalParameter,
                                JACL::PhysicalParameter,
                                JACL::PhysicalParameter,
                                JACL::PhysicalParameter,
                                JACL::PhysicalParameter,
                                JACL::PhysicalParameter>;

    SS ss_;

//    JACL::Simulator* sim_;

};
