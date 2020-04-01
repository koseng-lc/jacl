/**
*   @author : koseng (Lintang)
*   @brief : GUI Header
*/

#pragma once

#include <boost/thread.hpp>

#include <QMainWindow>
#include <QGroupBox>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QGridLayout>
#include <QPushButton>
#include <QImage>
#include <QImageReader>
#include <QDial>
#include <QFile>
#include <QComboBox>
#include <QMenu>
#include <QDialog>

#include <jacl/jacl.hpp>
#include "controller_dialog.h"
#include "controller_sim.h"
#include "system_sim.h"
#include "observer_sim.h"
#include "posctrl_sim.h"
#include "spdctrl_sim.h"

//-- TODO : Add GroupBox for Control Input

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
private:
    void setupWidgets();
    void setupControllerDialog();
    void setupMenus();
    void setupActions();

    void closedLoopProcess();
    double angularSpeed2Voltage(double _speed, double _torque);
private:
    //-- GUI stuff

    Ui::MainWindow *ui;

    QWidget* main_widget_;
    QGridLayout* main_layout_;

    //-- Params
    QGroupBox* params_gb_;
    QGridLayout* params_gl_;
    QLabel* params_label_[6];
    QDoubleSpinBox* params_dsb_[6];

    //-- Command
    QGroupBox* command_gb_;
    QGridLayout* command_gl_;
    QPushButton* perturb_pb_;
    QPushButton* reset_pb_;
    QPushButton* simulate_pb_;
    QPushButton* set_input_pb_;

    //-- Input
    QGroupBox* input_gb_;
    QGridLayout* input_gl_;
    QLabel* torque_in_label_;
    QDoubleSpinBox* torque_in_dsb_;
    QLabel* voltage_in_label_;
    QDoubleSpinBox* voltage_in_dsb_;

    //-- Fault
    QGroupBox* fault_gb_;
    QGridLayout* fault_gl_;
    QLabel* target_label_;
    QComboBox* target_cb_;
    QPushButton* details_pb_;
    QDial* bias_dial_;
    QDoubleSpinBox* bias_dsb_;
    QLabel* bias_label_;
    QDial* scale_dial_;
    QDoubleSpinBox* scale_dsb_;
    QLabel* scale_label_;
    QDial* dead_zone_dial_;
    QDoubleSpinBox* dead_zone_dsb_;
    QLabel* dead_zone_label_;

    //-- Watermark
    QWidget* watermark_widget_;
    QGridLayout* watermark_gl_;
    QImage* logo_image_;
    QLabel* logo_image_label_;
    QLabel* title_label_;
    QLabel* author_label_;

    //-- Menu stuff
    QMenu* tools_menu_;
    QAction* controller_act_;

    //-- Controller Dialog
    jacl::ControllerDialog* controller_dialog_;
//    QDialog* controller_dialog_;
    QGridLayout* ctrl_gl_;
    QGroupBox* ctrl_in_gb_;    

private:
    //-- Control stuff

    //-- Closed-loop process
    boost::thread cl_thread_;
    boost::mutex cl_mtx_;
    std::atomic<bool> cl_status_;
    const double TIME_STEP{.02};

    enum ParamIndex{
        iBm, // Viscous Friction (N m s / rad)
        iJm, // Rotor Intertia (kg m^2)
        iKi, // Torque Constant (N m / A)
        iLa, // Armature Inductance (H)
        iKb, // Back EMF Constant (V s / rad)
        iRa  // Armature Resistance (Ohm)
    };
    //--
    jacl::PhysicalParameter bm;
    jacl::PhysicalParameter jm;
    jacl::PhysicalParameter ki;
    jacl::PhysicalParameter la;
    jacl::PhysicalParameter kb;
    jacl::PhysicalParameter ra;

    using LinearStateSpace =
        jacl::LinearStateSpace<3, 2, 3,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter>;

    LinearStateSpace ss_;
    jacl::SystemSim<LinearStateSpace> system_sim_;
    jacl::ObserverSim<LinearStateSpace> observer_sim_;
    //--
    using SIMO = jacl::LinearStateSpace<3, 1, 3,
                                    jacl::PhysicalParameter,
                                    jacl::PhysicalParameter,
                                    jacl::PhysicalParameter,
                                    jacl::PhysicalParameter,
                                    jacl::PhysicalParameter,
                                    jacl::PhysicalParameter>;
    SIMO simo_;

    //--
    using GRealization =
        jacl::LinearStateSpace<3, 8, 9,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter>;

    GRealization G_;

    using PRealization = jacl::LinearStateSpace<9, 2, 3>;
    PRealization P_;

    using InterConnMat = jacl::LinearStateSpace<9, 10, 8>;
    InterConnMat ICM_;

    using Controller = jacl::LinearStateSpace<9, 3, 2>;
    Controller K_;
    jacl::ControllerSim<Controller> controller_sim_;

    using HInf = jacl::synthesis::HInf<InterConnMat, 5, 8>;
    HInf* hinf_;

    //-- H-infinity position control of dc motor
    using SISOPos = jacl::LinearStateSpace<3,1,1>;
    SISOPos siso_pos_;
    using InterConnMatPos = jacl::LinearStateSpace<6, 4, 3>;
    InterConnMatPos icm_pos_;
    using PosCtrl = jacl::LinearStateSpace<6, 1, 1>;
    PosCtrl k_pos_;
    using HInfPC = jacl::synthesis::HInf<InterConnMatPos, 2, 3>;
    HInfPC* hinf_pc_;
    jacl::PosCtrlSim<PosCtrl> posctrl_sim_;
    void setupPositionController();

    //-- H-infinity speed control of dc motor
    using InterConnMatSpd = jacl::LinearStateSpace<5,4,3>;
    InterConnMatSpd icm_spd_;
    using SpdCtrl = jacl::LinearStateSpace<5, 1, 1>;
    SpdCtrl k_spd_;
    using HInfSC = jacl::synthesis::HInf<InterConnMatSpd, 2, 3>;
    HInfSC* hinf_sc_;
    jacl::SpdCtrlSim<SpdCtrl> spdctrl_sim_;
    void setupSpeedController();

    arma::mat ref_;
    int control_mode_;        

private Q_SLOTS:
    void perturbAct();
    void resetAct();
    void simulateAct();
    void setInputAct();
    void biasDialConv(double _val);
    void scaleDialConv(double _val);
    void deadZoneDialConv(double _val);
    void biasDSBConv(int _val);
    void scaleDSBConv(int _val);
    void deadZoneDSBConv(int _val);

    void openControllerDialog();
    void refAct();
    void modeAct(int _val);
};