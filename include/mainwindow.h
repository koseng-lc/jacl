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

#include <jacl>

#include "controller_dialog.h"

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

    //-- Closed-loop process
    boost::thread cl_thread_;
    boost::mutex cl_mtx_;
    bool cl_status_;
    void closedLoopProcess();
    const double TIME_STEP{.02};

    void setupWidgets();
    void setupControllerDialog();
    void setupMenus();
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

    jacl::PhysicalParameter bm;
    jacl::PhysicalParameter jm;
    jacl::PhysicalParameter ki;
    jacl::PhysicalParameter la;
    jacl::PhysicalParameter kb;
    jacl::PhysicalParameter ra;

    using StateSpace =
        jacl::StateSpace<3, 2, 3,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter>;

    StateSpace ss_;
//    jacl::Simulator sim_;
    jacl::SystemSim<StateSpace> system_sim_;    

    arma::mat::fixed<3, 3> obs_gain_;
    jacl::ObserverSim<StateSpace> observer_sim_;

    using GRealization =
        jacl::StateSpace<3, 8, 9,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter>;

    GRealization G_;

    using PRealization = jacl::StateSpace<9, 2, 3>;
    PRealization P_;

    using InterConnMat = jacl::StateSpace<9, 10, 8>;
    InterConnMat ICM_;

    using Controller = jacl::StateSpace<9, 3, 2>;
    Controller K_;
    jacl::ControllerSim<Controller> controller_sim_;

    using HInf = jacl::synthesis::HInf<InterConnMat, 5, 8>;
    HInf* h_inf_;


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
};
