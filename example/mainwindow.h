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

//-- TODO : Add GroupBox for Control Input

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT
    jacl::py_stuff::PyEmbedHandler py_handler_;
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
    QPushButton* fault_pb_;
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
    const double SAMPLING_PERIOD{.001};

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
    arma::vec::fixed<6> weight_;

    using LinearStateSpace =
        jacl::state_space::Linear<double,3, 2, 3,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter>;

    LinearStateSpace ss_;
    jacl::system::Continuous<LinearStateSpace> csys_;
    jacl::Plotter<jacl::system::Continuous<LinearStateSpace> > csys_plt_;    
    jacl::system::ContinuousObserver<LinearStateSpace> cobserver_;
    jacl::Plotter<jacl::system::ContinuousObserver<LinearStateSpace> > cobserver_plt_;

    //-- Discrete SIMO DC motor
    using SIMO = jacl::state_space::Linear<double,3, 1, 3,
                                    jacl::PhysicalParameter,
                                    jacl::PhysicalParameter,
                                    jacl::PhysicalParameter,
                                    jacl::PhysicalParameter,
                                    jacl::PhysicalParameter,
                                    jacl::PhysicalParameter>;
    SIMO simo_;
    jacl::state_space::Linear<double,3,1,3> dsimo_;
    jacl::system::Discrete<jacl::state_space::Linear<double,3,1,3>> dsys_simo_;
    jacl::system::DiscreteObserver<jacl::state_space::Linear<double,3,1,3>> dobserver_simo_;
    jacl::Plotter<jacl::system::Discrete<jacl::state_space::Linear<double,3,1,3>>> dsys_simo_plt_;
    jacl::Plotter<jacl::system::DiscreteObserver<jacl::state_space::Linear<double,3,1,3>>> dobserver_simo_plt_;
    arma::mat din_;    
    jacl::diagnosis::IFD<jacl::system::Discrete<jacl::state_space::Linear<double,3,1,3>>> ifd_;
    using SIFD = jacl::diagnosis::SIFD<jacl::system::Discrete<jacl::state_space::Linear<double,3,1,3>>,0>;
    SIFD sifd_;
    void makeFault(arma::vec* _out);
    void setupSIMODCMotor();
    //-- Discrete H-infinity controller
    using MReal = jacl::state_space::Linear<double,9,1,1,
                                         jacl::PhysicalParameter,
                                         jacl::PhysicalParameter,
                                         jacl::PhysicalParameter,
                                         jacl::PhysicalParameter,
                                         jacl::PhysicalParameter,
                                         jacl::PhysicalParameter>;    
    using MICM = jacl::state_space::Linear<double,MReal::n_states, MReal::n_inputs+3, MReal::n_outputs+2>;
    using MSys = jacl::system::Discrete<MICM>;
    MReal m_real_;
    MICM m_icm_;
    MSys m_sys_;
    using PosReal = jacl::state_space::Linear<double,3, 1, 1>;
    using PosICM = jacl::state_space::Linear<double,PosReal::n_states, PosReal::n_inputs+3, PosReal::n_outputs+2>;
    using PosSys = jacl::system::Discrete<PosICM>;
    using PosDHinf = jacl::synthesis::DHinf<PosSys,2,3>;
    using PosDCtrl = jacl::state_space::Linear<double,PosReal::n_states, PosReal::n_outputs, PosReal::n_inputs>;
    PosReal pos_real_;
    PosICM pos_icm_;
    PosSys pos_sys_;
    PosDHinf* pos_dhinf_;
    PosDCtrl pos_dctrl_;
    jacl::system::Discrete<PosDCtrl> pos_dctrl_sys_;
    jacl::Plotter<jacl::system::Discrete<PosDCtrl>> pos_dctrl_plt_;

    //-- Another stuff
    using GRealization =
        jacl::state_space::Linear<double,3, 8, 9,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter,
                         jacl::PhysicalParameter>;
    GRealization G_;
    using PRealization = jacl::state_space::Linear<double,9, 2, 3>;
    PRealization P_;
    using InterConnMat = jacl::state_space::Linear<double,9, 10, 8>;
    InterConnMat ICM_;
    // using HInf = jacl::synthesis::Hinf<InterConnMat, 5, 8>;
    // HInf* hinf_;

    //-- Continuous H-infinity position control of dc motor
    using SISOPos = jacl::state_space::Linear<double,3,1,1>;
    SISOPos siso_pos_;
    using InterConnMatPos = jacl::state_space::Linear<double,6, 4, 3>;
    InterConnMatPos icm_pos_;
    using PosCtrl = jacl::state_space::Linear<double,6, 1, 1>;
    PosCtrl k_pos_;
    using HInfPC = jacl::synthesis::Hinf<jacl::system::Continuous<InterConnMatPos>,2,3>;
    HInfPC* hinf_pc_;
    jacl::system::Continuous<PosCtrl> posctrl_sys_;
    jacl::Plotter<jacl::system::Continuous<PosCtrl>> posctrl_plt_;
    void setupPositionController();

    //-- H-infinity speed control of dc motor
    using InterConnMatSpd = jacl::state_space::Linear<double,5,4,3>;
    InterConnMatSpd icm_spd_;
    using SpdCtrl = jacl::state_space::Linear<double,5,1,1>;
    SpdCtrl k_spd_;
    using HInfSC = jacl::synthesis::Hinf<jacl::system::Continuous<InterConnMatSpd>,2,3>;
    HInfSC* hinf_sc_;
    jacl::system::Continuous<SpdCtrl> spdctrl_sys_;
    jacl::Plotter<jacl::system::Continuous<SpdCtrl>> spdctrl_plt_;
    void setupSpeedController();    

    //-- Non-Linear Pendulum Model
    enum {
        ix1,
        ix2,
        iu1,
        ipg,
        ipl,
        ipk,
        ipm
    };
    jacl::PhysicalParameter pg_;
    jacl::PhysicalParameter pl_;
    jacl::PhysicalParameter pk_;
    jacl::PhysicalParameter pm_;
    using NLPendulum = jacl::state_space::NonLinear<double,2,1,1,
                            jacl::PhysicalParameter,
                            jacl::PhysicalParameter,
                            jacl::PhysicalParameter,
                            jacl::PhysicalParameter>;
    NLPendulum nlp_;
    using NLPSys = jacl::system::Continuous<NLPendulum>;
    NLPSys nlp_sys_;
    void setupNLP();                        

    arma::mat ref_;
    int control_mode_;        

Q_SIGNALS:
    void setDataMonitor(QVector<double> _data);

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
