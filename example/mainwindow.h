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
    explicit MainWindow(double _bias, double _scale, double _dead_zone,
                        double _gamma_pos, double _gamma_spd, QWidget *parent = 0);
    double bias_, scale_, dead_zone_;
    double gamma_pos_, gamma_spd_;
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
    const double SAMPLING_PERIOD{.0001};
    //-- speed controller sampling period
    static constexpr double SPD_SP{.0001};

    enum ParamIndex{
        iBm, // Viscous Friction (N m s / rad)
        iJm, // Rotor Intertia (kg m^2)
        iKi, // Torque Constant (N m / A)
        iLa, // Armature Inductance (H)
        iKb, // Back EMF Constant (V s / rad)
        iRa  // Armature Resistance (Ohm)
    };
    //--
    jacl::PhysicalParameter<> bm;
    jacl::PhysicalParameter<> jm;
    jacl::PhysicalParameter<> ki;
    jacl::PhysicalParameter<> la;
    jacl::PhysicalParameter<> kb;
    jacl::PhysicalParameter<> ra;
    arma::vec::fixed<6> weight_;

    using LinearStateSpace =
        jacl::state_space::Linear<double,3, 2, 3,
                         jacl::PhysicalParameter<>,
                         jacl::PhysicalParameter<>,
                         jacl::PhysicalParameter<>,
                         jacl::PhysicalParameter<>,
                         jacl::PhysicalParameter<>,
                         jacl::PhysicalParameter<>>;
    LinearStateSpace ss_;
    jacl::system::Continuous<LinearStateSpace> csys_;
    jacl::Plotter<decltype(csys_)> csys_plt_;    
    jacl::system::ContinuousObserver<LinearStateSpace> cobserver_;
    jacl::Plotter<decltype(cobserver_)> cobserver_plt_;

    //-- Discrete SIMO DC motor
    using SIMO = jacl::state_space::Linear<double,3, 1, 3,
                                    jacl::PhysicalParameter<>,
                                    jacl::PhysicalParameter<>,
                                    jacl::PhysicalParameter<>,
                                    jacl::PhysicalParameter<>,
                                    jacl::PhysicalParameter<>,
                                    jacl::PhysicalParameter<>>;
    SIMO simo_;
    jacl::state_space::Linear<double,3,1,3> dsimo_;
    jacl::system::Discrete<decltype(dsimo_)> dsys_simo_;
    jacl::system::DiscreteObserver<decltype(dsimo_)> dobserver_simo_;
    jacl::Plotter<decltype(dsys_simo_)> dsys_simo_plt_;
    jacl::Plotter<decltype(dobserver_simo_)> dobserver_simo_plt_;
    arma::mat din_;    
    jacl::diagnosis::IFD<decltype(dsys_simo_)> ifd_;
    using SIFD = jacl::diagnosis::SIFD<decltype(dsys_simo_),0>;
    SIFD sifd_;
    void makeFault(arma::vec* _out);
    void setupSIMODCMotor();

    //-- Discrete H-infinity controller
    using MReal = jacl::state_space::Linear<double,9,1,1,
                                         jacl::PhysicalParameter<>,
                                         jacl::PhysicalParameter<>,
                                         jacl::PhysicalParameter<>,
                                         jacl::PhysicalParameter<>,
                                         jacl::PhysicalParameter<>,
                                         jacl::PhysicalParameter<>>;
    MReal m_real_;
    jacl::state_space::Linear<double, MReal::n_states, MReal::n_inputs+3, MReal::n_outputs+2> m_icm_;
    //-- Position controller
    using PosReal = jacl::state_space::Linear<double, 3, 1, 1>;
    PosReal pos_real_;
    jacl::state_space::Linear<PosReal::scalar_t ,PosReal::n_states, PosReal::n_inputs+3, PosReal::n_outputs+2> pos_icm_;
    jacl::system::Discrete<decltype(pos_icm_)> pos_sys_;    
    jacl::state_space::Linear<PosReal::scalar_t, PosReal::n_states, PosReal::n_outputs, PosReal::n_inputs> pos_dctrl_;
    jacl::system::Discrete<decltype(pos_dctrl_)> pos_dctrl_sys_;
    jacl::Plotter<jacl::system::Discrete<decltype(pos_dctrl_)>> pos_dctrl_plt_;    
    //-- Speed controller
    using SpdReal = jacl::state_space::Linear<double, 2, 1, 1>;
    SpdReal spd_real_;
    jacl::state_space::Linear<SpdReal::scalar_t, SpdReal::n_states, SpdReal::n_inputs+3, SpdReal::n_outputs+2> spd_icm_;
    jacl::system::Discrete<decltype(spd_icm_)> spd_sys_;    
    jacl::state_space::Linear<SpdReal::scalar_t, SpdReal::n_states, SpdReal::n_outputs, SpdReal::n_inputs> spd_dctrl_;    
    jacl::system::Discrete<decltype(spd_dctrl_)> spd_dctrl_sys_;
    jacl::Plotter<jacl::system::Discrete<decltype(spd_dctrl_)>> spd_dctrl_plt_;
    void setupDiscreteController();

    //-- Continuous H-infinity controller
    using SISOPos = jacl::state_space::Linear<double,3,1,1>;
    SISOPos siso_pos_;
    //-- Position controller
    using CPosReal = jacl::state_space::Linear<double, 3, 1, 1>;
    CPosReal c_pos_real_;
    jacl::state_space::Linear<CPosReal::scalar_t, CPosReal::n_states, CPosReal::n_inputs+3, CPosReal::n_outputs+2> c_pos_icm_;
    jacl::system::Continuous<decltype(c_pos_icm_)> c_pos_sys_;    
    jacl::state_space::Linear<CPosReal::scalar_t, CPosReal::n_states, CPosReal::n_inputs, CPosReal::n_outputs> c_pos_ctrl_;
    jacl::system::Continuous<decltype(c_pos_ctrl_)> c_pos_ctrl_sys_;
    jacl::Plotter<decltype(c_pos_ctrl_sys_)> c_pos_ctrl_plt_;
    //-- Speed controller
    using CSpdReal = jacl::state_space::Linear<double, 2, 1, 1>;
    CSpdReal c_spd_real_;
    jacl::state_space::Linear<CSpdReal::scalar_t, CSpdReal::n_states, CSpdReal::n_inputs+3, CSpdReal::n_outputs+2> c_spd_icm_;
    jacl::system::Continuous<decltype(c_spd_icm_)> c_spd_sys_;    
    jacl::state_space::Linear<CSpdReal::scalar_t, CSpdReal::n_states, CSpdReal::n_inputs, CSpdReal::n_outputs> c_spd_ctrl_;
    jacl::system::Continuous<decltype(c_spd_ctrl_)> c_spd_ctrl_sys_;
    jacl::Plotter<decltype(c_spd_ctrl_sys_)> c_spd_ctrl_plt_;
    void setupContinuousController();    

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
    jacl::PhysicalParameter<> pg_;
    jacl::PhysicalParameter<> pl_;
    jacl::PhysicalParameter<> pk_;
    jacl::PhysicalParameter<> pm_;
    using NLPendulum = jacl::state_space::NonLinear<double,2,1,1,
                            jacl::PhysicalParameter<>,
                            jacl::PhysicalParameter<>,
                            jacl::PhysicalParameter<>,
                            jacl::PhysicalParameter<>>;
    NLPendulum nlp_;
    using NLPSys = jacl::system::Continuous<NLPendulum>;
    NLPSys nlp_sys_;
    void setupNLP();                        

    arma::mat ref_;
    int control_mode_;        

    static auto test(){
        using cek_t = jacl::linear_algebra::Matrix<3,3>;
        using cek2_t = jacl::linear_algebra::Matrix<3,4>;
        constexpr cek_t cek(1.,2.,3.,4.,5.,6.,7.,8.,9.);
        constexpr cek2_t cek2(1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.);
        constexpr cek2_t cek3(cek*cek2);
        if constexpr(cek(0)){
            
        }
        for(int i(0); i < cek2_t::n_elems; i++){
            std::cout << cek3(i) << "\t" << std::endl;
        }
        for(int i(0); i < cek_t::n_elems; i++){
            std::cout << cek(i) << "\t" << std::endl;
        }     
        return 0;
    }

    auto test2(){
        auto ret = test();
    }

Q_SIGNALS:
    void setDataMonitor(arma::vec _data);
    void setSensorStatus(arma::Col<uint8_t> _status);

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
