/**
*   @author : koseng (Lintang)
*   @brief : GUI Source
*/

#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(double _bias,
                       double _scale,
                       double _dead_zone,
                       double _gamma_pos,
                       double _gamma_spd,
                       QWidget *parent)
    : QMainWindow(parent)
    , bias_(_bias), scale_(_scale), dead_zone_(_dead_zone)
    , gamma_pos_(_gamma_pos), gamma_spd_(_gamma_spd)
    , ui(new Ui::MainWindow)
    //-- DC motor parameters - Maxon 2322.983-11.225-200
    //-- Max speed : 994.8376724999999 rad/s (9500 rpm)
    //-- Nominal voltage : 24 V
    //-- No load speed : 881.7403371 (8420 rpm)
    , bm(0.0436e-3), jm(5.71e-7), ki(30.7e-3), la(1.97e-3), kb(30.7e-3), ra(21.6)
    , weight_({0.1, 0.1, 0.1, 0.1, 0.1, 0.1})
    //-- MIMO DC motor
    , ss_(&bm, &jm, &ki, &la, &kb, &ra)    
    , csys_(&ss_)
    , csys_plt_(&csys_,{1,2,3,4,5,6,7,8})
    , cobserver_(&ss_, arma::zeros<arma::mat>(3,3))
    , cobserver_plt_(&cobserver_, {1,2,3,6,7,8})
    //-- SIMO DC motor
    , simo_(&bm, &jm, &ki, &la, &kb, &ra)
    , dsys_simo_(&dsimo_, SAMPLING_PERIOD)
    , dobserver_simo_(&dsimo_, arma::zeros<arma::mat>(3,3), SAMPLING_PERIOD)
    , dsys_simo_plt_(&dsys_simo_, {1,2,3,4}, SAMPLING_PERIOD)
    , dobserver_simo_plt_(&dobserver_simo_, {1,2,3}, SAMPLING_PERIOD)
    , ifd_(&dsys_simo_, {10.,9.,8.})
    , sifd_(&dsys_simo_, {10., .05, .07})
    , m_real_(&bm,&jm,&ki,&la,&kb,&ra)
    , m_sys_(&m_icm_, SAMPLING_PERIOD)
    //-- discrete position controller
    , pos_sys_(&pos_icm_, SAMPLING_PERIOD)
    , pos_dctrl_sys_(&pos_dctrl_, SAMPLING_PERIOD)
    , pos_dctrl_plt_(&pos_dctrl_sys_, {4,5}, SAMPLING_PERIOD)
    //-- discrete speed controller
    , spd_sys_(&spd_icm_, SPD_SP)
    , spd_dctrl_sys_(&spd_dctrl_, SPD_SP)
    , spd_dctrl_plt_(&spd_dctrl_sys_, {4,5}, SPD_SP)
    , ref_(3, 1, arma::fill::zeros)
    , cl_status_(true)
    , posctrl_sys_(&k_pos_)
    , posctrl_plt_(&posctrl_sys_,{7,8})
    , spdctrl_sys_(&k_spd_)
    , spdctrl_plt_(&spdctrl_sys_,{6,7})
    , pg_(9.8), pl_(1.0), pk_(0.1), pm_(0.5)
    , nlp_(&pg_, &pl_, &pk_, &pm_)
    , nlp_sys_(&nlp_){
    //-- GUI    
    ui->setupUi(this);
    this->removeToolBar(ui->mainToolBar);
    QFile style_file("../example/gui/dark_style_lintang.qss");
    style_file.open(QFile::ReadOnly);
    QString style_sheet(style_file.readAll());
    this->setStyleSheet(style_sheet);

    //-- DC Motor Open-Loop
    std::cout << "Preparing system ..." << std::endl;
    {
        LinearStateSpace::formula_t A11 = JC(ss_,.0); LinearStateSpace::formula_t A12 = JC(ss_,1.0); LinearStateSpace::formula_t A13 = JC(ss_,.0);
        LinearStateSpace::formula_t A21 = JC(ss_,.0); LinearStateSpace::formula_t A22 = JE(ss_,-_p(iBm)/_p(iJm)); LinearStateSpace::formula_t A23 = JE(ss_,_p(iKi)/_p(iJm));
        LinearStateSpace::formula_t A31 = JC(ss_,.0); LinearStateSpace::formula_t A32 = JE(ss_,-_p(iKb)/_p(iLa)); LinearStateSpace::formula_t A33 = JE(ss_,-_p(iRa)/_p(iLa));

        LinearStateSpace::formulas_t fA{
            A11,A12,A13,
            A21,A22,A23,
            A31,A32,A33
        };

        LinearStateSpace::formula_t B11 = JC(ss_, .0); LinearStateSpace::formula_t B12 = JC(ss_, .0);
        LinearStateSpace::formula_t B21 = JC(ss_, .0); LinearStateSpace::formula_t B22 = JE(ss_, -1.0/_p(iJm));
        LinearStateSpace::formula_t B31 = JE(ss_, 1.0/_p(iLa)); LinearStateSpace::formula_t B32 = JC(ss_, .0);

        LinearStateSpace::formulas_t fB{
            B11, B12,
            B21, B22,
            B31, B32
        };

        LinearStateSpace::formula_t C11 = JC(ss_, 1.0); LinearStateSpace::formula_t C12 = JC(ss_, .0); LinearStateSpace::formula_t C13 = JC(ss_, .0);
        LinearStateSpace::formula_t C21 = JC(ss_, .0); LinearStateSpace::formula_t C22 = JC(ss_, 1.0); LinearStateSpace::formula_t C23 = JC(ss_, .0);
        LinearStateSpace::formula_t C31 = JC(ss_, .0); LinearStateSpace::formula_t C32 = JC(ss_, .0); LinearStateSpace::formula_t C33 = JC(ss_, 1.0);

        LinearStateSpace::formulas_t fC{
            C11, C12, C13,
            C21, C22, C23,
            C31, C32, C33
        };

        LinearStateSpace::formula_t D11 = JC(ss_, .0); LinearStateSpace::formula_t D12 = JC(ss_, .0);
        LinearStateSpace::formula_t D21 = JC(ss_, .0); LinearStateSpace::formula_t D22 = JC(ss_, .0);
        LinearStateSpace::formula_t D31 = JC(ss_, .0); LinearStateSpace::formula_t D32 = JC(ss_, .0);

        LinearStateSpace::formulas_t fD{
            D11, D12,
            D21, D22,
            D31, D32,
        };

        ss_.setA(fA);
        ss_.setB(fB);
        ss_.setC(fC);
        ss_.setD(fD);
    }

    ss_.A().print("A : ");
    ss_.B().print("B : ");
    ss_.C().print("C : ");
    ss_.D().print("D : ");

    std::cout << "Controllable : " << jacl::lti_common::controllable(ss_.A(), ss_.B()) << std::endl;
    std::cout << "Observable : " << jacl::lti_common::observable(ss_.A(), ss_.C()) << std::endl;

    //-- Observer for MIMO DC motor
    arma::mat observer_K;
    arma::mat cpoles{-50,-51,-52};
    jacl::pole_placement::KautskyNichols(&ss_, cpoles, &observer_K, jacl::pole_placement::PolePlacementType::Observer);
    cobserver_.setGain(observer_K.t());

    setupSIMODCMotor();
    setupDiscreteController();
    setupPositionController();
    setupSpeedController();    

    //-- plotter
    csys_plt_.init();
    csys_plt_.setTitle("DC Motor");
    csys_plt_.setDelay() = SAMPLING_PERIOD;
    csys_plt_.setPlotName({"Angular Position", "Angular Velocity", "Current",
                       "Voltage In", "Torque In",
                       "Angular Position", "Angular Velocity", "Current"});

    cobserver_plt_.init();
    cobserver_plt_.setTitle("Full-order Luenberger Observer of DC Motor");
    cobserver_plt_.setDelay() = SAMPLING_PERIOD;
    cobserver_plt_.setPlotName({"Est. Position", "Est. Velocity", "Est. Current"
                              ,"Est. Out Position", "Est. Out Velocity", "Est. Out Current"});

    setupWidgets();
    setupControllerDialog();
    setupMenus();
    setupActions();

    //-- Transient reponse
    {
        arma::vec out(3,1,arma::fill::zeros),in(1,1,arma::fill::zeros);
        arma::vec err(in),est(out);
        din_ = in;
        control_mode_ = jacl::traits::toUType(jacl::ControllerDialog::ControlMode::Position);
        sifd_.init({{.8,.5}}, "SIFD", {"Est. Curr.","Est. Spd.","Est. Pos."});
        constexpr std::size_t ND(1000);
        std::vector<double> resp_err1(ND,.0), resp_err2(ND,.0), resp_err3(ND,.0);        
        auto min = std::numeric_limits<double>::min();
        auto max = std::numeric_limits<double>::max();
        std::pair<double, double> peak1(min,max), peak2(min,max), peak3(min,max);
        jacl::Random<double,3,std::normal_distribution<double>> noise_gen(
            {0,0,0},
            {{10.,0,0},
             {0,10.,0},
             {0,0,5.}}
        );
        arma::vec::fixed<3> noise;

        for(std::size_t i(0); i < ND; i++){

            if(i > ND/64){
                ref_(control_mode_) = 1.693;
            }

            //-- Error between reference and output
            err(control_mode_) = ref_(control_mode_) - out(control_mode_);

            //-- Discrete
            din_ = pos_dctrl_sys_.convolve(err);

            //-- saturation
            if(din_(0) > 24.)
                din_(0) = 24.;
            else if(din_(0) < -24.)
                din_(0) = -24.;
            
            out = dsys_simo_.convolve(din_);            
            noise = noise_gen();
            // out(1) += noise(1);
            // out(2) += noise(2);

            if(i > ND/64){
                out(control_mode_) += bias_;
                if(std::fabs(out(control_mode_)) < dead_zone_)
                    out(control_mode_) = .0;
                out(control_mode_) *= scale_;
            }           

            SIFD::diag_pack_t diag_pack = sifd_.detect(din_, out);            
            resp_err1[i] = std::get<1>(diag_pack[0])(0);            
            resp_err2[i] = std::get<1>(diag_pack[1])(0);
            resp_err3[i] = std::get<1>(diag_pack[2])(0);            

            if(resp_err1[i] > std::get<0>(peak1))
                std::get<0>(peak1) = resp_err1[i];
            if(resp_err1[i] < std::get<1>(peak1))
                std::get<1>(peak1) = resp_err1[i];

            if(resp_err2[i] > std::get<0>(peak2))
                std::get<0>(peak2) = resp_err2[i];
            if(resp_err2[i] < std::get<1>(peak2))
                std::get<1>(peak2) = resp_err2[i];

            if(resp_err3[i] > std::get<0>(peak3))
                std::get<0>(peak3) = resp_err3[i];
            if(resp_err3[i] < std::get<1>(peak3))
                std::get<1>(peak3) = resp_err3[i];
        }
        // std::cout << "Position error max : " << std::get<0>(peak1)
                //   << " ; min : " << std::get<1>(peak1) << std::endl;
        // jacl::plot(resp_err1, SAMPLING_PERIOD, "Position error", {"Position sensor error"});
        std::cout << "Velocity error max : " << std::get<0>(peak2)
                  << " ; min : " << std::get<1>(peak2) << std::endl;
        std::transform(resp_err2.begin(), resp_err2.end(),
                    resp_err2.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, 60./(2.*M_PI)));
        jacl::plot(resp_err2, SAMPLING_PERIOD, "Error speed response sim - r = 97(deg)", {"Error speed response sim (rpm)"});
        std::cout << "Current error max : " << std::get<0>(peak3)
                  << " ; min : " << std::get<1>(peak3) << std::endl;
        std::transform(resp_err3.begin(), resp_err3.end(),
                    resp_err3.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, 1000));          
        jacl::plot(resp_err3, SAMPLING_PERIOD, "Error current response sim - r = 97(deg)", {"Error current response sim (mA)"});
    }

    std::vector<double> position_response,speed_response,current_response;    
    jacl::parser::readArray(&position_response, "../sample/position_sample.txt");
    jacl::analysis::transient_data_t transient_data1 = jacl::analysis::transient(position_response, {90}, SAMPLING_PERIOD, "Position response real - r = 90(deg)");
    std::cout << "[Real Pos] Rise time : " << jacl::analysis::getRiseTime(transient_data1) << std::endl;
    std::cout << "[Real Pos] Peak time : " << jacl::analysis::getPeakTime(transient_data1) << std::endl;
    std::cout << "[Real Pos] Overshoot : " << jacl::analysis::getOvershoot(transient_data1) << std::endl;
    std::cout << "[Real Pos] Settling time : " << jacl::analysis::getSettlingTime(transient_data1) << std::endl;
    jacl::parser::readArray(&speed_response, "../sample/speed_sample.txt");
    jacl::analysis::transient_data_t transient_data2 = jacl::analysis::transient(speed_response, {1260}, SAMPLING_PERIOD, "Speed response real - r = 1260(rpm)");
    std::cout << "[Real Spd] Rise time : " << jacl::analysis::getRiseTime(transient_data2) << std::endl;
    std::cout << "[Real Spd] Peak time : " << jacl::analysis::getPeakTime(transient_data2) << std::endl;
    std::cout << "[Real Spd] Overshoot : " << jacl::analysis::getOvershoot(transient_data2) << std::endl;
    std::cout << "[Real Spd] Settling time : " << jacl::analysis::getSettlingTime(transient_data2) << std::endl;
    jacl::parser::readArray(&current_response, "../sample/current_sample.txt");
    jacl::plot(current_response, SAMPLING_PERIOD, "Current response real - r = 149(deg)", {"Current response real (mA)"});

    std::vector<double> err_speed_response,err_current_response;
    jacl::parser::readArray(&err_speed_response, "../sample/error_speed_sample.txt");
    jacl::plot(err_speed_response, SAMPLING_PERIOD, "Error speed response real - r = 97(deg)", {"Error speed response real (rpm)"});
    jacl::parser::readArray(&err_current_response, "../sample/error_current_sample.txt");
    jacl::plot(err_current_response, SAMPLING_PERIOD, "Error current response real - r = 97(deg)", {"Error current response real (mA)"});

    cl_thread_ = boost::thread{boost::bind(&MainWindow::closedLoopProcess, this)};
}

MainWindow::~MainWindow(){
    cl_status_ = false;
    cl_thread_.join();
    
    delete ui;
}

void MainWindow::setupWidgets(){

    //-- Params
    params_gb_ = new QGroupBox;
    params_gl_ = new QGridLayout;
    params_label_[iBm] = new QLabel;
    params_label_[iBm]->setText(tr("Viscous Friction : "));
    params_dsb_[iBm] = new QDoubleSpinBox;
    params_dsb_[iBm]->setSingleStep(1e-6);
    params_dsb_[iBm]->setDecimals(6);
    params_dsb_[iBm]->setMaximum(1.0);
    params_dsb_[iBm]->setMinimum(.0);
    params_dsb_[iBm]->setValue(bm.nominal);
    params_dsb_[iBm]->adjustSize();
//    params_dsb_[iBm]->resize(params_dsb_[iBm]->minimumSizeHint().width(),
//                params_dsb_[iBm]->minimumSizeHint().height());

    params_label_[iJm] = new QLabel;
    params_label_[iJm]->setText(tr("Rotor Inertia : "));
    params_dsb_[iJm] = new QDoubleSpinBox;
    params_dsb_[iJm]->setSingleStep(1e-7);
    params_dsb_[iJm]->setDecimals(7);
    params_dsb_[iJm]->setMaximum(1.0);
    params_dsb_[iJm]->setMinimum(.0);
    params_dsb_[iJm]->setValue(jm.nominal);
    params_dsb_[iJm]->adjustSize();

    params_label_[iKi] = new QLabel;
    params_label_[iKi]->setText(tr("Torque Constant : "));
    params_dsb_[iKi] = new QDoubleSpinBox;
    params_dsb_[iKi]->setSingleStep(1e-3);
    params_dsb_[iKi]->setDecimals(3);
    params_dsb_[iKi]->setMaximum(10.0);
    params_dsb_[iKi]->setMinimum(.0);
    params_dsb_[iKi]->setValue(ki.nominal);
    params_dsb_[iKi]->adjustSize();

    params_label_[iLa] = new QLabel;
    params_label_[iLa]->setText(tr("Armature Inductance : "));
    params_dsb_[iLa] = new QDoubleSpinBox;
    params_dsb_[iLa]->setSingleStep(1e-5);
    params_dsb_[iLa]->setDecimals(5);
    params_dsb_[iLa]->setMaximum(10.0);
    params_dsb_[iLa]->setMinimum(.0);
    params_dsb_[iLa]->setValue(la.nominal);
    params_dsb_[iLa]->adjustSize();

    params_label_[iKb] = new QLabel;
    params_label_[iKb]->setText(tr("Back EMF : "));
    params_dsb_[iKb] = new QDoubleSpinBox;
    params_dsb_[iKb]->setSingleStep(1e-3);
    params_dsb_[iKb]->setDecimals(3);
    params_dsb_[iKb]->setMaximum(10.0);
    params_dsb_[iKb]->setMinimum(.0);
    params_dsb_[iKb]->setValue(kb.nominal);
    params_dsb_[iKb]->adjustSize();

    params_label_[iRa] = new QLabel;
    params_label_[iRa]->setText(tr("Armature Resistance : "));
    params_dsb_[iRa] = new QDoubleSpinBox;
    params_dsb_[iRa]->setSingleStep(1e-2);
    params_dsb_[iRa]->setDecimals(2);
    params_dsb_[iRa]->setMaximum(50.0);
    params_dsb_[iRa]->setMinimum(.0);
    params_dsb_[iRa]->setValue(ra.nominal);
    params_dsb_[iRa]->adjustSize();

    params_gl_->addWidget(params_label_[iBm],0,0,1,1);
    params_gl_->addWidget(params_label_[iJm],1,0,1,1);
    params_gl_->addWidget(params_label_[iKi],2,0,1,1);
    params_gl_->addWidget(params_label_[iLa],3,0,1,1);
    params_gl_->addWidget(params_label_[iKb],4,0,1,1);
    params_gl_->addWidget(params_label_[iRa],5,0,1,1);
    params_gl_->addWidget(params_dsb_[iBm],  0,1,1,1);
    params_gl_->addWidget(params_dsb_[iJm],  1,1,1,1);
    params_gl_->addWidget(params_dsb_[iKi],  2,1,1,1);
    params_gl_->addWidget(params_dsb_[iLa],  3,1,1,1);
    params_gl_->addWidget(params_dsb_[iKb],  4,1,1,1);
    params_gl_->addWidget(params_dsb_[iRa],  5,1,1,1);

    params_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),6,2);

    params_gb_->setLayout(params_gl_);
    params_gb_->setTitle(tr("Parameters"));

    //-- Command

    command_gb_ = new QGroupBox;
    command_gl_ = new QGridLayout;

    perturb_pb_ = new QPushButton;
    perturb_pb_->setText(tr("Perturb !"));

    reset_pb_ = new QPushButton;
    reset_pb_->setText(tr("Reset"));

    simulate_pb_ = new QPushButton;
    simulate_pb_->setText(tr("Simulate"));

    set_input_pb_ = new QPushButton;
    set_input_pb_->setText(tr("Disturb !"));

    command_gl_->addWidget(perturb_pb_  , 0,0,1,1);
    command_gl_->addWidget(reset_pb_    , 0,1,1,1);
    command_gl_->addWidget(simulate_pb_ , 1,0,1,1);
    command_gl_->addWidget(set_input_pb_, 1,1,1,1);
    command_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Ignored,QSizePolicy::Expanding),2,2);

    command_gb_->setLayout(command_gl_);
    command_gb_->setTitle(tr("Command"));

    //-- Input
    input_gb_ = new QGroupBox;
    input_gl_ = new QGridLayout;

    voltage_in_label_ = new QLabel;
    voltage_in_label_->setText(tr("Voltage : "));
    voltage_in_dsb_ = new QDoubleSpinBox;
    voltage_in_dsb_->setSingleStep(1e-3);
    voltage_in_dsb_->setDecimals(3);
    voltage_in_dsb_->setMinimum(-24.0);
    voltage_in_dsb_->setMaximum(24.0); // 24 volt
    voltage_in_dsb_->setValue(.0);
    voltage_in_dsb_->adjustSize();

    torque_in_label_ = new QLabel;
    torque_in_label_->setText(tr("Torque : "));
    torque_in_dsb_ = new QDoubleSpinBox;
    torque_in_dsb_->setSingleStep(1e-3);
    torque_in_dsb_->setDecimals(3);
    torque_in_dsb_->setMinimum(.0);
    torque_in_dsb_->setMaximum(10.0);
    torque_in_dsb_->setValue(.0);
    torque_in_dsb_->adjustSize();

    input_gl_->addWidget(voltage_in_label_, 0,0,1,1);
    input_gl_->addWidget(voltage_in_dsb_,   0,1,1,1);
    input_gl_->addWidget(torque_in_label_,  1,0,1,1);
    input_gl_->addWidget(torque_in_dsb_,    1,1,1,1);
    input_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Ignored,QSizePolicy::Expanding),2,2);

    input_gb_->setLayout(input_gl_);
    input_gb_->setTitle(tr("Input Disturbance"));

    //-- Fault
    fault_gb_ = new QGroupBox;
    fault_gl_ = new QGridLayout;

    target_label_ = new QLabel;
    target_label_->setText(tr("Sensor : "));

    target_cb_ = new QComboBox;
    target_cb_->addItem(tr("Position"));
    target_cb_->addItem(tr("Velocity"));
    target_cb_->addItem(tr("Current"));

    fault_pb_ = new QPushButton;
    fault_pb_->setText(tr("Fault !"));

    bias_label_ = new QLabel;
    bias_label_->setText(tr("Bias : "));
    bias_dsb_ =  new QDoubleSpinBox;
    // Settings DSB params
    bias_dsb_->setValue(.0);
    bias_dial_ = new QDial;

    scale_label_ = new QLabel;
    scale_label_->setText(tr("Scale : "));
    scale_dsb_ =  new QDoubleSpinBox;
    // Settings DSB params    
    scale_dsb_->setMaximum(2.);
    scale_dsb_->setMinimum(.0);
    scale_dsb_->setValue(1.);
    scale_dsb_->setDecimals(3);
    scale_dsb_->setSingleStep(1e-3);
    scale_dial_ = new QDial;    
    scale_dial_->setMaximum(2000);
    scale_dial_->setMinimum(0);
    scale_dial_->setValue(1000);
    scale_dial_->setSingleStep(1);

    dead_zone_label_ = new QLabel;
    dead_zone_label_->setText(tr("Dead Zone : "));
    dead_zone_dsb_ =  new QDoubleSpinBox;
    // Settings DSB params
    dead_zone_dsb_->setValue(.0);
    dead_zone_dial_ = new QDial;

    fault_gl_->addWidget(target_label_,    0,0,1,1);
    fault_gl_->addWidget(target_cb_,       0,1,1,1);
    fault_gl_->addWidget(fault_pb_,      0,2,1,1);
    fault_gl_->addWidget(bias_label_,      1,0,1,1);
    fault_gl_->addWidget(bias_dsb_,        1,1,1,1);
    fault_gl_->addWidget(bias_dial_,       1,2,1,1);
    fault_gl_->addWidget(scale_label_,     2,0,1,1);
    fault_gl_->addWidget(scale_dsb_,       2,1,1,1);
    fault_gl_->addWidget(scale_dial_,      2,2,1,1);
    fault_gl_->addWidget(dead_zone_label_, 3,0,1,1);
    fault_gl_->addWidget(dead_zone_dsb_,   3,1,1,1);
    fault_gl_->addWidget(dead_zone_dial_,  3,2,1,1);
    fault_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),4,3);

    fault_gb_->setLayout(fault_gl_);
    fault_gb_->setTitle(tr("Fault"));

    //-- Watermark
    watermark_widget_ = new QWidget;
    watermark_gl_ = new QGridLayout;

    QImageReader image_reader("../example/gui/Logo_Universitas_Gadjah_Mada.png");
    logo_image_ = new QImage;
    *logo_image_ = image_reader.read().scaled(QSize(200,200), Qt::KeepAspectRatio);
    logo_image_label_ = new QLabel;
    logo_image_label_->setPixmap(QPixmap::fromImage(*logo_image_));

    title_label_ = new QLabel;
    title_label_->setText(tr("Control and Sensor Fault Isolation\n"
                             "   Design on DC Motor along with\n"
                             "           Remote Monitoring"));
    title_label_->setStyleSheet("font: bold italic");

    watermark_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding), 0,0);
    watermark_gl_->addWidget(logo_image_label_, 1,0,1,1);
    watermark_gl_->setAlignment(logo_image_label_, Qt::AlignCenter);
    watermark_gl_->addWidget(title_label_,      2,0,1,1);
    watermark_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding), 3,1);
    watermark_widget_->setLayout(watermark_gl_);

    //-- Main
    main_widget_ = new QWidget;
    main_layout_ = new QGridLayout;
    main_layout_->addWidget(params_gb_,        0,0,2,1);
    main_layout_->addWidget(input_gb_,         0,1,1,1);
    main_layout_->addWidget(command_gb_,       1,1,1,1);
    main_layout_->addWidget(fault_gb_,         2,0,1,1);
    main_layout_->addWidget(watermark_widget_, 2,1,1,1);

//    main_layout_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),1,2);
    main_layout_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),3,2);
    main_widget_->setLayout(main_layout_);

    this->setCentralWidget(main_widget_);
    this->adjustSize();
    this->layout()->setSizeConstraint(QLayout::SetFixedSize);
    this->setWindowTitle(tr("Control Simulator"));

}

void MainWindow::setupControllerDialog(){
    controller_dialog_ = new jacl::ControllerDialog(this);
}

void MainWindow::setupMenus(){
    controller_act_ = new QAction(tr("Controller"), this);
    controller_act_->setStatusTip(tr("Open controller dialog"));
    connect(controller_act_, SIGNAL(triggered()), this, SLOT(openControllerDialog()));
    tools_menu_ = this->menuBar()->addMenu(tr("Tools"));
    tools_menu_->addAction(controller_act_);
}

void MainWindow::setupActions(){
    connect(perturb_pb_, SIGNAL(clicked()), this, SLOT(perturbAct()));
    connect(reset_pb_, SIGNAL(clicked()), this, SLOT(resetAct()));
    connect(simulate_pb_, SIGNAL(clicked()), this, SLOT(simulateAct()));
    connect(set_input_pb_, SIGNAL(clicked()), this, SLOT(setInputAct()));

    connect(bias_dsb_, SIGNAL(valueChanged(double)), this, SLOT(biasDialConv(double)));
    connect(scale_dsb_, SIGNAL(valueChanged(double)), this, SLOT(scaleDialConv(double)));
    connect(dead_zone_dsb_, SIGNAL(valueChanged(double)), this, SLOT(deadZoneDialConv(double)));

    connect(bias_dial_, SIGNAL(valueChanged(int)), this, SLOT(biasDSBConv(int)));
    connect(scale_dial_, SIGNAL(valueChanged(int)), this, SLOT(scaleDSBConv(int)));
    connect(dead_zone_dial_, SIGNAL(valueChanged(int)), this, SLOT(deadZoneDSBConv(int)));

    connect(controller_dialog_, SIGNAL(setRefSig()), this, SLOT(refAct()));
    connect(controller_dialog_, SIGNAL(setModeSig(int)), this, SLOT(modeAct(int)));
    qRegisterMetaType<arma::vec>("arma::vec");
    connect(this, SIGNAL(setDataMonitor(arma::vec)), controller_dialog_, SLOT(setDataMonitor(arma::vec)));
    qRegisterMetaType<arma::Col<uint8_t>>("arma::Col<uint8_t>");
    connect(this, SIGNAL(setSensorStatus(arma::Col<uint8_t>)), controller_dialog_, SLOT(setSensorStatus(arma::Col<uint8_t>)));
}

void MainWindow::perturbAct(){
    for(auto &pair:{std::pair<std::decay<jacl::PhysicalParameter<>>::type*, QDoubleSpinBox*>
                        {&bm,params_dsb_[iBm]},
                        {&jm,params_dsb_[iJm]},
                        {&ki,params_dsb_[iKi]},
                        {&la,params_dsb_[iLa]},
                        {&kb,params_dsb_[iKb]},
                        {&ra,params_dsb_[iRa]}})
        pair.first->perturbed = pair.second->value();

    // cobserver_plt_.updateVar();
//    sim_.setLinearStateSpace(ss_.A(), ss_.B(), ss_.C(), ss_.D());
}

void MainWindow::resetAct(){

    params_dsb_[iBm]->setValue(bm.nominal);
    params_dsb_[iJm]->setValue(jm.nominal);
    params_dsb_[iKi]->setValue(ki.nominal);
    params_dsb_[iLa]->setValue(la.nominal);
    params_dsb_[iKb]->setValue(kb.nominal);
    params_dsb_[iRa]->setValue(ra.nominal);

    perturbAct();
}

void MainWindow::simulateAct(){
    //-- Continuous
    // csys_plt_.start();
    // cobserver_plt_.start();
    // posctrl_plt_.start();
    // spdctrl_plt_.start();
    // cobserver_plt_.simulate();

    //-- Discrete
    dsys_simo_plt_.start();
    pos_dctrl_plt_.start();
    // spd_dctrl_plt_.start();
    // dobserver_simo_plt_.start();
    sifd_.viewSignals();
}

void MainWindow::setInputAct(){
    arma::mat in(2, 1);
    in(0) += voltage_in_dsb_->value();
    in(1) += torque_in_dsb_->value();
    din_ += in.submat(0,0,0,0);
}

void MainWindow::biasDialConv(double _val){
    bias_dial_->setValue((int)_val);
}

void MainWindow::scaleDialConv(double _val){
    scale_dial_->setValue(_val * 1000);
}

void MainWindow::deadZoneDialConv(double _val){
    dead_zone_dial_->setValue((int)_val);
}

void MainWindow::biasDSBConv(int _val){
    bias_dsb_->setValue((double)_val);
}

void MainWindow::scaleDSBConv(int _val){
    scale_dsb_->setValue((double)_val * .001);
}

void MainWindow::deadZoneDSBConv(int _val){
    dead_zone_dsb_->setValue((double)_val);
}

void MainWindow::openControllerDialog(){
    controller_dialog_->show();
}

void MainWindow::refAct(){
    ref_(0) = controller_dialog_->getPosition();
    ref_(1) = controller_dialog_->getVelocity();
    ref_(2) = controller_dialog_->getCurrent();
}

void MainWindow::modeAct(int _val){
    control_mode_ = _val;
}

void MainWindow::makeFault(arma::vec* _out){
    int target_idx = target_cb_->currentIndex();
    (*_out)(target_idx) += bias_dsb_->value();
    if(std::fabs((*_out)(target_idx)) < dead_zone_dsb_->value())
        (*_out)(target_idx) = 0;
    (*_out)(target_idx) *= scale_dsb_->value();
}

void MainWindow::setupSIMODCMotor(){
    //-- SIMO DC motor open-loop
    {
        SIMO::formulas_t fA{
            JC(simo_,.0),JC(simo_,1.0),JC(simo_,.0),
            JC(simo_,.0),JE(simo_,-_p(iBm)/_p(iJm)),JE(simo_,_p(iKi)/_p(iJm)),
            JC(simo_,.0),JE(simo_,-_p(iKb)/_p(iLa)),JE(simo_,-_p(iRa)/_p(iLa))
        };

        SIMO::formulas_t fB{
            JC(simo_, .0),
            JC(simo_, .0),
            JE(simo_, 1.0/_p(iLa))
        };

        SIMO::formulas_t fC{
            JC(simo_, 1.0), JC(simo_, .0), JC(simo_, .0),
            JC(simo_, .0), JC(simo_, 1.0), JC(simo_, .0),
            JC(simo_, .0), JC(simo_, .0), JC(simo_, 1.0)
        };

        SIMO::formulas_t fD{
            JC(simo_, .0),
            JC(simo_, .0),
            JC(simo_, .0)
        };
        simo_.setA(fA);
        simo_.setB(fB);
        simo_.setC(fC);
        simo_.setD(fD);
    }    
    jacl::lti_common::StateSpacePack dsimo = jacl::lti_common::discretize(simo_, SAMPLING_PERIOD);    
    dsimo_.setA(std::get<0>(dsimo)); dsimo_.setB(std::get<1>(dsimo));
    dsimo_.setC(std::get<2>(dsimo)); dsimo_.setD(std::get<3>(dsimo));
    dsimo_.A().print("Ad : "); dsimo_.B().print("Bd : ");
    dsimo_.C().print("Cd : "); dsimo_.D().print("Dd : ");
    jacl::parser::saveStateSpace(dsimo_, "motor_dc_dsimo.jacl");

    arma::mat dobsv_gain;
    arma::mat dpoles{-0.8,-0.85,0.83};
    jacl::pole_placement::KautskyNichols(&dsimo_, dpoles, &dobsv_gain, jacl::pole_placement::PolePlacementType::Observer);

    dobserver_simo_.setGain(dobsv_gain.t());
    dobsv_gain.print("\nDiscrete Observer Gain : ");
    jacl::parser::saveGain(dobsv_gain, "discrete_observer_gain.jacl");
//    jacl::parser::readGain(&obsv_gain, "discrete_observer_gain.jacl");

    dsys_simo_plt_.init();
    dsys_simo_plt_.setTitle("SIMO DC motor");
    dsys_simo_plt_.setPlotName({"Position", "Velocity", "Current", "Voltage"});

    dobserver_simo_plt_.init();
    dobserver_simo_plt_.setTitle("Observer SIMO DC motor");
    dobserver_simo_plt_.setPlotName({"Est. Position", "Est. Velocity", "Est. Current"});

    {
        MReal::formulas_t fA{
            //-- Row-1
            JC(m_real_, .0), JC(m_real_, 1.), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            //-- Row-2
            JC(m_real_, .0), JE(m_real_, -_p(iBm)/_p(iJm)), JE(m_real_, _p(iKi)/_p(iJm)),
            JE(m_real_, weight_(iBm)*-1./_p(iJm)), JE(m_real_, weight_(iJm)*-1./_p(iJm)), JE(m_real_, weight_(iKi)*-1./_p(iJm)), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            // JE(m_real_, -1./_p(iJm)), JE(m_real_, -1./_p(iJm)), JE(m_real_, -1./_p(iJm)), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            //-- Row-3
            JC(m_real_, .0), JE(m_real_, -_p(iKb)/_p(iLa)), JE(m_real_, -_p(iRa)/_p(iLa)),
            JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JE(m_real_, weight_(iLa)*-1./_p(iLa)), JE(m_real_, weight_(iKb)*-1./_p(iLa)), JE(m_real_, weight_(iRa)*-1./_p(iLa)),
            // JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JE(m_real_, -1./_p(iLa)), JE(m_real_, -1./_p(iLa)), JE(m_real_, -1./_p(iLa)),
            //-- Row-4
            JC(m_real_, .0), JC(m_real_, 1.), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            //-- Row-5
            JC(m_real_, .0), JE(m_real_, -_p(iBm)/_p(iJm)), JE(m_real_, _p(iKi)/_p(iJm)),
            JE(m_real_, weight_(iBm)*-1./_p(iJm)), JE(m_real_, weight_(iJm)*-1./_p(iJm)), JE(m_real_, weight_(iKi)*-1./_p(iJm)), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            // JE(m_real_, -1./_p(iJm)), JE(m_real_, -1./_p(iJm)), JE(m_real_, -1./_p(iJm)), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            //-- Row-6
            JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, 1.), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            //-- Row-7
            JC(m_real_, .0), JE(m_real_, -_p(iKb)/_p(iLa)), JE(m_real_, -_p(iRa)/_p(iLa)),
            JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JE(m_real_, weight_(iLa)*-1./_p(iLa)), JE(m_real_, weight_(iKb)*-1./_p(iLa)), JE(m_real_, weight_(iRa)*-1./_p(iLa)),
            // JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JE(m_real_, -1./_p(iLa)), JE(m_real_, -1./_p(iLa)), JE(m_real_, -1./_p(iLa)),
            //-- Row-8
            JC(m_real_, .0), JC(m_real_, 1.), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            //-- Row-9
            JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, 1.), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0)
        };

        MReal::formulas_t fB{
            JC(m_real_, .0),        
            JC(m_real_, .0),
            JE(m_real_, -1./_p(iLa)),
            JC(m_real_, .0),
            JC(m_real_, .0),
            JC(m_real_, .0),
            JE(m_real_, -1./_p(iLa)),
            JC(m_real_, .0),
            JC(m_real_, .0)
        };

        MReal::formulas_t fC{
            JC(m_real_, 1.), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0)
            // JC(m_real_, .0), JC(m_real_, 1.), JC(m_real_, .0),
            // JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, 1.)
        };

        MReal::formulas_t fD{
            JC(m_real_, .0)
            // JC(m_real_, .0),
            // JC(m_real_, .0)
        };
        m_real_.setA(fA);
        m_real_.setB(fB);
        m_real_.setC(fC);
        m_real_.setD(fD);
    }
    jacl::lti_common::StateSpacePack dm_real_ = jacl::lti_common::discretize(m_real_, SAMPLING_PERIOD);
    {
        arma::mat temp1, temp2, temp3;
        arma::mat zeros9x1(9,1,arma::fill::zeros);
        arma::mat zeros1x9(1,9,arma::fill::zeros);
        arma::mat zeros1x1(1,1,arma::fill::zeros);
        arma::mat eye1x1(1,1,arma::fill::eye);

        temp1 = arma::join_horiz(zeros9x1, std::get<1>(dm_real_));
        temp2 = arma::join_horiz(zeros9x1, temp1);
        arma::mat B = arma::join_horiz(std::get<1>(dm_real_), temp2);

        temp1 = arma::join_vert(zeros1x9, -std::get<2>(dm_real_));
        arma::mat C = arma::join_vert(-std::get<2>(dm_real_), temp1);

        temp1 = arma::join_horiz(zeros1x1, eye1x1);
        temp2 = arma::join_horiz(zeros1x1, temp1);
        temp3 = arma::join_horiz(zeros1x1, arma::join_horiz(zeros1x1, zeros1x1));
        arma::mat D11 = arma::join_vert(temp2, temp3);

        arma::mat D12 = arma::join_vert(zeros1x1, eye1x1);

        temp1 = arma::join_horiz(-eye1x1,eye1x1);
        arma::mat D21 = arma::join_horiz(zeros1x1,temp1);
        
        arma::mat D22 = zeros1x1;

        temp1 = arma::join_horiz(D11, D12);
        temp2 = arma::join_horiz(D21, D22);
        arma::mat D = arma::join_vert(temp1, temp2);

        m_icm_.setA(std::get<0>(dm_real_));
        m_icm_.setB(B);
        m_icm_.setC(C);
        m_icm_.setD(D);
    }    
}

void MainWindow::setupDiscreteController(){

    //-- position controller
    {
        auto weight1(.35);
        auto weight2(1.);
        pos_real_.setA(dsimo_.A());
        pos_real_.setB(dsimo_.B());
        pos_real_.setC({1.,0.,0.});
        pos_real_.setD(arma::zeros(1,1));

        arma::mat temp1, temp2, temp3;
        arma::mat zeros3x1(3,1,arma::fill::zeros);
        arma::mat zeros1x3(1,3,arma::fill::zeros);
        arma::mat zeros1x1(1,1,arma::fill::zeros);
        arma::mat eye1x1(1,1,arma::fill::eye);

        temp1 = arma::join_horiz(zeros3x1, pos_real_.B());
        temp2 = arma::join_horiz(zeros3x1, temp1);        
        arma::mat B = arma::join_horiz(pos_real_.B(), temp2);
        arma::mat B1 = B.head_cols(3);
        arma::mat B2 = B.tail_cols(1);

        temp1 = arma::join_vert(zeros1x3, -pos_real_.C() * weight1);
        arma::mat C = arma::join_vert(-pos_real_.C(), temp1);
        arma::mat C1 = C.head_rows(2);
        arma::mat C2 = C.tail_rows(1);

        temp1 = arma::join_horiz(zeros1x1, eye1x1);
        temp2 = arma::join_horiz(zeros1x1, temp1);
        temp3 = arma::join_horiz(zeros1x1, arma::join_horiz(zeros1x1, zeros1x1));
        arma::mat D11 = arma::join_vert(temp2, temp3);

        arma::mat D12 = arma::join_vert(zeros1x1, eye1x1 * weight2);

        temp1 = arma::join_horiz(-eye1x1 * weight1, eye1x1 * weight1);
        arma::mat D21 = arma::join_horiz(zeros1x1,temp1);
        
        arma::mat D22 = zeros1x1;

        temp1 = arma::join_horiz(D11, D12);
        temp2 = arma::join_horiz(D21, D22);
        arma::mat D = arma::join_vert(temp1, temp2);

        pos_icm_.setA(pos_real_.A());
        pos_icm_.setB(B);
        pos_icm_.setC(C);
        pos_icm_.setD(D);

        pos_icm_.A().print("[Pos] A_icm : ");
        pos_icm_.B().print("[Pos] B_icm : ");
        pos_icm_.C().print("[Pos] C_icm : ");
        pos_icm_.D().print("[Pos] D_icm : ");
        pos_dhinf_ = new PosDHinf(&pos_sys_, gamma_pos_, {1e-4,1e-4,1e-4}, {1e-4,1e-4,1e-4});
        auto K( pos_dhinf_->solve() );
        pos_dctrl_.setA(std::get<0>(K)); pos_dctrl_.setB(std::get<1>(K));
        pos_dctrl_.setC(std::get<2>(K)); pos_dctrl_.setD(std::get<3>(K));
        jacl::parser::saveStateSpace(pos_dctrl_, "position_controller.jacl");

        //-- closed-loop
        jacl::state_space::Linear<double, 6,1,1> cl_ss;
        jacl::system::Discrete<jacl::state_space::Linear<double, 6,1,1>> cl_sys(&cl_ss, SAMPLING_PERIOD);
        {
            arma::mat temp1;

            temp1 = pos_real_.B()*pos_dctrl_.D()*pos_real_.C();
            arma::mat A = arma::join_cols(
                arma::join_rows(pos_real_.A() - temp1, pos_real_.B()*pos_dctrl_.C()),
                arma::join_rows(-pos_dctrl_.B()*pos_real_.C(), pos_dctrl_.A())
            );
            
            arma::mat B = arma::join_cols(
                pos_real_.B()*pos_dctrl_.D(),
                pos_dctrl_.B()
            );

            cl_ss.setA(A);
            cl_ss.setB(B);
            cl_ss.setC(arma::join_rows(pos_real_.C(), arma::zeros(1,3)));
            cl_ss.setD(pos_real_.D());
        }

        cl_ss.A().print("[Pos] A_cl : ");
        cl_ss.B().print("[Pos] B_cl : ");
        cl_ss.C().print("[Pos] C_cl : ");
        cl_ss.D().print("[Pos] D_cl : ");
        jacl::analysis::transient_data_t transient_data =
            jacl::analysis::transient<jacl::system::Discrete<jacl::state_space::Linear<double,6,1,1>>,5000>(cl_sys, {1.57},
                "Closed-loop position control simulation - r = 1.57 (rad or 90 degree)");

        std::cout << "[Pos] Rise time : " << jacl::analysis::getRiseTime(transient_data) << std::endl;
        std::cout << "[Pos] Peak time : " << jacl::analysis::getPeakTime(transient_data) << std::endl;
        std::cout << "[Pos] Overshoot : " << jacl::analysis::getOvershoot(transient_data) << std::endl;
        std::cout << "[Pos] Settling time : " << jacl::analysis::getSettlingTime(transient_data) << std::endl;

        //-- Nominal analysis

        jacl::state_space::Linear<double, 6, 3, 2> llft;
        arma::mat R = eye1x1 - D22*pos_dctrl_.D();
        arma::mat R_inv = arma::inv(R);
        arma::mat S = eye1x1 - pos_dctrl_.D()*D22;
        arma::mat S_inv = arma::inv(S);
        temp1 = B2*S_inv*pos_dctrl_.D()*C2;
        temp2 = pos_dctrl_.B()*R_inv*D22*pos_dctrl_.C();
        llft.setA(
            arma::join_cols(
                arma::join_rows(pos_icm_.A() + temp1, B2*S_inv*pos_dctrl_.C()),
                arma::join_rows(pos_dctrl_.B()*R_inv*C2, pos_dctrl_.A() + temp2)
            )
        );

        llft.setB(
            arma::join_cols(
                B1 + B2*pos_dctrl_.D()*D21,
                pos_dctrl_.B()*R_inv*D21
            )
        );

        temp1 = D12*S_inv*pos_dctrl_.D()*C2;
        llft.setC(arma::join_rows(C1 + temp1, D12*S_inv*pos_dctrl_.C()));

        temp1 = D12*S_inv*pos_dctrl_.D()*D21;
        llft.setD(D11 + temp1);

        llft.A().print("[Pos] LLFT A : ");
        llft.B().print("[Pos] LLFT B : ");
        llft.C().print("[Pos] LLFT C : ");
        llft.D().print("[Pos] LLFT D : ");
        std::cout << "[Pos] Nominal stability : " << std::boolalpha << (bool)jacl::analysis::nominalStability(pos_real_, pos_dctrl_) << std::endl;
        std::cout << "[Pos] Nominal performance : " << std::boolalpha << (bool)jacl::analysis::nominalPerformance(llft, gamma_pos_) << std::endl;

        pos_dctrl_plt_.init();
        pos_dctrl_plt_.setTitle("Discrete Position Controller");
        pos_dctrl_plt_.setPlotName({"Error","Input"});
    }

    //-- speed controller
    jacl::system::Discrete<jacl::state_space::Linear<double,2,1,1>> spd_sys(&spd_real_, SPD_SP);
    jacl::state_space::Linear<double,2,1,1> spd;
    {        
        spd.setA(simo_.A().submat(1,1,2,2));
        spd.setB(simo_.B().tail_rows(2));
        spd.setC({1.,0.});
        spd.setD(arma::zeros(1,1));

        jacl::lti_common::StateSpacePack dspd = jacl::lti_common::discretize(spd, SPD_SP);

        spd_real_.setA(std::get<0>(dspd));
        spd_real_.setB(std::get<1>(dspd));
        spd_real_.setC(std::get<2>(dspd));
        spd_real_.setD(std::get<3>(dspd));

        arma::mat temp1, temp2, temp3;
        arma::mat zeros2x1(2,1,arma::fill::zeros);
        arma::mat zeros1x2(1,2,arma::fill::zeros);
        arma::mat zeros1x1(1,1,arma::fill::zeros);
        arma::mat eye1x1(1,1,arma::fill::eye);
        
        temp1 = arma::join_horiz(zeros2x1, spd_real_.B());
        temp2 = arma::join_horiz(zeros2x1, temp1);
        arma::mat B = arma::join_horiz(spd_real_.B(), temp2);
        arma::mat B1 = B.head_cols(3);
        arma::mat B2 = B.tail_cols(1);

        temp1 = arma::join_vert(zeros1x2, -spd_real_.C());
        arma::mat C = arma::join_vert(-spd_real_.C(), temp1);
        arma::mat C1 = C.head_rows(2);
        arma::mat C2 = C.tail_rows(1);

        temp1 = arma::join_horiz(zeros1x1, eye1x1);
        temp2 = arma::join_horiz(zeros1x1, temp1);
        temp3 = arma::join_horiz(zeros1x1, arma::join_horiz(zeros1x1, zeros1x1));
        arma::mat D11 = arma::join_vert(temp2, temp3);

        arma::mat D12 = arma::join_vert(zeros1x1, eye1x1);

        temp1 = arma::join_horiz(-eye1x1,eye1x1);
        arma::mat D21 = arma::join_horiz(zeros1x1,temp1);
        
        arma::mat D22 = zeros1x1;

        temp1 = arma::join_horiz(D11, D12);
        temp2 = arma::join_horiz(D21, D22);
        arma::mat D = arma::join_vert(temp1, temp2);

        spd_icm_.setA(spd_real_.A());
        spd_icm_.setB(B);
        spd_icm_.setC(C);
        spd_icm_.setD(D);

        spd_icm_.A().print("[Spd] A_icm : ");
        spd_icm_.B().print("[Spd] B_icm : ");
        spd_icm_.C().print("[Spd] C_icm : ");
        spd_icm_.D().print("[Spd] D_icm : ");
        spd_dhinf_ = new SpdDHinf(&spd_sys_, gamma_spd_, {1e-4,1e-4}, {1e-4,1e-4});
        auto K( spd_dhinf_->solve() );
        spd_dctrl_.setA(std::get<0>(K)); spd_dctrl_.setB(std::get<1>(K));
        spd_dctrl_.setC(std::get<2>(K)); spd_dctrl_.setD(std::get<3>(K));
        jacl::parser::saveStateSpace(spd_dctrl_, "speed_controller.jacl");

        //-- closed-loop
        jacl::state_space::Linear<double, 4,1,1> cl_ss;
        jacl::system::Discrete<jacl::state_space::Linear<double, 4,1,1>> cl_sys(&cl_ss, SPD_SP);
        {
            arma::mat temp1;

            temp1 = spd_real_.B()*spd_dctrl_.D()*spd_real_.C();
            arma::mat A = arma::join_cols(
                arma::join_rows(spd_real_.A() - temp1, spd_real_.B()*spd_dctrl_.C()),
                arma::join_rows(-spd_dctrl_.B()*spd_real_.C(), spd_dctrl_.A())
            );
            
            arma::mat B = arma::join_cols(
                spd_real_.B()*spd_dctrl_.D(),
                spd_dctrl_.B()
            );

            cl_ss.setA(A);
            cl_ss.setB(B);
            cl_ss.setC(arma::join_rows(spd_real_.C(), arma::zeros(1,2)));
            cl_ss.setD(spd_real_.D());
        }
        cl_ss.A().print("[Spd] A_cl : ");
        cl_ss.B().print("[Spd] B_cl : ");
        cl_ss.C().print("[Spd] C_cl : ");
        cl_ss.D().print("[Spd] D_cl : ");
        jacl::analysis::transient_data_t transient_data =
            jacl::analysis::transient<jacl::system::Discrete<jacl::state_space::Linear<double, 4,1,1>>,500>(cl_sys, {131},
                "Closed-loop speed control simulation - r = 1260 (rpm) or 131 rad/s");

        std::cout << "[Spd] Rise time : " << jacl::analysis::getRiseTime(transient_data) << std::endl;
        std::cout << "[Spd] Peak time : " << jacl::analysis::getPeakTime(transient_data) << std::endl;
        std::cout << "[Spd] Overshoot : " << jacl::analysis::getOvershoot(transient_data) << std::endl;
        std::cout << "[Spd] Settling time : " << jacl::analysis::getSettlingTime(transient_data) << std::endl;

        //-- Nominal analysis
        jacl::state_space::Linear<double,4,3,2> llft;
        arma::mat R = eye1x1 - D22*spd_dctrl_.D();
        arma::mat R_inv = arma::inv(R);
        arma::mat S = eye1x1 - spd_dctrl_.D()*D22;
        arma::mat S_inv = arma::inv(S);
        temp1 = B2*S_inv*spd_dctrl_.D()*C2;
        temp2 = spd_dctrl_.B()*R_inv*D22*spd_dctrl_.C();
        llft.setA(
            arma::join_cols(
                arma::join_rows(spd_icm_.A() + temp1, B2*S_inv*spd_dctrl_.C()),
                arma::join_rows(spd_dctrl_.B()*R_inv*C2, spd_dctrl_.A() + temp2)
            )
        );

        llft.setB(
            arma::join_cols(
                B1 + B2*spd_dctrl_.D()*D21,
                spd_dctrl_.B()*R_inv*D21
            )
        );

        temp1 = D12*S_inv*spd_dctrl_.D()*C2;
        llft.setC(arma::join_rows(C1 + temp1, D12*S_inv*spd_dctrl_.C()));

        temp1 = D12*S_inv*spd_dctrl_.D()*D21;
        llft.setD(D11 + temp1);

        llft.A().print("[Spd] LLFT A : ");
        llft.B().print("[Spd] LLFT B : ");
        llft.C().print("[Spd] LLFT C : ");
        llft.D().print("[Spd] LLFT D : ");
        std::cout << "[Spd] Nominal stability : " << std::boolalpha << (bool)jacl::analysis::nominalStability(spd_real_, spd_dctrl_) << std::endl;
        std::cout << "[Spd] Nominal performance : " << std::boolalpha << (bool)jacl::analysis::nominalPerformance(llft, gamma_spd_) << std::endl;

        spd_dctrl_plt_.init();
        spd_dctrl_plt_.setTitle("Discrete Speed Controller");
        spd_dctrl_plt_.setPlotName({"Error","Input"});
    }

}

void MainWindow::setupPositionController(){

    siso_pos_.setA({{0.,1.,0.},
                    {0.,-.6667,722.2},
                    {0,-83.87,-3339}});
    siso_pos_.setB(arma::colvec({0.,
                    0.,
                    1613.}));
    siso_pos_.setC({1.,0.,0.});
    siso_pos_.setD(arma::colvec({0.}));

    InterConnMatPos::formulas_t fA = {
        JC(icm_pos_, -1.444e+04), JC(icm_pos_, -7780), JC(icm_pos_, -2912), JC(icm_pos_, -841.9), JC(icm_pos_,-269), JC(icm_pos_, 0),
        JC(icm_pos_, 8192), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 4096), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 1024), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 256), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 1), JC(icm_pos_, 0)
    };


    InterConnMatPos::formulas_t fB = {
        JC(icm_pos_, 64), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 64),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0)
    };

    InterConnMatPos::formulas_t fC = {
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 2.276e-06), JC(icm_pos_, 0.04105), JC(icm_pos_, 0.5249), JC(icm_pos_, 77.95),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 2.276e-06), JC(icm_pos_, 0.04105), JC(icm_pos_, 0.5249), JC(icm_pos_, 77.95)
    };

    InterConnMatPos::formulas_t fD = {
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 1), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 1),
        JC(icm_pos_, 0), JC(icm_pos_, -1), JC(icm_pos_, 1), JC(icm_pos_, 0)
    };    

    icm_pos_.setA(fA);
    icm_pos_.setB(fB);
    icm_pos_.setC(fC);
    icm_pos_.setD(fD);
    jacl::system::Continuous<InterConnMatPos> icm_pos_sys(&icm_pos_);
    hinf_pc_ = new HInfPC(&icm_pos_sys, 2.7882);
    jacl::lti_common::StateSpacePack K( hinf_pc_->solve() );
    k_pos_.setA(std::get<0>(K));
    k_pos_.setB(std::get<1>(K));
    k_pos_.setC(std::get<2>(K));
    k_pos_.setD(std::get<3>(K));    

    k_pos_.A().print("A : ");
    k_pos_.B().print("B : ");
    k_pos_.C().print("C : ");
    k_pos_.D().print("D : ");

    

    //-- Example read
    /*PosCtrl controller;
    jacl::parser::readStateSpace(&controller, "position_controller.jacl");
    controller.A().print("Controller A : ");
    controller.B().print("Controller B : ");
    controller.C().print("Controller C : ");
    controller.D().print("Controller D : ");*/

    posctrl_plt_.init();
    posctrl_plt_.setTitle("Position Controller");
    posctrl_plt_.setDelay() = SAMPLING_PERIOD;
    posctrl_plt_.setPlotName({"Position Error", "Voltage"});
}

void MainWindow::setupSpeedController(){
    InterConnMatSpd::formulas_t fA = {
        JC(icm_spd_, -1.444e+04), JC(icm_spd_, -7780), JC(icm_spd_, -2912), JC(icm_spd_, -841.9), JC(icm_spd_, -269),
        JC(icm_spd_, 8192), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 4096), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 1024), JC(icm_spd_, 0), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 256), JC(icm_spd_, 0)
    };

    InterConnMatSpd::formulas_t fB = {
        JC(icm_spd_, 128), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 128),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0)
    };

    InterConnMatSpd::formulas_t fC = {
        JC(icm_spd_, 0), JC(icm_spd_, 0.004661), JC(icm_spd_, 24.99), JC(icm_spd_, 80.17), JC(icm_spd_, 48.49),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 0.004661), JC(icm_spd_, 24.99), JC(icm_spd_, 80.17), JC(icm_spd_, 48.49)
    };

    InterConnMatSpd::formulas_t fD = {
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 1), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 1),
        JC(icm_spd_, 0), JC(icm_spd_, -1), JC(icm_spd_, 1), JC(icm_spd_, 0)
    };

    icm_spd_.setA(fA);
    icm_spd_.setB(fB);
    icm_spd_.setC(fC);
    icm_spd_.setD(fD);
    jacl::system::Continuous<InterConnMatSpd> icm_spd_sys(&icm_spd_);
    hinf_sc_ = new HInfSC(&icm_spd_sys, 3.1);
    jacl::lti_common::StateSpacePack K( hinf_sc_->solve() );
    k_spd_.setA(std::get<0>(K));
    k_spd_.setB(std::get<1>(K));
    k_spd_.setC(std::get<2>(K));
    k_spd_.setD(std::get<3>(K));    

    k_spd_.A().print("Kspd_A : ");
    k_spd_.B().print("Kspd_B : ");
    k_spd_.C().print("Kspd_C : ");
    k_spd_.D().print("Kspd_D : ");

    spdctrl_plt_.init();
    spdctrl_plt_.setTitle("Speed Controller");
    spdctrl_plt_.setDelay() = SAMPLING_PERIOD;
    spdctrl_plt_.setPlotName({"Velocity Error", "Voltage"});
}

void MainWindow::setupNLP(){
    {        
        NLPendulum::formulas_t state = {
            JE(nlp_, _p(ix2)),
            JE(nlp_, -1*(_p(ipg)/_p(ipl))*sin(_p(ix1)) - (_p(ipk)/_p(ipm))*_p(ix1)), 
        };
        NLPendulum::formulas_t output = {
            JE(nlp_, _p(ix1))
        };
    }
}

double MainWindow::angularSpeed2Voltage(double _speed, double _torque){
    return _torque*ra.perturbed/ki.perturbed + ki.perturbed*_speed;
}

void MainWindow::closedLoopProcess(){
    arma::vec out(3,1,arma::fill::zeros),in(1,1,arma::fill::zeros);
    arma::vec err(in),est(out);
    din_ = in;
    control_mode_ = jacl::traits::toUType(jacl::ControllerDialog::ControlMode::Position);
//    arma::mat last_err(err);
//    arma::mat diff(err);
//    auto Kp(10.), Kd(1.);
    // ifd_.init({{-.76,-.65}, {-.63,-.51}, {-.86,-.72}});
    // sifd_.init({{-.11,-.15}}, "SIFD", {"Est. Curr.","Est. Spd.","Est. Pos."});
    arma::cx_vec p = jacl::lti_common::poles(dsimo_);    
    p.print("Discrete DC motor poles : ");
    p = jacl::lti_common::poles(simo_);
    p.print("Continuous DC motor poles : ");        
    while(cl_status_){
        //-- PD Control
       /*err(1) = ref_(1) - out(1);
       diff(1) = err(1) - last_err(1);
       last_err(1) = err(1);
       err.print("Error : ");
       in(0) = angularSpeed2Voltage(Kp*err(1) + Kd*diff(1),.0);
       if(in(0) > 12.0){
           in(0) = 12.0;
       }
       if(in(0) < -12.0){
           in(0) = -12.0;
       }*/

        //-- Error between reference and output
        err(control_mode_) = ref_(control_mode_) - out(control_mode_);

        //-- H-infinity Controller
        // if(control_mode_ == jacl::traits::toUType(jacl::ControllerDialog::ControlMode::Position)){
        //     err(1) = 0;
        //     in.submat(0,0,0,0) = posctrl_sys_.convolve(err.submat(0,0,0,0));
        // }else if(control_mode_ == jacl::traits::toUType(jacl::ControllerDialog::ControlMode::Velocity)){
        //     err(0) = 0;
        //     in.submat(0,0,0,0) = spdctrl_sys_.convolve(err.submat(1,0,1,0));
        // }else{
        //     //-- do nothing
        // }

        //-- Continuous
        // out = csys_.convolve(in);
        // est = cobserver_.convolve(in, out);

        //-- Discrete
        din_ = pos_dctrl_sys_.convolve(err);

        //-- saturation
        if(din_(0) > 24.)
            din_(0) = 24.;
        else if(din_(0) < -24.)
            din_(0) = -24.;
        
        out = dsys_simo_.convolve(din_);
        makeFault(&out);
        SIFD::diag_pack_t diag_pack = sifd_.detect(din_, out);        
        Q_EMIT setDataMonitor(
            arma::join_cols(out,
                arma::join_cols(
                    std::get<0>(diag_pack[0]),
                    arma::join_cols(std::get<0>(diag_pack[1]),
                                    std::get<0>(diag_pack[2])
                    )
                )
            )
        );
        arma::Col<uint8_t> status{std::get<2>(diag_pack[0]), std::get<2>(diag_pack[1]), std::get<2>(diag_pack[2])};
        Q_EMIT setSensorStatus(status);
        boost::this_thread::sleep_for(boost::chrono::milliseconds(20));
    }
}
