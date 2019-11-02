/**
*   @author : koseng (Lintang)
*   @brief : GUI Source
*/

#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
    , bm(.000048)
    , jm(7.2e-6)
    , ki(.052)
    , la(.00062)
    , kb(.052)
    , ra(2.07)
    , ss_(bm, jm, ki, la, kb, ra)
//    , obs_gain_({{ 101.0,     1.0},
//                 {-311.9, -3143.4},
//                 { 139.8,  1367.6}})
    , obs_gain_({{ .0, .0},
                 { .0, .0},
                 { .0, .0}})
//    , sim_(3,2,2)
    , system_sim_(&ss_)
    , observer_sim_(&ss_, obs_gain_){

    ui->setupUi(this);

    this->removeToolBar(ui->mainToolBar);

    QFile style_file("../gui/dark_style_lintang.qss");
    style_file.open(QFile::ReadOnly);
    QString style_sheet(style_file.readAll());
    this->setStyleSheet(style_sheet);

    /*
     * Notice ! : ss keyword is preserved variable name
     */
    SS::Formulas fA{
        JACL_CONST_SS(.0),                    JACL_CONST_SS(1.0),                     JACL_CONST_SS(.0),
        JACL_CONST_SS(.0), JACL_SS(-ss.param(iBm)/ss.param(iJm)),  JACL_SS(ss.param(iKi)/ss.param(iJm)),
        JACL_CONST_SS(.0), JACL_SS(-ss.param(iKb)/ss.param(iLa)), JACL_SS(-ss.param(iRa)/ss.param(iLa))
    };

    SS::Formulas fB{
                 JACL_CONST_SS(.0),           JACL_CONST_SS(.0),
                 JACL_CONST_SS(.0), JACL_SS(-1.0/ss.param(iJm)),
        JACL_SS(1.0/ss.param(iLa)),           JACL_CONST_SS(.0)
    };

    SS::Formulas fC{
        JACL_CONST_SS(1.0),  JACL_CONST_SS(.0), JACL_CONST_SS(.0),
         JACL_CONST_SS(.0), JACL_CONST_SS(1.0), JACL_CONST_SS(.0)
    };

    SS::Formulas fD{
        JACL_CONST_SS(.0), JACL_CONST_SS(.0),
        JACL_CONST_SS(.0), JACL_CONST_SS(.0)
    };


    ss_.param(iBm) = bm.nominal;
    ss_.param(iJm) = jm.nominal;
    ss_.param(iKi) = ki.nominal;
    ss_.param(iLa) = la.nominal;
    ss_.param(iKb) = kb.nominal;
    ss_.param(iRa) = ra.nominal;

    ss_.setA(fA);
    ss_.setB(fB);
    ss_.setC(fC);
    ss_.setD(fD);
    ss_.formulaToMat();

    std::cout << "Preparing system ..." << std::endl;
    ss_.A().print("A : ");
    ss_.B().print("B : ");
    ss_.C().print("C : ");
    ss_.D().print("D : ");

    // Test KautskyNichols
    JACL::Mat observer_K;
//    JACL::Mat poles{-10,-9,-5};
    JACL::Mat poles{-50,-51,-52};
    JACL::PolePlacement::KautskyNichols(&ss_, poles, &observer_K);

    observer_K.print("Observer Gain : ");

    observer_sim_.setGain(observer_K.t());

    system_sim_.init();
    system_sim_.setTitle("DC Motor");
    system_sim_.setDelay() = .02;
    system_sim_.setPlotName({"Angular Position", "Angular Velocity", "Current",
                       "Voltage In", "Torque In",
                       "Angular Position", "Angular Velocity"});
    system_sim_.updateVariables();

    observer_sim_.init();
    observer_sim_.setTitle("Full-order Luenberger Observer of DC Motor");
    observer_sim_.setDelay() = .02;
    observer_sim_.setPlotName({"Est. Position", "Est. Velocity", "Est. Current"
                              ,"Est. Out Position", "Est. Out Velocity"});
    observer_sim_.updateVariables();

//    sim_.init();
//    sim_.setTitle("DC Motor Simulation");
//    sim_.setDelay() = .0;
//    sim_.setPlotName({"Angular Position", "Angular Velocity", "Current",
//                      "Torque In", "Voltage In",
//                      "Angular Position", "Angular Velocity"});
//    sim_.setStateSpace(ss_.A(), ss_.B(), ss_.C(), ss_.D());    

    setupWidgets();
    setupActions();

    JACL::StateSpace<12,4,4> another_ss;

}

MainWindow::~MainWindow(){

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
    set_input_pb_->setText(tr("Set Input"));

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
    input_gb_->setTitle(tr("Input"));

    //-- Fault

    fault_gb_ = new QGroupBox;

    fault_gl_ = new QGridLayout;

    target_label_ = new QLabel;
    target_label_->setText(tr("Sensor : "));

    target_cb_ = new QComboBox;
    target_cb_->addItem(tr("Position"));
    target_cb_->addItem(tr("Velocity"));
    target_cb_->addItem(tr("Current"));

    details_pb_ = new QPushButton;
    details_pb_->setText(tr("Details"));

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
    scale_dsb_->setValue(.0);
    scale_dial_ = new QDial;

    dead_zone_label_ = new QLabel;
    dead_zone_label_->setText(tr("Dead Zone : "));
    dead_zone_dsb_ =  new QDoubleSpinBox;
    // Settings DSB params
    dead_zone_dsb_->setValue(.0);
    dead_zone_dial_ = new QDial;

    fault_gl_->addWidget(target_label_,    0,0,1,1);
    fault_gl_->addWidget(target_cb_,       0,1,1,1);
    fault_gl_->addWidget(details_pb_,      0,2,1,1);
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

    QImageReader image_reader("../gui/Logo_Universitas_Gadjah_Mada.png");
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
}

void MainWindow::perturbAct(){

    for(int idx(iBm); idx <= iRa; idx++)
        ss_.param(idx) = params_dsb_[idx]->value();

    ss_.formulaToMat();
    system_sim_.updateVariables();
    observer_sim_.updateVariables();
//    sim_.setStateSpace(ss_.A(), ss_.B(), ss_.C(), ss_.D());
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

//    sim_.simulate();
    system_sim_.simulate();
    observer_sim_.simulate();
}

void MainWindow::setInputAct(){
    JACL::Mat in(2, 1);
    in(0) = voltage_in_dsb_->value();
    in(1) = torque_in_dsb_->value();
//    sim_.setInput(in);
    system_sim_.setInput(in);
    observer_sim_.setInput(in);
}

void MainWindow::biasDialConv(double _val){
    bias_dial_->setValue((int)_val);
}

void MainWindow::scaleDialConv(double _val){
    scale_dial_->setValue((int)_val);
}

void MainWindow::deadZoneDialConv(double _val){
    dead_zone_dial_->setValue((int)_val);
}

void MainWindow::biasDSBConv(int _val){
    bias_dsb_->setValue((double)_val);
}

void MainWindow::scaleDSBConv(int _val){
    scale_dsb_->setValue((double)_val);
}

void MainWindow::deadZoneDSBConv(int _val){
    dead_zone_dsb_->setValue((double)_val);
}
