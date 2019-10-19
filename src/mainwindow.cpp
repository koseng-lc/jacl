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
    , jm(.00000072)
    , ki(.052)
    , la(.00062)
    , kb(.052)
    , ra(2.07)
    , ss_(bm, jm, ki, la, kb, ra)
    , sim_(3,2,2){

    ui->setupUi(this);

    this->removeToolBar(ui->mainToolBar);

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

    sim_.setTitle("DC Motor Simulation");
    sim_.setDelay() = .0;
    sim_.init();    
    sim_.setStateSpace(ss_.A(), ss_.B(), ss_.C(), ss_.D());

    setupWidgets();
    setupActions();

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
    params_dsb_[iJm]->setSingleStep(1e-8);
    params_dsb_[iJm]->setDecimals(8);
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

    command_gl_->addWidget(perturb_pb_ , 0,0,1,1);
    command_gl_->addWidget(reset_pb_   , 1,0,1,1);
    command_gl_->addWidget(simulate_pb_, 2,0,1,1);

    command_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),3,1);

    command_gb_->setLayout(command_gl_);
    command_gb_->setTitle(tr("Command"));

    //-- Main

    main_widget_ = new QWidget;

    main_layout_ = new QGridLayout;

    main_layout_->addWidget(params_gb_,  0,0,2,1);
    main_layout_->addWidget(command_gb_, 0,1,1,1);

    main_layout_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),1,2);
    main_layout_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),2,2);

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
}

void MainWindow::perturbAct(){

    for(int idx(iBm); idx <= iRa; idx++)
        ss_.param(idx) = params_dsb_[idx]->value();

    ss_.formulaToMat();

    sim_.setStateSpace(ss_.A(), ss_.B(), ss_.C(), ss_.D());
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

    sim_.simulate();
}
