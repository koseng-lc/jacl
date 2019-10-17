/**
*   @author : koseng (Lintang)
*   @brief : GUI Source
*/

#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
    , ss_(bm, jm){

    ui->setupUi(this);

    this->removeToolBar(ui->mainToolBar);

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

    ss_.param(iBm) = .000048;
    ss_.param(iJm) = .00000072;
    ss_.param(iKi) = .052;
    ss_.param(iLa) = .00062;
    ss_.param(iKb) = .052;
    ss_.param(iRa) = 2.07;

    ss_.setA(fA);
    ss_.setB(fB);
    ss_.setC(fC);
    ss_.setD(fD);
    ss_.formulaToMat();

//    sim_ = new JACL::Simulator();

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
    params_dsb_[iBm]->setValue(.000048);

    params_label_[iJm] = new QLabel;
    params_label_[iJm]->setText(tr("Rotor Inertia : "));
    params_dsb_[iJm] = new QDoubleSpinBox;
    params_dsb_[iJm]->setValue(.0000072);

    params_label_[iKi] = new QLabel;
    params_label_[iKi]->setText(tr("Torque Constant : "));
    params_dsb_[iKi] = new QDoubleSpinBox;
    params_dsb_[iKi]->setValue(.052);

    params_label_[iLa] = new QLabel;
    params_label_[iLa]->setText(tr("Armature Inductance : "));
    params_dsb_[iLa] = new QDoubleSpinBox;
    params_dsb_[iLa]->setValue(.00062);

    params_label_[iKb] = new QLabel;
    params_label_[iKb]->setText(tr("Back EMF : "));
    params_dsb_[iKb] = new QDoubleSpinBox;
    params_dsb_[iKb]->setValue(.052);

    params_label_[iRa] = new QLabel;
    params_label_[iRa]->setText(tr("Armature Resistance : "));
    params_dsb_[iRa] = new QDoubleSpinBox;
    params_dsb_[iRa]->setValue(2.07);

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

    command_gl_->addWidget(perturb_pb_, 0,0,1,1);
    command_gl_->addWidget(reset_pb_  , 1,0,1,1);

    command_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),2,1);

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

}

void MainWindow::setupActions(){

}
