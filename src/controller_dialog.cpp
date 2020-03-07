#include "controller_dialog.h"

using namespace jacl;

ControllerDialog::ControllerDialog(QWidget* parent)
    : QDialog(parent){

    input_label_[Voltage] = new QLabel;
    input_label_[Voltage]->setText(tr("Voltage"));
    input_dsb_[Voltage] = new QDoubleSpinBox;
    input_dsb_[Voltage]->setDecimals(3);
    input_dsb_[Voltage]->setMinimum(-24.);
    input_dsb_[Voltage]->setMaximum(-24.);
    input_dsb_[Voltage]->setSingleStep(1e-3);
    input_dsb_[Voltage]->setValue(.0);

    input_label_[Torque] = new QLabel;
    input_label_[Torque]->setText(tr("Torque"));
    input_dsb_[Torque] = new QDoubleSpinBox;
    input_dsb_[Torque]->setDecimals(3);
    input_dsb_[Torque]->setMinimum(0.);
    input_dsb_[Torque]->setMaximum(10.);
    input_dsb_[Torque]->setSingleStep(1e-3);
    input_dsb_[Torque]->setValue(.0);

    set_input_pb_ = new QPushButton;
    set_input_pb_->setText(tr("Set Input"));

    input_gl_ = new QGridLayout;
    input_gl_->addWidget(input_label_[Voltage],0,0,1,1);
    input_gl_->addWidget(input_dsb_[Voltage],0,1,1,1);
    input_gl_->addWidget(input_label_[Torque],1,0,1,1);
    input_gl_->addWidget(input_dsb_[Torque],1,1,1,1);
    input_gl_->addWidget(set_input_pb_,2,0,1,2);
    input_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),3,3);

    input_gb_ = new QGroupBox;
    input_gb_->setTitle(tr("Input"));
    input_gb_->setLayout(input_gl_);

    main_layout_ = new QGridLayout;
    main_layout_->addWidget(input_gb_,0,0,1,1);
    main_layout_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),1,1);
    this->setLayout(main_layout_);
    this->setWindowTitle(tr("Controller"));

    //-- Actions
    connect(set_input_pb_,SIGNAL(clicked()), this, SLOT(setInputAct()));

}

ControllerDialog::~ControllerDialog(){

}

void ControllerDialog::setInputAct(){
    Q_EMIT setInputSig();
}
