#include "controller_dialog.h"

using namespace jacl;

ControllerDialog::ControllerDialog(QWidget* parent)
    : QDialog(parent){

    ref_label_[Position] = new QLabel;
    ref_label_[Position]->setText(tr("Position"));
    ref_dsb_[Position] = new QDoubleSpinBox;
    ref_dsb_[Position]->setDecimals(0);
    ref_dsb_[Position]->setMinimum(-180.);
    ref_dsb_[Position]->setMaximum(180.);
    ref_dsb_[Position]->setSingleStep(1);
    ref_dsb_[Position]->setValue(.0);

    ref_label_[Velocity] = new QLabel;
    ref_label_[Velocity]->setText(tr("Velocity"));
    ref_dsb_[Velocity] = new QDoubleSpinBox;
    ref_dsb_[Velocity]->setDecimals(3);
    ref_dsb_[Velocity]->setMinimum(-100.);
    ref_dsb_[Velocity]->setMaximum(100.);
    ref_dsb_[Velocity]->setSingleStep(1e-3);
    ref_dsb_[Velocity]->setValue(.0);

    ref_label_[Current] = new QLabel;
    ref_label_[Current]->setText(tr("Current"));
    ref_dsb_[Current] = new QDoubleSpinBox;
    ref_dsb_[Current]->setDecimals(3);
    ref_dsb_[Current]->setMinimum(-2.);
    ref_dsb_[Current]->setMaximum(2.);
    ref_dsb_[Current]->setSingleStep(1e-3);
    ref_dsb_[Current]->setValue(.0);

    set_ref_pb_ = new QPushButton;
    set_ref_pb_->setText(tr("Set Ref"));

    ref_gl_ = new QGridLayout;
    ref_gl_->addWidget(ref_label_[Position],0,0,1,1);
    ref_gl_->addWidget(ref_dsb_[Position],0,1,1,1);
    ref_gl_->addWidget(ref_label_[Velocity],1,0,1,1);
    ref_gl_->addWidget(ref_dsb_[Velocity],1,1,1,1);
    ref_gl_->addWidget(ref_label_[Current],2,0,1,1);
    ref_gl_->addWidget(ref_dsb_[Current],2,1,1,1);
    ref_gl_->addWidget(set_ref_pb_,3,0,1,3);
    ref_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),4,4);

    ref_gb_ = new QGroupBox;
    ref_gb_->setTitle(tr("Reference"));
    ref_gb_->setLayout(ref_gl_);

    main_layout_ = new QGridLayout;
    main_layout_->addWidget(ref_gb_,0,0,1,1);
    main_layout_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),1,1);
    this->setLayout(main_layout_);
    this->setWindowTitle(tr("Controller"));

    //-- Actions
    connect(set_ref_pb_,SIGNAL(clicked()), this, SLOT(setRefAct()));

}

ControllerDialog::~ControllerDialog(){

}

void ControllerDialog::setRefAct(){
    Q_EMIT setRefSig();
}
