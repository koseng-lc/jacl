#include "controller_dialog.h"

using namespace jacl;

ControllerDialog::ControllerDialog(QWidget* parent)
    : QDialog(parent){

    ref_label_[Position] = new QLabel;
    ref_label_[Position]->setText(tr("Position"));
    ref_dsb_[Position] = new QDoubleSpinBox;
    ref_dsb_[Position]->setDecimals(3);
    ref_dsb_[Position]->setMinimum(-180.);
    ref_dsb_[Position]->setMaximum(180.);
    ref_dsb_[Position]->setSingleStep(1e-3);
    ref_dsb_[Position]->setValue(.0);
    ref_dsb_[Position]->setEnabled(false);

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

    //--
    mode_rb_[PosControl] = new QRadioButton;
    mode_rb_[PosControl]->setText(tr("Position Mode"));

    mode_rb_[VelControl] = new QRadioButton;
    mode_rb_[VelControl]->setText(tr("Velocity Mode"));
    mode_rb_[VelControl]->setChecked(true);

    mode_gl_ = new QGridLayout;
    mode_gl_->addWidget(mode_rb_[PosControl],0,0,1,1);
    mode_gl_->addWidget(mode_rb_[VelControl],1,0,1,1);
    mode_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),2,1);

    mode_gb_ = new QGroupBox;
    mode_gb_->setLayout(mode_gl_);
    mode_gb_->setTitle(tr("Control Mode"));
    //--

    main_layout_ = new QGridLayout;
    main_layout_->addWidget(ref_gb_,0,0,1,1);
    main_layout_->addWidget(mode_gb_,1,0,1,1);
    main_layout_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),2,1);
    this->setLayout(main_layout_);
    this->setWindowTitle(tr("Controller"));

    //-- Actions
    connect(set_ref_pb_,SIGNAL(clicked()), this, SLOT(setRefAct()));
    connect(mode_rb_[PosControl], SIGNAL(clicked()), this, SLOT(setModeAct()));
    connect(mode_rb_[VelControl], SIGNAL(clicked()), this, SLOT(setModeAct()));
}

ControllerDialog::~ControllerDialog(){

}

void ControllerDialog::setRefAct(){
    Q_EMIT setRefSig();
}

void ControllerDialog::setModeAct(){
    if(sender() == mode_rb_[PosControl]){
        Q_EMIT setModeSig((int)PosControl);
        ref_dsb_[Position]->setEnabled(true);
        ref_dsb_[Velocity]->setValue(.0);
        ref_dsb_[Velocity]->setEnabled(false);
    }else{
        Q_EMIT setModeSig((int)VelControl);
        ref_dsb_[Position]->setValue(.0);
        ref_dsb_[Position]->setEnabled(false);
        ref_dsb_[Velocity]->setEnabled(true);
    }
}
