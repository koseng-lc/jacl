#include "controller_dialog.h"

using namespace jacl;

ControllerDialog::ControllerDialog(QWidget* parent)
    : QDialog(parent){

    ref_label_[traits::toUType(RefType::Position)] = new QLabel;
    ref_label_[traits::toUType(RefType::Position)]->setText(tr("RefType::Position"));
    ref_dsb_[traits::toUType(RefType::Position)] = new QDoubleSpinBox;
    ref_dsb_[traits::toUType(RefType::Position)]->setDecimals(3);
    ref_dsb_[traits::toUType(RefType::Position)]->setMinimum(-180.);
    ref_dsb_[traits::toUType(RefType::Position)]->setMaximum(180.);
    ref_dsb_[traits::toUType(RefType::Position)]->setSingleStep(1e-3);
    ref_dsb_[traits::toUType(RefType::Position)]->setValue(.0);
    ref_dsb_[traits::toUType(RefType::Position)]->setEnabled(false);

    ref_label_[traits::toUType(RefType::Velocity)] = new QLabel;
    ref_label_[traits::toUType(RefType::Velocity)]->setText(tr("RefType::Velocity"));
    ref_dsb_[traits::toUType(RefType::Velocity)] = new QDoubleSpinBox;
    ref_dsb_[traits::toUType(RefType::Velocity)]->setDecimals(3);
    ref_dsb_[traits::toUType(RefType::Velocity)]->setMinimum(-100.);
    ref_dsb_[traits::toUType(RefType::Velocity)]->setMaximum(100.);
    ref_dsb_[traits::toUType(RefType::Velocity)]->setSingleStep(1e-3);
    ref_dsb_[traits::toUType(RefType::Velocity)]->setValue(.0);

    ref_label_[traits::toUType(RefType::Current)] = new QLabel;
    ref_label_[traits::toUType(RefType::Current)]->setText(tr("RefType::Current"));
    ref_dsb_[traits::toUType(RefType::Current)] = new QDoubleSpinBox;
    ref_dsb_[traits::toUType(RefType::Current)]->setDecimals(3);
    ref_dsb_[traits::toUType(RefType::Current)]->setMinimum(-2.);
    ref_dsb_[traits::toUType(RefType::Current)]->setMaximum(2.);
    ref_dsb_[traits::toUType(RefType::Current)]->setSingleStep(1e-3);
    ref_dsb_[traits::toUType(RefType::Current)]->setValue(.0);

    set_ref_pb_ = new QPushButton;
    set_ref_pb_->setText(tr("Set Ref"));

    ref_gl_ = new QGridLayout;
    ref_gl_->addWidget(ref_label_[traits::toUType(RefType::Position)],0,0,1,1);
    ref_gl_->addWidget(ref_dsb_[traits::toUType(RefType::Position)],0,1,1,1);
    ref_gl_->addWidget(ref_label_[traits::toUType(RefType::Velocity)],1,0,1,1);
    ref_gl_->addWidget(ref_dsb_[traits::toUType(RefType::Velocity)],1,1,1,1);
    ref_gl_->addWidget(ref_label_[traits::toUType(RefType::Current)],2,0,1,1);
    ref_gl_->addWidget(ref_dsb_[traits::toUType(RefType::Current)],2,1,1,1);
    ref_gl_->addWidget(set_ref_pb_,3,0,1,3);
    ref_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),4,4);

    ref_gb_ = new QGroupBox;
    ref_gb_->setTitle(tr("Reference"));
    ref_gb_->setLayout(ref_gl_);

    //--
    mode_rb_[traits::toUType(ControlMode::Position)] = new QRadioButton;
    mode_rb_[traits::toUType(ControlMode::Position)]->setText(tr("Position Mode"));

    mode_rb_[traits::toUType(ControlMode::Velocity)] = new QRadioButton;
    mode_rb_[traits::toUType(ControlMode::Velocity)]->setText(tr("Velocity Mode"));
    mode_rb_[traits::toUType(ControlMode::Velocity)]->setChecked(true);

    mode_gl_ = new QGridLayout;
    mode_gl_->addWidget(mode_rb_[traits::toUType(ControlMode::Position)],0,0,1,1);
    mode_gl_->addWidget(mode_rb_[traits::toUType(ControlMode::Velocity)],1,0,1,1);
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
    connect(mode_rb_[traits::toUType(ControlMode::Position)], SIGNAL(clicked()), this, SLOT(setModeAct()));
    connect(mode_rb_[traits::toUType(ControlMode::Velocity)], SIGNAL(clicked()), this, SLOT(setModeAct()));
}

ControllerDialog::~ControllerDialog(){

}

void ControllerDialog::setRefAct(){
    Q_EMIT setRefSig();
}

void ControllerDialog::setModeAct(){
    if(sender() == mode_rb_[traits::toUType(ControlMode::Position)]){
        Q_EMIT setModeSig(traits::toUType(ControlMode::Position));
        ref_dsb_[traits::toUType(RefType::Position)]->setEnabled(true);
        ref_dsb_[traits::toUType(RefType::Velocity)]->setValue(.0);
        ref_dsb_[traits::toUType(RefType::Velocity)]->setEnabled(false);
    }else{
        Q_EMIT setModeSig(traits::toUType(ControlMode::Velocity));
        ref_dsb_[traits::toUType(RefType::Position)]->setValue(.0);
        ref_dsb_[traits::toUType(RefType::Position)]->setEnabled(false);
        ref_dsb_[traits::toUType(RefType::Velocity)]->setEnabled(true);
    }
}
