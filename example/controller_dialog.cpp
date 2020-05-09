#include "controller_dialog.h"

using namespace jacl;

ControllerDialog::ControllerDialog(QWidget* parent)
    : QDialog(parent){

    ref_label_[traits::toUType(RefType::Position)] = new QLabel;
    ref_label_[traits::toUType(RefType::Position)]->setText(tr("Position"));
    ref_dsb_[traits::toUType(RefType::Position)] = new QDoubleSpinBox;
    ref_dsb_[traits::toUType(RefType::Position)]->setDecimals(3);
    ref_dsb_[traits::toUType(RefType::Position)]->setMinimum(-180.);
    ref_dsb_[traits::toUType(RefType::Position)]->setMaximum(180.);
    ref_dsb_[traits::toUType(RefType::Position)]->setSingleStep(1e-3);
    ref_dsb_[traits::toUType(RefType::Position)]->setValue(.0);
    ref_dsb_[traits::toUType(RefType::Position)]->setEnabled(false);

    ref_label_[traits::toUType(RefType::Velocity)] = new QLabel;
    ref_label_[traits::toUType(RefType::Velocity)]->setText(tr("Velocity"));
    ref_dsb_[traits::toUType(RefType::Velocity)] = new QDoubleSpinBox;
    ref_dsb_[traits::toUType(RefType::Velocity)]->setDecimals(3);
    ref_dsb_[traits::toUType(RefType::Velocity)]->setMinimum(-100.);
    ref_dsb_[traits::toUType(RefType::Velocity)]->setMaximum(100.);
    ref_dsb_[traits::toUType(RefType::Velocity)]->setSingleStep(1e-3);
    ref_dsb_[traits::toUType(RefType::Velocity)]->setValue(.0);

    ref_label_[traits::toUType(RefType::Current)] = new QLabel;
    ref_label_[traits::toUType(RefType::Current)]->setText(tr("Current"));
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

    //-- Data Monitor
    dm_out_label_[traits::toUType(ControlMode::Position)] = new QLabel;
    dm_out_label_[traits::toUType(ControlMode::Position)]->setText(tr("Position : "));
    dm_out_lined_[traits::toUType(ControlMode::Position)] = new QLineEdit;
    dm_out_lined_[traits::toUType(ControlMode::Position)]->setText(QString::number(0.));
    dm_out_lined_[traits::toUType(ControlMode::Position)]->setReadOnly(true);

    dm_out_label_[traits::toUType(ControlMode::Velocity)] = new QLabel;
    dm_out_label_[traits::toUType(ControlMode::Velocity)]->setText(tr("Velocity : "));
    dm_out_lined_[traits::toUType(ControlMode::Velocity)] = new QLineEdit;
    dm_out_lined_[traits::toUType(ControlMode::Velocity)]->setText(QString::number(0.));
    dm_out_lined_[traits::toUType(ControlMode::Velocity)]->setReadOnly(true);

    dm_out_label_[traits::toUType(ControlMode::Current)] = new QLabel;
    dm_out_label_[traits::toUType(ControlMode::Current)]->setText(tr("Current : "));
    dm_out_lined_[traits::toUType(ControlMode::Current)] = new QLineEdit;
    dm_out_lined_[traits::toUType(ControlMode::Current)]->setText(QString::number(0.));
    dm_out_lined_[traits::toUType(ControlMode::Current)]->setReadOnly(true);

    dm_est_label_[traits::toUType(ControlMode::Position)] = new QLabel;
    dm_est_label_[traits::toUType(ControlMode::Position)]->setText(tr("Est. Position : "));
    dm_est_lined_[traits::toUType(ControlMode::Position)] = new QLineEdit;
    dm_est_lined_[traits::toUType(ControlMode::Position)]->setText(QString::number(0.));
    dm_est_lined_[traits::toUType(ControlMode::Position)]->setReadOnly(true);

    dm_est_label_[traits::toUType(ControlMode::Velocity)] = new QLabel;
    dm_est_label_[traits::toUType(ControlMode::Velocity)]->setText(tr("Est. Velocity : "));
    dm_est_lined_[traits::toUType(ControlMode::Velocity)] = new QLineEdit;
    dm_est_lined_[traits::toUType(ControlMode::Velocity)]->setText(QString::number(0.));
    dm_est_lined_[traits::toUType(ControlMode::Velocity)]->setReadOnly(true);

    dm_est_label_[traits::toUType(ControlMode::Current)] = new QLabel;
    dm_est_label_[traits::toUType(ControlMode::Current)]->setText(tr("Est. Current : "));
    dm_est_lined_[traits::toUType(ControlMode::Current)] = new QLineEdit;
    dm_est_lined_[traits::toUType(ControlMode::Current)]->setText(QString::number(0.));
    dm_est_lined_[traits::toUType(ControlMode::Current)]->setReadOnly(true);

    dm_gl_ = new QGridLayout;
    dm_gl_->addWidget(dm_out_label_[traits::toUType(ControlMode::Position)], 0,0,1,1);
    dm_gl_->addWidget(dm_out_lined_[traits::toUType(ControlMode::Position)], 0,1,1,1);
    dm_gl_->addWidget(dm_out_label_[traits::toUType(ControlMode::Velocity)], 1,0,1,1);
    dm_gl_->addWidget(dm_out_lined_[traits::toUType(ControlMode::Velocity)], 1,1,1,1);
    dm_gl_->addWidget(dm_out_label_[traits::toUType(ControlMode::Current)], 2,0,1,1);
    dm_gl_->addWidget(dm_out_lined_[traits::toUType(ControlMode::Current)], 2,1,1,1);
    dm_gl_->addWidget(dm_est_label_[traits::toUType(ControlMode::Position)], 3,0,1,1);
    dm_gl_->addWidget(dm_est_lined_[traits::toUType(ControlMode::Position)], 3,1,1,1);
    dm_gl_->addWidget(dm_est_label_[traits::toUType(ControlMode::Velocity)], 4,0,1,1);
    dm_gl_->addWidget(dm_est_lined_[traits::toUType(ControlMode::Velocity)], 4,1,1,1);
    dm_gl_->addWidget(dm_est_label_[traits::toUType(ControlMode::Current)], 5,0,1,1);
    dm_gl_->addWidget(dm_est_lined_[traits::toUType(ControlMode::Current)], 5,1,1,1);
    dm_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),6,2);

    dm_gb_ = new QGroupBox;
    dm_gb_->setTitle(tr("Data Monitor"));
    dm_gb_->setLayout(dm_gl_);

    main_layout_ = new QGridLayout;
    main_layout_->addWidget(ref_gb_,0,0,1,1);
    main_layout_->addWidget(mode_gb_,1,0,2,1);
    main_layout_->addWidget(dm_gb_,0,1,2,1);
    main_layout_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),3,1);
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

void ControllerDialog::setDataMonitor(QVector<double> _data){
    for(int i(0); i < 6; i++){
        if(i < 3)
            dm_out_lined_[i]->setText(QString::number(_data[i]));
        else
            dm_est_lined_[i-3]->setText(QString::number(_data[i]));
    }
}