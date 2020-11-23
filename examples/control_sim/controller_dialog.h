#pragma once

#include <QDialog>
#include <QGridLayout>
#include <QGroupBox>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QLabel>
#include <QRadioButton>
#include <QLineEdit>

#include <jacl/traits.hpp>

namespace jacl{

class ControllerDialog: public QDialog{
    Q_OBJECT
public:
    ControllerDialog(QWidget* parent=nullptr);
    ~ControllerDialog();
    enum class RefType: int{
        Position,
        Velocity,
        Current
    };
    enum class ControlMode: int{
        Position,
        Velocity,
        Current
    };

    inline auto getPosition() const -> double{
        return ref_dsb_[traits::toUType(RefType::Position)]->value();
    }
    inline auto getVelocity() const -> double{
        return ref_dsb_[traits::toUType(RefType::Velocity)]->value();
    }
    inline auto getCurrent() const -> double{
        return ref_dsb_[traits::toUType(RefType::Current)]->value();
    }    

private:
    QGridLayout* main_layout_;
    //-- Reference Part
    QGroupBox* ref_gb_;
    QLabel* ref_label_[3];
    QDoubleSpinBox* ref_dsb_[3];
    QPushButton* set_ref_pb_;
    QGridLayout* ref_gl_;
    //-- Mode Part
    QGroupBox* mode_gb_;
    QLabel* mode_label_[2];
    QRadioButton* mode_rb_[2];
    QGridLayout* mode_gl_;
    //-- Data Monitor
    QGroupBox* dm_gb_;
    QGridLayout* dm_gl_;
    QLabel* dm_out_label_[3];
    QLabel* dm_est_label_[3];
    QLineEdit* dm_out_lined_[3];
    QLineEdit* dm_est_lined_[3];
    //-- Fault Monitor
    QGroupBox* fm_gb_;
    QGridLayout* fm_gl_;
    QLabel* fm_label_[3];
    QLabel* fm_status_[3];

Q_SIGNALS:
    void setRefSig();
    void setModeSig(int);    

private Q_SLOTS:
    void setRefAct();
    void setModeAct();
    void setDataMonitor(arma::vec _data);
    void setSensorStatus(arma::Col<uint8_t> _status);
    
};

}
