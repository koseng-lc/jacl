#pragma once

#include <QDialog>
#include <QGridLayout>
#include <QGroupBox>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QLabel>

namespace jacl{

class ControllerDialog: public QDialog{
    Q_OBJECT
public:
    ControllerDialog(QWidget* parent=nullptr);
    ~ControllerDialog();
    inline double getPosition() const{
        return ref_dsb_[Position]->value();
    }
    inline double getVelocity() const{
        return ref_dsb_[Velocity]->value();
    }
    inline double getCurrent() const{
        return ref_dsb_[Current]->value();
    }
private:
    enum RefType{
        Position,
        Velocity,
        Current
    };
    QGridLayout* main_layout_;
    QGroupBox* ref_gb_;
    QLabel* ref_label_[3];
    QDoubleSpinBox* ref_dsb_[3];
    QPushButton* set_ref_pb_;
    QGridLayout* ref_gl_;
Q_SIGNALS:
    void setRefSig();
private Q_SLOTS:
    void setRefAct();
};

}
