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
    inline double getVoltage() const{
        return input_dsb_[Voltage]->value();
    }
    inline double getTorque() const{
        return input_dsb_[Torque]->value();
    }
private:
    enum InputType{
        Voltage,
        Torque
    };
    QGridLayout* main_layout_;
    QGroupBox* input_gb_;
    QLabel* input_label_[2];
    QDoubleSpinBox* input_dsb_[2];
    QPushButton* set_input_pb_;
    QGridLayout* input_gl_;
Q_SIGNALS:
    void setInputSig();
private Q_SLOTS:
    void setInputAct();
};

}
