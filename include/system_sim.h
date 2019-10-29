/**
*   @author : koseng (Lintang)
*   @brief : Open-Loop System Simulator
*/

#pragma once

#include "sim.h"
#include "state_space.h"

namespace JACL{

template <class SSpace>
class SystemSim: public Sim<SSpace>{
public:

    SystemSim(SSpace* _ss, double _time_step = 1e-4);
    ~SystemSim();

    void setInput(const Mat& _in);
    void updateVariables();

protected:

    Mat signalCalc();

private:

    SSpace* ss_;

    Mat u_;
    Mat prev_u_;

    //-- System Variable
    Mat state_;
    Mat prev_state_;
    Mat state_trans_;
    Mat output_;

    void calcState();
    void calcOutput();

    double dt_;

};

template <class SSpace>
SystemSim<SSpace>::SystemSim(SSpace* _ss, double _time_step)
    : Sim<SSpace >(_ss->A().n_rows +_ss->B().n_cols + _ss->C().n_rows)
    , ss_(_ss)
    , u_(_ss->B().n_cols, 1, arma::fill::zeros)
    , prev_u_(u_)
    , state_(_ss->A().n_rows, 1, arma::fill::zeros)
    , prev_state_(state_)
    , dt_(_time_step){

    updateVariables();

}

template <class SSpace>
SystemSim<SSpace>::~SystemSim(){
    ss_ = nullptr;
}

template <class SSpace>
void SystemSim<SSpace>::setInput(const Mat& _in){

    u_ = _in;
}

template <class SSpace>
void SystemSim<SSpace>::updateVariables(){

    state_trans_ = arma::expmat(ss_->A() * dt_);
}

template <class SSpace>
Mat SystemSim<SSpace>::signalCalc(){
    static Mat term1, term2, term3, term4, term5;

    term1 = state_trans_ * prev_state_;
    term2 = ss_->B() * u_;
    term3 = state_trans_ * ss_->B() * prev_u_;
    term4 = (term2 + term3) * (dt_ * .5);

    state_ = term1 + term4;

    term5 = ss_->C() * prev_state_;
    output_ = term5 + (ss_->D() * u_);

    prev_state_ = state_;
    prev_u_ = u_;

    return arma::join_cols(state_, arma::join_cols(u_, output_));
}

}
