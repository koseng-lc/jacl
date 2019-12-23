/**
*   @author : koseng (Lintang)
*   @brief : Open-Loop System Simulator
*/

#pragma once

#include "sim.h"
#include "state_space.h"

namespace jacl{

template <class _StateSpace>
class SystemSim: public Sim<_StateSpace>{
public:
    SystemSim(_StateSpace* _ss, double _time_step = 1e-4);
    ~SystemSim();

    void setInput(const arma::mat& _in);
    void updateVariables();

protected:
    arma::mat signalCalc();

private:
    _StateSpace* ss_;

    arma::mat u_;
    arma::mat prev_u_;

    arma::mat Ad_;
    arma::mat Bd_;
    arma::mat Cd_;
    arma::mat Dd_;

    //-- System Variable
    arma::mat state_;
    arma::mat prev_state_;
    arma::mat state_trans_;
    arma::mat output_;

    void calcState();
    void calcOutput();

    double dt_;

};

template <class _StateSpace>
SystemSim<_StateSpace>::SystemSim(_StateSpace* _ss, double _time_step)
    : Sim<_StateSpace >(_ss->A().n_rows +_ss->B().n_cols + _ss->C().n_rows)
    , ss_(_ss)
    , u_(_ss->B().n_cols, 1, arma::fill::zeros)
    , prev_u_(u_)
    , state_(_ss->A().n_rows, 1, arma::fill::zeros)
    , prev_state_(state_)
    , dt_(_time_step){

    // updateVariables();

}

template <class _StateSpace>
SystemSim<_StateSpace>::~SystemSim(){
    ss_ = nullptr;
}

template <class _StateSpace>
void SystemSim<_StateSpace>::setInput(const arma::mat& _in){
    u_ = _in;
}

template <class _StateSpace>
void SystemSim<_StateSpace>::updateVariables(){
    state_trans_ = arma::expmat(ss_->A() * dt_);
}

template <class _StateSpace>
arma::mat SystemSim<_StateSpace>::signalCalc(){

    static arma::mat term1, term2, term3, term4, term5;

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
