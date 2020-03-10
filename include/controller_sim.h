#pragma once

#include "sim.h"
#include "controller.h"

namespace jacl{

template <class _StateSpace>
class ControllerSim: public Sim<_StateSpace>{
public:
    ControllerSim(_StateSpace* _ss, double _time_step=1e-4);
    ~ControllerSim();
    auto updateVariables() -> void;
    auto setInput(const arma::mat& _in) -> void;
    auto getOutputSig() const -> arma::mat;

protected:
    arma::mat signalCalc();

private:
    _StateSpace* ss_;
    arma::mat u_;
    arma::mat prev_u_;
    arma::mat state_;
    arma::mat prev_state_;
    arma::mat state_trans_;
    arma::mat output_;

    double dt_;
};

template <class _StateSpace>
ControllerSim<_StateSpace>::ControllerSim(_StateSpace* _ss, double _time_step)
    : Sim<_StateSpace>(_StateSpace::n_inputs + _StateSpace::n_outputs)
    , ss_(_ss)
    , u_(_StateSpace::n_inputs, 1, arma::fill::zeros)
    , prev_u_(u_)
    , state_(_StateSpace::n_states, 1, arma::fill::zeros)
    , prev_state_(state_)
    , output_(_StateSpace::n_outputs, 1, arma::fill::zeros)
    , dt_(_time_step){

}

template <class _StateSpace>
ControllerSim<_StateSpace>::~ControllerSim(){
    ss_ = nullptr;
}

template <class _StateSpace>
auto ControllerSim<_StateSpace>::setInput(const arma::mat& _in) -> void{
    u_ = _in;
}

template <class _StateSpace>
auto ControllerSim<_StateSpace>::updateVariables() -> void{
    state_trans_ = arma::expmat(ss_->A() * dt_);
}

template <class _StateSpace>
auto ControllerSim<_StateSpace>::signalCalc() -> arma::mat{
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

    return arma::join_cols(u_, output_);
}

template <class _StateSpace>
auto ControllerSim<_StateSpace>::getOutputSig() const -> arma::mat{
    return output_;
}

}
