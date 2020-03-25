#pragma once

#include "sim.h"

namespace jacl{

template <class _Controller>
class PosCtrlSim final:public Sim<_Controller>{
public:
public:
    PosCtrlSim(_Controller* _k, double _time_step=1e-4);
    ~PosCtrlSim();
    auto updateVariables() -> void;
    auto setInput(const arma::mat& _in) -> void;
    auto getOutputSig() const -> arma::mat;
    auto propagate(const arma::mat& _in) -> arma::mat;

protected:
    auto signalCalc() -> arma::mat override;
    auto getSig() -> arma::mat override;

private:
    _Controller* k_;
    arma::mat u_;
    arma::mat prev_u_;
    arma::mat state_;
    arma::mat prev_state_;
    arma::mat state_trans_;
    arma::mat output_;

    arma::mat plotted_sig_;

    double dt_;
};

template <class _Controller>
PosCtrlSim<_Controller>::PosCtrlSim(_Controller* _k, double _time_step)
    : Sim<_Controller>(_Controller::n_inputs + _Controller::n_outputs)
    , k_(_k)
    , u_(_Controller::n_inputs, 1, arma::fill::zeros)
    , prev_u_(u_)
    , state_(_Controller::n_states, 1, arma::fill::zeros)
    , prev_state_(state_)
    , output_(_Controller::n_outputs, 1, arma::fill::zeros)
    , dt_(_time_step){

}

template <class _Controller>
PosCtrlSim<_Controller>::~PosCtrlSim(){
    k_ = nullptr;
}

template <class _Controller>
auto PosCtrlSim<_Controller>::setInput(const arma::mat& _in) -> void{
    u_ = _in;
}

template <class _Controller>
auto PosCtrlSim<_Controller>::updateVariables() -> void{
    state_trans_ = arma::expmat(k_->A() * dt_);
}

template <class _Controller>
auto PosCtrlSim<_Controller>::signalCalc() -> arma::mat{
    boost::mutex::scoped_lock lk(this->sig_mtx_);
    static arma::mat term1, term2, term3, term4, term5;

    term1 = state_trans_ * prev_state_;
    term2 = k_->B() * u_;
    term3 = state_trans_ * k_->B() * prev_u_;
    term4 = (term2 + term3) * (dt_ * .5);

    state_ = term1 + term4;

    term5 = k_->C() * prev_state_;
    output_ = term5 + (k_->D() * u_);

    prev_state_ = state_;
    prev_u_ = u_;

    return arma::join_cols(u_, output_);
}

template <class _Controller>
auto PosCtrlSim<_Controller>::getSig() -> arma::mat{
    boost::mutex::scoped_lock lk(this->sig_mtx_);
    return plotted_sig_;
}

template <class _Controller>
auto PosCtrlSim<_Controller>::getOutputSig() const -> arma::mat{
    boost::mutex::scoped_lock lk(this->sig_mtx_);
    return output_;
}

template <class _Controller>
auto PosCtrlSim<_Controller>::propagate(const arma::mat &_in) -> arma::mat{
    setInput(_in);
    plotted_sig_ = signalCalc();
    return output_;
}

}
