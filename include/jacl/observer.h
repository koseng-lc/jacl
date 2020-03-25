/**
*   @author : koseng (Lintang)
*   @brief : Simple Luenberger Full-Order Observer
*/

#pragma once

#include "state_space.h"

namespace jacl{

namespace{
    namespace linalg = linear_algebra;
}

template <class _StateSpace>
class Observer{
public:
    Observer(_StateSpace* _ss, const arma::mat& _K, double _time_step = 1e-4);
    ~Observer();

    arma::mat calcStateEst();
    void calcOutputEst();

    void setInput(const arma::mat& _u){
        u_ = _u;
    }

    void setTimeStep(double _time_step){
        dt_ = _time_step;
    }

    void setGain(const arma::mat& _K){
        assert(_K.n_rows == ss_->A().n_rows && _K.n_cols == ss_->C().n_rows);
        K_ = _K;
        updateVariables();
    }

    void updateVariables();


private:
    _StateSpace* ss_;

    // Observer Gain
    arma::mat K_;
    arma::mat H_;

    // Stability matrix of observer
    arma::mat A_hat_;

    arma::mat u_;
    arma::mat prev_u_;

    //-- System Variable
    arma::mat state_;
    arma::mat prev_state_;
    arma::mat state_trans_;
    arma::mat output_;
    arma::mat prev_output_;

    //-- Estimation Variable
    // State
    arma::mat state_est_;
    arma::mat output_est_;
    // Prev state
    arma::mat prev_state_est_;
    // State transition
    arma::mat state_trans_est_;

    void calcState();
    void calcOutput();

    double dt_;

};

template <class _StateSpace>
Observer<_StateSpace>::Observer(_StateSpace* _ss, const arma::mat& _K, double _time_step)
    : ss_(_ss)
    , u_(_ss->B().n_cols, 1, arma::fill::zeros)
    , prev_u_(u_)
    , state_(_ss->A().n_rows, 1, arma::fill::zeros)
    , prev_state_(state_)
    , output_(_ss->C().n_rows, 1, arma::fill::zeros)
    , prev_output_(output_)
    , state_est_(state_)
    , prev_state_est_(state_)
    , dt_(_time_step){

    // there are no other member initialization due this one
    assert(_K.n_rows == _ss->A().n_rows && _K.n_cols == _ss->C().n_rows);

    // for Full-order observer
    K_ = _K;

    // updateVariables();

}

template <class _StateSpace>
Observer<_StateSpace>::~Observer(){

    ss_ = nullptr;
}

template <class _StateSpace>
void Observer<_StateSpace>::calcState(){

    static arma::mat term1, term2, term3, term4;

    term1 = state_trans_ * prev_state_;
    term2 = ss_->B() * u_;
    term3 = state_trans_ * ss_->B() * prev_u_;
    term4 = (term2 + term3) * (dt_ * .5);

    state_ = term1 + term4;

    calcOutput();    

}

template <class _StateSpace>
void Observer<_StateSpace>::calcOutput(){

    arma::mat term1 = ss_->C() * prev_state_;
    output_ = term1 + (ss_->D() * u_);
}

template <class _StateSpace>
arma::mat Observer<_StateSpace>::calcStateEst(){

    calcState();

    static arma::mat term1, term2, term3, term4, term5;
    term1 = state_trans_est_ * prev_state_est_;
    term2 = K_ * output_ + H_ * u_;
    term3 = K_ * prev_output_ + H_ * prev_u_;
    term4 = state_trans_ * term3;
    term5 = (term2 + term4) * (dt_ * .5);

    state_est_ = term1 + term5;

    calcOutputEst();

    prev_state_est_ = state_est_;
    prev_state_ = state_;
    prev_u_ = u_;
    prev_output_ = output_;

    return arma::join_cols(state_est_, output_est_);

}

template <class _StateSpace>
void Observer<_StateSpace>::calcOutputEst(){

    output_est_ = ss_->C() * prev_state_est_;
}

template <class _StateSpace>
void Observer<_StateSpace>::updateVariables(){

    // system
    state_trans_ = arma::expmat(ss_->A() * dt_);

    // for Full-order observer
//    K_.print("K : ");
    A_hat_ = ss_->A() - (K_ * ss_->C());
    H_ = ss_->B();
    state_trans_est_ = arma::expmat(A_hat_ * dt_);

}

}
