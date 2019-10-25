// TODO : Check the SSpace template argument are real StateSpace or not

#pragma once

#include "defs.h"
#include "state_space.h"

namespace JACL{

template <class SSpace>
class Observer{
public:

    Observer(SSpace* _ss, const Mat& _K, double _time_step = 1e-4);
    ~Observer();

    Mat calcStateEst();
    void calcOutputEst();

    void setInput(const Mat& _u){
        u_ = _u;
    }

    void setTimeStep(double _time_step){
        dt_ = _time_step;
    }

    void updateVariables();

private:

    SSpace* ss_;

    // Observer Gain
    Mat K_;
    Mat H_;

    // Stability matrix of observer
    Mat A_hat_;

    Mat u_;
    Mat prev_u_;

    //-- System Variable
    Mat state_;
    Mat prev_state_;
    Mat state_trans_;
    Mat st_I_;
    Mat output_;

    //-- Estimation Variable
    // State
    Mat state_est_;
    Mat output_est_;
    // Prev state
    Mat prev_state_est_;
    // State transition
    Mat state_trans_est_;
    // State transition + Identity
    Mat st_est_I_;

    void calcState();
    void calcOutput();

    double dt_;
};

template <class SSpace>
Observer<SSpace>::Observer(SSpace* _ss, const Mat& _K, double _time_step)
    : ss_(_ss)
    , u_(_ss->B().n_cols, 1, arma::fill::zeros)
    , prev_u_(u_)
    , state_(_ss->A().n_rows, 1, arma::fill::zeros)
    , prev_state_(state_)
    , state_est_(state_)
    , prev_state_est_(state_)
    , dt_(_time_step){

    // there are no other member initialization due this one
    assert(_K.n_rows == _ss->A().n_rows && _K.n_cols == _ss->C().n_rows);

    // for Full-order observer
    K_ = _K;

    updateVariables();

}

template <class SSpace>
Observer<SSpace>::~Observer(){

    ss_ = nullptr;
}

template <class SSpace>
void Observer<SSpace>::calcState(){

    Mat term1, term2, term3;
    term1 = state_trans_ * prev_state_;
    term2 = ss_->B() * (u_ + prev_u_);
    term3 = st_I_ * (dt_ * .5);

    state_ = term1 + (term3 * term2);

    calcOutput();

    prev_state_ = state_;
    prev_u_ = u_;

}

template <class SSpace>
void Observer<SSpace>::calcOutput(){
    Mat term1 = ss_->C() * prev_state_;
    output_ = term1 + (ss_->D() * u_);
}

template <class SSpace>
Mat Observer<SSpace>::calcStateEst(){

    calcState();

    Mat term1, term2, term3, term4, term5, term6;
    term1 = state_trans_ * prev_state_;
    term2 = K_ * output_;
    term3 = H_ * u_;
    term4 = term2 + term3;
    term5 = st_est_I_ * (dt_ * .5);
    term6 = term1 + (term5 * term4);

    state_est_ = term1 + term6;

    calcOutputEst();

    prev_state_est_ = state_est_;

    return arma::join_cols(state_est_, output_est_);

}

template <class SSpace>
void Observer<SSpace>::calcOutputEst(){

    output_est_ = ss_->C() * prev_state_est_;
}

template <class SSpace>
void Observer<SSpace>::updateVariables(){

    // system
    state_trans_ = arma::expmat(ss_->A() * dt_);
    st_I_ = state_trans_ + Mat(arma::size(ss_->A()), arma::fill::eye);

    // for Full-order observer
    A_hat_ = ss_->A() - (K_ * ss_->C());
    H_ = ss_->B();
    state_trans_est_ = arma::expmat(A_hat_ * dt_);
    st_est_I_ = state_trans_est_ + Mat(arma::size(ss_->A()), arma::fill::eye);

}

}
