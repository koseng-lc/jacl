// TODO : Check the SSpace template argument are real StateSpace or not

#pragma once

#include "defs.h"
#include "state_space.h"

#include <armadillo>

namespace JACL{

template <class SSpace>
class Observer{
public:
    Observer(SSpace* _ss, Mat _K);
    ~Observer();

private:
//    using SS = SSpace;
    SSpace* ss_;

    // Observer Gain
    Mat K_;
    Mat H_;

    // Stability matrix of observer
    Mat A_hat_;
};

template <class SSpace>
Observer<SSpace>::Observer(SSpace* _ss, Mat _K)
    : ss_(_ss){

    assert(_K.n_rows == _ss->A().n_rows && _K.n_cols == _ss->C().n_rows);

    // for Full-order observer
    K_ = _K;
    A_hat_ = ss_->A() - K_ * ss_->C();
    H_ = ss_->B();
}

template <class SSpace>
Observer<SSpace>::~Observer(){

}

}
