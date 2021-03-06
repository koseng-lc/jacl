/**
*   @author : koseng (Lintang)
*   @brief : Simple Discrete Luenberger Full-Order Observer
*/

#pragma once

#include <jacl/system/base_system.hpp>

namespace jacl::system{

template <class _StateSpace>
class DiscreteObserver:public Observer<_StateSpace>{
public:
    DiscreteObserver(_StateSpace* _ss, const arma::mat& _K, double _time_step = 1e-4)
        : Observer<_StateSpace>(_ss, _K, _time_step){}
    ~DiscreteObserver(){} 

protected:
    auto dstate() -> typename Observer<_StateSpace>::state_t{
        static arma::mat term1, term2;
        term1 = this->A_hat_ * this->prev_state_;
        term2 = this->K_ * this->meas_ + this->H_ * this->in_;
        return  term1 + term2;
    }
    void updateVar() override{
        this->A_hat_ = this->ss_->A() - (this->K_ * this->ss_->C());
        this->H_ = this->ss_->B();
    }

};

} // namespace jacl::system
