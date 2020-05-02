/**
*   @author : koseng (Lintang)
*   @brief : Simple Continuous Luenberger Full-Order Observer
*/

#pragma once

#include <jacl/system/base_system.hpp>

namespace jacl{ namespace system{

template <class _StateSpace>
class ContinuousObserver:public Observer<_StateSpace>{
public:
    ContinuousObserver(_StateSpace* _ss, const arma::mat& _K, double _time_step = 1e-4)
        : Observer<_StateSpace>(_ss, _K, _time_step){}
    ~ContinuousObserver(){}
        
protected:
    auto dstate() -> typename Observer<_StateSpace>::state_t{
        static arma::mat term1, term2, term3, term4, term5;

        term1 = this->state_trans_ * this->prev_state_;
        term2 = this->K_ * this->meas_ + this->H_ * this->in_;
        term3 = this->K_ * this->prev_meas_ + this->H_ * this->prev_in_;
        term4 = this->state_trans_ * term3;
        term5 = (term2 + term4) * (this->dt_ * .5);

        this->prev_in_ = this->in_;

        return term1 + term5;
    }
    auto updateVar() -> void override{
        this->A_hat_ = this->ss_->A() - (this->K_ * this->ss_->C());
        this->H_ = this->ss_->B();
        this->state_trans_ = arma::expmat(this->A_hat_ * this->dt_);
    }

};

} } // namespace jacl::system
