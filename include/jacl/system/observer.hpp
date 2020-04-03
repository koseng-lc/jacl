/**
*   @author : koseng (Lintang)
*   @brief : Simple Luenberger Full-Order Observer
*/

#pragma once

#include <jacl/system/base_system.hpp>

namespace jacl{

template <class _StateSpace>
class Observer:public BaseSystem<_StateSpace>{
public:
    Observer(_StateSpace* _ss, const arma::mat& _K, double _time_step = 1e-4)
        : BaseSystem<_StateSpace>(_ss, _time_step)
        , K_(_K)
        , meas_(this->out_)
        , prev_meas_(meas_){}
    ~Observer(){}
    
    void updateVar() override{
        A_hat_ = this->ss_->A() - (K_ * this->ss_->C());
        H_ = this->ss_->B();
        this->state_trans_ = arma::expmat(A_hat_ * this->dt_);
    }
    //-- overload
    auto convolve(const arma::vec& _in, const arma::vec& _meas) -> arma::vec{
        meas_ = _meas;
        arma::vec est( convolve(_in) );
        prev_meas_ = meas_;
        return est;
    }
    auto convolve(const arma::vec& _in) -> arma::vec override{
        setIn(_in);
        this->state_ = dstate();        
        this->out_ = output();
        this->prev_state_ = this->state_;
        return this->out_;
    }
    void setGain(const arma::mat& _K){
        assert(_K.n_rows == ss_->A().n_rows && _K.n_cols == ss_->C().n_rows);
        K_ = _K;
        updateVar();
    }

protected:
    auto setIn(const arma::vec& _in) -> void{
        this->in_ = _in;
    }
    auto dstate() -> arma::vec{
        static arma::mat term1, term2, term3, term4, term5;

        term1 = this->state_trans_ * this->prev_state_;
        term2 = K_ * meas_ + H_ * this->in_;
        term3 = K_ * prev_meas_ + H_ * this->prev_in_;
        term4 = this->state_trans_ * term3;
        term5 = (term2 + term4) * (this->dt_ * .5);

        this->prev_in_ = this->in_;        

        return term1 + term5;
    }
    auto output() -> arma::vec{
        return this->ss_->C() * this->prev_state_;
    }

private:
    arma::mat K_;
    arma::mat H_;
    arma::mat A_hat_;
    arma::vec meas_;
    arma::vec prev_meas_;

};

}
