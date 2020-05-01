/**
*   @author : koseng (Lintang)
*   @brief : Simple Luenberger Full-Order Observer
*/

#pragma once

#include <jacl/system/base_system.hpp>

namespace jacl{ namespace system{

template <class _StateSpace>
class Observer:public BaseSystem<_StateSpace>{
public:
    Observer(_StateSpace* _ss, const arma::mat& _K, double _time_step = 1e-4)
        : BaseSystem<_StateSpace>(_ss, _time_step)
        , K_(_K)
        , meas_(this->out_)
        , prev_meas_(meas_){}
    ~Observer(){}
        
    //-- overload
    auto convolve(const arma::vec& _in, const arma::vec& _meas) -> arma::vec{
        meas_ = _meas;
        arma::vec est( convolve(_in) );
        prev_meas_ = meas_;
        return est;
    }    
    void setGain(const arma::mat& _K){
        assert(_K.n_rows == _StateSpace::n_states && _K.n_cols == _StateSpace::n_outputs);
        K_ = _K;
        updateVar();
    }

protected:
    auto convolve(const arma::vec& _in)
        -> typename BaseSystem<_StateSpace>::output_t override{
        setIn(_in);
        this->state_ = dstate();        
        this->out_ = output();
        this->prev_state_ = this->state_;
        return this->out_;
    }
    auto setIn(const arma::vec& _in) -> void{
        this->in_ = _in;
    }
    virtual auto dstate() -> arma::vec = 0;
    auto output() -> arma::vec override{
        return this->ss_->C() * this->prev_state_;
    }
    virtual auto updateVar() -> void = 0;

protected:
    arma::mat K_;
    arma::mat H_;
    arma::mat A_hat_;
    arma::vec meas_;
    arma::vec prev_meas_;

};

} } // namespace jacl::system
