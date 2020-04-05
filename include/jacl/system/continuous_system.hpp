#pragma once

#include <jacl/system/base_system.hpp>

namespace jacl{ namespace system {

template <class _StateSpace>
class ContinuousSystem:public BaseSystem<_StateSpace>{
public:
    ContinuousSystem(_StateSpace* _ss, double _time_step = 1e-4)
        : BaseSystem<_StateSpace>(_ss, _time_step){}
    ~ContinuousSystem(){}

protected:
    auto setIn(const arma::vec& _in)-> void override{
        this->in_ = _in;
    }
    auto dstate() -> arma::vec override{
        static arma::mat term1, term2, term3, term4;        
        term1 = this->state_trans_ * this->prev_state_;
        term2 = this->ss_->B() * this->in_;
        term3 = this->state_trans_ * this->ss_->B() * this->prev_in_;
        term4 = (term2 + term3) * (this->dt_ * .5);        
        this->prev_in_ = this->in_;
        return term1 + term4;
    }    
    auto output() -> arma::vec override{
        static arma::mat term1;
        term1 = this->ss_->C() * this->prev_state_;
        return term1 + (this->ss_->D() * this->in_);
    }
    auto updateVar() -> void override{
        this->state_trans_ = arma::expmat(this->ss_->A() * this->dt_);
    }
};

template <std::size_t ns,
          std::size_t ni,
          std::size_t no,
          class PhysicalParam,
          class ...Rest>
class ContinuousSystem<NonLinearStateSpace<ns,ni,no,PhysicalParam,Rest...> >
        :public BaseSystem<NonLinearStateSpace<ns,ni,no,PhysicalParam,Rest...> >{
private:
    typedef NonLinearStateSpace<ns,ni,no,PhysicalParam,Rest...> _StateSpace;

public:
    ContinuousSystem(_StateSpace* _ss, double _time_step = 1e-4)
        : BaseSystem<_StateSpace>(_ss, _time_step){}
    ~ContinuousSystem(){}    

protected:
    auto setIn(const arma::vec& _in)
        -> void override{
        detail::NonLinearStateSpaceClient<_StateSpace>::setSig(this->ss_, arma::join_cols(this->prev_state_, _in));
    }
    auto dstate()
        -> arma::vec override{
        return detail::NonLinearStateSpaceClient<_StateSpace>::dstate(this->ss_) * this->dt_;
    }
    auto output()
        -> arma::vec override{
        return detail::NonLinearStateSpaceClient<_StateSpace>::output(this->ss_);
    }
    auto updateVar() -> void override{
        
    }
};

} } // namespace jacl::system