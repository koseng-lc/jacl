#pragma once

#include <jacl/system/base_system.hpp>

namespace jacl{

template <class _StateSpace>
class DiscreteSystem:public BaseSystem<_StateSpace>{
public:
    DiscreteSystem(_StateSpace* _ss, double _time_step = 1e-4)
        : BaseSystem<_StateSpace>(_ss, _time_step){}
    ~DiscreteSystem(){}

protected:    
    auto setIn(const arma::vec& _in) -> void override{
        this->in_ = _in;
    }
    auto dstate() -> arma::vec override{
        static arma::mat term1, term2;
        term1 = this->ss_->A() * this->prev_state_;
        term2 = this->ss_->B() * this->in_;
        return term1 + term2;
    }
    auto output() -> arma::vec override{
        static arma::mat term1;
        term1 = this->ss_->C() * this->prev_state_;
        return term1 + (this->ss_->D() * this->in_);
    }
};

template <std::size_t ns,
          std::size_t ni,
          std::size_t no,
          class PhysicalParam,
          class ...Rest>
class DiscreteSystem<NonLinearStateSpace<ns,ni,no,PhysicalParam,Rest...> >
        :public BaseSystem<NonLinearStateSpace<ns,ni,no,PhysicalParam,Rest...> >{
private:
    typedef NonLinearStateSpace<ns,ni,no,PhysicalParam,Rest...> _StateSpace;

public:
    DiscreteSystem(_StateSpace* _ss, double _time_step = 1e-4)
        : BaseSystem<_StateSpace>(_ss, _time_step){}
    ~DiscreteSystem(){}

 protected:
    auto setIn(const arma::vec& _in)
        -> void override{
        detail::NonLinearStateSpaceClient<_StateSpace>::setSig(this->ss_, arma::join_cols(this->prev_state_, _in));
    }
    auto dstate()
        -> arma::vec override{
        return detail::NonLinearStateSpaceClient<_StateSpace>::dstate(this->ss_);
    }
    auto output()
        -> arma::vec override{
        return detail::NonLinearStateSpaceClient<_StateSpace>::output(this->ss_);
    }
};

}