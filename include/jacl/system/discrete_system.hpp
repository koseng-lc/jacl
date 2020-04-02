#pragma once

#include <jacl/system/base_system.hpp>

namespace jacl{

template <class _StateSpace>
class DiscreteSystem:public BaseSystem<_StateSpace>{
public:
    DiscreteSystem(_StateSpace* _ss, double _time_step = 1e-4)
        : BaseSystem<_StateSpace>(_ss, _time_step){}
    ~DiscreteSystem(){}
    inline auto setIn(const arma::vec& _in) -> void override{
        this->in_ = _in;
    }
    inline auto dstate() -> arma::vec override{
        static arma::mat term1, term2;
        term1 = this->ss_->A() * this->prev_state_;
        term2 = this->ss_->B() * this->in_;
        this->state_ = term1 + term2;
        this->prev_state_ = this->state_;
        return this->state_;
    }
    inline auto output() -> arma::vec override{
        static arma::mat term1;
        term1 = this->ss_->C() * this->prev_state_;
        return term1 + (this->ss_->D() * this->in_);
    }
    inline auto evaluate() -> arma::vec override{
        return arma::join_cols(dstate(), output());
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
 
    inline auto setIn(const arma::vec& _in)
        -> void override{
        detail::NonLinearStateSpaceClient<_StateSpace>::setSig(this->ss_, arma::join_cols(this->prev_state_, _in));
    }
    virtual inline auto dstate()
        -> arma::vec override{
        this->state_ = detail::NonLinearStateSpaceClient<_StateSpace>::dstate(this->ss_);
        this->prev_state_ = this->state_;
        return this->state_;
    }
    inline auto output()
        -> arma::vec override{
        return detail::NonLinearStateSpaceClient<_StateSpace>::output(this->ss_);
    }
    inline auto evaluate() -> arma::vec override{
        return arma::join_cols(dstate(), output());
    }
};

}