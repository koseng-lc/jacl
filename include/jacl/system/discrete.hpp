/**
*   @author : koseng (Lintang)
*   @brief : Discrete system
*/

#pragma once

#include <jacl/system/base_system.hpp>

namespace jacl{ namespace system{

template <class _StateSpace>
class Discrete:public BaseSystem<_StateSpace>{
public:
    Discrete(_StateSpace* _ss, double _time_step)
        : BaseSystem<_StateSpace>(_ss, _time_step){}
    ~Discrete(){}
    auto samplingPeriod(){
        return this->dt_;
    }
protected:    
    void setIn(const typename BaseSystem<_StateSpace>::input_t& _in) override{
        this->in_ = _in;
    }
    auto dstate() -> typename BaseSystem<_StateSpace>::state_t override{
        static arma::mat term1, term2;
        term1 = this->ss_->A() * this->prev_state_;
        term2 = this->ss_->B() * this->in_;
        return term1 + term2;
    }
    auto output() -> typename BaseSystem<_StateSpace>::output_t override{
        static arma::mat term1;
        term1 = this->ss_->C() * this->prev_state_;
        return term1 + (this->ss_->D() * this->in_);
    }
    void updateVar() override{
        
    }
};

template <typename Scalar,
          std::size_t ns,
          std::size_t ni,
          std::size_t no,
          class PhysicalParam,
          class ...Rest>
class Discrete<::jacl::state_space::NonLinear<Scalar,ns,ni,no,PhysicalParam,Rest...> >
        :public BaseSystem<::jacl::state_space::NonLinear<Scalar,ns,ni,no,PhysicalParam,Rest...> >{
private:
    typedef ::jacl::state_space::NonLinear<Scalar,ns,ni,no,PhysicalParam,Rest...> _StateSpace;

public:
    Discrete(_StateSpace* _ss, double _time_step = 1e-4)
        : BaseSystem<_StateSpace>(_ss, _time_step){}
    ~Discrete(){}

 protected:
    void setIn(const typename BaseSystem<_StateSpace>::input_t& _in) override{
        ::jacl::state_space::detail::NonLinearStateSpaceClient<_StateSpace>::setSig(this->ss_, arma::join_cols(this->prev_state_, _in));
    }
    auto dstate() -> typename BaseSystem<_StateSpace>::state_t override{
        return ::jacl::state_space::detail::NonLinearStateSpaceClient<_StateSpace>::dstate(this->ss_);
    }
    auto output() -> typename BaseSystem<_StateSpace>::output_t override{
        return ::jacl::state_space::detail::NonLinearStateSpaceClient<_StateSpace>::output(this->ss_);
    }
    void updateVar() override{
        
    }
};

} } // namespace jacl::system