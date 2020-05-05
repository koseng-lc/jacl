/**
*   @author : koseng (Lintang)
*   @brief : Continuous system
*/

#pragma once

#include <jacl/system/base_system.hpp>

namespace jacl{ namespace system{

template <class _StateSpace>
class Continuous:public BaseSystem<_StateSpace>{
public:
    Continuous(_StateSpace* _ss, double _time_step = 1e-4)
        : BaseSystem<_StateSpace>(_ss, _time_step){}
    ~Continuous(){}

protected:
    void setIn(const typename BaseSystem<_StateSpace>::input_t& _in) override{
        this->in_ = _in;
    }
    auto dstate() -> typename BaseSystem<_StateSpace>::state_t override{
        static arma::mat term1, term2, term3, term4;        
        term1 = this->state_trans_ * this->prev_state_;
        term2 = this->ss_->B() * this->in_;
        term3 = this->state_trans_ * this->ss_->B() * this->prev_in_;
        term4 = (term2 + term3) * (this->dt_ * .5);        
        this->prev_in_ = this->in_;
        return term1 + term4;
    }    
    auto output() -> typename BaseSystem<_StateSpace>::output_t override{
        static arma::mat term1;
        term1 = this->ss_->C() * this->prev_state_;
        return term1 + (this->ss_->D() * this->in_);
    }
    void updateVar() override{
        this->state_trans_ = arma::expmat(this->ss_->A() * this->dt_);
    }
};

template <typename Scalar,
          std::size_t ns,
          std::size_t ni,
          std::size_t no,
          class PhysicalParam,
          class ...Rest>
class Continuous<::jacl::state_space::NonLinear<Scalar,ns,ni,no,PhysicalParam,Rest...> >
        :public BaseSystem<::jacl::state_space::NonLinear<Scalar,ns,ni,no,PhysicalParam,Rest...> >{
private:
    typedef ::jacl::state_space::NonLinear<Scalar,ns,ni,no,PhysicalParam,Rest...> _StateSpace;

public:
    Continuous(_StateSpace* _ss, double _time_step = 1e-4)
        : BaseSystem<_StateSpace>(_ss, _time_step){}
    ~Continuous(){}    

protected:
    void setIn(const typename BaseSystem<_StateSpace>::input_t& _in) override{
        ::jacl::state_space::detail::NonLinearStateSpaceClient<_StateSpace>::setSig(this->ss_, arma::join_cols(this->prev_state_, _in));
    }
    auto dstate() -> typename BaseSystem<_StateSpace>::state_t override{
        return ::jacl::state_space::detail::NonLinearStateSpaceClient<_StateSpace>::dstate(this->ss_) * this->dt_;
    }
    auto output() -> typename BaseSystem<_StateSpace>::output_t override{
        return ::jacl::state_space::detail::NonLinearStateSpaceClient<_StateSpace>::output(this->ss_);
    }
    void updateVar() override{
        
    }
};

} } // namespace jacl::system