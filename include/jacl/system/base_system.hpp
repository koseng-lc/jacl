/**
*   @author : koseng (Lintang)
*   @brief : Interface class for System
*/

#pragma once

#include <jacl/pattern/observer.hpp>
#include <jacl/state_space/linear.hpp>
#include <jacl/state_space/nonlinear.hpp>
#include <jacl/medium.hpp>

namespace jacl{ namespace system{

template <class _StateSpace>
class BaseSystem:public pattern::Observer
                ,public ::jacl::state_space::detail::NonLinearStateSpaceClient<_StateSpace>{
public:
    static constexpr std::size_t n_states{_StateSpace::n_states};
    static constexpr std::size_t n_inputs{_StateSpace::n_inputs};
    static constexpr std::size_t n_outputs{_StateSpace::n_outputs};

    using state_space_t = _StateSpace;
    using base_t = BaseSystem;
    using scalar_t = typename _StateSpace::scalar_t;
    using state_t = typename arma::Col<scalar_t>::template fixed<n_states>;
    using input_t = typename arma::Col<scalar_t>::template fixed<n_inputs>;
    using output_t = typename arma::Col<scalar_t>::template fixed<n_outputs>;

    BaseSystem(_StateSpace* _ss, double _time_step = 1e-4)
        : pattern::Observer(_ss)
        , state_(arma::fill::zeros)
        , prev_state_(state_)
        , state_trans_(arma::fill::zeros)
        , in_(arma::fill::zeros)
        , prev_in_(in_)
        , out_(arma::fill::zeros)
        , ss_(_ss)
        , dt_(_time_step){

    }
    ~BaseSystem(){}    

    inline auto recapitulate()
        -> typename arma::Col<scalar_t>::template fixed<n_states+n_inputs+n_outputs>{
        return arma::join_cols(prev_state_, arma::join_cols(in_, out_));
    }
    virtual auto convolve(const input_t& _in)
        -> output_t{
        setIn(_in);        
        this->state_ = dstate();        
        this->out_ = output();
        this->prev_state_ = this->state_;
        return this->out_;
    }  
    
    auto reset(){
        state_.fill(.0);
        prev_state_ = state_;
        in_.fill(.0);
        prev_in_ = in_;
        out_.fill(.0);
    }
    inline auto dt() const{ return dt_; }
    
protected:
    virtual void setIn(const input_t& _in) = 0;
    virtual auto dstate() -> state_t = 0;
    virtual auto output() -> output_t = 0;
    virtual void updateVar(){}
    void update() override{
        updateVar();
    }
    auto setSubject(_StateSpace* _ss){
        if(_ss){
            ss_ = _ss;
            this->s_ = ss_;
            this->s_->attach(this);
        }
    }
    auto A() -> typename state_space_t::state_matrix_t { return this->ss_->A(); }
    auto B() -> typename state_space_t::input_matrix_t { return this->ss_->B(); }
    auto C() -> typename state_space_t::output_matrix_t { return this->ss_->C(); }
    auto D() -> typename state_space_t::feedforward_matrix_t { return this->ss_->D(); }

protected:
    state_t state_;
    state_t prev_state_;
    typename state_space_t::state_matrix_t state_trans_;
    input_t in_;
    input_t prev_in_;
    output_t out_;
    double dt_;
    _StateSpace* ss_;
    
private:
    friend class ::jacl::system::detail::BaseSystemClient<BaseSystem>;
};

} } // namespace jacl::system

