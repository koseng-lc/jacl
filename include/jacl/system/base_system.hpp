#pragma once

#include <jacl/pattern/observer.hpp>
#include <jacl/linear_state_space.hpp>
#include <jacl/nonlinear_state_space.hpp>
#include <jacl/medium.hpp>

namespace jacl{ namespace system{

template <class _StateSpace>
class BaseSystem:public pattern::Observer, public ::jacl::state_space::detail::NonLinearStateSpaceClient<_StateSpace>{
public:
    static constexpr std::size_t n_states{_StateSpace::n_states};
    static constexpr std::size_t n_inputs{_StateSpace::n_inputs};
    static constexpr std::size_t n_outputs{_StateSpace::n_outputs};
    using StateSpace = _StateSpace;
    using base = BaseSystem;

    BaseSystem(_StateSpace* _ss, double _time_step = 1e-4)
        : pattern::Observer(_ss)
        , state_(_StateSpace::n_states,1, arma::fill::zeros)
        , prev_state_(state_)
        , state_trans_(_StateSpace::n_states, _StateSpace::n_states, arma::fill::zeros)
        , in_(_StateSpace::n_inputs, 1, arma::fill::zeros)
        , prev_in_(in_)
        , out_(_StateSpace::n_outputs, 1, arma::fill::zeros)
        , ss_(_ss)
        , dt_(_time_step){

    }
    ~BaseSystem(){}    

    inline auto recapitulate() -> arma::vec{
        return arma::join_cols(prev_state_, arma::join_cols(in_, out_));
    }
    virtual auto convolve(const arma::vec& _in) -> arma::vec{
        setIn(_in);        
        this->state_ = dstate();
        this->prev_state_ = this->state_;
        this->out_ = output();
        return this->out_;
    }
    auto A() -> arma::mat { return this->ss_->A(); }
    auto B() -> arma::mat { return this->ss_->B(); }
    auto C() -> arma::mat { return this->ss_->C(); }
    auto D() -> arma::mat { return this->ss_->D(); }
    
    auto reset(){
        state_ = arma::zeros(_StateSpace::n_states, 1);
        prev_state_ = state_;
        in_ = arma::zeros(_StateSpace::n_inputs, 1);
        prev_in_ = in_;
        out_ = state_ = arma::zeros(_StateSpace::n_outputs, 1);
    }

protected:
    virtual void setIn(const arma::vec& _in) = 0;
    virtual auto dstate() -> arma::vec = 0;
    virtual auto output() -> arma::vec = 0;
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

protected:
    arma::vec state_;
    arma::vec prev_state_;
    arma::mat state_trans_;
    arma::vec in_;
    arma::vec prev_in_;
    arma::vec out_;
    double dt_;
    _StateSpace* ss_;
    
private:
    friend class ::jacl::system::detail::BaseSystemClient<BaseSystem>;
};

} } // namespace jacl::system

