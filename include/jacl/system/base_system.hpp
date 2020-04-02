#pragma once

#include <jacl/linear_state_space.hpp>
#include <jacl/nonlinear_state_space.hpp>

namespace jacl{
template <class _StateSpace>
class BaseSystem;
namespace detail {
template<class _StateSpace>
    class NonLinearStateSpaceClient{
    public:
        inline static auto setSig(_StateSpace* _ss, const arma::vec& _sig) -> void{
            _ss->sig_ = arma::conv_to<decltype(_ss->sig_)>::from(_sig);
        }
        inline static auto dstate(_StateSpace* _ss) -> arma::vec {
            arma::vec xdot(_StateSpace::n_states, 1, arma::fill::zeros);
            int i(0);
            for(const auto& fn:_ss->state_fn_){
                xdot(i) = fn(*_ss);
                ++i;
            }
            return xdot;
        }
        inline static auto output(_StateSpace* _ss) -> arma::vec {
            arma::vec out(_StateSpace::n_outputs, 1, arma::fill::zeros);
            int i(0);
            for(const auto& fn:_ss->output_fn_){
                out(i) = fn(*_ss);
                ++i;
            }
            return out;
        }
        friend class BaseSystem<_StateSpace>;
    };
} // namespace detail

template <class _StateSpace>
class BaseSystem{
public:
    BaseSystem(_StateSpace* _ss, double _time_step = 1e-4)
        : ss_(_ss)
        , dt_(_time_step){};
    ~BaseSystem(){}

protected:
    virtual inline auto setIn(const arma::vec& _in) -> void = 0;
    virtual inline auto dstate() -> arma::vec = 0;
    virtual inline auto output() -> arma::vec = 0;
    virtual inline auto evaluate() -> arma::vec = 0;

protected:
    arma::vec state_;
    arma::vec prev_state_;
    arma::mat state_trans_;
    arma::vec in_;
    arma::vec prev_in_;
    double dt_;
    _StateSpace* ss_;
};

} // namespace jacl

