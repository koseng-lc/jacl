#pragma once

#include <jacl/linear_state_space.hpp>
#include <jacl/nonlinear_state_space.hpp>

namespace jacl{

template<class _StateSpace>
class NonLinearStateSpaceClient{
    public:
    static auto setSig(_StateSpace* _ss, const arma::vec& _sig) -> arma::vec{
        _ss->sig_ = arma::conv_to<decltype(_ss->sig_)>(_sig);
    }
    static auto dstate(_StateSpace* _ss) -> arma::vec {
        arma::vec xdot(_StateSpace::n_states, 1, arma::fill::zeros);
        for(auto& x:xdot)
            x = _ss->state_fn(*_ss);
        return xdot;
    }
    static auto output(_StateSpace* _ss) -> arma::vec {
        arma::vec out(_StateSpace::n_outputs, 1, arma::fill::zeros);
        for(auto& o:out)
            out = _ss->state_fn(*_ss);
        return out;
    }
    friend class BaseSystem<_StateSpace>;
};

template <class _StateSpace>
class BaseSystem{
    public:
    BaseSystem(_StateSpace* _ss):ss_(_ss){};
    ~BaseSystem(){}
    protected:
    auto setSig(const arma::vec& _sig) -> void{
        NonLinearStateSpaceClient<_StateSpace>::setSig(ss_, _sig);
    }    
    protected:
    _StateSpace* ss_;
};

}

