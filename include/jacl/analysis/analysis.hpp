#pragma once

#include <jacl/lti_common.hpp>

namespace jacl{ namespace analysis{

template <typename _Plant, typename _Controller>
auto nominalStability(const _Plant& _p, const _Controller& _k){
    typename arma::Mat<typename _Plant::scalar_t>::template
        fixed<_Plant::n_states + _Controller::n_states,
         _Controller::n_inputs> temp2 = arma::join_cols(
        _p.B()*_k.D(),
        _k.B()
    );

    typename arma::Mat<typename _Plant::scalar_t>::template
        fixed<_Plant::n_outputs,
         _Plant::n_states + _Controller::n_states> temp3 = arma::join_rows(
        _p.C(),
        _p.D()*_k.C()
    );

    typename arma::Mat<typename _Plant::scalar_t>::template
        fixed<_Controller::n_inputs,
         _Plant::n_outputs> temp4 = 
        arma::eye(_Controller::n_inputs, _Plant::n_outputs)
        - _p.D()*_k.D();
    temp4 = arma::inv(temp4);

    typename arma::Mat<typename _Plant::scalar_t>::template
        fixed<_Plant::n_states + _Controller::n_states,
         _Plant::n_states + _Controller::n_states> temp1 = arma::join_cols(
        arma::join_rows(_p.A(), _p.B()*_k.C()),
        arma::join_rows(arma::zeros(_Controller::n_states, _Plant::n_states), _k.A())
    );

    typename arma::Mat<typename _Plant::scalar_t>::template
        fixed<_Plant::n_states + _Controller::n_states,
         _Plant::n_states + _Controller::n_states> A_bar = temp1 + temp2*temp4*temp3;
         
    return lti_common::isStable(A_bar,false)
            & lti_common::stabilizable(temp1, temp2)
            & lti_common::detectability(temp1, temp3);
}

template <typename _StateSpace>
auto nominalPerformance(const _StateSpace& _ss, double _obj){
    return lti_common::approxInfNorm(_ss, _obj*10., _obj*0.1) < _obj;
}

template <typename _StateSpace>
auto ssv(const _StateSpace& _ss){
    auto ubound(.0);
    auto lbound(.0);
    return std::make_pair(ubound, lbound);
}

} } // namespace::analysis