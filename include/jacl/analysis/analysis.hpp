/**
*   @author : koseng (Lintang)
*   @brief : jacl analysis stuff
*/

#pragma once

#include <jacl/lti_common.hpp>
#include <jacl/traits.hpp>

namespace jacl::analysis{

template <typename _Plant, typename _Controller,
          typename = std::enable_if_t<::jacl::traits::is_system_v<_Plant>>,
          typename = std::enable_if_t<::jacl::traits::is_system_v<_Controller>>>
auto nominalStability(_Plant _p, _Controller _k){
    
    constexpr auto continuous = ::jacl::traits::is_continuous_system_v<_Plant>;

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

    auto a = lti_common::isStable<continuous>(A_bar);
    auto b = lti_common::stabilizable<continuous>(temp1, temp2);
    auto c = lti_common::detectability<continuous>(temp1, temp3);

    return a & b & c;
}

template <typename _System,
          typename = std::enable_if_t<::jacl::traits::is_system_v<_System>>>
auto nominalPerformance(_System _sys, double _obj){
    
    return  lti_common::approxInfNorm(_sys) < _obj;
}

template <typename _StateSpace>
auto ssv(const _StateSpace& _ss){
    auto ubound(.0);
    auto lbound(.0);
    return std::make_pair(ubound, lbound);
}

} // namespace jacl::analysis