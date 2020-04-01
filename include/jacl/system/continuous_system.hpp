#pragma once

#include <jacl/system/base_system.hpp>

namespace jacl{
template <class _StateSpace>
class ContinuousSystem:public BaseSystem<_StateSpace>{
public:
    ContinuousSystem(){};
    ~ContinuousSystem(){}

    auto dstate(const arma::vec& _state, const arma::vec& _in) -> decltype(_state){
        arma::vec xdot;
    }
};
}