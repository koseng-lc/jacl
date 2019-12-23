#pragma once

#include "state_space.h"

namespace jacl{
template <class _StateSpace>
class LFT{
public:
    LFT(_StateSpace* _ss);
    ~LFT();

protected:
    //-- Realization of LFT
    _StateSpace* ss_;

};

template <class _StateSpace>
LFT<_StateSpace>::LFT(_StateSpace* _ss)
    : ss_(_ss){

}

template <class _StateSpace>
LFT<_StateSpace>::~LFT(){

}

}
