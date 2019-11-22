#pragma once

#include "state_space.h"

namespace jacl{
template <class SSpace>
class LFT{
public:
    LFT(SSpace* _ss);
    ~LFT();

protected:
    //-- Realization of LFT
    //-- This one consist A, B2, C2, D22
    SSpace* ss_;

};

template <class SSpace>
LFT<SSpace>::LFT(SSpace* _ss)
    : ss_(_ss){

}

template <class SSpace>
LFT<SSpace>::~LFT(){

}

}
