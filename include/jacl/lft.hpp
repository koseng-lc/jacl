#pragma once

#include "linear_state_space.hpp"

namespace jacl{
template <class _StateSpace,
          std::size_t z1_size,
          std::size_t w1_size,
          std::size_t z2_size,
          std::size_t w2_size,
          std::size_t y_size,
          std::size_t u_size>
class LFT{
public:
    LFT(_StateSpace* _ss);
    ~LFT();

    static constexpr auto n_states{ _StateSpace::n_states };
    static constexpr auto z1_sz { z1_size };
    static constexpr auto w1_sz { w1_size };
    static constexpr auto z2_sz { z2_size };
    static constexpr auto w2_sz { w2_size };
    static constexpr auto y_sz { y_size };
    static constexpr auto u_sz { u_size };

protected:
    //-- Realization of LFT
    _StateSpace* ss_;

};

template <class _StateSpace,
          std::size_t z1_size,
          std::size_t w1_size,
          std::size_t z2_size,
          std::size_t w2_size,
          std::size_t y_size,
          std::size_t u_size>
LFT<_StateSpace, z1_size, w1_size, z2_size, w2_size, y_size, u_size>::LFT(_StateSpace* _ss)
    : ss_(_ss){

}

template <class _StateSpace,
          std::size_t z1_size,
          std::size_t w1_size,
          std::size_t z2_size,
          std::size_t w2_size,
          std::size_t y_size,
          std::size_t u_size>
LFT<_StateSpace, z1_size, w1_size, z2_size, w2_size, y_size, u_size>::~LFT(){

}

}
