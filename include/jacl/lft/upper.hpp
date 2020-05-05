/**
*   @author : koseng (Lintang)
*   @brief : Upper Linear Fractional Transformation
*/

#pragma once

#include <jacl/state_space/linear.hpp>
#include <jacl/lft/lft.hpp>

namespace jacl{ namespace lft{

template <class _StateSpace, std::size_t z1_size, std::size_t w1_size, std::size_t y1_size, std::size_t u1_size>
class UpperLFT:public LFT<_StateSpace, z1_size, w1_size, 0, 0, y1_size, u1_size>{
public:
    UpperLFT(_StateSpace* _ss);
    UpperLFT(const arma::mat& _A, const arma::mat& _B1, const arma::mat& _B2,
             const arma::mat& _C1, const arma::mat& _D11, const arma::mat& _D12,
             const arma::mat& _C2, const arma::mat& _D21, const arma::mat& _D22);
    ~UpperLFT(){}
};

template <class _StateSpace, std::size_t z1_size, std::size_t w1_size, std::size_t y1_size, std::size_t u1_size>
UpperLFT<_StateSpace, z1_size, w1_size, y1_size, u1_size>::UpperLFT(_StateSpace* _ss)
    : LFT<_StateSpace, z1_size, w1_size, 0, 0, y1_size, u1_size>(_ss,
      _ss->A(), _ss->B().head_cols(w1_size), _ss->B().tail_cols(u1_size),
      _ss->C().head_rows(z1_size), _ss->D().submat(0, 0, z1_size - 1, w1_size - 1), _ss->D().submat(0, w1_size, z1_size - 1, w1_size + u1_size - 1),
      _ss->C().tail_rows(y1_size), _ss->D().submat(z1_size, 0, z1_size + y1_size - 1, w1_size - 1), _ss->D().submat(z1_size, w1_size, z1_size + y1_size - 1, w1_size + u1_size - 1)){      
}

template <class _StateSpace, std::size_t z1_size, std::size_t w1_size, std::size_t y1_size, std::size_t u1_size>
UpperLFT<_StateSpace, z1_size, w1_size, y1_size, u1_size>::UpperLFT(
            const arma::mat& _A, const arma::mat& _B1, const arma::mat& _B2,
            const arma::mat& _C1, const arma::mat& _D11, const arma::mat& _D12,
            const arma::mat& _C2, const arma::mat& _D21, const arma::mat& _D22)
    : LFT<_StateSpace, z1_size, w1_size, 0, 0, y1_size, u1_size>(nullptr,
        _A, _B1, _B2,
        _C1, _D11, _D12,
        _C2, _D12, _D22,){
}

} } // namespace jacl::lft
