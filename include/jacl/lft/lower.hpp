/**
*   @author : koseng (Lintang)
*   @brief : Lower Linear Fractional Transformation
*/

#pragma once

#include <jacl/state_space/linear.hpp>
#include <jacl/lft/lft.hpp>

namespace jacl{ namespace lft{

template <class _StateSpace, std::size_t z2_size, std::size_t w2_size, std::size_t y2_size, std::size_t u2_size>
class LowerLFT:public LFT<_StateSpace, 0, 0, z2_size, w2_size, y2_size, u2_size>{
public:
    LowerLFT(_StateSpace* _ss);
    LowerLFT(const arma::mat& _A, const arma::mat& _B1, const arma::mat& _B2,
             const arma::mat& _C1, const arma::mat& _D11, const arma::mat& _D12,
             const arma::mat& _C2, const arma::mat& _D21, const arma::mat& _D22);
    ~LowerLFT(){}
};

template <class _StateSpace, std::size_t z2_size, std::size_t w2_size, std::size_t y2_size, std::size_t u2_size>
LowerLFT<_StateSpace, z2_size, w2_size, y2_size, u2_size>::LowerLFT(_StateSpace* _ss)
    : LFT<_StateSpace, 0, 0, z2_size, w2_size, y2_size, u2_size>(_ss,
      _ss->A(), _ss->B().head_cols(w2_size), _ss->B().tail_cols(u2_size),
      _ss->C().head_rows(z2_size), _ss->D().submat(0, 0, z2_size - 1, w2_size - 1), _ss->D().submat(0, w2_size, z2_size - 1, w2_size + u2_size - 1),
      _ss->C().tail_rows(y2_size), _ss->D().submat(z2_size, 0, z2_size + y2_size - 1, w2_size - 1), _ss->D().submat(z2_size, w2_size, z2_size + y2_size - 1, w2_size + u2_size - 1)){
}

template <class _StateSpace, std::size_t z2_size, std::size_t w2_size, std::size_t y2_size, std::size_t u2_size>
LowerLFT<_StateSpace, z2_size, w2_size, y2_size, u2_size>::LowerLFT(
            const arma::mat& _A, const arma::mat& _B1, const arma::mat& _B2,
            const arma::mat& _C1, const arma::mat& _D11, const arma::mat& _D12,
            const arma::mat& _C2, const arma::mat& _D21, const arma::mat& _D22)
    : LFT<_StateSpace, 0, 0, z2_size, w2_size, y2_size, u2_size>(nullptr,
        _A, _B1, _B2,
        _C1, _D11, _D12,
        _C2, _D12, _D22){
}

} } // namespace jacl::lft
