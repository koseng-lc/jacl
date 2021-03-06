/**
*   @author : koseng (Lintang)
*   @brief : Upper Linear Fractional Transformation
*/

#pragma once

#include <jacl/state_space/linear.hpp>
#include <jacl/lft/lft.hpp>

namespace jacl::lft{

template <typename _StateSpace, std::size_t z_size, std::size_t w_size, std::size_t y_size, std::size_t u_size>
class UpperLFT:public LFT<_StateSpace, y_size, u_size, z_size, w_size>{
public:
    UpperLFT(_StateSpace* _ss);
    UpperLFT(const arma::mat& _A, const arma::mat& _B1, const arma::mat& _B2,
             const arma::mat& _C1, const arma::mat& _D11, const arma::mat& _D12,
             const arma::mat& _C2, const arma::mat& _D21, const arma::mat& _D22);
    ~UpperLFT(){}
};

template <typename _StateSpace, std::size_t z_size, std::size_t w_size, std::size_t y_size, std::size_t u_size>
UpperLFT<_StateSpace, z_size, w_size, y_size, u_size>::UpperLFT(_StateSpace* _ss)
    : LFT<_StateSpace, y_size, u_size, z_size, w_size>(_ss,
      _ss->A(), _ss->B().head_cols(w_size), _ss->B().tail_cols(u_size),
      _ss->C().head_rows(z_size), _ss->D().submat(0, 0, z_size - 1, w_size - 1), _ss->D().submat(0, w_size, z_size - 1, w_size + u_size - 1),
      _ss->C().tail_rows(y_size), _ss->D().submat(z_size, 0, z_size + y_size - 1, w_size - 1), _ss->D().submat(z_size, w_size, z_size + y_size - 1, w_size + u_size - 1)){      
}

template <typename _StateSpace, std::size_t z_size, std::size_t w_size, std::size_t y_size, std::size_t u_size>
UpperLFT<_StateSpace, z_size, w_size, y_size, u_size>::UpperLFT(
            const arma::mat& _A, const arma::mat& _B1, const arma::mat& _B2,
            const arma::mat& _C1, const arma::mat& _D11, const arma::mat& _D12,
            const arma::mat& _C2, const arma::mat& _D21, const arma::mat& _D22)
    : LFT<_StateSpace, y_size, u_size, z_size, w_size>(nullptr,
        _A, _B1, _B2,
        _C1, _D11, _D12,
        _C2, _D12, _D22){
}

} // namespace jacl::lft
