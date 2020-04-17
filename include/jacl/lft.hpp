#pragma once

#include <jacl/linear_state_space.hpp>

namespace jacl{ namespace lft{
template <class _StateSpace,
          std::size_t z1_size,
          std::size_t w1_size,
          std::size_t z2_size,
          std::size_t w2_size,
          std::size_t y_size,
          std::size_t u_size>
class LFT{
public:
    LFT(_StateSpace* _ss,
        const arma::mat& _A, const arma::mat& _B1, const arma::mat& _B2,
        const arma::mat& _C1, const arma::mat& _D11, const arma::mat& _D12,
        const arma::mat& _C2, const arma::mat& _D21, const arma::mat& _D22);
    ~LFT(){}

    //-- setters and getters
    inline const auto& A() const{ return A_; }

    inline const auto& B1() const{ return B1_; }
    inline auto& B1(){ return B1_; }
    inline const auto& B2() const{ return B2_; }
    inline auto& B2(){ return B2_; }

    inline const auto& C1() const{ return C1_; }
    inline auto& C1(){ return C1_; }
    inline const auto& C2() const{ return C2_; }
    inline auto& C2(){ return C2_; }

    inline const auto& D11() const{ return D11_; }
    inline auto& D11(){ return D11_; }
    inline const auto& D12() const{ return D12_; }
    inline auto& D12(){ return D12_; }
    inline const auto& D21() const{ return D21_; }
    inline auto& D21(){ return D21_; }
    inline const auto& D22() const{ return D22_; }
    inline auto& D22(){ return D22_; }

    inline const auto& delta() const{ return delta_; }
    inline auto& delta(){ return delta_; }

public:
    static constexpr auto n_states{ _StateSpace::n_states };
    static constexpr auto z1_sz { z1_size };
    static constexpr auto w1_sz { w1_size };
    static constexpr auto z2_sz { z2_size };
    static constexpr auto w2_sz { w2_size };
    static constexpr auto y_sz { y_size };
    static constexpr auto u_sz { u_size };

protected:
     //-- Realization
    const arma::mat& A_;
    arma::mat B1_;
    arma::mat B2_;
    arma::mat C1_;
    arma::mat C2_;
    arma::mat D11_;
    arma::mat D12_;
    arma::mat D21_;
    arma::mat D22_;

    //-- Delta Matrix
    arma::mat delta_;

    _StateSpace* ss_;
};
template <class _StateSpace,
          std::size_t z1_size,
          std::size_t w1_size,
          std::size_t z2_size,
          std::size_t w2_size,
          std::size_t y_size,
          std::size_t u_size>
LFT<_StateSpace, z1_size, w1_size, z2_size, w2_size, y_size, u_size>::LFT(_StateSpace* _ss,
            const arma::mat& _A, const arma::mat& _B1, const arma::mat& _B2,
            const arma::mat& _C1, const arma::mat& _D11, const arma::mat& _D12,
            const arma::mat& _C2, const arma::mat& _D21, const arma::mat& _D22)
    : ss_(_ss)
    , A_(_A)
    , B1_(_B1)
    , B2_(_B2)
    , C1_(_C1)
    , C2_(_C2)
    , D11_(_D11)
    , D12_(_D12)
    , D21_(_D21)
    , D22_(_D22){
}

} } // namespace jacl::lft
