/**
*   @author : koseng (Lintang)
*   @brief : jacl Linear Fractional Transformation base class
*/

#pragma once

#include <jacl/state_space/linear.hpp>

namespace jacl{ namespace lft{
template <typename _StateSpace,
          std::size_t z1_size,
          std::size_t w1_size,
          std::size_t z2_size,
          std::size_t w2_size>
class LFT{
public:

    using state_space_t = _StateSpace;

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

    template <typename _LFT>
    auto operator*(const _LFT& _lft){
        using ss_t = ::jacl::state_space::Linear<typename state_space_t::scalar_t,
                                    state_space_t::n_states + _LFT::state_space_t::n_states,
                                    w1_sz + _LFT::w2_sz,
                                    z1_sz + _LFT::z2_sz>;
        ss_t ss;

        using R_t = typename arma::Mat<typename state_space_t::scalar_t>::template fixed<z2_size, _LFT::w1_sz>;
        using Rh_t = typename arma::Mat<typename state_space_t::scalar_t>::template fixed<_LFT::z1_sz, w2_sz>;

        R_t R( arma::eye(w2_size, _LFT::w1_sz) - D22()*_lft.D11() );
        R_t R_inv( arma::inv(R));
        Rh_t Rh( arma::eye(_LFT::z1_sz, w2_sz) - _lft.D11()*D22() );
        Rh_t Rh_inv( arma::inv(Rh) );

        arma::mat temp1, temp2, temp3, temp4;

        temp1 = B2()*Rh_inv;
        temp2 = temp1*_lft.D11();
        temp3 = _lft.B1()*R_inv;
        temp4 = temp3*D22();
        typename arma::Mat<typename state_space_t::scalar_t
                 >::template fixed<state_space_t::n_states + _LFT::state_space_t::n_states,
                                   state_space_t::n_states + _LFT::state_space_t::n_states> A( arma::join_cols(
            arma::join_rows(A() + temp2*C2(), temp1*_lft.C1()),
            arma::join_rows(temp3*C2(), _lft.A() + temp4*_lft.C1()) )
        );
        typename arma::Mat<typename state_space_t::scalar_t
                 >::template fixed<state_space_t::n_states + _LFT::state_space_t::n_states,
                                   w1_sz + _LFT::w2_sz> B( arma::join_cols(
            arma::join_rows(B1() + temp2*D21(), temp1*_lft.D12()),
            arma::join_rows(temp3*D21(), _lft.B2() + temp4*_lft.D12()) )
        );
        temp1 = D12()*_lft.D11()*R_inv;
        temp2 = D12()*Rh_inv;
        temp3 = _lft.D21()*R_inv;
        temp4 = temp3*D22();
        typename arma::Mat<typename state_space_t::scalar_t
                 >::template fixed<z1_sz + _LFT::z2_sz,
                                   state_space_t::n_states + _LFT::state_space_t::n_states> C( arma::join_cols(
            arma::join_rows(C1() + temp1*C2(), temp2*_lft.C1()),
            arma::join_rows(temp3*C2(), _lft.C2() + temp4*_lft.C1()) )
        );
        typename arma::Mat<typename state_space_t::scalar_t
                 >::template fixed<z1_sz + _LFT::z2_sz,
                                   w1_sz + _LFT::w2_sz> D( arma::join_cols(
            arma::join_rows(D11() + temp1*D21(), temp2*_lft.D12()),
            arma::join_rows(temp3*D21(), _lft.D22() + temp4*_lft.D12()) )
        );
        LFT<ss_t, z1_sz, w1_sz, _LFT::z2_sz, _LFT::w2_sz> lft(A,B,C,D);
        
        return lft;
    }

public:
    static constexpr auto n_states{ _StateSpace::n_states };
    static constexpr auto z1_sz { z1_size };
    static constexpr auto w1_sz { w1_size };
    static constexpr auto z2_sz { z2_size };
    static constexpr auto w2_sz { w2_size };

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
template <typename _StateSpace,
          std::size_t z1_size,
          std::size_t w1_size,
          std::size_t z2_size,
          std::size_t w2_size>
LFT<_StateSpace, z1_size, w1_size, z2_size, w2_size>::LFT(_StateSpace* _ss,
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
