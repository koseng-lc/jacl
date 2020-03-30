/**
*   @author : koseng (Lintang)
*/

#pragma once

#include "state_space.hpp"
#include "lft.hpp"

namespace jacl{

template <class _StateSpace, std::size_t z1_size, std::size_t w1_size, std::size_t y1_size, std::size_t u1_size>
class UpperLFT:public LFT<_StateSpace, z1_size, w1_size, 0, 0, y1_size, u1_size>{
public:
    UpperLFT(_StateSpace* _ss);
    ~UpperLFT();

    //-- setters and getters
    inline auto A() const -> arma::mat const&{ return this->ss_->A(); }

    inline auto B1() const -> arma::mat const & { return B1_; }
    inline auto B1() -> arma::mat& { return B1_; }
    inline auto B2() const -> arma::mat const & { return B2_; }
    inline auto B2() -> arma::mat& { return B2_; }

    inline auto C1() const -> arma::mat const & { return C1_; }
    inline auto C1() -> arma::mat& { return C1_; }
    inline auto C2() const -> arma::mat const & { return C2_; }
    inline auto C2() -> arma::mat& { return C2_; }

    inline auto D11() const -> arma::mat const & { return D11_; }
    inline auto D11() -> arma::mat& { return D11_; }
    inline auto D12() const -> arma::mat const & { return D12_; }
    inline auto D12() -> arma::mat& { return D12_; }
    inline auto D21() const -> arma::mat const & { return D21_; }
    inline auto D21() -> arma::mat& { return D21_; }
    inline auto D22() const -> arma::mat const & { return D22_; }
    inline auto D22() -> arma::mat& { return D22_; }

    inline auto delta() const -> arma::mat const & { return delta_; }
    auto delta() -> arma::mat& { return delta_; }

private:
    //-- Realization
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

};

template <class _StateSpace, std::size_t z1_size, std::size_t w1_size, std::size_t y1_size, std::size_t u1_size>
UpperLFT<_StateSpace, z1_size, w1_size, y1_size, u1_size>::UpperLFT(_StateSpace* _ss)
    : LFT<_StateSpace,z1_size, w1_size, 0, 0, y1_size, u1_size>(_ss){

    B1_ = this->ss_->B().head_cols(w1_size);
    B2_ = this->ss_->B().tail_cols(u1_size);

    C1_ = this->ss_->C().head_rows(z1_size);
    C2_ = this->ss_->C().tail_rows(y1_size);

    D11_ = this->ss_->D().submat(0, 0, z1_size - 1, w1_size - 1);
    D12_ = this->ss_->D().submat(0, w1_size,
                                 z1_size - 1, w1_size + u1_size - 1);
    D21_ = this->ss_->D().submat(z1_size, 0,
                                 z1_size + y1_size - 1, w1_size - 1);
    D22_ = this->ss_->D().submat(z1_size, w1_size,
                                 z1_size + y1_size - 1, w1_size + u1_size - 1);

}

template <class _StateSpace, std::size_t z1_size, std::size_t w1_size, std::size_t y1_size, std::size_t u1_size>
UpperLFT<_StateSpace, z1_size, w1_size, y1_size, u1_size>::~UpperLFT(){

}

}
