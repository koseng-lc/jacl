/**
*   @author : koseng (Lintang)
*   @brief : Lower Linear Fractional Transformation
*/

#pragma once

#include "state_space.h"
#include "lft.h"

namespace jacl{

template <class _StateSpace, std::size_t z2_size, std::size_t w2_size, std::size_t y2_size, std::size_t u2_size>
class LowerLFT:public LFT<_StateSpace>{
public:
    LowerLFT(_StateSpace* _ss);
    ~LowerLFT();

    //-- setters and getters
    inline arma::mat const& A() const{ return this->ss_->A(); }

    inline arma::mat const& B1() const{ return B1_; }
    inline arma::mat& B1(){ return B1_; }
    inline arma::mat const& B2() const{ return B2_; }
    inline arma::mat& B2(){ return B2_; }

    inline arma::mat const& C1() const{ return C1_; }
    inline arma::mat& C1(){ return C1_; }
    inline arma::mat const& C2() const{ return C2_; }
    inline arma::mat& C2(){ return C2_; }

    inline arma::mat const& D11() const{ return D11_; }
    inline arma::mat& D11(){ return D11_; }
    inline arma::mat const& D12() const{ return D12_; }
    inline arma::mat& D12(){ return D12_; }
    inline arma::mat const& D21() const{ return D21_; }
    inline arma::mat& D21(){ return D21_; }
    inline arma::mat const& D22() const{ return D22_; }
    inline arma::mat& D22(){ return D22_; }

    inline arma::mat const& delta() const{ return delta_; }
    arma::mat& delta(){ return delta_; }

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

template <class _StateSpace, std::size_t z2_size, std::size_t w2_size, std::size_t y2_size, std::size_t u2_size>
LowerLFT<_StateSpace, z2_size, w2_size, y2_size, u2_size>::LowerLFT(_StateSpace* _ss)
    : LFT<_StateSpace>(_ss){

    B1_ = this->ss_->B().head_cols(w2_size);
    B2_ = this->ss_->B().tail_cols(u2_size);

    C1_ = this->ss_->C().head_rows(z2_size);
    C2_ = this->ss_->C().tail_rows(y2_size);

    D11_ = this->ss_->D().submat(0, 0, z2_size - 1, w2_size - 1);
    D12_ = this->ss_->D().submat(0, w2_size,
                                  z2_size - 1, w2_size + u2_size - 1);
    D21_ = this->ss_->D().submat(z2_size, 0,
                                  z2_size + y2_size - 1, w2_size - 1);
    D22_ = this->ss_->D().submat(z2_size, w2_size,
                                 z2_size + y2_size - 1, w2_size + u2_size - 1);
}

template <class _StateSpace, std::size_t z2_size, std::size_t w2_size, std::size_t y2_size, std::size_t u2_size>
LowerLFT<_StateSpace, z2_size, w2_size, y2_size, u2_size>::~LowerLFT(){

}

}
