/**
*   @author : koseng (Lintang)
*/

#pragma once

#include "state_space.h"
#include "lft.h"

namespace jacl{

template <class _StateSpace, int num_uncertainty>
class UpperLFT:public LFT<_StateSpace>{
public:
    UpperLFT(_StateSpace* _ss);
    ~UpperLFT();

    //-- setters and getters
    inline arma::mat const& A() const{ return this->ss_->A(); }

    inline arma::mat const& B1() const{ return B1_; }
    inline arma::mat const& B2() const{ return B1_; }

    inline arma::mat const& C1() const{ return C1_; }
    inline arma::mat const& C2() const{ return C2_; }

    inline arma::mat const& D11() const{ return D11_; }

    inline arma::mat const& D12() const{ return D12_; }

    inline arma::mat const& D21() const{ return D21_; }

    inline arma::mat const& D22() const{ return D22_; }

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

template <class _StateSpace, int num_uncertainty>
UpperLFT<_StateSpace, num_uncertainty>::UpperLFT(_StateSpace* _ss)
    : LFT<_StateSpace>(_ss){

    B1_ = this->ss_->B().head_cols(num_uncertainty);
    B2_ = this->ss_->B().tail_cols(this->ss_->numInputs() - num_uncertainty);

    C1_ = this->ss_->C().head_rows(num_uncertainty);
    C2_ = this->ss_->C().tail_rows(this->ss_->numOutputs() - num_uncertainty);

    D11_ = this->ss_->D().submat(0, 0, num_uncertainty - 1, num_uncertainty - 1);
    D12_ = this->ss_->D().submat(0, num_uncertainty,
                                  num_uncertainty - 1, this->ss_->numInputs() - 1);
    D21_ = this->ss_->D().submat(num_uncertainty, 0,
                                  this->ss_->numOutputs() - 1, num_uncertainty - 1);
    D22_ = this->ss_->D().submat(num_uncertainty, num_uncertainty,
                                  this->ss_->numOutputs() - 1, this->ss_->numInputs() - 1);

}

template <class _StateSpace, int num_uncertainty>
UpperLFT<_StateSpace, num_uncertainty>::~UpperLFT(){

}

}
