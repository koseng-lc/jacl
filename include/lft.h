#pragma once

#include "state_space.h"

namespace jacl{
template <int num_perturbation = 0, class SSpace>
class ULFT{
public:
    ULFT(SSpace* _ss);
    ~ULFT();

    arma::mat const& B1() const{
        return B1_;
    }

    arma::mat const& C1() const{
        return C1_;
    }

    arma::mat const& D11() const{
        return D11_;
    }

    arma::mat const& D12() const{
        return D12_;
    }

    arma::mat const& D21() const{
        return D21_;
    }

private:

    //-- Realization of Upper LFT
    //-- This one consist A, B2, C2, D22
    SSpace* ss_;

    //-- And the remains B1, C1, D11, D12, D21
    arma::mat B1_;
    arma::mat C1_;
    arma::mat D11_;
    arma::mat D12_;
    arma::mat D21_;

    //-- Delta Matrix
    arma::mat delta_;

};

template <int num_perturbation, class SSpace>
ULFT<num_perturbation, SSpace>::ULFT(SSpace* _ss)
    : ss_(_ss){

}

template <int num_perturbation, class SSpace>
ULFT<num_perturbation, SSpace>::~ULFT(){

}

}
