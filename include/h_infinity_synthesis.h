#pragma once

#include "linear_algebra.h"
#include "are.h"

namespace jacl{

namespace synthesis{

namespace{
    namespace linalg = linear_algebra;
}

template <class SSpace, int num_perturbation = 0>
class Hinf{
public:
    Hinf(SSpace* _ss);
    ~Hinf();

    void init();

private:
    SSpace* ss_;

    ARE<SSpace> are_solver1_;
    ARE<SSpace> are_solver2_;

    //-- Partitioning Realization of Interconnection Matrix
    arma::mat A;
    arma::mat B1;
    arma::mat B2;
    arma::mat C1;
    arma::mat C2;
    arma::mat D11;
    arma::mat D12;
    arma::mat D21;
    arma::mat D22;

};

template <class SSpace, int num_perturbation>
Hinf<SSpace, num_perturbation>::Hinf(SSpace* _ss)
    : ss_(_ss)
    , are_solver1_(ss_)
    , are_solver2_(ss_){

}

template <class SSpace, int num_perturbation>
Hinf<SSpace, num_perturbation>::~Hinf(){

}

template <class SSpace, int num_perturbation>
void Hinf<SSpace, num_perturbation>::init(){
    A = ss_->A();

    B1 = ss_->B().head_cols(num_perturbation);
    B2 = ss_->B().tail_cols(ss_->numStates());
    C1 = ss_->C().head_rows(num_perturbation);
    C2 = ss_->C().tail_rows(ss_->numStates());
}

}

}
