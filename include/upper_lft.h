/**
*   @author : koseng (Lintang)
*/

#pragma once

#include "state_space.h"
#include "lft.h"

namespace jacl{

template <int num_states,
          int num_inputs,
          int num_outputs,
          int num_uncertainty = 0,
          class SSpace>
class UpperLFT:public LFT<SSpace>{
public:
    UpperLFT(SSpace* _ss);
    ~UpperLFT();

    //-- setters and getters
    inline arma::mat const& B1() const{ return B1_; }
    arma::mat& B1(){ return B1_; }

    inline arma::mat const& C1() const{ return C1_; }
    arma::mat& C1(){ return C1_; }

    inline arma::mat const& D11() const{ return D11_; }
    arma::mat& D11(){ return D11_; }

    inline arma::mat const& D12() const{ return D12_; }
    arma::mat& D12(){ return D12_; }

    inline arma::mat const& D21() const{ return D21_; }
    arma::mat& D21(){ return D21_; }

    inline arma::mat const& delta() const{ return delta_; }
    arma::mat& delta(){ return delta_; }

    bool checkAssumption1();
    bool checkAssumption2();
    bool checkAssumption3();
    bool checkAssumption4();
    bool checkAllAssumption();

private:
    StateSpace<num_states,
               num_inputs + num_uncertainty,
               num_outputs + num_uncertainty,
               PhysicalParameter,
               PhysicalParameter,
               PhysicalParameter,
               PhysicalParameter,
               PhysicalParameter,
               PhysicalParameter> realization_;

    //-- And the remains B1, C1, D11, D12, D21
    arma::mat B1_;
    arma::mat C1_;
    arma::mat D11_;
    arma::mat D12_;
    arma::mat D21_;

    //-- Delta Matrix
    arma::mat delta_;

};

template <int num_states,
          int num_inputs,
          int num_outputs,
          int num_uncertainty,
          class SSpace>
UpperLFT<num_states, num_inputs, num_outputs, num_uncertainty, SSpace>::UpperLFT(SSpace* _ss)
    : LFT(_ss){

}

template <int num_states,
          int num_inputs,
          int num_outputs,
          int num_uncertainty,
          class SSpace>
UpperLFT<num_states, num_inputs, num_outputs, num_uncertainty, SSpace>::~UpperLFT(){

}

template <int num_states,
          int num_inputs,
          int num_outputs,
          int num_uncertainty,
          class SSpace>
bool UpperLFT<num_states, num_inputs, num_outputs, num_uncertainty, SSpace>::checkAssumption1(){

}

template <int num_states,
          int num_inputs,
          int num_outputs,
          int num_uncertainty,
          class SSpace>
bool UpperLFT<num_states, num_inputs, num_outputs, num_uncertainty, SSpace>::checkAssumption2(){

}

template <int num_states,
          int num_inputs,
          int num_outputs,
          int num_uncertainty,
          class SSpace>
bool UpperLFT<num_states, num_inputs, num_outputs, num_uncertainty, SSpace>::checkAssumption3(){

}

template <int num_states,
          int num_inputs,
          int num_outputs,
          int num_uncertainty,
          class SSpace>
bool UpperLFT<num_states, num_inputs, num_outputs, num_uncertainty, SSpace>::checkAssumption4(){

}

template <int num_states,
          int num_inputs,
          int num_outputs,
          int num_uncertainty,
          class SSpace>
bool UpperLFT<num_states, num_inputs, num_outputs, num_uncertainty, SSpace>::checkAllAssumption(){

}

}
