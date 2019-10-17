/**
*   @author : koseng (Lintang)
*   @brief : JACL State Space module
*/

#pragma once

#include <functional>
#include <vector>

#include <boost/range/combine.hpp>
#include <boost/tuple/tuple.hpp>

#include <armadillo>

#include "physical_parameter.h"

namespace JACL{

using Mat = arma::mat;

template<int num_states, int num_inputs, int num_outputs, class T, class ...Args>
class StateSpace{
public:
    using Formula = std::function<double(StateSpace)>;
    using Formulas = std::vector<Formula>;

public:
    StateSpace(T _param, Args... _rest);
    ~StateSpace();

    Mat A() const{
        return A_;
    }

    Mat B() const{
        return B_;
    }

    Mat C() const{
        return C_;
    }

    Mat D() const{
        return D_;
    }

    bool controlable();
    bool observable();

    double param(int _index) const{
        return params_[_index].perturbed;
    }

    double& param(int _index){
        return params_[_index].perturbed;
    }

    void setA(Formulas f){
        fA_ = f;
    }

    void setB(Formulas f){
        fB_ = f;
    }

    void setC(Formulas f){
        fC_ = f;
    }

    void setD(Formulas f){
        fD_ = f;
    }

    void formulaToMat();

private:

    Formulas fA_;
    Formulas fB_;
    Formulas fC_;
    Formulas fD_;

    std::vector<PhysicalParameter > params_;

    Mat A_;
    Mat B_;
    Mat C_;
    Mat D_;

    void push(PhysicalParameter _param){
        params_.push_back(_param);
    }

    void push(PhysicalParameter _param, Args... _rest){
        push(_param);
        push(_rest...);
    }

};

template<int num_states, int num_inputs, int num_outputs, class T, class ...Args>
StateSpace<num_states, num_inputs, num_outputs, T, Args...>::StateSpace(T _param, Args... _rest){

    push(_param, _rest...);

    A_.resize(num_states, num_states);
    B_.resize(num_states, num_inputs);
    C_.resize(num_outputs, num_states);
    D_.resize(num_outputs, num_inputs);

}

template<int num_states, int num_inputs, int num_outputs, class T, class ...Args>
StateSpace<num_states, num_inputs, num_outputs, T, Args...>::~StateSpace(){

}

template<int num_states, int num_inputs, int num_outputs, class T, class ...Args>
bool StateSpace<num_states, num_inputs, num_outputs, T, Args...>::controlable(){

}

template<int num_states, int num_inputs, int num_outputs, class T, class ...Args>
bool StateSpace<num_states, num_inputs, num_outputs, T, Args...>::observable(){

}

template<int num_states, int num_inputs, int num_outputs, class T, class ...Args>
void StateSpace<num_states, num_inputs, num_outputs, T, Args...>::formulaToMat(){

    assert(fA_.size () > 0 && fB_.size() > 0 && fC_.size() > 0 && fD_.size() > 0);

    for(int i(0); i < num_states; i++){
        for(int j(0); j < num_states; j++){
            A_(i, j) = fA_[i * num_states + j](*this);
        }
    }

    for(int i(0); i < num_states; i++){
        for(int j(0); j < num_inputs; j++){
            B_(i, j) = fB_[i * num_inputs + j](*this);
        }
    }

    for(int i(0); i < num_outputs; i++){
        for(int j(0); j < num_states; j++){
            C_(i, j) = fC_[i * num_states + j](*this);
        }
    }

    for(int i(0); i < num_outputs; i++){
        for(int j(0); j < num_inputs; j++){
            D_(i, j) = fD_[i * num_inputs + j](*this);
        }
    }

}

}
