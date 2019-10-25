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

    void setA(const Formulas& f){
        fA_.clear();
        fA_.insert(fA_.end(), f.begin(), f.end());
    }

    void setB(const Formulas& f){
        fB_.clear();
        fB_.insert(fB_.end(), f.begin(), f.end());
    }

    void setC(const Formulas& f){
        fC_.clear();
        fC_.insert(fC_.end(), f.begin(), f.end());
    }

    void setD(const Formulas& f){
        fD_.clear();
        fD_.insert(fD_.end(), f.begin(), f.end());
    }

    void formulaToMat();

private:

    Formulas fA_;
    Formulas fB_;
    Formulas fC_;
    Formulas fD_;

    // TODO : change the arguement to pointer
    std::vector<PhysicalParameter > params_;

    Mat A_;
    Mat B_;
    Mat C_;
    Mat D_;

    void push(T _param){
        params_.push_back(_param);
    }

    template <class ...Rest>
    void push(T _param, Rest... _rest){
        push(_param);
        push(_rest...);
    }

};

template<int num_states, int num_inputs, int num_outputs, class T, class ...Args>
StateSpace<num_states, num_inputs, num_outputs, T, Args...>::StateSpace(T _param, Args... _rest){

//    std::cout << sizeof...(Args) << std::endl;

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
