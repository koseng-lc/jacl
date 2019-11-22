/**
*   @author : koseng (Lintang)
*   @brief : jacl State Space module
*/

#pragma once

#include <functional>
#include <vector>
#include <type_traits>

#include <boost/range/combine.hpp>
#include <boost/tuple/tuple.hpp>

#include <armadillo>

#include "physical_parameter.h"
#include "linear_algebra.h"

namespace jacl{

namespace{
    namespace linalg = linear_algebra;
}

template<int num_states, int num_inputs, int num_outputs, class PhysicalParam = int, class ...Rest>
class StateSpace{
public:
    using Formula = std::function<double(StateSpace)>;
    using Formulas = std::vector<Formula>;

public:
    template <typename _PhysicalParam = PhysicalParam,
              typename std::enable_if<
                  std::is_same<typename std::decay<_PhysicalParam>::type,
                               PhysicalParameter>::value, int>::type* = nullptr>
    StateSpace(_PhysicalParam* _param, Rest*... _rest)
        : A_(num_states, num_states)
        , B_(num_states, num_inputs)
        , C_(num_outputs, num_states)
        , D_(num_outputs, num_inputs){

    //    std::cout << sizeof...(Rest) << std::endl;

        push(_param, _rest...);
    }

    StateSpace();
    ~StateSpace();

    arma::mat const& A() const{
        return A_;
    }

    arma::mat const& B() const{
        return B_;
    }

    arma::mat const& C() const{
        return C_;
    }

    arma::mat const& D() const{
        return D_;
    }

    double const& param(int _index) const{
        return params_[_index]->perturbed;
    }

    double& param(int _index){
        return params_[_index]->perturbed;
    }

    auto setA(const Formulas& f) -> void{
        fA_.clear();
        fA_.insert(fA_.end(), f.begin(), f.end());
    }

    auto setB(const Formulas& f) -> void{
        fB_.clear();
        fB_.insert(fB_.end(), f.begin(), f.end());
    }

    auto setC(const Formulas& f) -> void{
        fC_.clear();
        fC_.insert(fC_.end(), f.begin(), f.end());
    }

    auto setD(const Formulas& f) -> void{
        fD_.clear();
        fD_.insert(fD_.end(), f.begin(), f.end());
    }

    void formulaToMat();

    inline int numStates() const{
        return num_states;
    }

    inline int numInputs() const{
        return num_inputs;
    }

    inline int numOutputs() const{
        return num_outputs;
    }

    inline double const& operator () (int _index) const{
        return param(_index);
    }

    inline double operator () (double _constant){
        return _constant;
    }

    inline double& operator () (int _index){
        return param(_index);
    }

//    struct Forml{

//    };

private:

    Formulas fA_;
    Formulas fB_;
    Formulas fC_;
    Formulas fD_;

    // TODO : change the arguement to pointer to prevent object slicing
    std::vector<PhysicalParameter* > params_;

    arma::mat A_;
    arma::mat B_;
    arma::mat C_;
    arma::mat D_;

    template <typename _PhysicalParam = PhysicalParam>
    auto push(_PhysicalParam* _param)
    -> typename std::enable_if<
    std::is_same<typename std::decay<_PhysicalParam>::type,
                 PhysicalParameter>::value, void>::type{
        params_.push_back(_param);
    }

    template <typename _PhysicalParam = PhysicalParam, class ...PhysicalParamRest>
    auto push(_PhysicalParam* _param, PhysicalParamRest*... _rest)
    -> typename std::enable_if<
    std::is_same<typename std::decay<_PhysicalParam>::type,
                 PhysicalParameter>::value, void>::type{
        push(_param);
        push(_rest...);
    }

};

//template <int num_states, int num_inputs, int num_outputs, class PhysicalParam, class ...Rest>
//StateSpace<num_states, num_inputs, num_outputs, PhysicalParam, Rest...>::StateSpace(PhysicalParam _param, Rest... _rest){


//}

template <int num_states, int num_inputs, int num_outputs, class PhysicalParam, class ...Rest>
StateSpace<num_states, num_inputs, num_outputs, PhysicalParam, Rest...>::StateSpace()
    : A_(num_states, num_states)
    , B_(num_states, num_inputs)
    , C_(num_outputs, num_states)
    , D_(num_outputs, num_inputs){
}

template<int num_states, int num_inputs, int num_outputs, class PhysicalParam, class ...Rest>
StateSpace<num_states, num_inputs, num_outputs, PhysicalParam, Rest...>::~StateSpace(){

}

template<int num_states, int num_inputs, int num_outputs, class PhysicalParam, class ...Rest>
void StateSpace<num_states, num_inputs, num_outputs, PhysicalParam, Rest...>::formulaToMat(){

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
