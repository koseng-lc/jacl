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

//-- force to have fixed size
template<std::size_t num_states, std::size_t num_inputs, std::size_t num_outputs, class PhysicalParam = int, class ...Rest>
class StateSpace{
public:
    typedef std::function<double(StateSpace)> Formula;
    typedef std::vector<Formula> Formulas;

public:
    static constexpr auto n_states{ num_states };
    static constexpr auto n_inputs{ num_inputs };
    static constexpr auto n_outputs{ num_outputs };

    template <typename _PhysicalParam = PhysicalParam,
              typename std::enable_if<
                  std::is_same<typename std::decay<_PhysicalParam>::type,
                               PhysicalParameter>::value, int>::type* = nullptr>
    StateSpace(_PhysicalParam* _param, Rest*... _rest)
        : A_(num_states, num_states)
        , B_(num_states, num_inputs)
        , C_(num_outputs, num_states)
        , D_(num_outputs, num_inputs){

        push(_param, _rest...);
    }

    StateSpace();
    ~StateSpace();

    arma::mat const& A() const{ return A_; }
    arma::mat const& B() const{ return B_; }
    arma::mat const& C() const{ return C_; }
    arma::mat const& D() const{ return D_; }

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

    auto setA(const arma::mat& _A) -> void{
        A_ = _A;
    }

    auto setB(const Formulas& f) -> void{
        fB_.clear();
        fB_.insert(fB_.end(), f.begin(), f.end());
    }

    auto setB(const arma::mat& _B) -> void{
        B_ = _B;
    }

    auto setC(const Formulas& f) -> void{
        fC_.clear();
        fC_.insert(fC_.end(), f.begin(), f.end());
    }

    auto setC(const arma::mat& _C) -> void{
        C_ = _C;
    }

    auto setD(const Formulas& f) -> void{
        fD_.clear();
        fD_.insert(fD_.end(), f.begin(), f.end());
    }

    auto setD(const arma::mat& _D) -> void{
        D_ = _D;
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

protected:
    Formulas fA_;
    Formulas fB_;
    Formulas fC_;
    Formulas fD_;

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

template <std::size_t num_states, std::size_t num_inputs, std::size_t num_outputs, class PhysicalParam, class ...Rest>
StateSpace<num_states, num_inputs, num_outputs, PhysicalParam, Rest...>::StateSpace()
    : A_(num_states, num_states)
    , B_(num_states, num_inputs)
    , C_(num_outputs, num_states)
    , D_(num_outputs, num_inputs){
}

template<std::size_t num_states, std::size_t num_inputs, std::size_t num_outputs, class PhysicalParam, class ...Rest>
StateSpace<num_states, num_inputs, num_outputs, PhysicalParam, Rest...>::~StateSpace(){

}

template<std::size_t num_states, std::size_t num_inputs, std::size_t num_outputs, class PhysicalParam, class ...Rest>
void StateSpace<num_states, num_inputs, num_outputs, PhysicalParam, Rest...>::formulaToMat(){

    assert(fA_.size () > 0 && fB_.size() > 0 && fC_.size() > 0 && fD_.size() > 0);

    for(std::size_t i(0); i < num_states; i++){
        for(std::size_t j(0); j < num_states; j++){
            A_(i, j) = fA_[i * num_states + j](*this);
        }
    }

    for(std::size_t i(0); i < num_states; i++){
        for(std::size_t j(0); j < num_inputs; j++){
            B_(i, j) = fB_[i * num_inputs + j](*this);
        }
    }

    for(std::size_t i(0); i < num_outputs; i++){
        for(std::size_t j(0); j < num_states; j++){
            C_(i, j) = fC_[i * num_states + j](*this);
        }
    }

    for(std::size_t i(0); i < num_outputs; i++){
        for(std::size_t j(0); j < num_inputs; j++){
            D_(i, j) = fD_[i * num_inputs + j](*this);
        }
    }

}

}
