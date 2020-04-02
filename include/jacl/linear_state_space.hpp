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

#include <jacl/physical_parameter.hpp>

namespace jacl{

//-- force to have fixed size
template<std::size_t num_states, std::size_t num_inputs, std::size_t num_outputs, class PhysicalParam = int, class ...Rest>
class LinearStateSpace{
public:
    typedef std::function<double(LinearStateSpace)> Formula;
    typedef std::vector<Formula> Formulas;

public:    
    template <typename _PhysicalParam = PhysicalParam,
              typename std::enable_if<
                  std::is_same<typename std::decay<_PhysicalParam>::type,
                               PhysicalParameter>::value, int>::type* = nullptr>
    LinearStateSpace(_PhysicalParam* _param, Rest*... _rest)
        : A_(num_states, num_states)
        , B_(num_states, num_inputs)
        , C_(num_outputs, num_states)
        , D_(num_outputs, num_inputs){

        push(_param, _rest...);
    }

    LinearStateSpace();
    ~LinearStateSpace();

    auto A() const -> arma::mat const&{ return A_; }
    auto B() const -> arma::mat const&{ return B_; }
    auto C() const -> arma::mat const&{ return C_; }
    auto D() const -> arma::mat const&{ return D_; }

    auto param(int _index) const -> decltype(std::declval<PhysicalParameter>().perturbed) const&{
        return params_[_index]->perturbed;
    }

    auto param(int _index) -> decltype(std::declval<PhysicalParameter>().perturbed)&{
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

    auto formulaToMat() -> void;

    inline auto numStates() const -> int{
        return num_states;
    }

    inline auto numInputs() const -> int{
        return num_inputs;
    }

    inline auto numOutputs() const -> int{
        return num_outputs;
    }

    inline auto operator () (int _index) const -> double const&{
        return param(_index);
    }

    inline auto operator () (double _constant) -> double const&{
        return _constant;
    }

    inline auto operator () (int _index) -> double&{
        return param(_index);
    }

protected:
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

public:
    static constexpr auto n_states{ num_states };
    static constexpr auto n_inputs{ num_inputs };
    static constexpr auto n_outputs{ num_outputs };

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
};

template <std::size_t num_states, std::size_t num_inputs, std::size_t num_outputs, class PhysicalParam, class ...Rest>
LinearStateSpace<num_states, num_inputs, num_outputs, PhysicalParam, Rest...>::LinearStateSpace()
    : A_(num_states, num_states)
    , B_(num_states, num_inputs)
    , C_(num_outputs, num_states)
    , D_(num_outputs, num_inputs){
}

template<std::size_t num_states, std::size_t num_inputs, std::size_t num_outputs, class PhysicalParam, class ...Rest>
LinearStateSpace<num_states, num_inputs, num_outputs, PhysicalParam, Rest...>::~LinearStateSpace(){

}

template<std::size_t num_states, std::size_t num_inputs, std::size_t num_outputs, class PhysicalParam, class ...Rest>
auto LinearStateSpace<num_states, num_inputs, num_outputs, PhysicalParam, Rest...>::formulaToMat() -> void{

    assert(fA_.size () > 0 && fB_.size() > 0 && fC_.size() > 0 && fD_.size() > 0);

    for(std::size_t i(0); i < num_states; i++){
        for(std::size_t j(0); j < num_states; j++){
            A_(i, j) = fA_[i * num_states + j](*this);
        }

        for(std::size_t j(0); j < num_inputs; j++){
            B_(i, j) = fB_[i * num_inputs + j](*this);
        }
    }

    for(std::size_t i(0); i < num_outputs; i++){
        for(std::size_t j(0); j < num_states; j++){
            C_(i, j) = fC_[i * num_states + j](*this);
        }

        for(std::size_t j(0); j < num_inputs; j++){
            D_(i, j) = fD_[i * num_inputs + j](*this);
        }
    }
}

}
