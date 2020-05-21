/**
*   @author : koseng (Lintang)
*   @brief : jacl state space module
*/

#pragma once

#include <functional>
#include <vector>
#include <type_traits>
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>

#include <jacl/pattern/observer.hpp>
#include <jacl/physical_parameter.hpp>

namespace jacl{ namespace state_space{

//-- force to have fixed size
template<typename Scalar,
         std::size_t num_states,
         std::size_t num_inputs,
         std::size_t num_outputs,
         typename PhysicalParam = int,
         typename ...Rest>
class Linear:public pattern::Subject{
public:
    using scalar_t = Scalar;
    using formula_t = typename std::function<scalar_t(Linear)>;
    using formulas_t = std::vector<formula_t>;        
    using state_matrix_t = typename arma::Mat<scalar_t>::template fixed<num_states,num_states>;
    using input_matrix_t = typename arma::Mat<scalar_t>::template fixed<num_states,num_inputs>;
    using output_matrix_t = typename arma::Mat<scalar_t>::template fixed<num_outputs,num_states>;
    using feedforward_matrix_t = typename arma::Mat<scalar_t>::template fixed<num_outputs,num_inputs>;

public:    
    template <typename _PhysicalParam = PhysicalParam,
              typename std::enable_if_t<
                  std::is_same<typename std::decay_t<_PhysicalParam>,
                               PhysicalParameter<scalar_t>>::value, int>* = nullptr>
    explicit Linear(_PhysicalParam* _param, Rest*... _rest){
        push(_param, _rest...);
    }
    Linear(const state_matrix_t& _A,
           const input_matrix_t& _B,
           const output_matrix_t& _C,
           const feedforward_matrix_t& _D)
        : A_(_A)
        , B_(_B)
        , C_(_C)
        , D_(_D){

    }
    Linear(){}
    ~Linear(){}

    const auto& A() const{ return A_; }
    const auto& B() const{ return B_; }
    const auto& C() const{ return C_; }
    const auto& D() const{ return D_; }    

    enum class FormulaToMatMode:std::uint8_t{
        A = 0b00000001,
        B = 0b00000010,
        C = 0b00000100,
        D = 0b00001000,
        All = 0b00001111,
    };

    auto setA(const formulas_t& f){
        fA_.clear();
        fA_.insert(fA_.end(), f.begin(), f.end());
        formulaToMat(FormulaToMatMode::A);        
    }

    auto setA(const state_matrix_t& _A){
        A_ = _A;
        this->notify();
    }

    auto setB(const formulas_t& f){
        fB_.clear();
        fB_.insert(fB_.end(), f.begin(), f.end());
        formulaToMat(FormulaToMatMode::B);
    }

    auto setB(const input_matrix_t& _B){
        B_ = _B;
        this->notify();
    }

    auto setC(const formulas_t& f){
        fC_.clear();
        fC_.insert(fC_.end(), f.begin(), f.end());
        formulaToMat(FormulaToMatMode::C);
    }

    auto setC(const output_matrix_t& _C){
        C_ = _C;
        this->notify();
    }

    auto setD(const formulas_t& f){
        fD_.clear();
        fD_.insert(fD_.end(), f.begin(), f.end());
        formulaToMat(FormulaToMatMode::D);
    }

    auto setD(const feedforward_matrix_t& _D){
        D_ = _D;
        this->notify();
    }

    auto formulaToMat(FormulaToMatMode _mode = FormulaToMatMode::All);

protected:
    template <typename _PhysicalParam = PhysicalParam>
    auto push(_PhysicalParam* _param)
    -> typename std::enable_if_t<
    std::is_same<typename std::decay_t<_PhysicalParam>,
                 PhysicalParameter<scalar_t>>::value, void>{
        params_.push_back(_param);
    }

    template <typename _PhysicalParam = PhysicalParam, typename ...PhysicalParamRest>
    auto push(_PhysicalParam* _param, PhysicalParamRest*... _rest)
    -> typename std::enable_if_t<
    std::is_same<typename std::decay_t<_PhysicalParam>,
                 PhysicalParameter<scalar_t>>::value, void>{
        push(_param);
        push(_rest...);
    }

    auto param(int _index) const
        -> decltype(std::declval<PhysicalParameter<scalar_t>>().perturbed) const&{
        return params_[_index]->perturbed;
    }

    auto param(int _index)
        -> decltype(std::declval<PhysicalParameter<scalar_t>>().perturbed)&{
        return params_[_index]->perturbed;
    }

public:
    inline auto operator () (int _index) const -> decltype(param(_index)) const&{
        return param(_index);
    }

    inline const auto& operator () (double _constant){
        return _constant;
    }

    inline auto operator () (int _index) -> decltype(param(_index))&{        
        return param(_index);
    }

public:
    static constexpr auto n_states{ num_states };
    static constexpr auto n_inputs{ num_inputs };
    static constexpr auto n_outputs{ num_outputs };

protected:
    formulas_t fA_;
    formulas_t fB_;
    formulas_t fC_;
    formulas_t fD_;

    std::vector<PhysicalParameter<scalar_t>* > params_;

    state_matrix_t A_;
    input_matrix_t B_;
    output_matrix_t C_;
    feedforward_matrix_t D_;

};

template<typename Scalar,
         std::size_t num_states,
         std::size_t num_inputs,
         std::size_t num_outputs,
         typename PhysicalParam,
         typename ...Rest>
auto Linear<Scalar,
                      num_states,
                      num_inputs,
                      num_outputs,
                      PhysicalParam,
                      Rest...>::formulaToMat(FormulaToMatMode _mode){
    assert(fA_.size () > 0 && fB_.size() > 0 && fC_.size() > 0 && fD_.size() > 0);
    uint8_t mode = static_cast<typename std::underlying_type<FormulaToMatMode>::type>(_mode);
    bool change_A = static_cast<typename std::underlying_type<FormulaToMatMode>::type>(FormulaToMatMode::A) & mode;
    bool change_B = static_cast<typename std::underlying_type<FormulaToMatMode>::type>(FormulaToMatMode::B) & mode;
    bool change_C = static_cast<typename std::underlying_type<FormulaToMatMode>::type>(FormulaToMatMode::C) & mode;
    bool change_D = static_cast<typename std::underlying_type<FormulaToMatMode>::type>(FormulaToMatMode::D) & mode;
    for(std::size_t i(0); i < num_states; i++){
        if(change_A){
            for(std::size_t j(0); j < num_states; j++)
                A_(i, j) = fA_[i * num_states + j](*this);
        }        
        if(change_B){
            for(std::size_t j(0); j < num_inputs; j++)
                B_(i, j) = fB_[i * num_inputs + j](*this);
        }
    }

    for(std::size_t i(0); i < num_outputs; i++){
        if(change_C){
            for(std::size_t j(0); j < num_states; j++)
                C_(i, j) = fC_[i * num_states + j](*this);
        }
        if(change_D){
            for(std::size_t j(0); j < num_inputs; j++)
                D_(i, j) = fD_[i * num_inputs + j](*this);
        }
    }
    this->notify();
}

} }
