/**
*   @author : koseng (Lintang)
*   @brief : jacl nonlinear state space module
*/

#pragma once

#include <functional>
#include <vector>
#include <type_traits>
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>

#include <jacl/pattern/observer.hpp>
#include <jacl/physical_parameter.hpp>
#include <jacl/medium.hpp>

namespace jacl{ namespace state_space{

//-- force to have fixed size
template<typename Scalar,
         std::size_t num_states,
         std::size_t num_inputs,
         std::size_t num_outputs,
         typename PhysicalParam = int,
         typename ...Rest>
class NonLinear:public pattern::Subject{
public:
    using scalar_t = Scalar;
    using formula_t = typename std::function<scalar_t(NonLinear)>;
    using formulas_t = std::vector<formula_t>;
    //-- it's needed in BaseSystem
    using state_matrix_t = typename arma::Mat<scalar_t>::template fixed<num_states,num_states>;
    using input_matrix_t = typename arma::Mat<scalar_t>::template fixed<num_states,num_inputs>;
    using output_matrix_t = typename arma::Mat<scalar_t>::template fixed<num_outputs,num_states>;
    using feedforward_matrix_t = typename arma::Mat<scalar_t>::template fixed<num_outputs,num_inputs>;

public:
    template <typename _PhysicalParam = PhysicalParam,
              typename std::enable_if_t<
                  std::is_same<typename std::decay_t<_PhysicalParam>,
                               PhysicalParameter<scalar_t>>::value, int>* = nullptr>
    explicit NonLinear(_PhysicalParam* _param, Rest*... _rest)
        : state_fn_(n_states + n_inputs)
        , output_fn_(n_outputs){
        push(_param, _rest...);
    }

    NonLinear()
        : state_fn_(n_states + n_inputs)
        , output_fn_(n_outputs){
    }
    ~NonLinear(){}

    inline auto& stateFn(){
        return state_fn_;
    }

    inline const auto& stateFn() const{
        return state_fn_;
    }

    inline auto& outputFn(){
        return output_fn_;
    }

    inline const auto& outputFn() const{
        return output_fn_;
    }

private:
    template <typename _PhysicalParam = PhysicalParam>
    auto push(_PhysicalParam* _param)
        -> typename std::enable_if_t<
        std::is_same<typename std::decay_t<_PhysicalParam>,
                 PhysicalParameter<scalar_t>>::value, void>{
        params_.push_back(_param);
    }

    template <typename _PhysicalParam = PhysicalParam, class ...PhysicalParamRest>
    auto push(_PhysicalParam* _param, PhysicalParamRest*... _rest)
        -> typename std::enable_if_t<
        std::is_same<typename std::decay_t<_PhysicalParam>,
                 PhysicalParameter<scalar_t>>::value, void>{
        push(_param);
        push(_rest...);
    }

    auto param(int _index) const -> decltype(std::declval<PhysicalParameter<scalar_t>>().perturbed) const&{
        return params_[_index]->perturbed;
    }

    auto param(int _index) -> decltype(std::declval<PhysicalParameter<scalar_t>>().perturbed)&{
        return params_[_index]->perturbed;
    } 
    
public:
    inline const auto& operator () (int _index) const{        
        return _index < (n_states + n_inputs) ?
            sig_[_index] : param(_index - (n_states + n_inputs));
    }

    inline const auto& operator () (double _constant){
        return _constant;
    }

    inline double& operator () (int _index){
        return _index < (n_states + n_inputs) ?
            sig_[_index] : param(_index - (n_states + n_inputs));
    }

public:
    static constexpr auto n_states{ num_states };
    static constexpr auto n_inputs{ num_inputs };
    static constexpr auto n_outputs{ num_outputs };

private:
    formulas_t state_fn_;
    formulas_t output_fn_;
    std::vector<double> sig_;
    std::vector<PhysicalParameter<scalar_t>* > params_;
    
private:
    friend class ::jacl::state_space::detail::NonLinearStateSpaceClient<NonLinear>;

};

} }