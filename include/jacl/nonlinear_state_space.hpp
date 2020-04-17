/**
*   @author : koseng (Lintang)
*   @brief : jacl Nonlinear State Space module
*/

#pragma once

#include <functional>
#include <vector>
#include <type_traits>

#include <boost/range/combine.hpp>
#include <boost/tuple/tuple.hpp>

#include <armadillo>

#include <jacl/pattern/observer.hpp>
#include <jacl/physical_parameter.hpp>
#include <jacl/medium.hpp>

namespace jacl{

//-- force to have fixed size
template<std::size_t num_states, std::size_t num_inputs, std::size_t num_outputs, class PhysicalParam = int, class ...Rest>
class NonLinearStateSpace:public pattern::Subject{
public:
    typedef std::function<double(NonLinearStateSpace)> Formula;
    typedef std::vector<Formula> Formulas;

public:
    template <typename _PhysicalParam = PhysicalParam,
              typename std::enable_if_t<
                  std::is_same<typename std::decay_t<_PhysicalParam>,
                               PhysicalParameter>::value, int>* = nullptr>
    NonLinearStateSpace(_PhysicalParam* _param, Rest*... _rest)
        : state_fn_(n_states + n_inputs)
        , output_fn_(n_outputs){

        push(_param, _rest...);
    }

    NonLinearStateSpace();
    ~NonLinearStateSpace();           

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
                 PhysicalParameter>::value, void>{
        params_.push_back(_param);
    }

    template <typename _PhysicalParam = PhysicalParam, class ...PhysicalParamRest>
    auto push(_PhysicalParam* _param, PhysicalParamRest*... _rest)
    -> typename std::enable_if_t<
    std::is_same<typename std::decay_t<_PhysicalParam>,
                 PhysicalParameter>::value, void>{
        push(_param);
        push(_rest...);
    }

    auto param(int _index) const -> decltype(std::declval<PhysicalParameter>().perturbed) const&{
        return params_[_index]->perturbed;
    }

    auto param(int _index) -> decltype(std::declval<PhysicalParameter>().perturbed)&{
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
    Formulas state_fn_;
    Formulas output_fn_;
    std::vector<double> sig_;
    std::vector<PhysicalParameter* > params_;
    
private:
    friend class ::jacl::state_space::detail::NonLinearStateSpaceClient<NonLinearStateSpace>;

};

template <std::size_t num_states, std::size_t num_inputs, std::size_t num_outputs, class PhysicalParam, class ...Rest>
NonLinearStateSpace<num_states, num_inputs, num_outputs, PhysicalParam, Rest...>::NonLinearStateSpace()
    : state_fn_(n_states + n_inputs)
    , output_fn_(n_outputs){
}

template<std::size_t num_states, std::size_t num_inputs, std::size_t num_outputs, class PhysicalParam, class ...Rest>
NonLinearStateSpace<num_states, num_inputs, num_outputs, PhysicalParam, Rest...>::~NonLinearStateSpace(){

}

}
