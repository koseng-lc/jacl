#pragma once

#include <tuple>

#include "lti_common.h"
#include "are.h"
#include "upper_lft.h"
#include "lower_lft.h"
#include "state_space.h"
#include "traits.h"

namespace jacl{

namespace synthesis{

namespace{
    namespace linalg = linear_algebra;
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
class Hinf{
public:
    template<typename T = _StateSpace,
             typename std::enable_if<
               traits::is_state_space<T>::value>::type* = nullptr>
    Hinf(T* _ss)
        : ss_(_ss)
        , llft_(ss_)
        , are_solver1_(ss_)
        , are_solver2_(ss_){

    }    

    ~Hinf();

    void init();
    void solve();

private:
    _StateSpace* ss_;

    using LLFT = LowerLFT<_StateSpace,
                          performance_size,
                          perturbation_size,
                          _StateSpace::n_outputs,
                          _StateSpace::n_inputs>;

    LLFT llft_;

    ARE<_StateSpace> are_solver1_;
    ARE<_StateSpace> are_solver2_;

    bool checkAssumption1();
    bool checkAssumption2();
    bool checkAssumption3();
    bool checkAssumption4();
    bool checkAllAssumption();

};

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
Hinf<_StateSpace,
    performance_size,
    perturbation_size>::~Hinf(){

}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
void Hinf<_StateSpace,
    performance_size,
    perturbation_size>::init(){

}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool Hinf<_StateSpace,
    performance_size,
    perturbation_size>::checkAssumption1(){

    bool ctrb = common::controlable(llft_.A(), llft_.B2());
    bool obsv = common::observable(llft_.A(), llft_.C2());
    std::cout << "Assumption 1 : " << std::endl;
    std::cout << std::boolalpha << ctrb << " ; " << obsv << std::endl;
    return ctrb & obsv;
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool Hinf<_StateSpace,
    performance_size,
    perturbation_size>::checkAssumption2(){

    bool ok(false); //-- RVO
    arma::mat I_in(ss_->numInputs(), ss_->numInputs(), arma::fill::eye);
    arma::mat expected_D12(
                arma::join_vert(
                    arma::mat(performance_size - ss_->numInputs(), ss_->numInputs(), arma::fill::zeros),
                    I_in));

    if(performance_size > ss_->numInputs()){
        if(arma::approx_equal(llft_.D12(), expected_D12, "absdiff", .0)){
            ok = true;
        }else{
            //-- do normalization here using SVD
        }
    }else if(performance_size == ss_->numInputs()){
        if(arma::approx_equal(llft_.D12(), I_in, "absdiff", .0)){
            std::cout << "[Hinf] Special cases trigerred !!!" << std::endl;
            //-- do special cases here
            ok = true;
        }else{
            //-- do normalization here
        }
    }else{ //-- if not full column rank
        ok = false;
    }

    bool ok2(false); //-- RVO
    arma::mat I_out(ss_->numOutputs(), ss_->numOutputs(), arma::fill::eye);
    arma::mat expected_D21(
                arma::join_horiz(
                    arma::mat(ss_->numOutputs(), perturbation_size - ss_->numOutputs(), arma::fill::zeros),
                    I_out));

    if(perturbation_size > ss_->numOutputs()){
        if(arma::approx_equal(llft_.D12(), expected_D21, "absdiff", .0)){
            ok2 = true;
        }else{
            //-- do normalization here using SVD
        }
    }else if(perturbation_size == ss_->numOutputs()){
        if(arma::approx_equal(llft_.D12(), I_out, "absdiff", .0)){
            std::cout << "[Hinf] Special cases trigerred !!!" << std::endl;
            //-- do special cases here
            ok2 = true;
        }else{
            //-- do normalization here
            ok2 = false;
        }
    }else{ //-- if not full row rank
        ok2 = false;
    }

    return ok & ok2;
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool Hinf<_StateSpace,
    performance_size,
    perturbation_size>::checkAssumption3(){
    return true;
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool Hinf<_StateSpace,
    performance_size,
    perturbation_size>::checkAssumption4(){
    return true;
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool Hinf<_StateSpace,
    performance_size,
    perturbation_size>::checkAllAssumption(){
    return checkAssumption1()
            & checkAssumption2()
            & checkAssumption3()
            & checkAssumption4();
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
void Hinf<_StateSpace,
    performance_size,
    perturbation_size>::solve(){

    std::cout << std::boolalpha << checkAllAssumption() << std::endl;
}

}

}
