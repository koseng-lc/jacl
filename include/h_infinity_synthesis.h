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
             typename std::enable_if<traits::is_state_space<T>::value>::type* = nullptr>
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

    static constexpr auto INPUT_SIZE{_StateSpace::n_inputs - perturbation_size};
    static constexpr auto OUTPUT_SIZE{_StateSpace::n_outputs - performance_size};

    using LLFT = LowerLFT<_StateSpace,
                          performance_size,
                          perturbation_size,
                          OUTPUT_SIZE,
                          INPUT_SIZE>;

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

    bool ctrb = common::stabilizable(llft_.A(), llft_.B2());//common::controlable(llft_.A(), llft_.B2());
    bool obsv = common::detectability(llft_.A(), llft_.C2());//common::observable(llft_.A(), llft_.C2());
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

    arma::mat U1,R1;
    arma::mat U2,R2;

    bool ok(false); //-- RVO
    arma::mat I_in(INPUT_SIZE, INPUT_SIZE, arma::fill::eye);
    arma::mat expected_D12 = arma::join_vert(arma::mat(performance_size - INPUT_SIZE, INPUT_SIZE, arma::fill::zeros), I_in);

    if(performance_size > INPUT_SIZE){
        if(!arma::approx_equal(llft_.D12(), expected_D12, "absdiff", .0)){
            ok = true;
        }else{
            //-- do normalization here using SVD

            arma::mat zeros1(performance_size - INPUT_SIZE, INPUT_SIZE, arma::fill::zeros);
            arma::mat zeros2(INPUT_SIZE, performance_size - INPUT_SIZE, arma::fill::zeros);
            arma::mat eye(performance_size - INPUT_SIZE, performance_size - INPUT_SIZE, arma::fill::eye);

            arma::mat U, V;
            arma::vec s;
            arma::svd(U,s,V,llft_.D12());

            arma::mat S = arma::diagmat(s);
            arma::mat recip_S(arma::size(S), arma::fill::zeros);
            for(int i(0); i < S.n_rows; i++){
                recip_S(i, i) = 1./S(i,i);
            }

            arma::mat T = arma::join_vert(
                            arma::join_horiz(zeros1, eye),
                            arma::join_horiz(recip_S, zeros2)
                        );
        }
    }else if(performance_size == INPUT_SIZE){
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
    arma::mat I_out(OUTPUT_SIZE, OUTPUT_SIZE, arma::fill::eye);
    arma::mat expected_D21 = arma::join_horiz(arma::mat(OUTPUT_SIZE, perturbation_size - OUTPUT_SIZE, arma::fill::zeros), I_out);

    if(perturbation_size > OUTPUT_SIZE){
        if(arma::approx_equal(llft_.D12(), expected_D21, "absdiff", .0)){
            ok2 = true;
        }else{
            //-- do normalization here using SVD

            arma::mat zeros1(OUTPUT_SIZE, perturbation_size - OUTPUT_SIZE, arma::fill::zeros);
            arma::mat zeros2(perturbation_size - OUTPUT_SIZE, OUTPUT_SIZE, arma::fill::zeros);
            arma::mat eye(perturbation_size - OUTPUT_SIZE, perturbation_size - OUTPUT_SIZE, arma::fill::eye);

            arma::mat U, V;
            arma::vec s;
            arma::svd(U,s,V,llft_.D21());

            arma::mat S = arma::diagmat(s);
            arma::mat recip_S(arma::size(S), arma::fill::zeros);
            for(int i(0); i < S.n_rows; i++){
                recip_S(i, i) = 1./S(i,i);
            }

            arma::mat T = arma::join_vert(
                            arma::join_horiz(zeros1, recip_S),
                            arma::join_horiz(eye, zeros2)
                        );

//            arma::mat R = arma::inv(T) * V.t();
//            U.print("U : ");
//            V.print("V : ");
//            R.print("R : ");
//            arma::mat test = R.t() * R;test.print("TEST : ");
        }
    }else if(perturbation_size == OUTPUT_SIZE){
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

    /*
     * Based on Lemma 12.6 in Essential Robust Control
     * This condition is equivalent with ((I-D*inv(R)*conj(D))*C, A-B*inv(R)*conj(D)*C) has no
     * unobservable modes on the jw axis
     */

    arma::mat D12_conj( arma::conj(llft_.D12()) );
    arma::mat R( D12_conj*llft_.D12() );
    arma::mat R_inv( arma::inv(R) );
    arma::mat I(performance_size, performance_size, arma::fill::eye);

    arma::mat temp1, temp2;
    temp1 = llft_.D12()*R_inv*D12_conj;
    temp2 = I - temp1;
    arma::mat C = temp2*llft_.C1();

    temp1 = llft_.B2()*R_inv*D12_conj;
    temp2 = temp1*llft_.C1();
    arma::mat A = llft_.A() - temp2;

    //-- do PBH test here

    return true;
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool Hinf<_StateSpace,
    performance_size,
    perturbation_size>::checkAssumption4(){

    arma::mat D21_conj( arma::conj(llft_.D21()) );
    arma::mat R( llft_.D21()*D21_conj );
    arma::mat R_inv( arma::inv(R) );
    arma::mat I(perturbation_size, perturbation_size, arma::fill::eye);

    arma::mat temp1, temp2;
    temp1 = llft_.B1()*D21_conj*R_inv;
    temp2 = temp1*llft_.C2();
    arma::mat A = llft_.A() - temp2;

    temp1 = D21_conj*R_inv*llft_.D21();
    temp2 = llft_.B1()*temp1;
    arma::mat B = llft_.B1()*temp2;

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
