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
class HInf{
public:
    template<typename T = _StateSpace,
             typename std::enable_if<traits::is_state_space<T>::value>::type* = nullptr>
    HInf(T* _ss)
        : ss_(_ss)
        , llft_(ss_)
        , are_solver1_(ss_)
        , are_solver2_(ss_)
        , gam_(1.0){

    }    

    ~HInf();

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

    double gam_;

};

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
HInf<_StateSpace,
    performance_size,
    perturbation_size>::~HInf(){

}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
void HInf<_StateSpace,
    performance_size,
    perturbation_size>::init(){

}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkAssumption1(){

    bool ctrb = common::stabilizable(llft_.A(), llft_.B2());
    bool obsv = common::detectability(llft_.A(), llft_.C2());
    std::cout << "Assumption 1 : " << std::endl;
    std::cout << std::boolalpha << ctrb << " ; " << obsv << std::endl;
    return ctrb & obsv;
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkAssumption2(){

    bool ok(true); //-- RVO
    arma::mat I_in(INPUT_SIZE, INPUT_SIZE, arma::fill::eye);
    arma::mat expected_D12 = arma::join_vert(arma::mat(performance_size - INPUT_SIZE, INPUT_SIZE, arma::fill::zeros), I_in);

    if(performance_size > INPUT_SIZE){
        if(!arma::approx_equal(llft_.D12(), expected_D12, "absdiff", .0)){

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

            arma::mat U1 = U*arma::inv(T);
            arma::mat U1_t = U1.t();
            arma::mat R1 = V.t();
            arma::mat R1_inv = arma::inv(R1);

            llft_.B2() = llft_.B2()*R1_inv;
            llft_.C1() = U1_t*llft_.C1();
            llft_.D11() = U1_t*llft_.D11();
            llft_.D12() = U1_t*llft_.D12()*R1_inv;
            llft_.D22() = llft_.D22()*R1_inv;
        }
    }else if(performance_size == INPUT_SIZE){
        if(arma::approx_equal(llft_.D12(), I_in, "absdiff", .0)){
            std::cout << "[HInf] Special cases trigerred !!!" << std::endl;
        }else{
            //-- do normalization here
        }
    }else{ //-- if not full column rank
        ~ok;
    }

    bool ok2(true); //-- RVO
    arma::mat I_out(OUTPUT_SIZE, OUTPUT_SIZE, arma::fill::eye);
    arma::mat expected_D21 = arma::join_horiz(arma::mat(OUTPUT_SIZE, perturbation_size - OUTPUT_SIZE, arma::fill::zeros), I_out);

    if(perturbation_size > OUTPUT_SIZE){
        if(arma::approx_equal(llft_.D12(), expected_D21, "absdiff", .0)){

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

            arma::mat U2 = arma::inv(T)*V.t();
            arma::mat U2_t = U2.t();
            arma::mat R2 = U;
            arma::mat R2_inv = arma::inv(R2);

            llft_.B1() = llft_.B1()*U2_t;
            llft_.D11() = llft_.D11()*U2_t;
            llft_.C2() = R2_inv*llft_.C2();
            llft_.D21() = R2_inv*llft_.D21()*U2_t;
            llft_.D22() = R2_inv*llft_.D22();

            //-- TODO : the matrix that isn't unitary ?
//            arma::mat R = arma::inv(T) * V.t();
//            U.print("U : ");
//            V.print("V : ");
//            R.print("R : ");
//            arma::mat test = R.t() * R;test.print("TEST : ");
        }
    }else if(perturbation_size == OUTPUT_SIZE){
        if(arma::approx_equal(llft_.D12(), I_out, "absdiff", .0)){
            std::cout << "[HInf] Special cases trigerred !!!" << std::endl;
            //-- do special cases here
        }else{
            //-- do normalization here
        }
    }else{ //-- if not full row rank
        ~ok2;
    }

    std::cout << "Assumption 2 : " << std::endl;
    std::cout << std::boolalpha << ok << " ; " << ok2 << std::endl;

    return ok & ok2;
}

/*
 * Based on Lemma 12.6 in Essential Robust Control
 * This condition is equivalent with ((I-D*inv(R)*conj(D))*C, A-B*inv(R)*conj(D)*C) has no
 * unobservable modes on the jw axis
 */
template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkAssumption3(){   

    arma::mat D12_t( llft_.D12().t() );
    arma::mat R( D12_t*llft_.D12() );
    arma::mat R_inv( arma::inv(R) );
    arma::mat I(performance_size, performance_size, arma::fill::eye);

    arma::mat temp1, temp2;
    temp1 = llft_.D12()*R_inv*D12_t;
    temp2 = I - temp1;
    arma::mat C = temp2*llft_.C1();

    temp1 = llft_.B2()*R_inv*D12_t;
    temp2 = temp1*llft_.C1();
    arma::mat A = llft_.A() - temp2;    

    bool ok = !common::hasUnobservableModeInImAxis(A, C);

    std::cout << "Assumption 3 : " << std::endl;
    std::cout << std::boolalpha << ok << std::endl;

    return !common::hasUnobservableModeInImAxis(A, C);
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkAssumption4(){

    arma::mat D21_t( llft_.D21().t() );
    arma::mat R( llft_.D21()*D21_t );
    arma::mat R_inv( arma::inv(R) );
    arma::mat I(perturbation_size, perturbation_size, arma::fill::eye);

    arma::mat temp1, temp2;
    temp1 = llft_.B1()*D21_t*R_inv;
    temp2 = temp1*llft_.C2();
    arma::mat A = llft_.A() - temp2;

    temp1 = D21_t*R_inv*llft_.D21();
    temp2 = (I - temp1);
    arma::mat B = llft_.B1()*temp2;

    bool ok = !common::hasUncontrollableModeInImAxis(A, B);

    std::cout << "Assumption 4 : " << std::endl;
    std::cout << std::boolalpha << ok << std::endl;

    return ok;
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool HInf<_StateSpace,
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
void HInf<_StateSpace,
    performance_size,
    perturbation_size>::solve(){

    if(checkAllAssumption()){
        arma::mat temp1, temp2, temp3, temp4;

        arma::mat A_t = ss_->A().t();
        arma::mat B_t = ss_->B().t();
        arma::mat B1_t = llft_.B1().t();
        arma::mat C_t = ss_->C().t();
        arma::mat C1_t = llft_.C1().t();

        arma::mat D1_ = ss_->D().head_rows(performance_size);
        arma::mat D1__t = D1_.t();
        arma::mat D_1 = ss_->D().head_cols(perturbation_size);
        arma::mat D_1_t = D_1.t();

        temp1 = gam_*gam_*arma::eye(performance_size, performance_size);
        temp2 = arma::zeros<arma::mat>(D1_.n_cols,  D1_.n_cols);
        temp2.submat(0, 0, performance_size - 1, performance_size - 1) = temp1;
        arma::mat R1 = D1__t * D1_ - temp2;
        arma::mat R1_inv = arma::inv(R1);

        temp1 = gam_*gam_*arma::eye(perturbation_size, perturbation_size);
        temp2 = arma::zeros<arma::mat>(D_1.n_rows, D_1.n_rows);
        temp2.submat(0, 0, perturbation_size - 1, perturbation_size - 1) = temp1;
        arma::mat R2 = D_1 * D_1_t - temp2;
        arma::mat R2_inv = arma::inv(R2);

        //-- Ricatti Domain

        temp1 = arma::join_vert(
                    arma::join_horiz(ss_->A(), arma::zeros<arma::mat>(arma::size(ss_->A()))),
                    arma::join_horiz(-C1_t*llft_.C1(), -A_t)
                    );
        temp2 = arma::join_vert(
                    ss_->B(),
                    -C1_t*D1_
                    );
        temp3 = arma::join_horiz(D1__t*llft_.C1(), B_t);
        temp4 = temp2*R1_inv*temp3;
        arma::mat H_inf = temp1 - temp4;

        temp1 = arma::join_vert(
                    arma::join_horiz(A_t, arma::zeros<arma::mat>(arma::size(ss_->A()))),
                    arma::join_horiz(-llft_.B1()*B1_t, -ss_->A())
                    );
        temp2 = arma::join_vert(
                    C_t,
                    -llft_.B1()*D_1_t
                    );
        temp3 = arma::join_horiz(
                    D_1*B1_t,
                    ss_->C()
                    );
        temp4 = temp2*R2_inv*temp3;
        arma::mat J_inf = temp1 - temp4;

    }
}

}

}
