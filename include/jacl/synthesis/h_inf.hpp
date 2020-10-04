/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of H-infinity synthesis.
*            A lot of code here, is constructed by copy constructor because it's easier to read.
*/

#pragma once

#include <jacl/lti_common.hpp>
#include <jacl/are.hpp>
#include <jacl/lft/lower.hpp>
#include <jacl/state_space/linear.hpp>
#include <jacl/linear_algebra.hpp>
#include <jacl/traits.hpp>
#include <jacl/numerical_methods.hpp>
#include <jacl/medium.hpp>

#define HINF_VERBOSE

namespace jacl{ namespace synthesis{

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
class Hinf:public ::jacl::system::detail::BaseSystemClient<typename _System::base_t>{
public:
    template<typename __System = _System,
             typename std::enable_if_t<::jacl::traits::is_continuous_system_v<__System>, int>* = nullptr>
    Hinf(__System* _sys, double _gam)
        : ss_(::jacl::system::detail::BaseSystemClient<typename _System::base_t>::ss(_sys))
        , llft_(ss_)
        , gam_(_gam){}

    ~Hinf();

    auto solve() -> ::jacl::lti_common::StateSpacePack;

private:
    auto checkAssumption1();
    auto checkAssumption2();
    auto checkAssumption3();
    auto checkAssumption4();
    auto checkAllAssumption();

    auto checkCondition1();
    auto checkCondition2();
    auto checkCondition3();
    auto checkCondition4();
    auto checkAllCondition();

private:
    typename _System::state_space_t* ss_;
    static constexpr auto INPUT_SIZE{_System::n_inputs - perturbation_size};
    static constexpr auto OUTPUT_SIZE{_System::n_outputs - performance_size};
    lft::LowerLFT<typename _System::state_space_t,
             performance_size,
             perturbation_size,
             OUTPUT_SIZE,
             INPUT_SIZE> llft_;

    //-- normalization factors
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<INPUT_SIZE, INPUT_SIZE> R1_, R1_inv_;
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<OUTPUT_SIZE, OUTPUT_SIZE> R2_, R2_inv_;

    double gam_;
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<_System::n_states, _System::n_states> X_inf_;
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<_System::n_states, _System::n_states> Y_inf_;
};

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
Hinf<_System,
    performance_size,
    perturbation_size>::~Hinf(){
    ss_ = nullptr;
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto Hinf<_System,
    performance_size,
    perturbation_size>::checkAssumption1(){

    auto ctrb = lti_common::stabilizable(llft_.A(), llft_.B2());
    auto obsv = lti_common::detectability(llft_.A(), llft_.C2());
#ifdef HINF_VERBOSE
    std::cout << "[Hinf] Assumption 1 : " << std::boolalpha << ctrb << " ; " << obsv << std::endl;
#endif
    return ctrb & obsv;
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto Hinf<_System,
    performance_size,
    perturbation_size>::checkAssumption2(){

    auto ok(true); //-- RVO
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<INPUT_SIZE, INPUT_SIZE> I_in(arma::fill::eye);
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<performance_size, INPUT_SIZE> expected_D12 = arma::join_vert(
            arma::mat(performance_size - INPUT_SIZE, INPUT_SIZE, arma::fill::zeros),
            I_in
    );

    if(performance_size > INPUT_SIZE){
        if(!arma::approx_equal(llft_.D12(), expected_D12, "absdiff", .0)){
            R1_ = R1_inv_ = arma::eye(INPUT_SIZE, INPUT_SIZE);
        }else{
            //-- do normalization here using SVD

            typename arma::Mat<typename _System::scalar_t>::template
                fixed<performance_size, performance_size> U;
            typename arma::Mat<typename _System::scalar_t>::template
                fixed<INPUT_SIZE, INPUT_SIZE> V;
            typename arma::Col<typename _System::scalar_t>::template
                fixed<INPUT_SIZE> s;
            arma::svd(U,s,V,llft_.D12());

            typename arma::Mat<typename _System::scalar_t>::template
                fixed<INPUT_SIZE, INPUT_SIZE> S = arma::diagmat(s);
            // decltype(S) recip_S(arma::fill::zeros);
            // for(int i(0); i < S.n_rows; i++){
            //     recip_S(i, i) = 1./S(i,i);
            // }

            decltype(U) U_flip( arma::fliplr( U ) );
            decltype(S) S_flip( arma::flipud( S ) );

            decltype(U) U1 = U_flip;
            decltype(U1) U1_t = arma::trans( U1 );
            R1_ = S_flip * arma::trans( V );
            R1_inv_ = arma::inv( R1_ );

            llft_.B2() = llft_.B2()*R1_inv_;
            llft_.C1() = U1_t*llft_.C1();
            llft_.D11() = U1_t*llft_.D11();
            llft_.D12() = U1_t*llft_.D12()*R1_inv_;
            llft_.D22() = llft_.D22()*R1_inv_;
        }
    }else if(performance_size == INPUT_SIZE){
        if(arma::approx_equal(llft_.D12(), I_in, "absdiff", .0)){
            std::cout << "[Hinf] Special cases trigerred !!!" << std::endl;
        }else{
            //-- do normalization here
        }
    }else{ //-- if not full column rank
        ~ok;
    }

    auto ok2(true); //-- RVO
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<OUTPUT_SIZE, OUTPUT_SIZE> I_out(arma::fill::eye);
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<OUTPUT_SIZE, perturbation_size> expected_D21 = arma::join_horiz(
            arma::mat(OUTPUT_SIZE, perturbation_size - OUTPUT_SIZE, arma::fill::zeros),
            I_out
    );

    if(perturbation_size > OUTPUT_SIZE){
        if(arma::approx_equal(llft_.D12(), expected_D21, "absdiff", .0)){
            R2_ = R2_inv_ = arma::eye(OUTPUT_SIZE, OUTPUT_SIZE);
        }else{
            //-- do normalization here using SVD

            typename arma::Mat<typename _System::scalar_t>::template
                fixed<OUTPUT_SIZE, OUTPUT_SIZE> U;
            typename arma::Mat<typename _System::scalar_t>::template
                fixed<perturbation_size, perturbation_size> V;
            typename arma::Col<typename _System::scalar_t>::template
                fixed<OUTPUT_SIZE> s;
            arma::svd(U,s,V,llft_.D21());

            typename arma::Mat<typename _System::scalar_t>::template
                fixed<OUTPUT_SIZE, OUTPUT_SIZE> S = arma::diagmat(s);
            // arma::mat recip_S(arma::size(S), arma::fill::zeros);
            // for(int i(0); i < S.n_rows; i++){
            //     recip_S(i, i) = 1./S(i,i);
            // }

            decltype(V) V_flip( arma::flipud( arma::trans(V) ) );
            decltype(S) S_flip( arma::fliplr(S) );

            R2_ = U * S_flip;
            R2_inv_ = arma::inv( R2_ );
            decltype(V) U2 = V_flip;
            decltype(U2) U2_t = arma::trans( U2 );

            llft_.B1() = llft_.B1()*U2_t;
            llft_.D11() = llft_.D11()*U2_t;
            llft_.C2() = R2_inv_*llft_.C2();
            llft_.D21() = R2_inv_*llft_.D21()*U2_t;
            llft_.D22() = R2_inv_*llft_.D22();
        }
    }else if(perturbation_size == OUTPUT_SIZE){
        if(arma::approx_equal(llft_.D12(), I_out, "absdiff", .0)){
            std::cout << "[Hinf] Special cases trigerred !!!" << std::endl;
            //-- do special cases here
        }else{
            //-- do normalization here
        }
    }else{ //-- if not full row rank
        ~ok2;
    }
#ifdef HINF_VERBOSE
    std::cout << "[Hinf] Assumption 2 : " << std::boolalpha << ok << " ; " << ok2 << std::endl;
#endif
    return ok & ok2;
}

/*
 * Based on Lemma 12.6 in Essential Robust Control
 * This condition is equivalent with ((I-D*inv(R)*conj(D))*C, A-B*inv(R)*conj(D)*C) has no
 * unobservable modes on the jw axis
 */
template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto Hinf<_System,
    performance_size,
    perturbation_size>::checkAssumption3(){

    typename arma::Mat<typename _System::scalar_t>::template
        fixed<INPUT_SIZE, performance_size> D12_t( llft_.D12().t() );
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<INPUT_SIZE, INPUT_SIZE> R( D12_t*llft_.D12() );
    decltype(R) R_inv( arma::inv(R) );
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<performance_size, performance_size> I(arma::fill::eye);

    arma::mat temp1, temp2;
    temp1 = llft_.D12()*R_inv*D12_t;
    temp2 = I - temp1;
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<performance_size, _System::n_states> C = temp2*llft_.C1();

    temp1 = llft_.B2()*R_inv*D12_t;
    temp2 = temp1*llft_.C1();
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<_System::n_states, _System::n_states> A = llft_.A() - temp2;    

    bool ok = !lti_common::hasUnobservableModeInImAxis(A, C);
#ifdef HINF_VERBOSE
    std::cout << "[Hinf] Assumption 3 : " << std::boolalpha << ok << std::endl;
#endif
    return ok;
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto Hinf<_System,
    performance_size,
    perturbation_size>::checkAssumption4(){

    typename arma::Mat<typename _System::scalar_t>::template
        fixed<perturbation_size, OUTPUT_SIZE> D21_t( llft_.D21().t() );
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<OUTPUT_SIZE, OUTPUT_SIZE> R( llft_.D21()*D21_t );
    decltype(R) R_inv( arma::inv(R) );
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<perturbation_size, perturbation_size> I(arma::fill::eye);

    arma::mat temp1, temp2;
    temp1 = llft_.B1()*D21_t*R_inv;
    temp2 = temp1*llft_.C2();
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<_System::n_states, _System::n_states> A = llft_.A() - temp2;

    temp1 = D21_t*R_inv*llft_.D21();
    temp2 = (I - temp1);
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<_System::n_states, perturbation_size> B = llft_.B1()*temp2;

    bool ok = !lti_common::hasUncontrollableModeInImAxis(A, B);
#ifdef HINF_VERBOSE
    std::cout << "[Hinf] Assumption 4 : " << std::boolalpha << ok << std::endl;
#endif
    return ok;
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto Hinf<_System,
    performance_size,
    perturbation_size>::checkAllAssumption(){
    return checkAssumption1()
            & checkAssumption2()
            & checkAssumption3()
            & checkAssumption4();
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto Hinf<_System,
    performance_size,
    perturbation_size>::checkCondition1(){

    typename arma::Mat<typename _System::scalar_t>::template
        fixed<performance_size-INPUT_SIZE, perturbation_size> D11_ =
            llft_.D11().head_rows(performance_size - INPUT_SIZE);
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<performance_size, perturbation_size - OUTPUT_SIZE> D_11 =
            llft_.D11().head_cols(perturbation_size - OUTPUT_SIZE);

    return gam_ > std::max(
        linear_algebra::largestSV(D11_),
        linear_algebra::largestSV(D_11)
    );
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto Hinf<_System,
    performance_size,
    perturbation_size>::checkCondition2(){

    return linear_algebra::isPosSemiDefinite(X_inf_);
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto Hinf<_System,
    performance_size,
    perturbation_size>::checkCondition3(){

    return linear_algebra::isPosSemiDefinite(Y_inf_);
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto Hinf<_System,
    performance_size,
    perturbation_size>::checkCondition4(){

    return linear_algebra::spectralRadius(
        arma::Mat<std::complex<typename _System::scalar_t>>(X_inf_*Y_inf_)) < gam_*gam_;
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto Hinf<_System,
    performance_size,
    perturbation_size>::checkAllCondition(){
    auto cond1 = checkCondition1();
    auto cond2 = checkCondition2();
    auto cond3 = checkCondition3();
    auto cond4 = checkCondition4();
#ifdef HINF_VERBOSE
    std::cout << "[Hinf] Condition 1 : " << std::boolalpha << cond1 << std::endl;
    std::cout << "[Hinf] Condition 2 : " << std::boolalpha << cond2 << std::endl;
    std::cout << "[Hinf] Condition 3 : " << std::boolalpha << cond3 << std::endl;
    std::cout << "[Hinf] Condition 4 : " << std::boolalpha << cond4 << std::endl;
#endif
    return cond1 & cond2 & cond3 & cond4;
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto Hinf<_System,
    performance_size,
    perturbation_size>::solve() -> ::jacl::lti_common::StateSpacePack{
#ifdef HINF_VERBOSE
    std::cout << "[Hinf] Interconnection matrix assumptions : " << std::endl;
#endif
    bool check_assumption = checkAllAssumption();
    assert(check_assumption && "[Hinf] The assumption made for interconnection matrix is not fulfill !");

    llft_.D12().print("[Hinf] Normalized D12 : ");
    llft_.D21().print("[Hinf] Normalized D21 : ");

    using namespace ::jacl::linear_algebra;

    arma::mat temp1, temp2, temp3, temp4;
    arma::cx_mat ctemp1, ctemp2, ctemp3;

    const typename arma::Mat<typename _System::scalar_t>::template
        fixed<_System::n_states, _System::n_states>& A = llft_.A();
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<_System::n_states, _System::n_states> A_t = arma::trans(A);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<_System::n_states, _System::n_states> cxA = toCx(A);
    const typename arma::Mat<typename _System::scalar_t>::template
        fixed<_System::n_states, perturbation_size+INPUT_SIZE>& B = ss_->B();
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<perturbation_size+INPUT_SIZE, _System::n_states> B_t = arma::trans(B);
    const typename arma::Mat<typename _System::scalar_t>::template
        fixed<performance_size+OUTPUT_SIZE, _System::n_states>& C = ss_->C();
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<_System::n_states, performance_size+OUTPUT_SIZE> C_t = arma::trans(C);
    const typename arma::Mat<typename _System::scalar_t>::template
        fixed<performance_size+OUTPUT_SIZE, perturbation_size+INPUT_SIZE>& D = ss_->D();
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<perturbation_size+INPUT_SIZE, performance_size+OUTPUT_SIZE> D_t = arma::trans(D);
    const typename arma::Mat<typename _System::scalar_t>::template
        fixed<_System::n_states, perturbation_size>& B1 = llft_.B1();
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<_System::n_states, perturbation_size> cxB1 = toCx(B1);
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<perturbation_size, _System::n_states> B1_t = arma::trans(B1);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<perturbation_size, _System::n_states> cxB1_t = toCx(B1_t);
    const typename arma::Mat<typename _System::scalar_t>::template
        fixed<_System::n_states, INPUT_SIZE>& B2 = llft_.B2();
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<_System::n_states, INPUT_SIZE> cxB2 = toCx(B2);
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<INPUT_SIZE, _System::n_states> B2_t = arma::trans(B2);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<INPUT_SIZE, _System::n_states> cxB2_t = toCx(B2_t);
    const typename arma::Mat<typename _System::scalar_t>::template
        fixed<performance_size, _System::n_states>& C1 = llft_.C1();
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<performance_size, _System::n_states> cxC1 = toCx(C1);
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<_System::n_states, performance_size> C1_t = arma::trans(C1);
    const typename arma::Mat<typename _System::scalar_t>::template
        fixed<OUTPUT_SIZE, _System::n_states>& C2 = llft_.C2();
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<OUTPUT_SIZE, _System::n_states> cxC2 = toCx(C2);
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<_System::n_states, OUTPUT_SIZE> C2_t = arma::trans(C2);
    const typename arma::Mat<typename _System::scalar_t>::template
        fixed<performance_size, perturbation_size>& D11 = llft_.D11();
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<performance_size, perturbation_size> cxD11 = toCx(D11);
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<perturbation_size, performance_size> D11_t = arma::trans(D11);
    const typename arma::Mat<typename _System::scalar_t>::template
        fixed<performance_size, INPUT_SIZE>& D12 = llft_.D12();
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<performance_size, INPUT_SIZE> cxD12 = toCx(D12);
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<INPUT_SIZE, performance_size> D12_t = arma::trans(D12);
    const typename arma::Mat<typename _System::scalar_t>::template
        fixed<OUTPUT_SIZE, perturbation_size>& D21 = llft_.D21();
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<OUTPUT_SIZE, perturbation_size> cxD21 = toCx(D21);
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<perturbation_size, OUTPUT_SIZE> D21_t = arma::trans(D21);
    const typename arma::Mat<typename _System::scalar_t>::template
        fixed<OUTPUT_SIZE, INPUT_SIZE>& D22 = llft_.D22();
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<OUTPUT_SIZE, INPUT_SIZE> cxD22 = toCx(D22);
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<INPUT_SIZE, OUTPUT_SIZE> D22_t = arma::trans(D22);

    typename arma::Mat<typename _System::scalar_t>::template
        fixed<performance_size, perturbation_size+INPUT_SIZE> D1_ = ss_->D().head_rows(performance_size);
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<perturbation_size+INPUT_SIZE, performance_size> D1__t = D1_.t();
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<performance_size+OUTPUT_SIZE, perturbation_size> D_1 = ss_->D().head_cols(perturbation_size);
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<perturbation_size, performance_size+OUTPUT_SIZE> D_1_t = D_1.t();

    temp1 = (gam_*gam_)*arma::eye(perturbation_size, perturbation_size);
    temp2 = arma::zeros<arma::mat>(D1_.n_cols,  D1_.n_cols);
    temp2.submat(0, 0, perturbation_size-1, perturbation_size-1) = temp1;
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<perturbation_size+INPUT_SIZE, perturbation_size+INPUT_SIZE> R1 = (D1__t * D1_) - temp2;
//    R1.print("R1 : ");
    //-- create exception here for singular R1 and warn to change the gamma
    decltype(R1) R1_inv = arma::inv(R1);

    temp1 = gam_*gam_*arma::eye(performance_size, performance_size);
    temp2 = arma::zeros<arma::mat>(D_1.n_rows, D_1.n_rows);
    temp2.submat(0, 0, performance_size-1, performance_size-1) = temp1;
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<performance_size+OUTPUT_SIZE, performance_size+OUTPUT_SIZE> R2 = (D_1 * D_1_t) - temp2;
//    R2.print("R2 : ");
    //-- create exception here for singular R2 and warn to change the gamma
    decltype(R2) R2_inv = arma::inv(R2);

    //-- Ricatti Domain

    temp1 = arma::join_vert(
                //-- create a function to simplify the zeros
                arma::join_horiz(A, arma::zeros<arma::mat>(arma::size(A))),
                arma::join_horiz(-C1_t*C1, -A_t)
            );
    temp2 = arma::join_vert(
                B,
                -C1_t*D1_
            );
    temp3 = arma::join_horiz(D1__t*C1, B_t);
    temp4 = temp2*R1_inv*temp3;
    typename arma::Mat<typename _System::scalar_t>::template
        fixed<_System::n_states<<1, _System::n_states<<1> H_inf = temp1 - temp4;

    temp1 = arma::join_vert(
                arma::join_horiz(A_t, arma::zeros<arma::mat>(arma::size(A))),
                arma::join_horiz(-B1*B1_t, -A)
                );
    temp2 = arma::join_vert(
                C_t,
                -B1*D_1_t
                );
    temp3 = arma::join_horiz(
                D_1*B1_t,
                C
                );
    temp4 = temp2*R2_inv*temp3;
    decltype(H_inf) J_inf = temp1 - temp4;

    ARE<typename _System::state_space_t> solver1(ss_);
    ARE<typename _System::state_space_t> solver2(ss_);

#ifdef HINF_VERBOSE
    //-- check whether Hamiltonian or not
    {
        decltype(H_inf) J = arma::join_cols(
            arma::join_rows(arma::zeros(arma::size(A)), -arma::eye(arma::size(A))),
            arma::join_rows(arma::eye(arma::size(A)), arma::zeros(arma::size(A)))
        );

        decltype(H_inf) check_H_inf = -J*H_inf*J;
        check_H_inf = -arma::trans(check_H_inf);
        H_inf.print("[Hinf] H_inf : ");
        check_H_inf.print("[Hinf] Check H_inf : ");

        decltype(H_inf) check_J_inf = -J*J_inf*J;
        check_J_inf = -arma::trans(check_J_inf);
        J_inf.print("[Hinf] J_inf : ");
        check_J_inf.print("[Hinf] Check J_inf : ");
    }
#endif

    solver1.setHamiltonianMatrix(H_inf);
    solver2.setHamiltonianMatrix(J_inf);

    X_inf_ = solver1.solve();
    Y_inf_ = solver2.solve();

#ifdef HINF_VERBOSE
    //--check the ARE solution
    {
        decltype(X_inf_) H11 = toCx(H_inf.submat(0,0,_System::n_states-1,_System::n_states-1));
        decltype(X_inf_) H12 = toCx(H_inf.submat(0, _System::n_states, _System::n_states-1, 2*_System::n_states-1));
        decltype(X_inf_) H21 = toCx(H_inf.submat(_System::n_states, 0, 2*_System::n_states-1, _System::n_states-1));
        decltype(X_inf_) H22 = toCx(H_inf.submat(_System::n_states, _System::n_states,
                                            2*_System::n_states-1, 2*_System::n_states-1));
        ctemp1 = -H22*X_inf_ + X_inf_*H11;
        ctemp2 = X_inf_*H12*X_inf_ - H21;
        ctemp3 = ctemp1 + ctemp2;
        ctemp3.print("[Hinf] ARE 1 rhs : ");

        decltype(Y_inf_) J11 = toCx(J_inf.submat(0,0,_System::n_states-1,_System::n_states-1));
        decltype(Y_inf_) J12 = toCx(J_inf.submat(0, _System::n_states, _System::n_states-1, 2*_System::n_states-1));
        decltype(Y_inf_) J21 = toCx(J_inf.submat(_System::n_states, 0, 2*_System::n_states-1, _System::n_states-1));
        decltype(Y_inf_) J22 = toCx(J_inf.submat(_System::n_states, _System::n_states,
                                            2*_System::n_states-1, 2*_System::n_states-1));
        ctemp1 = -J22*Y_inf_ + Y_inf_*J11;
        ctemp2 = Y_inf_*J12*Y_inf_ - J21;
        ctemp3 = ctemp1 + ctemp2;
        ctemp3.print("[Hinf] ARE 2 rhs : ");
    }
#endif

#ifdef HINF_VERBOSE
    std::cout << "[Hinf] Solution of ARE : " << std::endl;
    X_inf_.print("[Hinf] X_inf : ");
    Y_inf_.print("[Hinf] Y_inf : ");
#endif

    ctemp1 = toCx(D1__t * C1);
    ctemp2 = toCx(B_t) * X_inf_;
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<perturbation_size+INPUT_SIZE, _System::n_states> F = -toCx(R1_inv)*(ctemp1 + ctemp2);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<perturbation_size, _System::n_states> F1_inf = F.head_rows(perturbation_size);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<INPUT_SIZE, _System::n_states> F2_inf = F.tail_rows(INPUT_SIZE);

    ctemp1 = toCx(B1 * D_1_t);
    ctemp2 = Y_inf_ * toCx(C_t);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<_System::n_states, performance_size+OUTPUT_SIZE> L = -(ctemp1 + ctemp2)*toCx(R2_inv);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<_System::n_states, performance_size> L1_inf = L.head_cols(performance_size);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<_System::n_states, OUTPUT_SIZE> L2_inf = L.tail_cols(OUTPUT_SIZE);

    //-- Partition of D, F1_inf and L1_inf
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<perturbation_size-OUTPUT_SIZE, _System::n_states> F11_inf = F1_inf.head_rows(perturbation_size - OUTPUT_SIZE);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<OUTPUT_SIZE, _System::n_states> F12_inf = F1_inf.tail_rows(OUTPUT_SIZE);

    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<_System::n_states, performance_size-INPUT_SIZE> L11_inf = L1_inf.head_cols(performance_size - INPUT_SIZE);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<_System::n_states, INPUT_SIZE> L12_inf = L1_inf.tail_cols(INPUT_SIZE);

    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<performance_size-INPUT_SIZE, perturbation_size-OUTPUT_SIZE> D1111 = toCx(D11).submat(0, 0,
                                          (performance_size - INPUT_SIZE) - 1, (perturbation_size - OUTPUT_SIZE) - 1);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<performance_size-INPUT_SIZE, OUTPUT_SIZE> D1112 = toCx(D11).submat(0, (perturbation_size - OUTPUT_SIZE),
                                          (performance_size - INPUT_SIZE) - 1, perturbation_size - 1);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<INPUT_SIZE, perturbation_size-OUTPUT_SIZE> D1121 = toCx(D11).submat((performance_size - INPUT_SIZE), 0,
                                            performance_size - 1, (perturbation_size - OUTPUT_SIZE) - 1);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<INPUT_SIZE, OUTPUT_SIZE> D1122 = toCx(D11).submat((performance_size - INPUT_SIZE), (perturbation_size - OUTPUT_SIZE),
                                            performance_size - 1,  perturbation_size - 1);

    bool check_cond = checkAllCondition();
#ifdef HINF_VERBOSE
    std::cout << "[Hinf] Condition : " << std::boolalpha << check_cond << std::endl;
#endif
    assert(check_cond && "[Hinf] Condition for existence controller that make the objective is less than gamma was failed !");

    //-- I assume the arbitrary system that have property less than to gamma is zero a.k.a central controller
    ctemp1 = (1./(gam_*gam_))*Y_inf_*X_inf_;
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<_System::n_states, _System::n_states> Z_inf = arma::inv(arma::eye<arma::cx_mat>(_System::n_states, _System::n_states) - ctemp1);    
    Z_inf.print("[Hinf] Z_inf : ");

    ctemp1 = -D1121*arma::trans(D1111);
    ctemp2 = (gam_*gam_)*arma::eye<arma::cx_mat>(performance_size - INPUT_SIZE, performance_size - INPUT_SIZE)
            - D1111*arma::trans(D1111);
    ctemp3 = ctemp1*arma::inv(ctemp2);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<INPUT_SIZE, OUTPUT_SIZE> D11_hat = ctemp3*D1112 - D1122;

    //-- calculate D12_hat and D21_hat
    ctemp1 = (gam_*gam_)*arma::eye(perturbation_size - OUTPUT_SIZE, perturbation_size - OUTPUT_SIZE)
            - arma::trans(D1111)*D1111;
    ctemp2 = arma::eye(INPUT_SIZE, INPUT_SIZE) - D1121*arma::inv(ctemp1)*arma::trans(D1121);
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<INPUT_SIZE, INPUT_SIZE> D12_hat = arma::chol(ctemp2, "lower");

    ctemp1 = (gam_*gam_)*arma::eye(performance_size - INPUT_SIZE, performance_size - INPUT_SIZE)
            - D1111*arma::trans(D1111);
    ctemp2 = arma::eye(OUTPUT_SIZE, OUTPUT_SIZE) - arma::trans(D1112)*arma::inv(ctemp1)*D1112;
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<OUTPUT_SIZE, OUTPUT_SIZE> D21_hat = arma::chol(ctemp2, "upper");

    // ctemp1 = B2 + L12_inf;
    // arma::cx_mat B2_hat = Z_inf*ctemp1*D12_hat;
    // arma::cx_mat C2_hat = -D21_hat*(C2 + F12_inf);
    // ctemp1 = -Z_inf*L2_inf;
    // ctemp2 = B2_hat*arma::inv(D12_hat)*D11_hat;
    // arma::cx_mat B1_hat = ctemp1 + ctemp2;
    // ctemp1 = D11_hat*arma::inv(D21_hat)*C2_hat;
    // arma::cx_mat C1_hat = F2_inf + ctemp1;
    // ctemp1 = B1_hat*arma::inv(D21_hat)*C2_hat;
    // arma::cx_mat A_hat = A + B*F + ctemp1;

    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<_System::n_states, INPUT_SIZE> B2_hat = (B2 + L12_inf)*D12_hat;
    ctemp1 = C2 + F12_inf;
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<OUTPUT_SIZE, _System::n_states> C2_hat = -D21_hat*ctemp1*Z_inf;
    ctemp1 = B2_hat*arma::inv(D12_hat)*D11_hat;
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<_System::n_states, OUTPUT_SIZE> B1_hat = -L2_inf + ctemp1;
    ctemp1 = F2_inf*Z_inf;
    ctemp2 = D11_hat*arma::inv(D21_hat)*C2_hat;
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<INPUT_SIZE, _System::n_states> C1_hat = ctemp1 + ctemp2;
    ctemp1 = B2_hat*arma::inv(D12_hat)*C1_hat;
    typename arma::Mat<std::complex<typename _System::scalar_t>>::template
        fixed<_System::n_states, _System::n_states> A_hat = A + L*C + ctemp1;

    //-- recover from normalization matrix
    B1_hat = B1_hat*toCx(R2_inv_);
    C1_hat = toCx(R1_inv_)*C1_hat;
    D11_hat = toCx(R1_inv_)*D11_hat*toCx(R2_inv_);
    D12_hat = toCx(R1_inv_)*D12_hat;
    D21_hat = D21_hat*toCx(R2_inv_);

#ifdef HINF_VERBOSE
    std::cout << "[Hinf] Controller Result : " << std::endl;
    A_hat.print("A_hat : ");
    B1_hat.print("B1_hat : ");
    B2_hat.print("B2_hat : ");
    C1_hat.print("C1_hat : ");
    C2_hat.print("C2_hat : ");
    D11_hat.print("D11_hat : ");
    D12_hat.print("D12_hat : ");
    D21_hat.print("D21_hat : ");
#endif

    return std::make_tuple(toReal(A_hat), toReal(B1_hat), toReal(C1_hat), toReal(D11_hat));
}

} } // namespace jacl::synthesis
