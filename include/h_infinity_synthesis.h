/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of H-infinity synthesis.
*            A lot of code here, is constructed by copy constructor because it's easier to read.
*/

#pragma once

#include "lti_common.h"
#include "are.h"
#include "upper_lft.h"
#include "lower_lft.h"
#include "state_space.h"
#include "traits.h"
#include "numerical_methods.h"

#define HINF_DEBUG

namespace jacl{

namespace synthesis{

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
class HInf{
public:
    template<typename T = _StateSpace,
             typename std::enable_if<traits::is_state_space<T>::value>::type* = nullptr>
    HInf(T* _ss, double _gam)
        : ss_(_ss)
        , llft_(ss_)
        , are_solver1_(ss_)
        , are_solver2_(ss_)
        , gam_(_gam){

    }    

    ~HInf();

    auto init() -> void;
    auto solve() -> common::StateSpacePack;

    template <typename SS1, typename SS2>
    auto starProducts(const SS1 &_ss1, const SS2 &_ss2) -> void{
        arma::cx_mat A_bar(SS1::n_states + SS2::n_states, SS1::n_states + SS2::n_states);
        arma::cx_mat B_bar(SS1::n_states + SS2::n_states, SS1::w1_size + SS2::u_sz);
        arma::cx_mat C_bar(SS1::z1_sz + SS2::y_sz, SS1::n_states + SS2::n_states);
        arma::cx_mat D_bar(SS1::z1_sz + SS2::y_sz, SS1::w1_size, + SS2::u_sz);
    }

private:
    auto checkAssumption1() -> bool;
    auto checkAssumption2() -> bool;
    auto checkAssumption3() -> bool;
    auto checkAssumption4() -> bool;
    auto checkAllAssumption() -> bool;

    auto checkCondition1() -> bool;
    auto checkCondition2() -> bool;
    auto checkCondition3() -> bool;
    auto checkCondition4() -> bool;
    auto checkAllCondition() -> bool;

    inline auto toCx(const arma::mat& _in) const -> arma::cx_mat{
        return arma::cx_mat(_in, arma::zeros<arma::mat>( arma::size(_in) ));
    }

    inline auto toReal(const arma::cx_mat& _in) const -> arma::mat{
        return arma::real(_in);
    }

    template <typename __StateSpace>
    auto isInfNormLessThan(double _gam, const __StateSpace& _ss) -> bool{
        arma::mat C_t = arma::trans( _ss.C() );
        arma::mat D_t = arma::trans( _ss.D() );
        arma::mat D_tD = D_t * _ss.D();
        arma::mat R = (_gam*_gam)*arma::eye(__StateSpace::n_inputs,
                                            __StateSpace::n_inputs)
                        - D_tD;
        arma::mat R_inv = arma::inv(R);

        //-- Hamiltonian matrix
        arma::mat H(__StateSpace::n_states << 1,
                    __StateSpace::n_states << 1);

        arma::mat temp1, temp2;

        //-- block 1,1
        temp1 = R_inv*D_t*_ss.C();
        temp2 = _ss.A() + temp1;
        H.submat(                         0,                          0,
                 __StateSpace::n_states - 1, __StateSpace::n_states - 1) = temp2;
        //-- block 2,2
        H.submat(__StateSpace::n_states, __StateSpace::n_states, __StateSpace::n_states * 2 - 1, __StateSpace::n_states * 2 - 1) = -arma::trans(temp2);
        //-- block 1,2
        H.submat(                         0,         __StateSpace::n_states,
                 __StateSpace::n_states - 1, __StateSpace::n_states * 2 - 1) = _ss.B()*R_inv*arma::trans(_ss.B());
        //-- block 2,1
        temp1 = _ss.D()*R_inv*D_t;
        temp2 = arma::eye(__StateSpace::n_outputs,
                          __StateSpace::n_outputs)
                 - temp1;
        H.submat(__StateSpace::n_states, 0, __StateSpace::n_states * 2 - 1, __StateSpace::n_states - 1) = -C_t*temp2*_ss.C();

        arma::cx_mat eigvec;
        arma::cx_vec eigval;
        arma::eig_gen(eigval, eigvec, H);
        arma::mat eigval_re = arma::real(eigval);
        auto ok(true);
        for(int i(0); i < eigval.n_rows; i++){
            if(std::fabs(eigval_re(i, 0)) < std::numeric_limits<double>::epsilon()){
                ~ok;
                break;
            }
        }
        return ok;
    }

private:
    _StateSpace* ss_;

    static constexpr auto INPUT_SIZE{_StateSpace::n_inputs - perturbation_size};
    static constexpr auto OUTPUT_SIZE{_StateSpace::n_outputs - performance_size};

    LowerLFT<_StateSpace,
             performance_size,
             perturbation_size,
             OUTPUT_SIZE,
             INPUT_SIZE> llft_;

    ARE<_StateSpace> are_solver1_;
    ARE<_StateSpace> are_solver2_;    

    double gam_;

    arma::cx_mat X_inf_;
    arma::cx_mat Y_inf_;

};

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
HInf<_StateSpace,
    performance_size,
    perturbation_size>::~HInf(){
    ss_ = nullptr;
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto HInf<_StateSpace,
    performance_size,
    perturbation_size>::init() -> void{

}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkAssumption1() -> bool{

    auto ctrb = common::stabilizable(llft_.A(), llft_.B2());
    auto obsv = common::detectability(llft_.A(), llft_.C2());
#ifdef HINF_DEBUG
    std::cout << "Assumption 1 : " << std::endl;
    std::cout << std::boolalpha << ctrb << " ; " << obsv << std::endl;
#endif
    return ctrb & obsv;
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkAssumption2() -> bool{

    auto ok(true); //-- RVO
    arma::mat I_in(INPUT_SIZE, INPUT_SIZE, arma::fill::eye);
    arma::mat expected_D12 = arma::join_vert(arma::mat(performance_size - INPUT_SIZE, INPUT_SIZE, arma::fill::zeros), I_in);

    if(performance_size > INPUT_SIZE){
        if(!arma::approx_equal(llft_.D12(), expected_D12, "absdiff", .0)){

        }else{
            //-- do normalization here using SVD

            arma::mat U, V;
            arma::vec s;
            arma::svd(U,s,V,llft_.D12());

            arma::mat S = arma::diagmat(s);
            arma::mat recip_S(arma::size(S), arma::fill::zeros);
            for(int i(0); i < S.n_rows; i++){
                recip_S(i, i) = 1./S(i,i);
            }

            arma::mat U_flip( arma::fliplr( U ) ), S_flip( arma::flipud( S ) );

            arma::mat U1 = U_flip;
            arma::mat U1_t = arma::trans( U1 );
            arma::mat R1 = S_flip * arma::trans( V );
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

    auto ok2(true); //-- RVO
    arma::mat I_out(OUTPUT_SIZE, OUTPUT_SIZE, arma::fill::eye);
    arma::mat expected_D21 = arma::join_horiz(arma::mat(OUTPUT_SIZE, perturbation_size - OUTPUT_SIZE, arma::fill::zeros), I_out);

    if(perturbation_size > OUTPUT_SIZE){
        if(arma::approx_equal(llft_.D12(), expected_D21, "absdiff", .0)){

        }else{
            //-- do normalization here using SVD

            arma::mat U, V;
            arma::vec s;
            arma::svd(U,s,V,llft_.D21());

            arma::mat S = arma::diagmat(s);
            arma::mat recip_S(arma::size(S), arma::fill::zeros);
            for(int i(0); i < S.n_rows; i++){
                recip_S(i, i) = 1./S(i,i);
            }

            arma::mat V_flip( arma::flipud( arma::trans(V) ) ), S_flip( arma::fliplr(S) );

            arma::mat R2 = U * S_flip;
            arma::mat R2_inv = arma::inv( R2 );
            arma::mat U2 = V_flip;
            arma::mat U2_t = arma::trans( U2 );

            llft_.B1() = llft_.B1()*U2_t;
            llft_.D11() = llft_.D11()*U2_t;
            llft_.C2() = R2_inv*llft_.C2();
            llft_.D21() = R2_inv*llft_.D21()*U2_t;
            llft_.D22() = R2_inv*llft_.D22();
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
#ifdef HINF_DEBUG
    std::cout << "Assumption 2 : " << std::endl;
    std::cout << std::boolalpha << ok << " ; " << ok2 << std::endl;
#endif
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
auto HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkAssumption3() -> bool{

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
#ifdef HINF_DEBUG
    std::cout << "Assumption 3 : " << std::endl;
    std::cout << std::boolalpha << ok << std::endl;
#endif
    return !common::hasUnobservableModeInImAxis(A, C);
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkAssumption4() -> bool{

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
#ifdef HINF_DEBUG
    std::cout << "Assumption 4 : " << std::endl;
    std::cout << std::boolalpha << ok << std::endl;
#endif
    return ok;
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkAllAssumption() -> bool{
    return checkAssumption1()
            & checkAssumption2()
            & checkAssumption3()
            & checkAssumption4();
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkCondition1() -> bool{

    arma::mat D11_ = llft_.D11().head_rows(performance_size - INPUT_SIZE);
    arma::mat D_11 = llft_.D11().head_cols(perturbation_size - OUTPUT_SIZE);

    arma::mat U, V;
    arma::vec s;
    arma::svd(U, s, V, D11_);
    auto sval1( s.max() );
    arma::svd(U, s, V, D_11);
    auto sval2( s.max() );

    return gam_ > std::max(sval1, sval2);
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkCondition2() -> bool{

    return linear_algebra::isPosSemiDefinite(X_inf_);
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkCondition3() -> bool{

    return linear_algebra::isPosSemiDefinite(Y_inf_);
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkCondition4() -> bool{
    arma::cx_mat temp( X_inf_*Y_inf_ );
    return linear_algebra::spectralRadius( std::move(temp) ) < gam_*gam_;
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkAllCondition() -> bool{
    return checkCondition1()
            & checkCondition2()
            & checkCondition3()
            & checkCondition4();
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto HInf<_StateSpace,
    performance_size,
    perturbation_size>::solve() -> common::StateSpacePack{
#ifdef HINF_DEBUG
    std::cout << "Interconnection matrix assumptions : " << std::endl;
#endif
    bool check_assumption = checkAllAssumption();
    assert(check_assumption && "The assumption made for interconnection matrix is not fulfill !");

    arma::mat temp1, temp2, temp3, temp4;
    arma::cx_mat ctemp1, ctemp2, ctemp3;

    arma::mat A_t = ss_->A().t();
    arma::mat B_t = ss_->B().t();
    arma::mat B1_t = llft_.B1().t();
    arma::mat C_t = ss_->C().t();
    arma::mat C1_t = llft_.C1().t();

    arma::mat D1_ = ss_->D().head_rows(performance_size);
    arma::mat D1__t = D1_.t();
    arma::mat D_1 = ss_->D().head_cols(perturbation_size);
    arma::mat D_1_t = D_1.t();

    temp1 = gam_*gam_*arma::eye(perturbation_size, perturbation_size);
    temp2 = arma::zeros<arma::mat>(D1_.n_cols,  D1_.n_cols);
    temp2.submat(0, 0, perturbation_size - 1, perturbation_size - 1) = temp1;
    arma::mat R1 = (D1__t * D1_) - temp2;
    //        R1.print("R1 : ");
    //-- create exception here for singular R1 and warn to change the gamma
    arma::mat R1_inv = arma::inv(R1);

    temp1 = gam_*gam_*arma::eye(performance_size, performance_size);
    temp2 = arma::zeros<arma::mat>(D_1.n_rows, D_1.n_rows);
    temp2.submat(0, 0, performance_size - 1, performance_size - 1) = temp1;
    arma::mat R2 = (D_1 * D_1_t) - temp2;
    //        R2.print("R2 : ");
    //-- create exception here for singular R2 and warn to change the gamma
    arma::mat R2_inv = arma::inv(R2);

    //-- Ricatti Domain

    temp1 = arma::join_vert(
                //-- create a function to simplify the zeros
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

    ARE<_StateSpace> solver1(ss_);
    ARE<_StateSpace> solver2(ss_);

    solver1.setHamiltonianMatrix(H_inf);
    solver2.setHamiltonianMatrix(J_inf);

    arma::cx_mat X_inf = solver1.solve();
    arma::cx_mat Y_inf = solver2.solve();

#ifdef HINF_DEBUG
    std::cout << ">>>>> Solution of ARE : " << std::endl;
    X_inf.print("X_inf : ");
    Y_inf.print("Y_inf : ");
#endif

    ctemp1 = toCx(D1__t * llft_.C1());
    ctemp2 = toCx(B_t) * X_inf;
    arma::cx_mat F = -toCx(R1_inv)*(ctemp1 + ctemp2);
    arma::cx_mat F1_inf = F.head_rows(perturbation_size);
    arma::cx_mat F2_inf = F.tail_rows(INPUT_SIZE);

    ctemp1 = toCx(llft_.B1() * D_1_t);
    ctemp2 = Y_inf * toCx(C_t);
    arma::cx_mat L = -(ctemp1 + ctemp2)*toCx(R2_inv);
    arma::cx_mat L1_inf = L.head_cols(performance_size);
    arma::cx_mat L2_inf = L.tail_cols(OUTPUT_SIZE);

    //-- Partition of D, F1_inf and L1_inf
    arma::cx_mat F11_inf = F1_inf.head_rows(perturbation_size - OUTPUT_SIZE);
    arma::cx_mat F12_inf = F1_inf.tail_rows(OUTPUT_SIZE);

    arma::cx_mat L11_inf = L1_inf.head_cols(performance_size - INPUT_SIZE);
    arma::cx_mat L12_inf = L1_inf.tail_cols(INPUT_SIZE);

    arma::cx_mat D1111 = toCx(llft_.D11()).submat(0, 0,
                                                  (performance_size - INPUT_SIZE) - 1, (perturbation_size - OUTPUT_SIZE) - 1);
    arma::cx_mat D1112 = toCx(llft_.D11()).submat(0, (perturbation_size - OUTPUT_SIZE),
                                                  (performance_size - INPUT_SIZE) - 1, perturbation_size - 1);
    arma::cx_mat D1121 = toCx(llft_.D11()).submat((performance_size - INPUT_SIZE), 0,
                                                  performance_size - 1, (perturbation_size - OUTPUT_SIZE) - 1);
    arma::cx_mat D1122 = toCx(llft_.D11()).submat((performance_size - INPUT_SIZE), (perturbation_size - OUTPUT_SIZE),
                                                  performance_size - 1,  perturbation_size - 1);

    X_inf_ = X_inf;
    Y_inf_ = Y_inf;

    bool check_cond = checkAllCondition();
#ifdef HINF_DEBUG
    std::cout << "Condition : " << std::boolalpha << check_cond << std::endl;
#endif
    assert(check_cond && "Condition for finding Controller that make the lower LFT is less than gamma was failed !");

    //-- I assume the arbitrary system that have property less than to gamma is zero
    ctemp1 = (gam_*gam_)*Y_inf*X_inf;
    arma::cx_mat Z_inf = arma::inv(arma::eye(_StateSpace::n_states, _StateSpace::n_states) - ctemp1);

    ctemp1 = -D1121*arma::trans(D1111);
    ctemp2 = (gam_*gam_)*arma::eye(performance_size - INPUT_SIZE, performance_size - INPUT_SIZE)
            - D1111*arma::trans(D1111);
    ctemp3 = ctemp1*arma::inv(ctemp2);
    arma::cx_mat D11_hat = ctemp3*D1112 - D1122;

    //-- calculate D12_hat and D21_hat
    ctemp1 = (gam_*gam_)*arma::eye(perturbation_size - OUTPUT_SIZE, perturbation_size - OUTPUT_SIZE)
            - arma::trans(D1111)*D1111;
    ctemp2 = arma::eye(INPUT_SIZE, INPUT_SIZE) - D1121*arma::inv(ctemp1)*arma::trans(D1121);
    arma::cx_mat D12_hat = arma::chol(ctemp2, "lower");

    ctemp1 = (gam_*gam_)*arma::eye(performance_size - INPUT_SIZE, performance_size - INPUT_SIZE)
            - D1111*arma::trans(D1111);
    ctemp2 = arma::eye(OUTPUT_SIZE, OUTPUT_SIZE) - arma::trans(D1112)*arma::inv(ctemp1)*D1112;
    arma::cx_mat D21_hat = arma::chol(ctemp2, "upper");

    ctemp1 = llft_.B2() + L12_inf;
    arma::cx_mat B2_hat = Z_inf*ctemp1*D12_hat;
    arma::cx_mat C2_hat = -D21_hat*(llft_.C2() + F12_inf);
    ctemp1 = -Z_inf*L2_inf;
    ctemp2 = B2_hat*arma::inv(D12_hat)*D11_hat;
    arma::cx_mat B1_hat = ctemp1 + ctemp2;
    ctemp1 = D11_hat*arma::inv(D21_hat)*C2_hat;
    arma::cx_mat C1_hat = F2_inf + ctemp1;
    ctemp1 = B1_hat*arma::inv(D21_hat)*C2_hat;
    arma::cx_mat A_hat = ss_->A() + ss_->B()*F + ctemp1;

#ifdef HINF_DEBUG
    std::cout << ">>>>> Controller Result : " << std::endl;
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

}

}
