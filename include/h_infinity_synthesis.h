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
        , gam_(15.){

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

    bool checkCondition1();
    bool checkCondition2();
    bool checkCondition3();
    bool checkCondition4();
    bool checkAllCondition();

    double gam_;

    arma::cx_mat X_inf_;
    arma::cx_mat Y_inf_;

    inline arma::cx_mat toCx(const arma::mat& _in) const{
        return arma::cx_mat(_in, arma::zeros<arma::mat>( arma::size(_in) ));
    }

    template <typename __StateSpace>
    bool isInfNormLessThan(double _gam, const __StateSpace& _ss){
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
        H.submat(__StateSpace::n_states,                     0, __StateSpace::n_states * 2 - 1,     __StateSpace::n_states - 1) = -C_t*temp2*_ss.C();

        arma::cx_mat eigvec;
        arma::cx_vec eigval;
        arma::eig_gen(eigval, eigvec, H);
        arma::mat eigval_re = arma::real(eigval);
        bool ok(true);
        for(int i(0); i < eigval.n_rows; i++){
            if(std::fabs(eigval_re(i, 0)) < std::numeric_limits<double>::epsilon()){
                ~ok;
                break;
            }
        }
        return ok;
    }

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

//            arma::mat zeros1(performance_size - INPUT_SIZE, INPUT_SIZE, arma::fill::zeros);
//            arma::mat zeros2(INPUT_SIZE, performance_size - INPUT_SIZE, arma::fill::zeros);
//            arma::mat eye(performance_size - INPUT_SIZE, performance_size - INPUT_SIZE, arma::fill::eye);

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

    bool ok2(true); //-- RVO
    arma::mat I_out(OUTPUT_SIZE, OUTPUT_SIZE, arma::fill::eye);
    arma::mat expected_D21 = arma::join_horiz(arma::mat(OUTPUT_SIZE, perturbation_size - OUTPUT_SIZE, arma::fill::zeros), I_out);

    if(perturbation_size > OUTPUT_SIZE){
        if(arma::approx_equal(llft_.D12(), expected_D21, "absdiff", .0)){

        }else{
            //-- do normalization here using SVD

//            arma::mat zeros1(OUTPUT_SIZE, perturbation_size - OUTPUT_SIZE, arma::fill::zeros);
//            arma::mat zeros2(perturbation_size - OUTPUT_SIZE, OUTPUT_SIZE, arma::fill::zeros);
//            arma::mat eye(perturbation_size - OUTPUT_SIZE, perturbation_size - OUTPUT_SIZE, arma::fill::eye);

            arma::mat U, V;
            arma::vec s;
            arma::svd(U,s,V,llft_.D21());

            arma::mat S = arma::diagmat(s);
            S.print("S : ");
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
bool HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkCondition1(){

    arma::mat D11_ = llft_.D11().head_rows(performance_size - INPUT_SIZE);
    arma::mat D_11 = llft_.D11().head_cols(perturbation_size - OUTPUT_SIZE);

    arma::mat U, V;
    arma::vec s;
    arma::svd(U, s, V, D11_);
    auto sval1( s.max() );
    arma::svd(U, s, V, D_11);
    auto sval2( s.max() );

//    std::cout << "VAL : " << sval1 << " : " << sval2 << std::endl;

    return gam_ > std::max(sval1, sval2);
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkCondition2(){
    return linear_algebra::isPosSemiDefinite(X_inf_);
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkCondition3(){

    return linear_algebra::isPosSemiDefinite(Y_inf_);
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkCondition4(){
    arma::cx_mat temp( X_inf_*Y_inf_ );
    return linear_algebra::spectralRadius( std::move(temp) ) < gam_*gam_;
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
bool HInf<_StateSpace,
    performance_size,
    perturbation_size>::checkAllCondition(){
    return checkCondition1()
            & checkCondition2()
            & checkCondition3()
            & checkCondition4();
}

template <class _StateSpace,
          std::size_t performance_size,
          std::size_t perturbation_size>
void HInf<_StateSpace,
    performance_size,
    perturbation_size>::solve(){

    if(checkAllAssumption()){
        arma::mat temp1, temp2, temp3, temp4;
        arma::cx_mat ctemp1, ctemp2, ctemp3, ctemp4;

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
        //-- create exception here for singular R1 and warn for change the gamma
        arma::mat R1_inv = arma::inv(R1);

        temp1 = gam_*gam_*arma::eye(performance_size, performance_size);
        temp2 = arma::zeros<arma::mat>(D_1.n_rows, D_1.n_rows);
        temp2.submat(0, 0, performance_size - 1, performance_size - 1) = temp1;
        arma::mat R2 = (D_1 * D_1_t) - temp2;
//        R2.print("R2 : ");
        //-- create exception here for singular R2 and warn for change the gamma
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
//        H_inf.print("H_inf : ");

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
//        J_inf.print("J_inf : ");

        ARE<_StateSpace> solver1(ss_);
        ARE<_StateSpace> solver2(ss_);

        solver1.setHamiltonianMatrix(H_inf);
        solver2.setHamiltonianMatrix(J_inf);

        arma::cx_mat X_inf = solver1.solve();
        arma::cx_mat Y_inf = solver2.solve();

        X_inf.print("X_inf : ");
        Y_inf.print("Y_inf : ");

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
        std::cout << "Condition : " << std::boolalpha << check_cond << std::endl;

    }
}

}

}
