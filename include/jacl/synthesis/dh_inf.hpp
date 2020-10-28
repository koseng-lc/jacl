/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of discrete H-infinity synthesis.
*            A lot of code here, is constructed by copy constructor because it's easier to read.
*/

#pragma once

#include <jacl/dare.hpp>
#include <jacl/traits.hpp>
#include <jacl/medium.hpp>
#include <jacl/lft/lower.hpp>
#include <jacl/lti_common.hpp>
#include <jacl/linear_algebra.hpp>
#include <jacl/numerical_methods.hpp>
#include <jacl/state_space/linear.hpp>

#define DHINF_VERBOSE

namespace jacl::synthesis{

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
class DHinf:public ::jacl::system::detail::BaseSystemClient<typename _System::base_t>{
public:
    template<typename __System = _System,
             typename std::enable_if_t<::jacl::traits::is_discrete_system_v<__System>, int>* = nullptr>
    DHinf(__System* _sys, double _gam,
          typename _System::state_t _reg1,
          typename _System::state_t _reg2)
        : ss_(::jacl::system::detail::BaseSystemClient<typename _System::base_t>::ss(_sys))
        , llft_(ss_)
        , gam_(_gam)
        , reg1_(_reg1)
        , reg2_(_reg2){}

    ~DHinf();

    auto solve() -> ::jacl::lti_common::StateSpacePack;

private:
    using scalar_t = typename _System::scalar_t;
    using cx_scalar_t = std::complex<scalar_t>;

private:
    auto checkAssumption1();
    auto checkAssumption2();
    auto checkAssumption3();
    auto checkAllAssumption();    

private:
    typename _System::state_space_t* ss_;
    static constexpr auto INPUT_SIZE{_System::n_inputs - perturbation_size};
    static constexpr auto OUTPUT_SIZE{_System::n_outputs - performance_size};
    lft::LowerLFT<typename _System::state_space_t,
             performance_size,
             perturbation_size,
             OUTPUT_SIZE,
             INPUT_SIZE> llft_;

    double gam_;   
    typename _System::state_t reg1_;
    typename _System::state_t reg2_;
};

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
DHinf<_System,
    performance_size,
    perturbation_size>::~DHinf(){
    ss_ = nullptr;
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto DHinf<_System,
    performance_size,
    perturbation_size>::checkAssumption1(){

    auto ctrb = ::jacl::lti_common::stabilizable<false>(llft_.A(), llft_.B2());
    auto obsv = ::jacl::lti_common::detectability<false>(llft_.A(), llft_.C2());
#ifdef DHINF_VERBOSE
    std::cout << "[DHinf] Assumption 1 : " << std::boolalpha << ctrb << " ; " << obsv << std::endl;
#endif
    return ctrb & obsv;
}

/*
 * This assumption is equivalent with Lemma 12.6 in Essential Robust Control
 * we only need to change the eigenvalue check on PBH test, to unit circle
 */
template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto DHinf<_System,
    performance_size,
    perturbation_size>::checkAssumption2(){
    
    typename arma::Mat<scalar_t>::template
        fixed<INPUT_SIZE, performance_size> D12_t( llft_.D12().t() );
    typename arma::Mat<scalar_t>::template
        fixed<INPUT_SIZE, INPUT_SIZE> R( D12_t*llft_.D12() );
    decltype(R) R_inv( arma::inv(R) );
    typename arma::Mat<scalar_t>::template
        fixed<performance_size, performance_size> I(arma::fill::eye);

    arma::mat temp1, temp2;
    temp1 = llft_.D12()*R_inv*D12_t;
    temp2 = I - temp1;
    typename arma::Mat<scalar_t>::template
        fixed<performance_size, _System::n_states> C = temp2*llft_.C1();

    temp1 = llft_.B2()*R_inv*D12_t;
    temp2 = temp1*llft_.C1();
    typename arma::Mat<scalar_t>::template
        fixed<_System::n_states, _System::n_states> A = llft_.A() - temp2;    

    auto ok = !::jacl::lti_common::hasUnobservableModeInUnitCircle(A, C);
   
#ifdef DHINF_VERBOSE
    std::cout << "[DHinf] Assumption 2 : " << std::boolalpha << ok << std::endl;
#endif
    return ok;
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto DHinf<_System,
    performance_size,
    perturbation_size>::checkAssumption3(){

    typename arma::Mat<scalar_t>::template
        fixed<perturbation_size, OUTPUT_SIZE> D21_t( llft_.D21().t() );
    typename arma::Mat<scalar_t>::template
        fixed<OUTPUT_SIZE, OUTPUT_SIZE> R( llft_.D21()*D21_t );
    decltype(R) R_inv( arma::inv(R) );
    typename arma::Mat<scalar_t>::template
        fixed<perturbation_size, perturbation_size> I(arma::fill::eye);

    arma::mat temp1, temp2;
    temp1 = llft_.B1()*D21_t*R_inv;
    temp2 = temp1*llft_.C2();
    typename arma::Mat<scalar_t>::template
        fixed<_System::n_states, _System::n_states> A = llft_.A() - temp2;

    temp1 = D21_t*R_inv*llft_.D21();
    temp2 = (I - temp1);
    typename arma::Mat<scalar_t>::template
        fixed<_System::n_states, perturbation_size> B = llft_.B1()*temp2;

    auto ok = !::jacl::lti_common::hasUncontrollableModeInUnitCircle(A, B);

#ifdef DHINF_VERBOSE
    std::cout << "[DHinf] Assumption 3 : " << std::boolalpha << ok << std::endl;
#endif
    return ok;
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto DHinf<_System,
    performance_size,
    perturbation_size>::checkAllAssumption(){
    return checkAssumption1()
            & checkAssumption2()
            & checkAssumption3();
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto DHinf<_System,
    performance_size,
    perturbation_size>::solve() -> ::jacl::lti_common::StateSpacePack{        

#ifdef DHINF_VERBOSE
    std::cout << "[DHinf] Interconnection matrix assumptions : " << std::endl;
#endif
    auto check_assumption = checkAllAssumption();
    assert(check_assumption && "[DHinf] The assumption made for interconnection matrix is not satisfied !");

    using namespace ::jacl::linear_algebra;

    arma::mat temp1, temp2, temp3, temp4, temp5;
    arma::cx_mat ctemp1, ctemp2, ctemp3, ctemp4, ctemp5, ctemp6;
    //-- Setup several frequently used matrix in calculations
    const typename arma::Mat<scalar_t>::template
        fixed<_System::n_states, _System::n_states>& A = llft_.A();
    typename arma::Mat<scalar_t>::template
        fixed<_System::n_states, _System::n_states> A_t = arma::trans(A);
    typename arma::Mat<cx_scalar_t>::template
        fixed<_System::n_states, _System::n_states> cxA = toCx(A);
    const typename arma::Mat<scalar_t>::template
        fixed<_System::n_states, perturbation_size+INPUT_SIZE>& B = ss_->B();
    const typename arma::Mat<scalar_t>::template
        fixed<_System::n_states, perturbation_size>& B1 = llft_.B1();
    typename arma::Mat<cx_scalar_t>::template
        fixed<_System::n_states, perturbation_size> cxB1 = toCx(B1);
    typename arma::Mat<scalar_t>::template
        fixed<perturbation_size, _System::n_states> B1_t = arma::trans(llft_.B1());
    typename arma::Mat<cx_scalar_t>::template
        fixed<perturbation_size, _System::n_states> cxB1_t = toCx(B1_t);
    const typename arma::Mat<scalar_t>::template
        fixed<_System::n_states, INPUT_SIZE>& B2 = llft_.B2();
    typename arma::Mat<cx_scalar_t>::template
        fixed<_System::n_states, INPUT_SIZE> cxB2 = toCx(B2);
    typename arma::Mat<scalar_t>::template
        fixed<INPUT_SIZE, _System::n_states> B2_t = arma::trans(llft_.B2());
    typename arma::Mat<cx_scalar_t>::template
        fixed<INPUT_SIZE, _System::n_states> cxB2_t = toCx(B2_t);
    const typename arma::Mat<scalar_t>::template
        fixed<performance_size+OUTPUT_SIZE, _System::n_states>& C = ss_->C();
    const typename arma::Mat<scalar_t>::template
        fixed<performance_size, _System::n_states>& C1 = llft_.C1();
    typename arma::Mat<cx_scalar_t>::template
        fixed<performance_size, _System::n_states> cxC1 = toCx(C1);
    typename arma::Mat<scalar_t>::template
        fixed<_System::n_states, performance_size> C1_t = arma::trans(llft_.C1());
    const typename arma::Mat<scalar_t>::template
        fixed<OUTPUT_SIZE, _System::n_states>& C2 = llft_.C2();
    typename arma::Mat<cx_scalar_t>::template
        fixed<OUTPUT_SIZE, _System::n_states> cxC2 = toCx(C2);
    typename arma::Mat<scalar_t>::template
        fixed<_System::n_states, OUTPUT_SIZE> C2_t = arma::trans(llft_.C2());
    const typename arma::Mat<scalar_t>::template
        fixed<performance_size, perturbation_size>& D11 = llft_.D11();
    typename arma::Mat<cx_scalar_t>::template
        fixed<performance_size, perturbation_size> cxD11 = toCx(D11);
    typename arma::Mat<scalar_t>::template
        fixed<perturbation_size, performance_size> D11_t = arma::trans(llft_.D11());
    const typename arma::Mat<scalar_t>::template
        fixed<performance_size, INPUT_SIZE>& D12 = llft_.D12();
    typename arma::Mat<cx_scalar_t>::template
        fixed<performance_size, INPUT_SIZE> cxD12 = toCx(D12);
    typename arma::Mat<scalar_t>::template
        fixed<INPUT_SIZE, performance_size> D12_t = arma::trans(llft_.D12());
    const typename arma::Mat<scalar_t>::template
        fixed<OUTPUT_SIZE, perturbation_size>& D21 = llft_.D21();
    typename arma::Mat<cx_scalar_t>::template
        fixed<OUTPUT_SIZE, perturbation_size> cxD21 = toCx(D21);
    typename arma::Mat<scalar_t>::template
        fixed<perturbation_size, OUTPUT_SIZE> D21_t = arma::trans(llft_.D21());
    const typename arma::Mat<scalar_t>::template
        fixed<OUTPUT_SIZE, INPUT_SIZE>& D22 = llft_.D22();
    typename arma::Mat<cx_scalar_t>::template
        fixed<OUTPUT_SIZE, INPUT_SIZE> cxD22 = toCx(D22);
    typename arma::Mat<scalar_t>::template
        fixed<INPUT_SIZE, OUTPUT_SIZE> D22_t = arma::trans(llft_.D22());

    typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, _System::n_states> Aclp;
    typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, _System::n_states> Ap;
    typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, perturbation_size> Ep;
    typename arma::Mat<cx_scalar_t>::template
            fixed<OUTPUT_SIZE, _System::n_states> C1p;
    typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE, _System::n_states> C2p;
    typename arma::Mat<cx_scalar_t>::template
            fixed<OUTPUT_SIZE, perturbation_size> D12p;
    typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE, INPUT_SIZE> D21p;
    typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE, perturbation_size> D22p;
    typename arma::Mat<cx_scalar_t>::template
            fixed<perturbation_size, _System::n_states> Z;
    {    
        typename arma::Mat<scalar_t>::template
            fixed<INPUT_SIZE+perturbation_size, INPUT_SIZE+perturbation_size> S =
            arma::join_cols(
                arma::join_rows(D12_t*D12, D12_t*D11),
                arma::join_rows(D11_t*D12, D11_t*D11 - (gam_*gam_)*arma::eye(perturbation_size,perturbation_size))
        );
        typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE+perturbation_size, INPUT_SIZE+perturbation_size> cxS = toCx(S);
        decltype(S) S_inv = arma::inv(S);

        typename arma::Mat<scalar_t>::template
            fixed<performance_size, INPUT_SIZE+perturbation_size> D12_D11 = arma::join_rows(D12, D11);
        typename arma::Mat<scalar_t>::template
            fixed<INPUT_SIZE+perturbation_size, performance_size> D12_D11_t = arma::trans(D12_D11);
        typename arma::Mat<scalar_t>::template
            fixed<_System::n_states, INPUT_SIZE+perturbation_size> B2_B1 = arma::join_rows(B2, B1);        

        typename arma::Mat<scalar_t>::template
            fixed<_System::n_states, _System::n_states> X = B2_B1*S_inv*arma::trans(B2_B1);
        decltype(X) X_t = arma::trans(X);
 
        temp1 = B2_B1*S_inv*D12_D11_t;
        typename arma::Mat<scalar_t>::template
            fixed<_System::n_states, _System::n_states> A_tilde = A - (temp1*C1);
        decltype(A_tilde) reg_A_tilde = A_tilde;
        if(arma::cond(A_tilde) > 1e6){
            //-- regularize the ill-conditioned matrix
            reg_A_tilde = A_tilde + arma::diagmat(reg1_)*arma::eye(arma::size(A));
            #ifdef DHINF_VERBOSE
            std::cout << "[DHinf] Regularization triggered" << std::endl;
            #endif
        }
        //-- more stable with solve() than with inv()
        decltype(A_tilde) A_tilde_tinv = arma::solve(reg_A_tilde.t(), arma::eye(arma::size(A)),
            arma::solve_opts::refine + arma::solve_opts::equilibrate + arma::solve_opts::allow_ugly);
    
        temp1 = D12_D11*S_inv*D12_D11_t;
        temp2 = C1_t*temp1*C1;
        typename arma::Mat<scalar_t>::template
            fixed<_System::n_states, _System::n_states> Q = C1_t*C1 - temp2;

        temp1 = A_tilde_tinv*Q;
        typename arma::Mat<scalar_t>::template
            fixed<_System::n_states<<1, _System::n_states<<1> H = arma::join_cols(
            arma::join_rows(A_tilde+X*temp1, -X*A_tilde_tinv),
            arma::join_rows(-temp1, A_tilde_tinv)
        );

        #ifdef DHINF_VERBOSE
        //-- check whether symplectic or not
        decltype(H) J = arma::join_cols(
            arma::join_rows(arma::zeros(arma::size(A)), -arma::eye(arma::size(A))),
            arma::join_rows(arma::eye(arma::size(A)), arma::zeros(arma::size(A)))
        );      
        decltype(H) check = H.t()*J*H;
        check.print("[DHinf] Check symplectic 1 : ");
        #endif

        typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, _System::n_states> P = ::jacl::dare::solve(H);
        decltype(P) P_t = arma::trans(P);
        
        typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, _System::n_states> cx_X = toCx(X);
        typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, _System::n_states> cxA_tilde = toCx(A_tilde);
        decltype(cxA_tilde) cxA_tilde_t = arma::trans(cxA_tilde);

        #ifdef DHINF_VERBOSE
        //-- check the DARE solution
        ctemp1 = cxA_tilde_t*P*cxA_tilde;
        ctemp2 = ctemp1 + toCx(C1_t*C1);
        ctemp3 = cxA_tilde_t*P*cx_X;
        ctemp4 = arma::inv(arma::eye(arma::size(P*cx_X)) + P*cx_X);
        ctemp5 = P*cxA_tilde;
        ctemp6 = ctemp3*ctemp4*ctemp5;
        decltype(cx_X) dare_rhs = ctemp2 - ctemp6 - P;
        dare_rhs.print("[DHinf] DARE RHS 1 : ");
        #endif

        //-- common term
        typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE, _System::n_states> cterm1 = cxB2_t*P*cxA + toCx(D12_t*C1);
        typename arma::Mat<cx_scalar_t>::template
            fixed<perturbation_size, _System::n_states> cterm2 = cxB1_t*P*cxA + toCx(D11_t*C1);
        typename arma::Mat<cx_scalar_t>::template
            fixed<perturbation_size, INPUT_SIZE> cterm3 = cxB1_t*P*cxB2 + toCx(D11_t*D12);
        typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE, perturbation_size> cterm4 = cxB2_t*P*cxB1 + toCx(D12_t*D11);
        typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE + perturbation_size, _System::n_states> cterm5 = arma::join_cols(cterm1, cterm2);

        ctemp1 = cxB2_t*P*cxB2;
        ctemp2 = toCx(D12_t*D12);
        typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE, INPUT_SIZE> V = ctemp1 + ctemp2;
        decltype(V) V_inv = arma::inv(V);

        ctemp1 = toCx(D11_t*D11);
        ctemp2 = cxB1_t*P*cxB1;
        ctemp3 = cterm3*V_inv*cterm4;
        typename arma::Mat<cx_scalar_t>::template
            fixed<perturbation_size, perturbation_size> R =
            (gam_*gam_)*arma::eye<arma::cx_mat>(perturbation_size, perturbation_size)
            - ctemp1 - ctemp2 + ctemp3;
        decltype(R) R_inv = arma::inv(R);

        ctemp1 = toCx(S);
        ctemp2 = toCx(B2_B1);
        typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE+perturbation_size, INPUT_SIZE+perturbation_size> G = ctemp1 + arma::trans(ctemp2)*P*ctemp2;

        decltype(G) G_inv = arma::inv(G);                

        ctemp1 = B2_B1*G_inv*cterm5;
        Aclp = cxA - ctemp1;
        assert(lti_common::isStable(Aclp, false) && "[DHinf] Aclp is not asymptotically stable !");
        
        decltype(R) R_mhp = arma::powmat((1/(gam_*gam_))*R, -.5);        
        decltype(V) V_mhp = arma::powmat(V, -.5);
        Z = cterm2 - cterm3*V_inv*cterm1;
        Ap = cxA + cxB1*R_inv*Z;
        Ep = cxB1*R_mhp;
        C1p = cxC2 + cxD21*R_inv*Z;
        ctemp1 = V_mhp*cterm4*R_inv*Z;
        C2p = V_mhp*cterm1 + ctemp1;
        D12p = cxD21*R_mhp;
        D21p = arma::powmat(V, .5); //-- another way to solve Cholesky factorization
        D22p = V_mhp*cterm4*R_mhp;
    }

    typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE, OUTPUT_SIZE> N;
    typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE, _System::n_states> M;
    typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, OUTPUT_SIZE> L;
    typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, _System::n_states> K;
    {
        typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, _System::n_states> Ap_t = arma::trans(Ap);
        typename arma::Mat<cx_scalar_t>::template
            fixed<perturbation_size, _System::n_states> Ep_t = arma::trans(Ep);
        typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, OUTPUT_SIZE> C1p_t = arma::trans(C1p);
        typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, INPUT_SIZE> C2p_t = arma::trans(C2p);
        typename arma::Mat<cx_scalar_t>::template
            fixed<perturbation_size, OUTPUT_SIZE> D12p_t = arma::trans(D12p);
        typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE, INPUT_SIZE> D21p_t = arma::trans(D21p);
        typename arma::Mat<cx_scalar_t>::template
            fixed<perturbation_size, INPUT_SIZE> D22p_t = arma::trans(D22p);        

        typename arma::Mat<cx_scalar_t>::template
            fixed<OUTPUT_SIZE+INPUT_SIZE, OUTPUT_SIZE+INPUT_SIZE> S = arma::join_cols(
            arma::join_rows(D12p*D12p_t, D12p*D22p_t),
            arma::join_rows(D22p*D12p_t, D22p*D22p_t - (gam_*gam_)*arma::eye<arma::cx_mat>(arma::size(D22p*D22p_t)))
        );
        decltype(S) S_inv = arma::inv(S);

        typename arma::Mat<cx_scalar_t>::template
            fixed<perturbation_size, OUTPUT_SIZE+INPUT_SIZE> D12_D22 = arma::join_rows(D12p_t, D22p_t);
        typename arma::Mat<cx_scalar_t>::template
            fixed<OUTPUT_SIZE+INPUT_SIZE, perturbation_size> D12_D22_t = arma::trans(D12_D22);
        typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, OUTPUT_SIZE+INPUT_SIZE> C1_C2 = arma::join_rows(C1p_t, C2p_t);

        typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, _System::n_states> X = C1_C2*S_inv*arma::trans(C1_C2);
        decltype(X) X_t = arma::trans(X);

        ctemp1 = C1_C2*S_inv*D12_D22_t;
        typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, _System::n_states> A_tilde = Ap_t - (ctemp1*Ep_t);
        decltype(A_tilde) reg_A_tilde = A_tilde;
        if(arma::cond(A_tilde) > 1e6){
            //-- regularize the ill-conditioned matrix
            reg_A_tilde = A_tilde + arma::diagmat(reg2_)*arma::eye(arma::size(Ap_t));
            #ifdef DHINF_VERBOSE
            std::cout << "[DHinf] Regularization triggered" << std::endl;
            #endif
        }
        //-- more stable with solve() than with inv()
        decltype(A_tilde) A_tilde_t = arma::trans(A_tilde);
        decltype(A_tilde) A_tilde_tinv = arma::solve(reg_A_tilde.t(), arma::eye<arma::cx_mat>(arma::size(A)),
            arma::solve_opts::refine + arma::solve_opts::equilibrate + arma::solve_opts::allow_ugly);

        ctemp1 = D12_D22*S_inv*D12_D22_t;
        ctemp2 = Ep*ctemp1*Ep_t;
        typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, _System::n_states> Q = Ep*Ep_t - ctemp2;

        ctemp1 = A_tilde_tinv*Q;
        typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states<<1, _System::n_states<<1> H = arma::join_cols(
            arma::join_rows(A_tilde+X*ctemp1, -X*A_tilde_tinv),
            arma::join_rows(-ctemp1, A_tilde_tinv)
        );

        #ifdef DHINF_VERBOSE
        //-- check whether symplectic or not
        typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states<<1, _System::n_states<<1> J = arma::join_cols(
            arma::join_rows(arma::zeros<arma::cx_mat>(arma::size(A)), -arma::eye<arma::cx_mat>(arma::size(A))),
            arma::join_rows(arma::eye<arma::cx_mat>(arma::size(A)), arma::zeros<arma::cx_mat>(arma::size(A)))
        );        
        decltype(J) check = H.t()*J*H;
        check.print("[DHihf] Check symplectic 2 : ");
        #endif

        typename arma::Mat<cx_scalar_t>::template
            fixed<_System::n_states, _System::n_states> Y = ::jacl::dare::solve(toReal(H));
        decltype(Y) Y_t = arma::trans(Y);
        
        #ifdef DHINF_VERBOSE
        //-- check the DARE solution
        ctemp1 = A_tilde_t*Y*A_tilde;
        ctemp2 = ctemp1 + Ep*Ep_t;
        ctemp3 = A_tilde_t*Y*X;
        ctemp4 = arma::inv(arma::eye<arma::cx_mat>(arma::size(Y*X)) + Y*X);
        ctemp5 = Y*A_tilde;
        ctemp6 = ctemp3*ctemp4*ctemp5;
        decltype(Y) dare_rhs = ctemp2 - ctemp6 - Y;
        dare_rhs.print("[DHinf] DARE RHS 2 : ");
        #endif

        //-- common term
        typename arma::Mat<cx_scalar_t>::template
            fixed<OUTPUT_SIZE, _System::n_states> cterm1 = C1p*Y*Ap_t + D12p*Ep_t;
        typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE, _System::n_states> cterm2 = C2p*Y*Ap_t + D22p*Ep_t;
        typename arma::Mat<cx_scalar_t>::template
            fixed<OUTPUT_SIZE, INPUT_SIZE> cterm3 = C1p*Y*C2p_t + D12p*D22p_t;
        typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE, OUTPUT_SIZE> cterm4 = C2p*Y*C1p_t + D22p*D12p_t;
        typename arma::Mat<cx_scalar_t>::template
            fixed<OUTPUT_SIZE+INPUT_SIZE, _System::n_states> cterm5 = arma::join_cols(cterm1, cterm2);

        typename arma::Mat<cx_scalar_t>::template
            fixed<OUTPUT_SIZE, OUTPUT_SIZE> V = D12p*D12p_t + C1p*Y*C1p_t;
        decltype(V) V_inv = arma::inv(V);

        ctemp1 = D22p*D22p_t + C2p*Y*C2p_t;
        ctemp2 = cterm4*V_inv*cterm3;
        typename arma::Mat<cx_scalar_t>::template
            fixed<INPUT_SIZE, INPUT_SIZE> R = (gam_*gam_)*arma::eye<arma::cx_mat>(INPUT_SIZE, INPUT_SIZE)
            - ctemp1 - ctemp2;

        decltype(S) G = S + arma::trans(C1_C2)*Y*C1_C2;
        decltype(G) G_inv = arma::inv(G);

        decltype(Ap) Acl = Ap - arma::trans(cterm5)*G_inv*arma::trans(C1_C2);
        assert(lti_common::isStable(Aclp, false) && "[DHinf] Acl is not asymptotically stable !");

        N = -arma::inv(D21p)*cterm4*V_inv;
        M = -(arma::inv(D21p)*C2p + N*C1p);
        L = B2*N + arma::trans(cterm1)*V_inv;
        K = Aclp - L*C1p;
    }

    /*
        Controller symbol adapted from A.A. Stoorvogel
        x(k+1) = K*x(k) + L*y(k)
        u(k) = M*x(k) + N*y(k)
    */
    #ifdef DHINF_VERBOSE
    std::cout << "[DHinf] Controller result : " << std::endl;
    K.print("Ak : "); L.print("Bk : ");
    M.print("Ck : "); N.print("Dk : ");
    #endif
    
    return std::make_tuple(toReal(K), toReal(L), toReal(M), toReal(N));
}

} // namespace jacl::synthesis
