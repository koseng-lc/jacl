/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of discrete H-infinity synthesis.
*            A lot of code here, is constructed by copy constructor because it's easier to read.
*/

#pragma once

#include <jacl/lti_common.hpp>
#include <jacl/dare.hpp>
#include <jacl/lower_lft.hpp>
#include <jacl/linear_state_space.hpp>
#include <jacl/linear_algebra.hpp>
#include <jacl/traits.hpp>
#include <jacl/numerical_methods.hpp>
#include <jacl/medium.hpp>

#define DHINF_VERBOSE

namespace jacl{ namespace synthesis{

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
class DHinf:public ::jacl::system::detail::BaseSystemClient<typename _System::base_t>{
public:
    template<typename __System = _System,
             typename std::enable_if_t<traits::is_discrete_system<__System>::value, int>* = nullptr>
    DHinf(__System* _sys, double _gam)
        : ss_(::jacl::system::detail::BaseSystemClient<typename _System::base_t>::ss(_sys))
        , llft_(ss_)
        , gam_(_gam){}

    ~DHinf();

    auto solve() -> ::jacl::common::StateSpacePack;

    template <typename SS1, typename SS2>
    auto starProducts(const SS1 &_ss1, const SS2 &_ss2){
        arma::cx_mat A_bar(SS1::n_states + SS2::n_states, SS1::n_states + SS2::n_states);
        arma::cx_mat B_bar(SS1::n_states + SS2::n_states, SS1::w1_size + SS2::u_sz);
        arma::cx_mat C_bar(SS1::z1_sz + SS2::y_sz, SS1::n_states + SS2::n_states);
        arma::cx_mat D_bar(SS1::z1_sz + SS2::y_sz, SS1::w1_size, + SS2::u_sz);
    }

private:
    auto checkAssumption1();
    auto checkAssumption2();
    auto checkAssumption3();
    auto checkAllAssumption();

    template <typename _StateSpace>
    auto isInfNormLessThan(double _gam, const _StateSpace& _ss){
        arma::mat C_t = arma::trans( _ss.C() );
        arma::mat D_t = arma::trans( _ss.D() );
        arma::mat D_tD = D_t * _ss.D();
        arma::mat R = (_gam*_gam)*arma::eye(_StateSpace::n_inputs,
                                            _StateSpace::n_inputs)
                        - D_tD;
        arma::mat R_inv = arma::inv(R);        

        arma::mat temp1, temp2, temp3, temp4;

        //-- block 1,1
        temp1 = R_inv*D_t*_ss.C();
        temp2 = _ss.A() + temp1;
        temp3 = _ss.D()*R_inv*D_t;
        temp4 = arma::eye(_StateSpace::n_outputs,
                          _StateSpace::n_outputs)
                 - temp3;
        //-- Hamiltonian matrix
        arma::mat H = arma::join_cols(
            arma::join_rows(temp2, _ss.B()*R_inv*arma::trans(_ss.B())),
            arma::join_rows(-C_t*temp4*_ss.C(), -arma::trans(temp2))
        );

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
    typename _System::state_space_t* ss_;
    static constexpr auto INPUT_SIZE{_System::n_inputs - perturbation_size};
    static constexpr auto OUTPUT_SIZE{_System::n_outputs - performance_size};
    lft::LowerLFT<typename _System::state_space_t,
             performance_size,
             perturbation_size,
             OUTPUT_SIZE,
             INPUT_SIZE> llft_;

    double gam_;
    arma::cx_mat X_inf_;
    arma::cx_mat Y_inf_;
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

    auto ctrb = common::stabilizable(llft_.A(), llft_.B2());
    auto obsv = common::detectability(llft_.A(), llft_.C2());
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
    perturbation_size>::solve() -> ::jacl::common::StateSpacePack{        

#ifdef DHINF_VERBOSE
    std::cout << "[DHinf] Interconnection matrix assumptions : " << std::endl;
#endif
    bool check_assumption = checkAllAssumption();
    assert(check_assumption && "[DHinf] The assumption made for interconnection matrix is not fulfill !");

    using namespace ::jacl::linear_algebra;

    arma::mat temp1, temp2, temp3, temp4, temp5;
    arma::cx_mat ctemp1, ctemp2, ctemp3, ctemp4, ctemp5, ctemp6;
    //-- Setup several frequently matrix operation
    const arma::mat& A = llft_.A();
    arma::mat A_t = arma::trans(A);
    arma::cx_mat cxA = toCx(A);
    const arma::mat& B = ss_->B();
    const arma::mat& B1 = llft_.B1();
    arma::cx_mat cxB1 = toCx(B1);
    arma::mat B1_t = arma::trans(llft_.B1());
    arma::cx_mat cxB1_t = toCx(B1_t);
    const arma::mat& B2 = llft_.B2();
    arma::cx_mat cxB2 = toCx(B2);
    arma::mat B2_t = arma::trans(llft_.B2());
    arma::cx_mat cxB2_t = toCx(B2_t);
    const arma::mat& C = ss_->C();
    const arma::mat& C1 = llft_.C1();
    arma::cx_mat cxC1 = toCx(C1);
    arma::mat C1_t = arma::trans(llft_.C1());
    const arma::mat& C2 = llft_.C2();
    arma::cx_mat cxC2 = toCx(C2);
    arma::mat C2_t = arma::trans(llft_.C2());
    const arma::mat& D11 = llft_.D11();
    arma::cx_mat cxD11 = toCx(D11);
    arma::mat D11_t = arma::trans(llft_.D11());
    const arma::mat& D12 = llft_.D12();
    arma::cx_mat cxD12 = toCx(D12);
    arma::mat D12_t = arma::trans(llft_.D12());
    const arma::mat& D21 = llft_.D21();
    arma::cx_mat cxD21 = toCx(D21);
    arma::mat D21_t = arma::trans(llft_.D21());
    const arma::mat& D22 = llft_.D22();
    arma::cx_mat cxD22 = toCx(D22);
    arma::mat D22_t = arma::trans(llft_.D22());

    arma::cx_mat Aclp;
    arma::cx_mat Z, Ap, Ep, C1p, C2p, D12p, D21p, D22p;
    {        
        arma::mat S = arma::join_cols(
            arma::join_rows(D12_t*D12, D12_t*D11),
            arma::join_rows(D11_t*D12, D11_t*D11 - (gam_*gam_)*arma::eye(perturbation_size,perturbation_size))
        );
        arma::cx_mat cxS = toCx(S);
        arma::mat S_inv = arma::inv(S);

        arma::mat D12_D11 = arma::join_rows(D12, D11);
        arma::mat D12_D11_t = arma::trans(D12_D11);
        arma::mat B2_B1 = arma::join_rows(B2, B1);        

        arma::mat X = B2_B1*S_inv*arma::trans(B2_B1);
        arma::mat X_t = arma::trans(X);
 
        temp1 = B2_B1*S_inv*D12_D11_t;
        arma::mat A_tilde = A - (temp1*C1);
        if(arma::cond(A_tilde) > 1e6){
            //-- regularize the ill-conditioned matrix
            A_tilde = A_tilde + .001*arma::eye(arma::size(A));
            #ifdef DHINF_VERBOSE
            std::cout << "[DHinf] Regularization triggered" << std::endl;
            #endif
        }
        //-- more stable with solve() than with inv()
        arma::mat A_tilde_tinv = arma::solve(A_tilde.t(), arma::eye(arma::size(A)),
            arma::solve_opts::refine);
    
        temp1 = D12_D11*S_inv*D12_D11_t;
        temp2 = C1_t*temp1*C1;
        arma::mat Q = C1_t*C1 - temp2;

        temp1 = A_tilde_tinv*Q;
        arma::mat H = arma::join_cols(
            arma::join_rows(A_tilde+X*temp1, -X*A_tilde_tinv),
            arma::join_rows(-temp1, A_tilde_tinv)
        );

        #ifdef DHINF_VERBOSE
        //-- check whether sympletic or not
        arma::mat J = arma::join_cols(
            arma::join_rows(arma::zeros(arma::size(A)), -arma::eye(arma::size(A))),
            arma::join_rows(arma::eye(arma::size(A)), arma::zeros(arma::size(A)))
        );      
        arma::mat check = H.t()*J*H;
        check.print("[DHinf] Check sympletic 1 : ");
        #endif

        DARE<typename _System::state_space_t> solver(ss_);
        solver.setSympleticMatrix(H);
        arma::cx_mat P = solver.solve();
        arma::cx_mat P_t = arma::trans(P);
        
        arma::cx_mat cx_X = toCx(X);
        arma::cx_mat cxA_tilde = toCx(A_tilde);
        arma::cx_mat cxA_tilde_t = arma::trans(cxA_tilde);

        #ifdef DHINF_VERBOSE
        //-- check the DARE solution
        ctemp1 = cxA_tilde_t*P*cxA_tilde;
        ctemp2 = ctemp1 + toCx(C1_t*C1);
        ctemp3 = cxA_tilde_t*P*cx_X;
        ctemp4 = arma::inv(arma::eye(arma::size(P*cx_X)) + P*cx_X);
        ctemp5 = P*cxA_tilde;
        ctemp6 = ctemp3*ctemp4*ctemp5;
        arma::cx_mat dare_rhs = ctemp2 - ctemp6 - P;
        dare_rhs.print("[DHinf] DARE RHS 1 : ");
        #endif

        //-- common term
        arma::cx_mat cterm1 = cxB2_t*P*cxA + toCx(D12_t*C1);
        arma::cx_mat cterm2 = cxB1_t*P*cxA + toCx(D11_t*C1);
        arma::cx_mat cterm3 = cxB1_t*P*cxB2 + toCx(D11_t*D12);
        arma::cx_mat cterm4 = cxB2_t*P*cxB1 + toCx(D12_t*D11);
        arma::cx_mat cterm5 = arma::join_cols(
            cterm1,
            cterm2
        );

        ctemp1 = cxB2_t*P*cxB2;
        ctemp2 = toCx(D12_t*D12);
        arma::cx_mat V = ctemp1 + ctemp2;
        arma::cx_mat V_inv = arma::inv(V);

        ctemp1 = toCx(D11_t*D11);
        ctemp2 = cxB1_t*P*cxB1;
        ctemp3 = cterm3*V_inv*cterm4;
        arma::cx_mat R = (gam_*gam_)*arma::eye<arma::cx_mat>(perturbation_size, perturbation_size)
            - ctemp1 - ctemp2 + ctemp3;
        arma::cx_mat R_inv = arma::inv(R);

        ctemp1 = toCx(S);
        ctemp2 = toCx(B2_B1);
        arma::cx_mat G = ctemp1 + arma::trans(ctemp2)*P*ctemp2;

        arma::cx_mat G_inv = arma::inv(G);

        // ctemp1 = toCx(A_t)*P*toCx(A);
        // ctemp2 = toCx(C1_t*C1);
        // arma::cx_mat D12_t_C1 = toCx(D12_t*C1);
        // arma::cx_mat D11_t_C1 = toCx(D11_t*C1);
        // ctemp3 = arma::join_cols(
        //     toCx(B2_t)*P_t*toCx(A) + D12_t_C1,
        //     toCx(B1_t)*P_t*toCx(A) + D11_t_C1
        // );
        // ctemp4 = arma::join_cols(
        //     toCx(B2_t)*P*toCx(A) + D12_t_C1,
        //     toCx(B1_t)*P*toCx(A) + D11_t_C1
        // );
        // ctemp5 = ctemp1 + ctemp2;
        // ctemp6 = arma::trans(ctemp3)*G_inv*ctemp4;
        // arma::cx_mat dare_rhs = ctemp5 - ctemp6 - P;                    

        //-- TODO : check whether is Acl asymptotically stable or not
        ctemp1 = B2_B1*G_inv*cterm5;
        Aclp = cxA - ctemp1;

        arma::cx_mat R_mhp = arma::powmat(R, -.5);
        arma::cx_mat V_mhp = arma::powmat(V, -.5);
        Z = cterm2 - cterm3*V_inv*cterm1;
        Ap = cxA + cxB1*R_inv*Z;
        Ep = cxB1*R_mhp;
        C1p = cxC2 + cxD21*R_inv*Z;
        ctemp1 = V_mhp*cterm4*R_inv*Z;
        C2p = V_mhp*cterm1 + ctemp1;
        D12p = cxD21*R_mhp;
        D21p = arma::powmat(V, .5);
        D22p = V_mhp*cterm4*R_mhp;
    }

    arma::cx_mat N, M, L, K;
    {
        arma::cx_mat Ap_t = arma::trans(Ap);
        arma::cx_mat Ep_t = arma::trans(Ep);
        arma::cx_mat C1p_t = arma::trans(C1p);
        arma::cx_mat C2p_t = arma::trans(C2p);
        arma::cx_mat D12p_t = arma::trans(D12p);
        arma::cx_mat D21p_t = arma::trans(D21p);
        arma::cx_mat D22p_t = arma::trans(D22p);        

        arma::cx_mat S = arma::join_cols(
            arma::join_rows(D12p*D12p_t, D12p*D22p_t),
            arma::join_rows(D22p*D12p_t, D22p*D22p_t - (gam_*gam_)*arma::eye<arma::cx_mat>(arma::size(D22p*D22p_t)))
        );
        arma::cx_mat S_inv = arma::inv(S);

        arma::cx_mat D12_D22 = arma::join_rows(D12p_t, D22p_t);
        arma::cx_mat D12_D22_t = arma::trans(D12_D22);
        arma::cx_mat C1_C2 = arma::join_rows(C1p_t, C2p_t);

        arma::cx_mat X = C1_C2*S_inv*arma::trans(C1_C2);
        arma::cx_mat X_t = arma::trans(X);

        ctemp1 = C1_C2*S_inv*D12_D22_t;
        arma::cx_mat A_tilde = Ap_t - (ctemp1*Ep_t);
        if(arma::cond(A_tilde) > 1e6){
            //-- regularize the ill-conditioned matrix
            A_tilde = A_tilde + .001*arma::eye(arma::size(Ap_t));
            #ifdef DHINF_VERBOSE
            std::cout << "[DHinf] Regularization triggered" << std::endl;
            #endif
        }
        //-- more stable with solve() than with inv()
        arma::cx_mat A_tilde_t = arma::trans(A_tilde);
        arma::cx_mat A_tilde_tinv = arma::solve(A_tilde_t, arma::eye<arma::cx_mat>(arma::size(A)),
            arma::solve_opts::refine);

        ctemp1 = D12_D22*S_inv*D12_D22_t;
        ctemp2 = Ep*ctemp1*Ep_t;
        arma::cx_mat Q = Ep*Ep_t - ctemp2;

        ctemp1 = A_tilde_tinv*Q;
        arma::cx_mat H = arma::join_cols(
            arma::join_rows(A_tilde+X*ctemp1, -X*A_tilde_tinv),
            arma::join_rows(-ctemp1, A_tilde_tinv)
        );

        #ifdef DHINF_VERBOSE
        //-- check whether sympletic or not
        arma::cx_mat J = arma::join_cols(
            arma::join_rows(arma::zeros<arma::cx_mat>(arma::size(A)), -arma::eye<arma::cx_mat>(arma::size(A))),
            arma::join_rows(arma::eye<arma::cx_mat>(arma::size(A)), arma::zeros<arma::cx_mat>(arma::size(A)))
        );        
        arma::cx_mat check = H.t()*J*H;
        check.print("[DHihf] Check sympletic 2 : ");
        #endif

        DARE<typename _System::state_space_t> solver(ss_);
        solver.setSympleticMatrix(toReal(H));
        arma::cx_mat Y = solver.solve();
        arma::cx_mat Y_t = arma::trans(Y);
        
        #ifdef DHINF_VERBOSE
        //-- check the DARE solution
        ctemp1 = A_tilde_t*Y*A_tilde;
        ctemp2 = ctemp1 + Ep*Ep_t;
        ctemp3 = A_tilde_t*Y*X;
        ctemp4 = arma::inv(arma::eye<arma::cx_mat>(arma::size(Y*X)) + Y*X);
        ctemp5 = Y*A_tilde;
        ctemp6 = ctemp3*ctemp4*ctemp5;
        arma::cx_mat dare_rhs = ctemp2 - ctemp6 - Y;
        dare_rhs.print("[DHinf] DARE RHS 2 : ");
        #endif

        //-- common term
        arma::cx_mat cterm1 = C1p*Y*Ap_t + D12p*Ep_t;
        arma::cx_mat cterm2 = C2p*Y*Ap_t + D22p*Ep_t;
        arma::cx_mat cterm3 = C1p*Y*C2p_t + D12p*D22p_t;
        arma::cx_mat cterm4 = C2p*Y*C1p_t + D22p*D12p_t;
        arma::cx_mat cterm5 = arma::join_cols(
            cterm1,
            cterm2
        );

        arma::cx_mat V = D12p*D12p_t + C1p*Y*C1p_t;
        arma::cx_mat V_inv = arma::inv(V);

        ctemp1 = D22p*D22p_t + C2p*Y*C2p_t;
        ctemp2 = cterm4*V_inv*cterm3;
        arma::cx_mat R = (gam_*gam_)*arma::eye<arma::cx_mat>(INPUT_SIZE, INPUT_SIZE)
            - ctemp1 - ctemp2;

        arma::cx_mat G = S + arma::trans(C1_C2)*Y*C1_C2;
        arma::cx_mat G_inv = arma::inv(G);

        //-- TODO : check whether is Acl asymptotically stable or not
        arma::cx_mat Acl = Ap - arma::trans(cterm5)*G_inv*arma::trans(C1_C2);

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

} } // namespace jacl::synthesis
