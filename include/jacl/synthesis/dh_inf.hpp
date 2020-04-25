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
class DHinf:public ::jacl::system::detail::BaseSystemClient<typename _System::base>{
public:
    template<typename __System = _System,
             typename std::enable_if_t<traits::is_discrete_system<__System>::value, int>* = nullptr>
    DHinf(__System* _sys, double _gam)
        : ss_(::jacl::system::detail::BaseSystemClient<typename _System::base>::ss(_sys))
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

    auto checkCondition1();
    auto checkCondition2();
    auto checkCondition3();
    auto checkCondition4();
    auto checkAllCondition();    

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
    typename _System::StateSpace* ss_;
    static constexpr auto INPUT_SIZE{_System::n_inputs - perturbation_size};
    static constexpr auto OUTPUT_SIZE{_System::n_outputs - performance_size};
    lft::LowerLFT<typename _System::StateSpace,
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
    std::cout << "Assumption 1 : " << std::boolalpha << ctrb << " ; " << obsv << std::endl;
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
    std::cout << "Assumption 2 : " << std::boolalpha << ok << std::endl;
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
    std::cout << "Assumption 3 : " << std::boolalpha << ok << std::endl;
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
    perturbation_size>::checkCondition1(){

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

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto DHinf<_System,
    performance_size,
    perturbation_size>::checkCondition2(){

    return linear_algebra::isPosSemiDefinite(X_inf_);
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto DHinf<_System,
    performance_size,
    perturbation_size>::checkCondition3(){

    return linear_algebra::isPosSemiDefinite(Y_inf_);
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto DHinf<_System,
    performance_size,
    perturbation_size>::checkCondition4(){
    arma::cx_mat temp( X_inf_*Y_inf_ );
    return linear_algebra::spectralRadius( std::move(temp) ) < gam_*gam_;
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto DHinf<_System,
    performance_size,
    perturbation_size>::checkAllCondition(){
    auto cond1 = checkCondition1();
    auto cond2 = checkCondition2();
    auto cond3 = checkCondition3();
    auto cond4 = checkCondition4();
#ifdef DHINF_VERBOSE
    std::cout << "Condition 1 : " << std::boolalpha << cond1 << std::endl;
    std::cout << "Condition 2 : " << std::boolalpha << cond2 << std::endl;
    std::cout << "Condition 3 : " << std::boolalpha << cond3 << std::endl;
    std::cout << "Condition 4 : " << std::boolalpha << cond4 << std::endl;
#endif
    return cond1 & cond2 & cond3 & cond4;
}

template <typename _System,
          std::size_t performance_size,
          std::size_t perturbation_size>
auto DHinf<_System,
    performance_size,
    perturbation_size>::solve() -> ::jacl::common::StateSpacePack{        

#ifdef DHINF_VERBOSE
    std::cout << ">>>>> Interconnection matrix assumptions : " << std::endl;
#endif
    bool check_assumption = checkAllAssumption();
    assert(check_assumption && "The assumption made for interconnection matrix is not fulfill !");

    arma::mat temp1, temp2, temp3, temp4, temp5;
    arma::cx_mat ctemp1, ctemp2, ctemp3, ctemp4, ctemp5, ctemp6;
    //-- Setup several frequently matrix operation
    const arma::mat& A = llft_.A();
    arma::mat A_t = arma::trans(A);
    const arma::mat& B = ss_->B();
    const arma::mat& B1 = llft_.B1();
    arma::mat B1_t = arma::trans(llft_.B1());
    const arma::mat& B2 = llft_.B2();
    arma::mat B2_t = arma::trans(llft_.B2());
    const arma::mat& C = ss_->C();
    const arma::mat& C1 = llft_.C1();
    arma::mat C1_t = arma::trans(llft_.C1());
    const arma::mat& C2 = llft_.C2();
    arma::mat C2_t = arma::trans(llft_.C2());
    const arma::mat& D11 = llft_.D11();
    arma::mat D11_t = arma::trans(llft_.D11());
    const arma::mat& D12 = llft_.D12();
    arma::mat D12_t = arma::trans(llft_.D12());

    {
        using namespace ::jacl::linear_algebra;
        arma::mat S = arma::join_cols(
            arma::join_rows(D12_t*D12, D12_t*D11),
            arma::join_rows(D11_t*D12, D11_t*D11 - (gam_*gam_)*arma::eye(perturbation_size,perturbation_size))
        );
        // S.print("S : ");
        arma::mat D12_D11 = arma::join_rows(D12, D11);
        arma::mat D12_D11_t = arma::trans(D12_D11);
        arma::mat B2_B1 = arma::join_rows(B2, B1);
        arma::mat S_inv = arma::inv(S);
        arma::mat X = B2_B1*S_inv*arma::trans(B2_B1);
        arma::mat X_t = arma::trans(X);
        std::cout << "1" << std::endl;
        std::cout << "2" << std::endl;
        temp1 = B2_B1*S_inv*D12_D11_t;
        arma::mat A_tilde = A - (temp1*C1);
        if(arma::cond(A_tilde) > 1e6){
            //-- regularize the ill-conditioned matrix
            A_tilde = A + .001*arma::eye(arma::size(A));
            std::cout << "[DHinf] Regularization triggered" << std::endl;
        }
        std::cout << "3" << std::endl;
        //-- more stable with solve() than with inv()
        arma::mat A_tilde_tinv = arma::solve(A_tilde.t(), arma::eye(arma::size(A)), arma::solve_opts::equilibrate);
        std::cout << "4" << std::endl;     
        temp1 = D12_D11*S_inv*D12_D11_t;
        temp2 = C1_t*temp1*C1;
        std::cout << "5" << std::endl;
        arma::mat Q = C1_t*C1 - temp2;
        temp1 = A_tilde_tinv*Q;
        arma::mat H = arma::join_cols(
            arma::join_rows(A_tilde+X*temp1, -X*A_tilde_tinv),
            arma::join_rows(-temp1, A_tilde_tinv)
        );
        //-- test
        arma::mat J = arma::join_cols(
            arma::join_rows(arma::zeros(arma::size(A)), -arma::eye(arma::size(A))),
            arma::join_rows(arma::eye(arma::size(A)), arma::zeros(arma::size(A)))
        );
        //--
        arma::mat check = H.t()*J*H;
        check.print("Check : ");
        DARE<typename _System::StateSpace> solver(ss_);
        solver.setSympleticMatrix(H);
        arma::cx_mat P = solver.solve();
        arma::cx_mat P_t = arma::trans(P);
        
        std::cout << "9" << std::endl;
        ctemp1 = toCx(B2_t)*P*toCx(B2);
        ctemp2 = toCx(D12_t)*toCx(D12);
        std::cout << "10" << std::endl;
        arma::cx_mat V = ctemp1 + ctemp2;
        arma::cx_mat V_inv = arma::inv(V);
        std::cout << "11" << std::endl;
        ctemp1 = toCx(D11_t)*toCx(D11);
        ctemp2 = toCx(B1_t)*P*toCx(B1);
        ctemp3 = toCx(B1_t)*P*toCx(B2) + toCx(D11_t)*toCx(D12);
        ctemp4 = toCx(B2_t)*P*toCx(B1) + toCx(D12_t)*toCx(D11);
        ctemp5 = ctemp3*V_inv*ctemp4;
        std::cout << "12" << std::endl;
        arma::cx_mat R = (gam_*gam_)*arma::eye<arma::cx_mat>(perturbation_size, perturbation_size)
            - ctemp1 - ctemp2 + ctemp5;

        ctemp1 = toCx(S);
        temp2 = arma::join_rows(
            B2, B1
        );
        ctemp2 = toCx(temp2);
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

        arma::cx_mat cx_X = toCx(X);
        arma::cx_mat cx_A = toCx(A_tilde);
        arma::cx_mat cx_A_t = arma::trans(cx_A);
        ctemp1 = cx_A_t*P_t*cx_A;
        ctemp2 = ctemp1 + toCx(C1_t*C1);
        ctemp3 = cx_A_t*P*cx_X;
        ctemp4 = arma::inv(arma::eye(arma::size(P*cx_X)) + P*cx_X);
        ctemp5 = P*cx_A;
        ctemp6 = ctemp3*ctemp4*ctemp5;
        arma::cx_mat dare_rhs = ctemp2 - ctemp6 - P;
        dare_rhs.print("Dare RHS : ");
    }

    {
        
    }

    return std::make_tuple(A, B1, C1, D11);
}

} } // namespace jacl::synthesis
