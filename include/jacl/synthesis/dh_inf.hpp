/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of discrete H-infinity synthesis.
*            A lot of code here, is constructed by copy constructor because it's easier to read.
*/

#pragma once

#include <jacl/lti_common.hpp>
#include <jacl/are.hpp>
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

    inline auto toCx(const arma::mat& _in) const
        -> decltype(arma::cx_mat(_in, arma::zeros<arma::mat>( arma::size(_in) ))){
        return arma::cx_mat(_in, arma::zeros<arma::mat>( arma::size(_in) ));
    }

    inline auto toReal(const arma::cx_mat& _in) const -> decltype(arma::real(_in)){
        return arma::real(_in);
    }

    template <typename _StateSpace>
    auto isInfNormLessThan(double _gam, const _StateSpace& _ss){
        arma::mat C_t = arma::trans( _ss.C() );
        arma::mat D_t = arma::trans( _ss.D() );
        arma::mat D_tD = D_t * _ss.D();
        arma::mat R = (_gam*_gam)*arma::eye(_StateSpace::n_inputs,
                                            _StateSpace::n_inputs)
                        - D_tD;
        arma::mat R_inv = arma::inv(R);

        //-- Hamiltonian matrix
        arma::mat H(_StateSpace::n_states << 1,
                    _StateSpace::n_states << 1);

        arma::mat temp1, temp2;

        //-- block 1,1
        temp1 = R_inv*D_t*_ss.C();
        temp2 = _ss.A() + temp1;
        H.submat(                         0,                          0,
                 _StateSpace::n_states - 1, _StateSpace::n_states - 1) = temp2;
        //-- block 2,2
        H.submat(_StateSpace::n_states, _StateSpace::n_states, _StateSpace::n_states * 2 - 1, _StateSpace::n_states * 2 - 1) = -arma::trans(temp2);
        //-- block 1,2
        H.submat(                         0,         _StateSpace::n_states,
                 _StateSpace::n_states - 1, _StateSpace::n_states * 2 - 1) = _ss.B()*R_inv*arma::trans(_ss.B());
        //-- block 2,1
        temp1 = _ss.D()*R_inv*D_t;
        temp2 = arma::eye(_StateSpace::n_outputs,
                          _StateSpace::n_outputs)
                 - temp1;
        H.submat(_StateSpace::n_states, 0, _StateSpace::n_states * 2 - 1, _StateSpace::n_states - 1) = -C_t*temp2*_ss.C();

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
    //-- Lemma 3.7 Essentials Robust Control
    auto invariantZeros(){

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
    arma::cx_mat ctemp1, ctemp2, ctemp3;
    //-- New code here
    //-- Setup several frequently matrix operation
    const arma::mat& A = llft_.A();
    const arma::mat& B1 = llft_.B1();
    arma::mat B1_t = arma::trans(llft_.B1());
    const arma::mat& B2 = llft_.B2();
    arma::mat B2_t = arma::trans(llft_.B2());
    const arma::mat& C1 = llft_.C1();
    arma::mat C1_t = arma::trans(llft_.C1());
    const arma::mat& C2 = llft_.C2();
    arma::mat C2_t = arma::trans(llft_.C2());
    const arma::mat& D11 = llft_.D11();
    arma::mat D11_t = arma::trans(llft_.D11());
    const arma::mat& D12 = llft_.D12();
    arma::mat D12_t = arma::trans(llft_.D12());
    //-- there exist symmetric matrix P that semi-positive definite
    arma::mat P(arma::size(A), arma::fill::zeros);
    temp1 = B2_t*P*B2;
    temp2 = D12_t*D12;
    arma::mat V = temp1 + temp2;
    arma::mat V_inv = arma::inv(V);
    temp1 = D11_t*D11;
    temp2 = B1_t*P*B1;
    temp3 = B1_t*P*B2 + D11_t*D12;
    temp4 = B2_t*P*B1 + D12_t*D11;
    temp5 = temp3*V_inv*temp4;
    arma::mat R = arma::eye(perturbation_size, perturbation_size)
        - temp1 - temp2 + temp5;
    //-- INPUT_SIZE + PERTURBATION_SIZE
    temp1 = arma::mat(INPUT_SIZE+perturbation_size, INPUT_SIZE+perturbation_size);

    arma::mat G_P(INPUT_SIZE+perturbation_size, INPUT_SIZE+perturbation_size);
    G_P(0,0,INPUT_SIZE-1,INPUT_SIZE-1) = D12_t*D12;
    G_P(0,INPUT_SIZE,INPUT_SIZE-1,INPUT_SIZE+perturbation_size-1) = D12_t*D11;
    G_P(INPUT_SIZE,0,INPUT_SIZE+perturbation_size-1,INPUT_SIZE-1) = D11_t*D12;
    G_P(INPUT_SIZE,INPUT_SIZE,
        INPUT_SIZE+perturbation_size-1,
        INPUT_SIZE+perturbation_size-1) = D11_t*D11 - arma::eye(perturbation_size,perturbation_size);


    return std::make_tuple(A, B1, C1, D11);
}

} } // namespace jacl::synthesis
