/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of pole placement methods
*/

#pragma once

#include <vector>
#include <iomanip>
#include <type_traits>

#include <boost/thread.hpp>

#include <jacl/linear_state_space.hpp>
#include <jacl/lti_common.hpp>

namespace jacl{ namespace pole_placement{

namespace detail{
    template <typename T>
    static auto VietaFormula(const arma::Col<T> _roots, arma::Col<T>* _coeff){
        *_coeff = arma::Col<T>(_roots.n_rows + 1, _roots.n_cols, arma::fill::zeros);
        (*_coeff)(0) = 1. + .0i;
        (*_coeff)(1) = _roots(0);
        std::complex<double> old_neighbor_val((*_coeff)(0));
        std::complex<double> temp(.0+.0i);
        for(int j(1); j < (_coeff->n_rows - 1); j++){
            for(int k(1); k <= j; k++){
                temp = (*_coeff)(k);
                (*_coeff)(k) = std::pow(-1, k)*(temp + old_neighbor_val*_roots(j));
                old_neighbor_val = temp;
            }
        }
        (*_coeff)(_coeff->n_rows - 1) = 1.;
        for(const auto& r:_roots)
            (*_coeff)(_coeff->n_rows - 1) *= r;
    }
}

enum class PolePlacementType{
    Controller,
    Observer
};

template <class _StateSpace>
static auto KautskyNichols(_StateSpace *_ss, const arma::mat& _poles, arma::mat* _K, PolePlacementType _type = PolePlacementType::Controller){
    assert(_poles.is_vec());

    arma::mat Q, R;
    arma::mat U0, U1, Z;    

    typename _StateSpace::state_matrix_t A( _ss->A() );
    arma::mat T; //type
    if(_type == PolePlacementType::Controller){
        T = _ss->B();
    }else if(_type == PolePlacementType::Observer){
        //-- tranpose won't change the eigenvalue
        A = arma::trans(A);
        T = arma::trans(_ss->C());
    }else{
        assert("Please choose appropriate type !");
    }

    int T_rows(T.n_rows);
    int T_cols(T.n_cols);

    //-- special treatment
    //-- due to the one of the span from orthogonal by QR of T is empty
    //-- with assumption that T is invertible
    //-- we use the original A, B or C from _ss
    if(T_rows == T_cols){
        arma::mat P(arma::diagmat(_poles));
        if(_type == PolePlacementType::Controller){
            *_K = arma::inv(_ss->B()) * (_ss->A() - P);
        }else if(_type == PolePlacementType::Observer){
            *_K = (_ss->A() - P) * arma::inv(_ss->C());
        }        
        return;
    }

    arma::qr(Q, R, T);         
    U0 = Q.head_cols(T_cols);
    U1 = Q.tail_cols( T_rows-T_cols);
    Z = R.head_rows(T_cols);

    //-- A Step
    std::vector<arma::mat> S_hat;
    std::vector<arma::mat> S;
    std::vector<arma::mat> another_R;

    arma::mat temp1, temp2, temp3, temp4;
    arma::mat temp5, temp6;


    for(auto p:_poles){
        temp1 = U1.t(); // place it in outer scope
        temp2 = arma::eye(arma::size(A)) * p;
        temp3 = A - temp2;
        temp4 = arma::trans(temp1 * temp3);

        int rows(temp4.n_rows);
        int cols(temp4.n_cols);

        arma::qr(temp5, temp6, temp4);

        S_hat.emplace_back(temp5.head_cols(cols));
        S.emplace_back(temp5.tail_cols(rows-cols));
        another_R.emplace_back(temp6.head_rows(cols));
    }

    //-- X Step

    // initialize X matrix
    arma::mat X(arma::size(A), arma::fill::randu);
    while(arma::rank(X) != X.n_rows)
        X = arma::randu<arma::mat>(arma::size(A));

    arma::mat X_rest;

    arma::mat Q_tilde;
    arma::mat y_tilde;
    arma::mat R_tilde;
    arma::mat Q1, R1;

    auto cond_number_tol(1e-3);
    auto norm(1e-6);
    auto cond(.0);
    auto prev_cond(cond);
    auto err(.0);

    int restX_span_cols = X.n_cols - 1;
    int restX_span_rows = X.n_rows;    

    do{
        for(int i(0); i < _poles.n_elem; i++){
            X_rest = X;
            X_rest.shed_col(i);

            arma::qr(Q1, R1, X_rest);

            Q_tilde = Q1.head_cols(restX_span_cols);
            y_tilde = Q1.tail_cols(restX_span_rows-restX_span_cols);
            R_tilde = R1.head_rows(restX_span_cols);

            // reuse the temporary variable
            temp1 = arma::trans(S[i]) * y_tilde;
            if(!arma::approx_equal(temp1, arma::zeros<arma::mat>(arma::size(temp1)),"absdiff",1e-6)){
                // calculate 2-norm
                norm = arma::norm(temp1, 2);    
                X.col(i) = (S[i] * temp1) / norm;
            }            
        }
        cond = arma::cond(X);
        err = std::fabs(cond - prev_cond);
        prev_cond = cond;
    }while(err > cond_number_tol); // iterate until well-conditioned

    //-- F Step

    arma::mat eigenvalue_diag(arma::diagmat(_poles));
    temp1 = arma::inv(X);
    temp2 = eigenvalue_diag * temp1;
    arma::mat M(X * temp2);

    temp3 = U0.t();
    temp4 = M - A;
    temp5 = temp3 * temp4;
    //-- because the original form on the paper is A + BK
    //-- but in general people will use form of A - BK
    *_K = -1 * arma::inv(Z) * temp5;
    if(_type == PolePlacementType::Observer)
        *_K = arma::trans(*_K);
}

template <typename _StateSpace>
static auto BassGura(_StateSpace* _ss, const arma::vec& _poles, arma::mat* _gain){
    const typename _StateSpace::state_matrix_t& A = _ss->A();
    const typename _StateSpace::output_matrix_t& C = _ss->C();
    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eig_gen(eigval, eigvec, A);
    arma::cx_vec a, a_hat;    
    detail::VietaFormula(eigval, &a);
    detail::VietaFormula(arma::conv_to<arma::cx_vec>::from(_poles), &a_hat);
    arma::vec a_real = arma::real(a); a_real.shed_row(0);
    arma::vec a_hat_real = arma::real(a_hat); a_hat_real.shed_row(0);
    arma::mat T(a_real.n_rows, a_real.n_rows, arma::fill::eye);
    for(int i(0),j(1); j < T.n_rows; j++,i++){
        T.submat(j,i,T.n_rows-1,i) = a_real.head_rows(T.n_rows - j);
    }
    arma::mat O;
    bool ok = common::observable(A,C,&O);
    if(!ok)
        std::cerr << "The system was not observable" << std::endl;
    arma::mat term1, term2;
    term1 = arma::inv(T*O);
    term2 = a_hat_real - a_real;
    *_gain = term1 * term2;  
}

} } // namespace jacl::pole_placement
