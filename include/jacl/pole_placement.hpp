/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of Kautsky-Nichols Robust pole placement
*/

#pragma once

#include <vector>
#include <iomanip>
#include <type_traits>

#include <boost/thread.hpp>

#include <jacl/linear_state_space.hpp>

namespace jacl{ namespace pole_placement{

enum class PolePlacementType{
    Controller,
    Observer
};

template <class _StateSpace>
auto KautskyNichols(_StateSpace *_ss, const arma::mat& _poles, arma::mat* _K, PolePlacementType _type = PolePlacementType::Controller) -> void{
    assert(_poles.is_vec());

    arma::mat Q, R;
    arma::mat U0, U1, Z;    

    arma::mat A( _ss->A() );
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

} } // namespace jacl::pole_placement
