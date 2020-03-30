/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of Kautsky-Nichols Robust pole placement
*/

#pragma once

#include <vector>
#include <iomanip>
#include <type_traits>

#include <boost/thread.hpp>

#include "state_space.hpp"

namespace jacl{

namespace pole_placement{

enum class PolePlacementType{
    Controller,
    Observer
};

template <class _StateSpace>
auto KautskyNichols(_StateSpace *_ss, const arma::mat& _poles, arma::mat* _K, PolePlacementType _type = PolePlacementType::Controller) -> void{

    assert(_poles.is_vec());

    arma::mat Q, R;
    arma::mat U0, U1, Z;    

    arma::mat T; //type
    if(_type == PolePlacementType::Controller)
        T = _ss->B();
    else if(_type == PolePlacementType::Observer)
        T = _ss->C();
    else
        assert("Please choose appropriate type !");

    arma::qr(Q, R, T);

    int T_rows(T.n_rows);
    int T_cols(T.n_cols);

    //-- special treatment
    //-- due to the one of the span from orthogonal by QR of T is empty
    //-- with assumption that T is invertible
    if(T_rows == T_cols){
        arma::mat inv_T(arma::inv(T));
        arma::mat P(arma::diagmat(_poles));
        *_K = (_ss->A() - P) * inv_T;
        return;
    }

    U0 = Q.submat(0, 0, T_rows-1, T_cols-1);
    U1 = Q.submat(0, T_cols, T_rows-1, T_rows-1);
    Z = R.submat(0, 0, T_cols-1, T_cols-1);


    //-- A Step

    std::vector<arma::mat> S_hat;
    std::vector<arma::mat> S;
    std::vector<arma::mat> another_R;

    arma::mat temp1, temp2, temp3, temp4;
    arma::mat temp5, temp6;


    for(auto p:_poles){

        temp1 = U1.t(); // place it in outer scope
        temp2 = arma::mat(arma::size(_ss->A()), arma::fill::eye) * p;
        temp3 = _ss->A() - temp2;
        temp4 = temp1 * temp3;
        temp4 = temp4.t();

        int rows(temp4.n_rows);
        int cols(temp4.n_cols);

        arma::qr(temp5, temp6, temp4);

        S_hat.push_back(temp5.submat(0, 0, rows-1, cols-1));

        S.push_back(temp5.submat(0, cols, rows-1, rows-1));

        another_R.push_back(temp6.submat(0, 0, cols-1, cols-1));
    }

    //-- X Step

    // initialize X matrix
    arma::mat X(arma::size(_ss->A()), arma::fill::eye);
    arma::mat X_rest;

    arma::mat Q_tilde;
    arma::mat y_tilde;
    arma::mat R_tilde;

    arma::mat Q1, R1;

    auto condition_number_tol(1e-3);
    auto norm(1e-6);
    auto cond(.0);
    auto prev_cond(cond);
    auto err(.0);

    int restX_span_cols = X.n_cols - 1;
    int restX_span_rows = X.n_rows;

    do{

        for(int i(0); i < _poles.n_elem; i++){

            // m x (n-1) size
            X_rest = X;
            X_rest.shed_col(i);

            arma::qr(Q1, R1, X_rest);

            Q_tilde = Q1.submat(0, 0, restX_span_rows-1, restX_span_cols-1);
            y_tilde = Q1.submat(0, restX_span_cols, restX_span_rows-1, restX_span_rows-1);
            R_tilde = R1.submat(0, 0, restX_span_cols-1, restX_span_cols-1);

            // reuse the temporary variable
            temp1 = S[i].t() * y_tilde;
            // calculate 2-norm
            norm = arma::norm(temp1, 2);            
            X.col(i) = (S[i] * temp1) / norm;
        }

        cond = arma::cond(X);
        err = std::fabs(cond - prev_cond);
        prev_cond = cond;

    }while(err > condition_number_tol); // iterate until well-conditioned

    //-- F Step

    arma::mat eigenvalue_diag(arma::diagmat(_poles));
    temp1 = arma::inv(X);
    temp2 = eigenvalue_diag * temp1;
    arma::mat M(X * temp2);

    temp3 = U0.t();
    temp4 = M - _ss->A();
    temp5 = temp3 * temp4;
    *_K = arma::inv(Z) * temp5;

}

}

}
