/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of Kautsky-Nichols Robust pole placement
*/

#pragma once

#include <vector>
#include <iomanip>

#include <boost/thread.hpp>

#include "state_space.h"

namespace JACL{

namespace PolePlacement{

// TODO : change to template of StateSpace
template <class SSpace>
void KautskyNichols(SSpace *_ss, const Mat& _poles, Mat* _K = nullptr){

    assert(_poles.is_vec());

    Mat Q, R;
    Mat U0, U1, Z;

    arma::qr(Q, R, _ss->B());

    int B_rows(_ss->B().n_rows);
    int B_cols(_ss->B().n_cols);

//    Q.print("Q : ");
//    R.print("R : ");

    U0 = Q.submat(0, 0, B_rows-1, B_cols-1);
    U1 = Q.submat(0, B_cols, B_rows-1, B_rows-1);
    Z = R.submat(0, 0, B_cols-1, B_cols-1);

//    U0.print("U0 : ");
//    U1.print("U1 : ");
//    Z.print("Z : ");

    //-- A Step

    std::vector<Mat> S_hat;
    std::vector<Mat> S;
    std::vector<Mat> another_R;

    Mat temp1, temp2, temp3, temp4;
    Mat temp5, temp6;

//    _poles.print("Poles : ");

    for(auto p:_poles){

        temp1 = U1.t();
        temp2 = Mat(arma::size(_ss->A()), arma::fill::eye) * p;
        temp3 = _ss->A() - temp2;
        temp4 = temp1 * temp3;
        temp4 = temp4.t();

//        temp4.print("Temp : ");

        int rows(temp4.n_rows);
        int cols(temp4.n_cols);

        arma::qr(temp5, temp6, temp4);
        S_hat.push_back(temp5.submat(0, 0, rows-1, cols-1));
//        S_hat.back().print("S_hat : ");
        S.push_back(temp5.submat(0, cols, rows-1, rows-1));
//        S.back().print("S : ");
        another_R.push_back(temp6.submat(0, 0, cols-1, cols-1));
    }

    //-- X Step

    // initialize X matrix
    Mat X(arma::size(_ss->A()), arma::fill::eye);
    Mat X_rest;

    Mat Q_tilde;
    Mat y_tilde;
    Mat R_tilde;

    Mat Q1, R1;

    auto condition_number_tol(1e-3);
    auto norm(1e-6);
    auto cond(.0);
    auto prev_cond(cond);
    auto err(.0);

//    X.print("X : ");

    do{

        int rows,cols;

        for(int i(0); i < _poles.n_elem; i++){

            // m x (n-1) size
            X_rest = X;
            X_rest.print("X before rest span : ");
            X_rest.shed_col(i);
            X_rest.print("X after rest span : ");
            rows = X_rest.n_rows;
            cols = X_rest.n_cols;
            arma::qr(Q1, R1, X_rest);
            Q1.print("Q1 : ");
            R1.print("R1 : ");
            Q_tilde = Q1.submat(0, 0, rows-1, cols-1);
            y_tilde = Q1.submat(0, cols, rows-1, rows-1);
            R_tilde = R1.submat(0, 0, cols-1, cols-1);
            Q_tilde.print("Q tilde : ");
            y_tilde.print("y tilde : ");
            R_tilde.print("R tilde : ");
            // reuse the temporary variable
//            S[i].print("S : ");
            temp1 = S[i].t() * y_tilde;
            // calculate 2-norm
            norm = arma::norm(temp1, 2);
            X.col(i) = (S[i] * temp1) / norm;
        }

        X.print("X : ");

        cond = arma::cond(X);
        err = std::fabs(cond - prev_cond);
        prev_cond = cond;
        std::cout << "Error : " << err << std::endl;
        std::cout << std::setprecision(6) << "Condition Number : " << cond << std::endl;

    }while(err > condition_number_tol);

    //-- F Step

    Mat eigenvalue_diag(arma::diagmat(_poles));
    temp1 = arma::inv(X);
    temp2 = eigenvalue_diag * temp1;
    Mat M(X * temp2);

    temp3 = U0.t();
    temp4 = M - _ss->A();
    temp5 = temp3 * temp4;
    *_K = arma::inv(Z) * temp5;

}

}

}
