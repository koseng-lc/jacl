#pragma once

#include <armadillo>

#include <vector>

namespace JACL{

namespace linear_algebra{

static void GramSchmidtProjection(int idx, const arma::mat& col, arma::mat* u, std::vector<arma::mat>* e){

    if(idx > 0){
        double norm = arma::dot(col, (*e)[idx-1]);
        *u -= (norm * (*e)[idx-1]);
        GramSchmidtProjection(idx - 1, col, u, e);
    }
}

static void QRDecomp(const arma::mat& in, arma::mat* Q, arma::mat* R){

    *Q = arma::mat(in.n_rows, in.n_cols, arma::fill::zeros);
    *R = arma::mat(arma::size(in), arma::fill::zeros);
    std::vector<arma::mat> e;

//    int n_rank = arma::rank(in);
    arma::mat u;

    //Gram-Schmidt Process
    auto norm(.0);
    arma::mat normalized_u;
    for(int i(0); i < in.n_cols; i++){
        u = in.col(i);
        GramSchmidtProjection(i, in.col(i), &u, &e);
        norm = arma::norm(u, 2);
        normalized_u = (norm == .0) ? u : u / norm;
        e.push_back(normalized_u);
        Q->col(i) = e[i];
    }

    for(int r(0); r < in.n_cols; r++){ // limit row iteration
        for(int c(r); c < in.n_cols; c++){
            (*R)(r, c) = arma::dot(in.col(c), Q->col(r));
        }
    }
}

static void QRAlgorithm(const arma::mat& in, arma::mat* T, arma::mat* U, int num_iter = 20){

    arma::mat updated_in(in);
    *U = arma::mat(arma::size(in), arma::fill::eye);
    arma::mat Q, R;
    for(int i(0); i < num_iter; i++){
//        std::cout << "TEST1" << std::endl;
//        updated_in.print("UIN : ");
        QRDecomp(updated_in, &Q, &R);
//        Q.print("Q : ");
//        R.print("R : ");
        updated_in = R * Q;
        (*U) = (*U) * Q;
//        std::cout << "TEST2" << std::endl;
    }

    *T = updated_in;
}

}

}
