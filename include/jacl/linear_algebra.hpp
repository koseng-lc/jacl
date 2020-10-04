/**
*   @author : koseng (Lintang)
*   @brief : Some linear algebra tools
*/

#pragma once

#include <vector>

#define ARMA_DONT_USE_WRAPPER
#include <armadillo>

namespace jacl::linear_algebra{

namespace detail{

    template <typename Type>
    static auto isPosSemiDefinite(const arma::Mat<Type>& in){
        arma::Mat<Type> hessian_in = .5*(in + arma::trans(in));
        arma::Mat<Type> R;
        return arma::chol(R, hessian_in);
    }

    template <typename Type>
    static auto isPosSemiDefinite(arma::Mat<Type>&& in){
        arma::Mat<Type> hessian_in = .5*(in + arma::trans(in));
        arma::Mat<Type> R;
        return arma::chol(R, hessian_in);
    }

    template <typename Type>
    static auto spectralRadius(const arma::Mat<Type>& in){
        arma::Mat<Type> eigvec;
        arma::Col<Type> eigval;
        arma::eig_gen(eigval, eigvec, in);
        auto max_mag_eigval(.0);
        auto mag_eigval(.0);
        for(int i(0); i < eigval.n_rows; i++){
            mag_eigval = std::abs(eigval(i));
            if(mag_eigval > max_mag_eigval)
                max_mag_eigval = mag_eigval;
        }
        return max_mag_eigval;
    }

    template <typename Type>
    static auto spectralRadius(arma::Mat<Type>&& in){
        arma::Mat<Type> eigvec;
        arma::Col<Type> eigval;
        arma::eig_gen(eigval, eigvec, in);
        auto max_mag_eigval(.0);
        auto mag_eigval(.0);
        for(int i(0); i < eigval.n_rows; i++){
            mag_eigval = std::abs(eigval(i));
            if(mag_eigval > max_mag_eigval)
                max_mag_eigval = mag_eigval;
        }
        return max_mag_eigval;
    }

    template <typename Type>
    static auto largestSV(const arma::Mat<Type>& in)
        -> Type{
        arma::Col<Type> s = arma::svd(in);
        return arma::max(s);
    }

    template <typename Type>
    static auto largestSV(arma::Mat<Type>&& in)
        -> Type{
        arma::Col<Type> s = arma::svd(in);
        return arma::max(s);
    }
}

// template <typename Type>
// static inline auto toCx(const arma::Mat<Type>& _in)
//     -> decltype(arma::Mat<std::complex<Type>>(_in, arma::zeros<arma::Mat<Type>>( arma::size(_in) ))){
//     return arma::Mat<std::complex<Type>>(_in, arma::zeros<arma::Mat<Type>>( arma::size(_in) ));
// }

static inline auto toCx(const arma::mat _in){
    return arma::cx_mat(_in, arma::zeros<arma::mat>( arma::size(_in) ));
}

template <typename Type>
static inline auto toReal(const arma::Mat<std::complex<Type>>& _in)
    -> decltype(arma::real(_in)){
    return arma::real(_in);
}

static void GramSchmidtProjection(int idx, const arma::mat& col, arma::mat* u, std::vector<arma::mat>* e){
    if(idx > 0){
        auto norm = arma::dot(col, (*e)[idx-1]);
        *u -= (norm * (*e)[idx-1]);
        GramSchmidtProjection(idx - 1, col, u, e);
    }
}

static auto QRDecomp(const arma::mat& in, arma::mat* Q, arma::mat* R){

    *Q = arma::mat(in.n_rows, in.n_cols, arma::fill::zeros);
    *R = arma::mat(arma::size(in), arma::fill::zeros);
    std::vector<arma::mat> e;

//    int n_rank = arma::rank(in);
    arma::mat u;

    //-- Gram-Schmidt Process
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

    for(int r(0); r < in.n_cols; r++){ //-- limit row iteration
        for(int c(r); c < in.n_cols; c++){
            (*R)(r, c) = arma::dot(in.col(c), Q->col(r));
        }
    }
}

static auto QRAlgorithm(const arma::mat& in, arma::mat* T, arma::mat* U, int num_iter = 20){
    arma::mat updated_in(in);
    *U = arma::mat(arma::size(in), arma::fill::eye);
    arma::mat Q, R;
    for(int i(0); i < num_iter; i++){
//        updated_in.print("UIN : ");
        QRDecomp(updated_in, &Q, &R);
//        Q.print("Q : ");
//        R.print("R : ");
        updated_in = R * Q;
        (*U) = (*U) * Q;
    }

    *T = updated_in;
}

//-- using universal reference

template <typename Matrix>
static auto isPosSemiDefinite(Matrix&& in){
    return detail::isPosSemiDefinite(std::forward<Matrix>(in));
}

template <typename Matrix>
static auto spectralRadius(Matrix&& in){
    return detail::spectralRadius(std::forward<Matrix>(in));
}

template <typename Matrix>
static auto largestSV(Matrix&& in){
    return detail::largestSV(std::forward<Matrix>(in));
}

} // namespace jacl::linear_algebra
