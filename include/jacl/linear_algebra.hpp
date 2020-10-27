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

    template <typename Matrix>
    static auto toCx(const Matrix& m){
        return arma::cx_mat(m, arma::zeros<arma::mat>( arma::size(m) ));
    }

    template <typename Matrix>
    static auto toCx(Matrix&& m){        
        return arma::cx_mat(m, arma::zeros<arma::mat>( arma::size(m) ));
    }

    template <typename Matrix>
    static inline auto toReal(const Matrix& m){
        return arma::real(m);
    }

    template <typename Matrix>
    static inline auto toReal(Matrix&& m){
        return arma::real(m);
    }

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

    template <std::size_t idx, typename C=arma::mat, typename U=arma::mat, typename E>
    static auto GramSchmidtProjection_iter(const C& col, U* u, const E& e){

        if constexpr(idx > 0){
            auto norm = arma::dot(col, e[idx-1]);
            *u -= (norm * e[idx-1]);
            GramSchmidtProjection_iter<idx-1>(col, u, e);
        }
    }

    template <std::size_t idx, typename I, typename Q, typename E>
    static auto GramSchmidtProjection(const I& in, Q* q, E* e){
        
        if constexpr(idx < I::n_cols){
            arma::mat u = in.col(idx);
            GramSchmidtProjection_iter<idx>(in.col(idx), &u, *e);
            auto norm = arma::norm(u, 2);
            arma::mat normalized_u = (norm == .0) ? u : u / norm;
            (*e)[idx] = normalized_u;
            q->col(idx) = (*e)[idx];
            GramSchmidtProjection<idx+1>(in, q, e);
        }
    }

    template <std::size_t col, std::size_t row, typename I, typename Q, typename R>
    static auto QRDecomp_col_iter(const I& in, const Q& q, R* r){

        if constexpr(col < I::n_cols){
            (*r)(row, col) = arma::dot(in.col(col), q.col(row));
            QRDecomp_col_iter<col+1,row>(in, q, r);
        }
    }

    template <std::size_t row, typename I, typename Q, typename R>
    static auto QRDecomp_row_iter(const I& in, const Q& q, R* r){

        if constexpr(row < I::n_cols){//-- set limit row iteration into cols size
            QRDecomp_col_iter<row,row>(in, q, r);
            QRDecomp_row_iter<row+1>(in, q, r);
        }
    }
}

template <typename Matrix>
static inline auto toCx(Matrix&& m){

    static_assert(arma::is_arma_type<Matrix>::value, "Please input arma_type only!");

    return detail::toCx(std::forward<Matrix>(m));
}

template <typename Matrix>
static inline auto toReal(Matrix&& m){

    static_assert(arma::is_arma_type<Matrix>::value, "Please input arma_type only!");
    
    return detail::toReal(std::forward<Matrix>(m));
}

template <typename I, typename Q, typename R>
static auto QRDecomp(const I& in, Q* q, R* r){

    static_assert(arma::is_arma_type<I>::value & arma::is_arma_type<Q>::value & arma::is_arma_type<R>::value, "Please input arma_type only!");

    *q = arma::mat::fixed<I::n_rows, I::n_cols>(arma::fill::zeros);
    *r = arma::mat::fixed<I::n_rows, I::n_cols>(arma::fill::zeros);
    std::array<arma::mat::fixed<I::n_rows,1>, I::n_cols> e;

    //-- Gram-Schmidt process
    detail::GramSchmidtProjection<0>(in, q, &e);

    detail::QRDecomp_row_iter<0>(in, *q, r);
}

//-- using universal reference

template <typename Matrix>
static auto isPosSemiDefinite(Matrix&& in){

    static_assert(arma::is_arma_type<Matrix>::value, "Please input arma_type only!");

    return detail::isPosSemiDefinite(std::forward<Matrix>(in));
}

template <typename Matrix>
static auto spectralRadius(Matrix&& in){

    static_assert(arma::is_arma_type<Matrix>::value, "Please input arma_type only!");

    return detail::spectralRadius(std::forward<Matrix>(in));
}

template <typename Matrix>
static auto largestSV(Matrix&& in){

    static_assert(arma::is_arma_type<Matrix>::value, "Please input arma_type only!");

    return detail::largestSV(std::forward<Matrix>(in));
}

} // namespace jacl::linear_algebra
