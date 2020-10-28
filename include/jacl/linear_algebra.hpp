/**
*   @author : koseng (Lintang)
*   @brief : Some linear algebra tools
*/

#pragma once

#include <vector>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <numpy/ndarrayobject.h>

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
        std::cout << "Move it !!!" << std::endl;
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

    static_assert(arma::is_arma_type<std::decay_t<Matrix>>::value, "Please input arma_type only!");

    return detail::toCx(std::forward<Matrix>(m));
}

template <typename Matrix>
static inline auto toReal(Matrix&& m){

    static_assert(arma::is_arma_type<std::decay_t<Matrix>>::value, "Please input arma_type only!");
    
    return detail::toReal(std::forward<Matrix>(m));
}

template <typename I, typename Q, typename R>
static auto QRDecomp(const I& in, Q* q, R* r){

    static_assert(arma::is_arma_type<I>::value
                  & arma::is_arma_type<Q>::value
                  & arma::is_arma_type<R>::value,
                  "Please input arma_type only!");

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

    static_assert(arma::is_arma_type<std::decay_t<Matrix>>::value, "Please input arma_type only!");

    return detail::isPosSemiDefinite(std::forward<Matrix>(in));
}

template <typename Matrix>
static auto spectralRadius(Matrix&& in){

    static_assert(arma::is_arma_type<std::decay_t<Matrix>>::value, "Please input arma_type only!");

    return detail::spectralRadius(std::forward<Matrix>(in));
}

template <typename Matrix>
static auto largestSV(Matrix&& in){

    static_assert(arma::is_arma_type<std::decay_t<Matrix>>::value, "Please input arma_type only!");

    return detail::largestSV(std::forward<Matrix>(in));
}

namespace{ //-- private namespace
        namespace py = boost::python;
        namespace np = boost::python::numpy;
}

static auto auxSchur(arma::mat h, std::tuple<arma::cx_mat, arma::cx_mat>* _TZ, std::string&& ord){
    // static_assert(arma::is_arma_type<H>::value, "Please input arma_type argument!");
    ::jacl::py_stuff::AcquireGIL lk;
    try{
        py::object sys = py::import("sys");
        sys.attr("path").attr("append")("../python");

        py::object aux_schur = py::import("aux_schur");
        Py_intptr_t h_shape[2] = {(int)h.n_rows, (int)h.n_cols};
        PyObject* np_h = PyArray_SimpleNewFromData(2, h_shape, NPY_FLOAT64, reinterpret_cast<void*>(h.memptr()));
        py::handle<> h_handle(np_h);
        py::object h_object(h_handle);

        py::object abc = aux_schur.attr("aux_schur")(h_object, ord);
        py::stl_input_iterator<py::object> begin(abc),end;
        std::vector<py::object> l1(begin, end);
        arma::mat dum[4];
        for(int i(0); i < l1.size(); i++){
            arma::mat temp(arma::size(h));
            std::vector<py::object> layer2(py::stl_input_iterator<py::object>(l1[i]), py::stl_input_iterator<py::object>());
            for(int j(0); j < layer2.size(); j++){
                std::vector<double> layer3(py::stl_input_iterator<double>(layer2[j]), py::stl_input_iterator<double>());
                for(int k(0); k < layer3.size(); k++)
                    temp(j,k) = layer3[k];
            }
            dum[i] = temp;
        }
        arma::cx_mat T(arma::size(h)),Z(arma::size(h));
        T.set_real(dum[0]);
        T.set_imag(dum[1]);
        Z.set_real(dum[2]);
        Z.set_imag(dum[3]);
        *_TZ = std::make_tuple(T,Z);
    }catch(py::error_already_set){
        PyErr_Print();
        return 1;
    }

    return 0;
}

} // namespace jacl::linear_algebra
