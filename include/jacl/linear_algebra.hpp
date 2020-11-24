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

using default_scalar_t = double;
using size_t = std::size_t;

template <size_t rows, size_t cols, typename Scalar = default_scalar_t>
class Matrix{
public:
    using scalar_t = Scalar;
    using ref_t = scalar_t&;
    using const_ref_t = const scalar_t&;
    using move_t = scalar_t&&;
    using transpose_t = Matrix<cols, rows, Scalar>;

    static constexpr auto n_rows{rows};
    static constexpr auto n_cols{cols};
    static constexpr auto n_elems{n_rows*n_cols};

    using data_t = std::array<scalar_t, n_rows*n_cols>;

public:

    template <typename ...Es>
    constexpr Matrix(Es... _elements)
        : data_{{std::forward<Es>(_elements)...}}{
        
    }

    constexpr Matrix()
        : data_{}{
        
    }

    //-- Does the copy and the move constructor are elided?
    /*template <typename _Matrix, typename Idx=std::make_index_sequence<n_elems>>
    constexpr Matrix(const _Matrix& _m){
        // std::cout << "[Matrix] Copy constructor !" << std::endl;
        static_assert(n_rows == _Matrix::n_rows
                        & n_cols == _Matrix::n_cols,
                        "[Matrix] Constructing must be in same dimension");
        copy_impl(_m, Idx{});        
    }*/

    /*template <typename _Matrix, typename Idx=std::make_index_sequence<n_elems>>
    constexpr Matrix(_Matrix&& _m){
        // std::cout << "[Matrix] Move constructor !" << std::endl;
        static_assert(n_rows == _Matrix::n_rows
                        & n_cols == _Matrix::n_cols,
                        "[Matrix] Constructing must be in same dimension");
        move_impl(_m, Idx{});        
    }*/

    template <typename _Matrix>
    constexpr auto operator *(_Matrix&& _m) const{        
        return mul_impl(std::forward<_Matrix>(_m));
    }

    template <typename _Matrix, typename Idx = std::make_index_sequence<n_elems>>
    constexpr auto operator +(_Matrix&& _m) const{
        return add_impl(std::forward<_Matrix>(_m), Idx{});
    }

    template <typename _Matrix, typename Idx = std::make_index_sequence<n_elems>>
    constexpr auto operator -(_Matrix&& _m) const{
        return sub_impl(std::forward<_Matrix>(_m), Idx{});
    }

    template <typename _Matrix, typename Idx = std::make_index_sequence<n_elems>>
    constexpr Matrix& operator =(const _Matrix& _m){
        static_assert(n_rows == _Matrix::n_rows
                        & n_cols == _Matrix::n_cols,
                        "[Matrix] Assignment must be in the same dimension");
        std::cout << "[Matrix] Copy assignment !" << std::endl;
        copy_impl(_m, Idx{});        
        return *this;
    }

    template <typename _Matrix, typename Idx = std::make_index_sequence<n_elems>>
    constexpr Matrix& operator =(_Matrix&& _m){
        static_assert(n_rows == _Matrix::n_rows
                        & n_cols == _Matrix::n_cols,
                        "[Matrix] Assignment must be in the same dimension");
        std::cout << "[Matrix] Move assignment !" << std::endl;
        move_impl(_m, Idx{});
        return *this;
    }

    constexpr ref_t operator ()(size_t _idx){
        return data_[_idx];
    }

    constexpr ref_t operator ()(size_t _row, size_t _col){
        return data_[_row*n_cols+_col];
    }

    constexpr const_ref_t operator ()(size_t _idx) const{
        return data_[_idx];
    }

    constexpr const_ref_t operator ()(size_t _row, size_t _col) const{
        return data_[_row*n_cols+_col];
    }

    template <size_t idx>
    constexpr move_t operator ()(size_t a=0){
        return std::get<idx>(std::move(data_));
    }

    template <typename Idx=std::make_index_sequence<n_elems>>
    constexpr auto zeros(){
        zeros_impl(Idx{});
    }

private:
    template <typename _Matrix, size_t... Is>
    constexpr auto copy_impl(const _Matrix& _m, std::index_sequence<Is...>){
        (((*this)(Is) = _m(Is)), ...);
    }

    template <typename _Matrix, size_t... Is>
    constexpr auto move_impl(_Matrix&& _m, std::index_sequence<Is...>){
        (((*this)(Is) = std::move(_m.template operator()<Is>(0))), ...);
    }

    template <size_t start, size_t end, typename _Scalar>
    constexpr auto fill(_Scalar s){
        if constexpr(start < end){
            return fill_impl<start+1, end>(s, s);
        }else{
            return {s};
        }
    }

    template <size_t start, size_t end, typename _Scalar, typename ...Rest>
    constexpr auto fill_impl(Scalar s, Rest... rest){
        if constexpr(start < end){
            return fill_impl<start+1, end>(s, s, rest...);
        }else{
            return {s,rest...};
        }
    }

    template <typename _Matrix, size_t... Is>
    constexpr auto add_impl(const _Matrix& _m, std::index_sequence<Is...>) const{
        Matrix<n_rows, n_cols> res{};
        res.zeros();
        ((res(Is) = data_[Is] + _m(Is)), ...);
        return res;
    }

    template <typename _Matrix, size_t... Is>
    constexpr auto sub_impl(const _Matrix& _m, std::index_sequence<Is...>) const{
        Matrix<n_rows, n_cols> res{};
        res.zeros();
        ((res(Is) = data_[Is] - _m(Is)), ...);
        return res;
    }

    template <typename _Matrix>
    constexpr auto mul_impl(const _Matrix& _m) const{
        // std::cout << "[Matrix] Copy !" << std::endl;
        Matrix<n_rows, _Matrix::n_cols> res{};
        res.zeros();
        for(size_t r(0); r < n_rows; r++){
            for(size_t c(0); c < _Matrix::n_cols; c++){
                for(size_t i(0); i < n_cols; i++){
                        res(r,c) += (*this)(r,i)*_m(i,c);
                }
            }
        }
        return res;       
    }

    template <typename _Matrix>
    constexpr auto mul_impl(_Matrix&& _m) const{
        // std::cout << "[Matrix] Move !" << std::endl;
        Matrix<n_rows, _Matrix::n_cols> res{};
        res.zeros();
        for(size_t r(0); r < n_rows; r++){
            for(size_t c(0); c < _Matrix::n_cols; c++){
                for(size_t i(0); i < n_cols; i++){
                        res(r,c) += (*this)(r,i)*_m(i,c);
                }
            }
        }
        return res;       
    }

    template <size_t... Is>
    constexpr auto zeros_impl(std::index_sequence<Is...>){
        ((data_[Is] = 0.), ...);
    }

private:
    data_t data_;
};

namespace detail{

    template <typename _Matrix>
    inline auto toCx(const _Matrix& m){
        return arma::cx_mat(m, arma::zeros<arma::mat>( arma::size(m) ));
    }

    template <typename _Matrix>
    inline auto toCx(_Matrix&& m){   
        return arma::cx_mat(m, arma::zeros<arma::mat>( arma::size(m) ));
    }

    template <typename _Matrix>
    inline auto toReal(const _Matrix& m){
        return arma::real(m);
    }

    template <typename _Matrix>
    inline auto toReal(_Matrix&& m){
        return arma::real(m);
    }

    template <std::size_t rows, std::size_t cols, typename Scalar>
    constexpr auto diagonalize(Scalar s) -> typename arma::Mat<Scalar>::template fixed<rows, cols>{
        return s*typename arma::Mat<Scalar>::template fixed<rows, cols>(arma::fill::eye);
    }

    template <typename T>
    constexpr auto abs(T t){
        if constexpr(traits::is_complex_v<T>){
            return std::sqrt(std::real(t)*std::real(t) + std::imag(t)*std::imag(t));
        }else{
            if constexpr(t<0){
                return -1*t;
            }else{
                return t;
            }
        }
    }

    template <typename Scalar = default_scalar_t>
    constexpr auto _sgn(Scalar s){
        if constexpr(s < Scalar(.0)) return Scalar(-1.0);
        else return Scalar(1.0);
    }

    template <typename Scalar = default_scalar_t>
    constexpr auto _round(Scalar s){
        if constexpr(s >= 0.) return int(s + .5);
        else return int(s - .5);
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

    template <typename _Matrix, std::size_t... Is>
    auto show_impl(const _Matrix& _m, std::index_sequence<Is...>){
        std::cout << std::setprecision(6);
        ((std::cout << (((_Matrix::n_cols- Is%_Matrix::n_cols) == 1)
                        ?(std::to_string(_m(Is))+"\n"):(std::to_string(_m(Is))+" "))), ...);
        std::cout << std::endl;
    }

    template <std::size_t _first_num, std::size_t _step, std::size_t _num_step, std::size_t... _Is>
    constexpr auto gen_specific_index_sequence(){
        if constexpr (_num_step == 0) return std::index_sequence<_Is...>{};
        else return gen_specific_index_sequence<_first_num-_step, _step, _num_step-1, _first_num-_step, _Is...>();
    }

    template <typename _Matrix, size_t... Is>
    constexpr auto getCol_impl(const _Matrix& _m, std::index_sequence<Is...>){
        return Matrix<_Matrix::n_rows, 1>(_m(Is)...);
    }

    template <typename _Matrix, size_t... Is>
    constexpr auto getRow_impl(const _Matrix& _m, std::index_sequence<Is...>){
        return Matrix<1, _Matrix::n_cols>(_m(Is)...);
    }

    template <typename _Matrix, std::size_t... Is>
    constexpr auto transpose_impl(const _Matrix& _m, std::index_sequence<Is...>){
        constexpr auto n_rows{_Matrix::n_rows};
        constexpr auto n_cols{_Matrix::n_cols};
        return typename _Matrix::transpose_t(_m((Is%n_rows)*n_cols + (Is/n_rows))...);
    }
}

template <std::size_t _first_num, std::size_t _step, std::size_t _num_step>
using make_specific_index_sequence = std::decay_t<
        decltype(detail::gen_specific_index_sequence<_first_num+_step, _step, _num_step>())>;

template <std::size_t col, typename _Matrix,
          typename Idx = make_specific_index_sequence<
                            _Matrix::n_elems-_Matrix::n_cols+(col%_Matrix::n_cols),
                            _Matrix::n_cols,_Matrix::n_rows>>
constexpr auto getCol(const _Matrix& _m){
    static_assert(col < _Matrix::n_cols, "[linear_algebra] The column is nothing!");
    return detail::getCol_impl(_m, Idx{});
}

template <std::size_t row, typename _Matrix,
          typename Idx = make_specific_index_sequence<
                            _Matrix::n_cols*(row+1)-1,
                            1,_Matrix::n_cols>>
constexpr auto getRow(const _Matrix& _m){
    static_assert(row < _Matrix::n_rows, "[linear_algebra] The row is nothing!");
    return detail::getRow_impl(_m, Idx{});
}

template <typename _Matrix, typename Idx = std::make_index_sequence<_Matrix::n_elems>>
constexpr auto transpose(const _Matrix& _m){
    return detail::transpose_impl(_m, Idx{});
}

template <typename _Matrix, typename Idx = std::make_index_sequence<_Matrix::n_elems>>
auto show(const _Matrix& _m, std::string&& _alias){
    std::cout << _alias << std::endl;
    detail::show_impl(_m, Idx{});
}

template <typename _Matrix>
inline auto toCx(_Matrix&& m){

    static_assert(arma::is_arma_type<std::decay_t<_Matrix>>::value, "Please input arma_type only!");

    return detail::toCx(std::forward<_Matrix>(m));
}

template <typename _Matrix>
inline auto toReal(_Matrix&& m){

    static_assert(arma::is_arma_type<std::decay_t<_Matrix>>::value, "Please input arma_type only!");
    
    return detail::toReal(std::forward<_Matrix>(m));
}

template <std::size_t rows, std::size_t cols=rows, typename Scalar>
constexpr auto diagonalize(Scalar s){
    return detail::diagonalize<rows,cols>(s);
}

template <typename T>
constexpr auto abs(T t){
    return detail::abs(t);
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

template <typename _Matrix>
static auto isPosSemiDefinite(_Matrix&& in){

    static_assert(arma::is_arma_type<std::decay_t<_Matrix>>::value, "Please input arma_type only!");

    return detail::isPosSemiDefinite(std::forward<_Matrix>(in));
}

template <typename _Matrix>
static auto spectralRadius(_Matrix&& in){

    static_assert(arma::is_arma_type<std::decay_t<_Matrix>>::value, "Please input arma_type only!");

    return detail::spectralRadius(std::forward<_Matrix>(in));
}

template <typename _Matrix>
static auto largestSV(_Matrix&& in){

    static_assert(arma::is_arma_type<std::decay_t<_Matrix>>::value, "Please input arma_type only!");

    return detail::largestSV(std::forward<_Matrix>(in));
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
