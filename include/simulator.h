/**
*   @author : koseng (Lintang)
*   @brief : Simulator
*/

#pragma once

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include <cassert>

#include <armadillo>

namespace JACL{

namespace py = boost::python;
namespace np = boost::python::numpy;

class Simulator{
public:
    Simulator(){

    }

    void retrieveStateSpace(arma::mat _A, arma::mat _B, arma::mat _C, arma::mat _D){

        //-- A must be square matrix and retain it's first shape
        assert(_A.is_square());
        assert(_A.n_rows == n_states_);

        //-- B shape must be equal to it's first shape
        assert(_B.n_rows == n_states_);
        assert(_B.n_cols == n_inputs_);

        //-- C shape must be equal to it's first shape
        assert(_C.n_rows == n_outputs_);
        assert(_C.n_cols == n_states_);

        //-- D shape must be equal to it's first shape
        assert(_D.n_rows == n_outputs_);
        assert(_D.n_cols == n_inputs_);

        //-- A must be square matrix
        assert(_A.is_square());

        //-- num of B rows must be equal to num of states
        assert(_B.n_rows == _A.n_rows);

        //-- num of C cols must be equal to num of states
        assert(_C.n_cols == _A.n_rows);

        //-- num of B cols must be equal to num of D cols
        assert(_B.n_cols == _D.n_cols);

        //-- num of C rows must be equal to num of D rows
        assert(_C.n_rows == _D.n_rows);

        n_states_ = _A.n_rows;
        n_inputs_ = _B.n_cols;
        n_outputs_ = _C.n_rows;

        A_.resize(_A.n_elem);
        B_.resize(_B.n_elem);
        C_.resize(_C.n_elem);
        D_.resize(_D.n_elem);

        A_.clear();
        B_.clear();
        C_.clear();
        D_.clear();

        for(int i(0); i < _A.n_rows; i++){
            std::vector<double> temp;
            temp = arma::conv_to<std::vector<double> >::from(_A.row(i));
            A_.insert(A_.end(), temp.begin(), temp.end());
        }

        for(int i(0); i < _B.n_rows; i++){
            std::vector<double> temp;
            temp = arma::conv_to<std::vector<double> >::from(_B.row(i));
            B_.insert(B_.end(), temp.begin(), temp.end());
        }

        for(int i(0); i < _C.n_rows; i++){
            std::vector<double> temp;
            temp = arma::conv_to<std::vector<double> >::from(_C.row(i));
            C_.insert(C_.end(), temp.begin(), temp.end());
        }

        for(int i(0); i < _D.n_rows; i++){
            std::vector<double> temp;
            temp = arma::conv_to<std::vector<double> >::from(_D.row(i));
            D_.insert(D_.end(), temp.begin(), temp.end());
        }
    }

    np::ndarray getA(){
        Py_intptr_t shape[2] = {3,3};
        np::ndarray res = np::zeros(2, shape, np::dtype::get_builtin<double>());
        std::copy(A_.begin(), A_.end(), reinterpret_cast<double*>(res.get_data()));
        return res;
    }

    np::ndarray getB(){
        Py_intptr_t shape[2] = {3,3};
        np::ndarray res = np::zeros(2, shape, np::dtype::get_builtin<double>());
        std::copy(B_.begin(), B_.end(), reinterpret_cast<double*>(res.get_data()));
        return res;
    }

    np::ndarray getC(){
        Py_intptr_t shape[2] = {3,3};
        np::ndarray res = np::zeros(2, shape, np::dtype::get_builtin<double>());
        std::copy(C_.begin(), C_.end(), reinterpret_cast<double*>(res.get_data()));
        return res;
    }

    np::ndarray getD(){
        Py_intptr_t shape[2] = {3,3};
        np::ndarray res = np::zeros(2, shape, np::dtype::get_builtin<double>());
        std::copy(D_.begin(), D_.end(), reinterpret_cast<double*>(res.get_data()));
        return res;
    }


private:

    int n_states_;
    int n_inputs_;
    int n_outputs_;

    std::vector<double> A_;
    std::vector<double> B_;
    std::vector<double> C_;
    std::vector<double> D_;

};

//BOOST_PYTHON_MODULE(simulator){

//    py::class_<Simulator>("Simulator")
//            .def("getA", &Simulator::getA)
//            .def("getB", &Simulator::getB)
//            .def("getC", &Simulator::getC)
//            .def("getD", &Simulator::getD)
//            ;
//}

}
