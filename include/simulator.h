/**
*   @author : koseng (Lintang)
*   @brief : Simulator
*/

#pragma once

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include <numpy/ndarrayobject.h>

#include <cassert>

#include <armadillo>

namespace JACL{

namespace py = boost::python;
namespace np = boost::python::numpy;

class Simulator{
public:
    Simulator(int _n_states, int _n_inputs, int _n_outputs)
        : n_states_(_n_states)
        , n_inputs_(_n_inputs)
        , n_outputs_(_n_outputs){

    }

    int init(){
        Py_Initialize();
        import_array();
        try{
            py::object sys = py::import("sys");
            sys.attr("path").attr("append")("../python");

            py::object sim_module = py::import("simulator");

            sim_ = sim_module.attr("Simulator")(n_states_, n_inputs_, n_outputs_);

            sim_.attr("setTitle")("Sim Test");

        }catch(py::error_already_set){
            PyErr_Print();
            return 1;
        }

        return 0;
    }

    void setStateSpace(const arma::mat& _A, const arma::mat& _B, const arma::mat& _C, const arma::mat& _D){

        //-- A must be square matrix and remain the same
        assert(_A.is_square());
        assert(_A.n_rows == n_states_);

        //-- B shape must be remain the same
        assert(_B.n_rows == n_states_);
        assert(_B.n_cols == n_inputs_);

        //-- C shape must be remain the same
        assert(_C.n_rows == n_outputs_);
        assert(_C.n_cols == n_states_);

        //-- D shape must be remain the same
        assert(_D.n_rows == n_outputs_);
        assert(_D.n_cols == n_inputs_);

        A_.clear();
        B_.clear();
        C_.clear();
        D_.clear();

        for(int i(0); i < n_states_; i++){
            std::vector<double> temp;
            temp = arma::conv_to<std::vector<double> >::from(_A.row(i));
            A_.insert(A_.end(), temp.begin(), temp.end());
        }

        for(int i(0); i < n_states_; i++){
            std::vector<double> temp;
            temp = arma::conv_to<std::vector<double> >::from(_B.row(i));
            B_.insert(B_.end(), temp.begin(), temp.end());
        }

        for(int i(0); i < n_outputs_; i++){
            std::vector<double> temp;
            temp = arma::conv_to<std::vector<double> >::from(_C.row(i));
            C_.insert(C_.end(), temp.begin(), temp.end());
        }

        for(int i(0); i < n_outputs_; i++){
            std::vector<double> temp;
            temp = arma::conv_to<std::vector<double> >::from(_D.row(i));
            D_.insert(D_.end(), temp.begin(), temp.end());
        }

        try{
            Py_intptr_t shape_A[2] = {n_states_, n_states_};
            //-- simple but deprecated
            PyObject* np_A = PyArray_SimpleNewFromData(2, shape_A, NPY_FLOAT64, reinterpret_cast<void*>(A_.data()));
            py::handle<> np_A_handle(np_A);
            py::object A_obj(np_A_handle);

            Py_intptr_t shape_B[2] = {n_states_, n_inputs_};
            //-- simple but deprecated
            PyObject* np_B = PyArray_SimpleNewFromData(2, shape_B, NPY_FLOAT64, reinterpret_cast<void*>(B_.data()));
            py::handle<> np_B_handle(np_B);
            py::object B_obj(np_B_handle);

            Py_intptr_t shape_C[2] = {n_outputs_, n_states_};
            //-- simple but deprecated
            PyObject* np_C = PyArray_SimpleNewFromData(2, shape_C, NPY_FLOAT64, reinterpret_cast<void*>(C_.data()));
            py::handle<> np_C_handle(np_C);
            py::object C_obj(np_C_handle);

            Py_intptr_t shape_D[2] = {n_outputs_, n_inputs_};
            //-- simple but deprecated
            PyObject* np_D = PyArray_SimpleNewFromData(2, shape_D, NPY_FLOAT64, reinterpret_cast<void*>(D_.data()));
            py::handle<> np_D_handle(np_D);
            py::object D_obj(np_D_handle);

            sim_.attr("setStateSpace")(A_obj, B_obj, C_obj, D_obj);
        }catch(py::error_already_set){
            PyErr_Print();
        }
    }    

    void simulate(){

        assert(!sim_.is_none());
        assert(A_.size() == n_states_ * n_states_);
        assert(B_.size() == n_states_ * n_inputs_);
        assert(C_.size() == n_outputs_ * n_states_);
        assert(D_.size() == n_outputs_ * n_inputs_);

        try{

            sim_.attr("setDelay")(delay_);
            sim_.attr("beginSimulation")();

        }catch(py::error_already_set){

            PyErr_Print();

        }

    }

    inline double& setDelay(){
        return delay_;
    }

private:

    py::object sim_;

    int n_states_;
    int n_inputs_;
    int n_outputs_;

    std::vector<double> A_;
    std::vector<double> B_;
    std::vector<double> C_;
    std::vector<double> D_;

    double delay_;


};

}
