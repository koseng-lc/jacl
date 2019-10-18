/*
 * There's something problem with Boost.Numpy
 *
 */

#include <boost/python.hpp>
#include <numpy/arrayobject.h>

#include <iostream>

namespace py = boost::python;

int main(int argc, char** argv){
    Py_Initialize();
    import_array();
    try{
        py::object sys = py::import("sys");
        sys.attr("path").attr("append")("../python");

        py::object sim_module = py::import("simulator");

        std::vector<double > A{1, 2, 3,
                               4, 5, 6,
                               7, 8, 9};

        Py_intptr_t shape_A[2] = {3,3};
        //-- simple but deprecated
        PyObject* npA_ = PyArray_SimpleNewFromData(2, shape_A, NPY_FLOAT64, reinterpret_cast<void*>(A.data()));
        py::handle<> npA_handle(npA_);
        py::object npA(npA_handle);

        std::vector<double > B{5, 1,
                               2, 3,
                               5, 9};

        Py_intptr_t shape_B[2] = {3,2};
        //-- simple but deprecated
        PyObject* npB_ = PyArray_SimpleNewFromData(2, shape_B, NPY_FLOAT64, reinterpret_cast<void*>(B.data()));
        py::handle<> npB_handle(npB_);
        py::object npB(npB_handle);

        std::vector<double > C{1, 3, 5,
                               9, 2, 6};

        Py_intptr_t shape_C[2] = {2,3};
        //-- simple but deprecated
        PyObject* npC_ = PyArray_SimpleNewFromData(2, shape_C, NPY_FLOAT64, reinterpret_cast<void*>(C.data()));
        py::handle<> npC_handle(npC_);
        py::object npC(npC_handle);

        std::vector<double > D{5, 1,
                               3, 11};

        Py_intptr_t shape_D[2] = {2,2};
        //-- simple but deprecated
        PyObject* npD_ = PyArray_SimpleNewFromData(2, shape_D, NPY_FLOAT64, reinterpret_cast<void*>(D.data()));
        py::handle<> npD_handle(npD_);
        py::object npD(npD_handle);


        py::object sim = sim_module.attr("Simulator")(3,2,2);

        sim.attr("setStateSpace")(npA, npB, npC, npD);
        sim.attr("setDelay")(0.0);
        sim.attr("setTitle")("Sim Test");
        sim.attr("beginSimulation")();

    }catch(py::error_already_set){
        PyErr_Print();
    }
    return 0;
}
