/**
*   @author : koseng (Lintang)
*   @brief : jacl Simulator / Plotter
*/

#pragma once

#include <boost/thread.hpp>
#include <boost/chrono.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/circular_buffer.hpp>

#include <numpy/ndarrayobject.h>

#include <cassert>
#include <algorithm>
#include <map>
#include <initializer_list>
#include <sstream>
#include <thread>
#include <queue>

#include "state_space.h"

namespace jacl{

namespace{
    namespace py = boost::python;
    namespace np = boost::python::numpy;
}

template <class _StateSpace>
class Sim{
public:    
    Sim(std::size_t _n_signals);
    ~Sim();

    auto init() -> int;

    auto simulate() -> void;

    auto setDelay() -> double&{
        return delay_;
    }

    auto setTitle(std::string&& _title) -> void{
        title_ = _title;
    }

    inline auto timeStep() const -> double{
        return d_time_;
    }

    auto setPlotName(std::initializer_list<std::string> _plot_name) -> void;

protected:    
    auto process() -> void;
    virtual auto signalCalc() -> arma::mat = 0;
    virtual auto getSig() -> arma::mat = 0;

    // more safe with RAII style
    class AcquireGIL{
    public:
        AcquireGIL():state(PyGILState_Ensure()){}
        ~AcquireGIL(){PyGILState_Release(state);}

    private:
        PyGILState_STATE state;

    };

protected:
    boost::thread process_thread_;
    mutable boost::mutex sig_mtx_;
    std::atomic<bool> running_;

    py::object sim_;

    std::size_t n_signals_;

    double d_time_;
    double view_interval_;
    std::size_t max_data_;

    boost::circular_buffer<double> signal_data_;
    boost::circular_buffer<double> time_data_;

    double delay_;
    std::string title_;
    std::vector<std::string> plot_name_;
    // to save the released thread state
    PyThreadState *py_state_;
};

template <class _StateSpace>
Sim<_StateSpace>::Sim(std::size_t _n_signals)
    : running_(true)
    , n_signals_(_n_signals)
    , d_time_(1e-4)
    , view_interval_(0.01)
    , max_data_(static_cast<std::size_t>(view_interval_/d_time_))
    , signal_data_(n_signals_ * max_data_)
    , time_data_(max_data_)
    , delay_(.02)
    , title_("Simulation"){

    for(std::size_t i(0); i < max_data_; i++){
        time_data_.push_back(-1.0);
        for(std::size_t j(0); j < n_signals_; j++){
            signal_data_.push_back(.0);
        }
    }
}

template <class _StateSpace>
Sim<_StateSpace>::~Sim(){

    running_ = false;    
    process_thread_.join();

    if(PyGILState_GetThisThreadState() == py_state_)
        PyThreadState_Swap(py_state_);
}

template <class _StateSpace>
auto Sim<_StateSpace>::init() -> int{

    if(!Py_IsInitialized())
        Py_Initialize();

    // idk why its automatically locked, contrary with common article

    if(PyGILState_Check()){
        PyEval_InitThreads();
        py_state_ = PyEval_SaveThread();
    }

    AcquireGIL lk;

    int ret;
    import_array1(ret);

    try{

        py::object sys = py::import("sys");
        sys.attr("path").attr("append")("../python");

        py::object sim_module = py::import("plotter");

        sim_ = sim_module.attr("Plotter")(n_signals_, d_time_, view_interval_);

    }catch(py::error_already_set){
        PyErr_Print();
        return 1;
    }

    return 0;
}

template <class _StateSpace>
auto Sim<_StateSpace>::process() -> void{

    arma::mat pres_signal;
    auto t(.0);    

    std::vector<double> temp_signal_data;
    std::vector<double> temp_time_data;

    AcquireGIL lk;

    while(running_){

//        pres_signal = signalCalc();
        pres_signal = getSig();

        for(const auto s:pres_signal)
            signal_data_.push_back(s);

        time_data_.push_back(t);

        try{

            temp_signal_data.clear();
            temp_time_data.clear();

            for(std::size_t i(0); i < signal_data_.size(); i++)
                temp_signal_data.push_back(signal_data_[i]);

            for(std::size_t i(0); i < time_data_.size(); i++)
                temp_time_data.push_back(time_data_[i]);

            Py_intptr_t signal_shape[2] = {(int)max_data_, (int)n_signals_};
            PyObject* np_signals = PyArray_SimpleNewFromData(2, signal_shape, NPY_FLOAT64, reinterpret_cast<void*>(temp_signal_data.data()));
            py::handle<> signals_handle(np_signals);
            py::object signals_obj(signals_handle);

            Py_intptr_t time_shape[1] = {(int)max_data_};
            PyObject* np_time = PyArray_SimpleNewFromData(1, time_shape, NPY_FLOAT64, reinterpret_cast<void*>(temp_time_data.data()));
            py::handle<> time_handle(np_time);
            py::object time_obj(time_handle);

            sim_.attr("setData")(signals_obj, time_obj);

        }catch(py::error_already_set){
            PyErr_Print();
        }

        boost::this_thread::sleep_for(boost::chrono::milliseconds((int)(delay_ * 1000.))); // 60 - Hz

        t += d_time_;

    }

}

template <class _StateSpace>
auto Sim<_StateSpace>::simulate() -> void{

    AcquireGIL lk;
    try{

        process_thread_ = boost::thread(boost::bind(&Sim<_StateSpace>::process, this));

        sim_.attr("setTitle")(title_.c_str());
        sim_.attr("setDelay")(delay_);
        sim_.attr("beginSimulation")();

    }catch(py::error_already_set){

        PyErr_Print();

    }

}

template <class _StateSpace>
auto Sim<_StateSpace>::setPlotName(std::initializer_list<std::string> _plot_name) -> void{

    plot_name_.clear();
    plot_name_.insert(plot_name_.end(), _plot_name.begin(), _plot_name.end());

    AcquireGIL lk;
    try{

        py::dict plot_name_dict;
        for(std::size_t i(0); i < plot_name_.size(); i++){
            std::stringstream ss;
            ss << "signal" << i;
            plot_name_dict[ss.str().c_str()] = plot_name_[i];

            sim_.attr("setPlotName")(plot_name_dict);

        }
    }catch(py::error_already_set){

        PyErr_Print();
    }
}

}
