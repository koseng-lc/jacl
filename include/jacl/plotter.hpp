/**
*   @author : koseng (Lintang)
*   @brief : jacl Simulator / Plotter
*/

#pragma once

#include <algorithm>
#include <map>
#include <initializer_list>
#include <sstream>
#include <thread>
#include <queue>
#include <cassert>

#include <boost/thread.hpp>
#include <boost/chrono.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/circular_buffer.hpp>

#include <numpy/ndarrayobject.h>

#include "linear_state_space.hpp"

namespace jacl{

namespace{
    namespace py = boost::python;
    namespace np = boost::python::numpy;
}

template <class _System>
class Plotter{
public:    
    Plotter(_System* _sys, std::initializer_list<std::size_t> _selected_sig);
    ~Plotter();

    auto init() -> int;
    auto start() -> void;
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
    std::vector<std::size_t> selected_sig_;

    double d_time_;
    double view_interval_;
    std::size_t max_data_;

    boost::circular_buffer<double> signal_data_;
    boost::circular_buffer<double> time_data_;

    double delay_;
    std::string title_;
    std::vector<std::string> plot_name_;
    // to save the released thread state
    PyThreadState* py_state_;

    _System* sys_;
};

template <class _System>
Plotter<_System>::Plotter(_System* _sys, std::initializer_list<std::size_t> _selected_sig)
    : running_(true)
    , n_signals_(_selected_sig.size())
    , d_time_(1e-4)
    , view_interval_(0.01)
    , max_data_(static_cast<std::size_t>(view_interval_/d_time_))
    , signal_data_(n_signals_ * max_data_)
    , time_data_(max_data_)
    , delay_(.02)
    , title_("Simulation")
    , sys_(_sys){

    for(std::size_t i(0); i < max_data_; i++){
        time_data_.push_back(-1.0);
        for(std::size_t j(0); j < n_signals_; j++){
            signal_data_.push_back(.0);
        }
    }
    // _selected_sig.size();
    selected_sig_.insert(selected_sig_.end(), _selected_sig.begin(), _selected_sig.end());
}

template <class _System>
Plotter<_System>::~Plotter(){

    running_ = false;    
    process_thread_.join();

    if(PyGILState_GetThisThreadState() == py_state_)
        PyThreadState_Swap(py_state_);
}

template <class _System>
auto Plotter<_System>::init() -> int{

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

template <class _System>
auto Plotter<_System>::process() -> void{

    arma::mat pres_signal;
    auto t(.0);    

    std::vector<double> temp_signal_data;
    std::vector<double> temp_time_data;

    AcquireGIL lk;

    while(running_){
        pres_signal = sys_->recapitulate();

        for(const auto& i:selected_sig_){
            signal_data_.push_back(pres_signal(i-1));
        }

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

template <class _System>
auto Plotter<_System>::start() -> void{

    AcquireGIL lk;
    try{

        process_thread_ = boost::thread(boost::bind(&Plotter<_System>::process, this));

        sim_.attr("setTitle")(title_.c_str());
        sim_.attr("setDelay")(delay_);
        sim_.attr("beginSimulation")();

    }catch(py::error_already_set){

        PyErr_Print();

    }

}

template <class _System>
auto Plotter<_System>::setPlotName(std::initializer_list<std::string> _plot_name) -> void{

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

} // namespace jacl