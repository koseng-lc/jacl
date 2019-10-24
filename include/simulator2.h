/**
*   @author : koseng (Lintang)
*   @brief : JACL Simulator2
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

namespace JACL{

namespace py = boost::python;
namespace np = boost::python::numpy;

template <class SSpace>
class Simulator2{
public:
    Simulator2(SSpace *_ss);
    ~Simulator2();

    int init();

    void setInput(const Mat& _input);

    void simulate();

    double& setDelay(){
        return delay_;
    }

    void setTitle(std::string&& _title){
        title_ = _title;
    }

    void setPlotName(std::initializer_list<std::string> _plot_name);

    void updateVariables();
private:

    boost::thread process_thread_;
    boost::mutex data_mtx_;
    void process();
    boost::mutex running_mtx_;
    bool running_;

    void setData();

    py::object sim_;

    SSpace* ss_;

    int n_states_;
    int n_inputs_;
    int n_outputs_;
    int n_signals_;

    Mat u_;

    double d_time_;
    double view_interval_;
    std::size_t max_data_;

    boost::circular_buffer<double> signal_data_;
    boost::circular_buffer<double> time_data_;

    double delay_;
    std::string title_;
    std::vector<std::string> plot_name_;

    Mat state_transition_;

    // State Transition + Identity
    Mat st_I_;

    class AcquireGIL{
    public:
        AcquireGIL():state(PyGILState_Ensure()){}
        ~AcquireGIL(){PyGILState_Release(state);}
    private:
        PyGILState_STATE state;

    };

    PyThreadState *py_state_;

};

template <class SSpace>
Simulator2<SSpace>::Simulator2(SSpace* _ss)
    : running_(true)
    , ss_(_ss)
    , n_states_(ss_->A().n_rows)
    , n_inputs_(ss_->B().n_cols)
    , n_outputs_(ss_->C().n_rows)
    , n_signals_(n_states_ + n_inputs_ + n_outputs_)
    , u_(n_inputs_, 1, arma::fill::zeros)
    , d_time_(1e-4)
    , view_interval_(0.01)
    , max_data_(static_cast<std::size_t>(view_interval_/d_time_))
    , signal_data_(n_signals_ * max_data_)
    , time_data_(max_data_)
    , delay_(.0)
    , title_("Simulation"){

    for(std::size_t i(0); i < max_data_; i++){
        time_data_.push_back(-1.0);
        for(std::size_t j(0); j < n_signals_; j++){
            signal_data_.push_back(.0);
        }
    }
}

template <class SSpace>
Simulator2<SSpace>::~Simulator2(){

    boost::mutex::scoped_lock lk(running_mtx_);
    running_ = false;    
    process_thread_.join();

    PyThreadState_Swap(py_state_);

}

template <class SSpace>
int Simulator2<SSpace>::init(){

    Py_Initialize();

    // idk why its automatically locked, contrary with common article
    py_state_ = PyEval_SaveThread();

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

template <class SSpace>
void Simulator2<SSpace>::updateVariables(){

    state_transition_ = arma::expmat(ss_->A());
    Mat I(n_states_, n_states_, arma::fill::eye);
    st_I_ = I + state_transition_;
}

template <class SSpace>
void Simulator2<SSpace>::setInput(const Mat& _input){

    setData();
    u_ = _input;
}

template <class SSpace>
void Simulator2<SSpace>::setData(){

}

template <class SSpace>
void Simulator2<SSpace>::process(){

    Mat term1, term2, term3, term4;
    Mat prev_u(u_);
    Mat state(n_states_, 1, arma::fill::zeros);
    Mat prev_state(state);
    Mat output;
    Mat pres_signal;
    auto t(.0);    

    AcquireGIL lk;

    while(running_){

        term1 = state_transition_ * prev_state;
        term2 = ss_->B() * (u_ + prev_u);
        term3 = st_I_ * (d_time_ * .5);

        state = term1 + (term3 * term2);

        term4 = ss_->C() * state;
        output = term4 + (ss_->D() * u_);

        pres_signal = arma::join_cols(arma::join_cols(state, u_), output);
//        pres_signal.print("Sig : ");
        for(const auto s:pres_signal)
            signal_data_.push_back(s);

        time_data_.push_back(t);

        prev_state = state;
        prev_u = u_;        
        std::cout << "Process : " << boost::this_thread::get_id() << std::endl;
//        setData();

//        AcquireGIL lk;

        try{

            std::vector<double> temp_signal_data;
            std::vector<double> temp_time_data;

            for(std::size_t i(0); i < signal_data_.size(); i++)
                temp_signal_data.push_back(signal_data_[i]);

            for(std::size_t i(0); i < time_data_.size(); i++)
                temp_time_data.push_back(time_data_[i]);

            Py_intptr_t signal_shape[2] = {(int)max_data_, n_signals_};
            PyObject* np_signals = PyArray_SimpleNewFromData(2, signal_shape, NPY_FLOAT64, reinterpret_cast<void*>(temp_signal_data.data()));
            py::handle<> signals_handle(np_signals);
            py::object signals_obj(signals_handle);

            Py_intptr_t time_shape[1] = {(int)max_data_};
            PyObject* np_time = PyArray_SimpleNewFromData(1, time_shape, NPY_FLOAT64, reinterpret_cast<void*>(temp_time_data.data()));
            py::handle<> time_handle(np_time);
            py::object time_obj(time_handle);

            sim_.attr("setData")(signals_obj, time_obj);
            std::cout << "setData : " << boost::this_thread::get_id() << std::endl;

        }catch(py::error_already_set){
            PyErr_Print();
        }
        boost::this_thread::sleep_for(boost::chrono::milliseconds{100});

        t += d_time_;

    }

}

template <class SSpace>
void Simulator2<SSpace>::simulate(){    
    AcquireGIL lk;
    try{

        process_thread_ = boost::thread(boost::bind(&Simulator2<SSpace>::process, this));

        sim_.attr("setDelay")(delay_);
        sim_.attr("beginSimulation")();

    }catch(py::error_already_set){

        PyErr_Print();

    }

}

template <class SSpace>
void Simulator2<SSpace>::setPlotName(std::initializer_list<std::string> _plot_name){
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
