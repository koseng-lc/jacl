/**
*   @author : koseng (Lintang)
*   @brief : Data of transient performance
*/

#pragma once

#include <jacl/traits.hpp>
#include <jacl/plotter.hpp>

namespace jacl{ namespace analysis{

using transient_data_t = std::tuple<double,double,double,double>;

static auto getRiseTime(const transient_data_t& _data){
    return std::get<0>(_data);
}

static auto getPeakTime(const transient_data_t& _data){
    return std::get<1>(_data);
}

static auto getOvershoot(const transient_data_t& _data){
    return std::get<2>(_data);
}

static auto getSettlingTime(const transient_data_t& _data){
    return std::get<3>(_data);
}

//-- we copy it because we call reset()
template <typename _System, std::size_t num_sample>
static auto transient(_System _sys, typename _System::input_t _input,
                      std::string&& _plot_name,
                      double _ts_threshold=.02)
    -> typename std::enable_if_t<
        ::jacl::traits::is_siso_v<typename _System::state_space_t>,transient_data_t>{
    _sys.reset();
    constexpr auto NUM_SAMPLE_DATA(num_sample);
    const typename _System::input_t in(_input);
    typename _System::output_t out(arma::fill::zeros);
    std::array<typename _System::scalar_t, NUM_SAMPLE_DATA> response;

    //-- rise time stuff
    auto tr_start(-1.);
    auto tr_end(-1.);
    auto tr(.0);

    //-- peak time stuff
    auto tp(.0);

    //-- settling time
    auto ts(.0);

    //-- overshoot
    auto overshoot(.0);
    auto peak_val(.0);

    for(auto i(0); i < NUM_SAMPLE_DATA; i++){
        out = _sys.convolve(in);                   
                
        if(out(0) >= 0.1*in(0) && tr_start < 0)
            tr_start = i;
        if(out(0) >= 0.9*in(0) && tr_end < 0)
            tr_end = i;

        if(out(0) > peak_val){
            peak_val = out(0);
            tp = i;
        }

        response[i] = out(0);
    }
    tr = (tr_end - tr_start)*_sys.dt();
    tp *= _sys.dt();
    if(tr_end < 0){ //-- overdamped
        overshoot = 0.;
    }else{
        overshoot = std::fabs(1. - peak_val/in(0)) * 100.;
    }    

    //-- wayback to find when the response reach 2% of error
    for(auto i(NUM_SAMPLE_DATA - 1); i >= 0; i--){
        if(std::fabs(1. - response[i]) >= _ts_threshold){
            ts = (decltype(ts))i*_sys.dt();
            break;
        }
    }
    std::vector<double> resp(response.begin(), response.end());
    jacl::plot(resp, _sys.dt(), std::move(_plot_name), {_plot_name});
    return std::make_tuple(tr,tp,overshoot,ts);
}

template <typename Scalar>
static auto transient(std::vector<Scalar>& _response,
                      const arma::Col<Scalar>& _in,
                      double _sampling_period,
                      std::string&& _plot_name,
                      double _ts_threshold=.02)
    -> transient_data_t{
     //-- rise time stuff
    auto tr_start(-1.);
    auto tr_end(-1.);
    auto tr(.0);

    //-- peak time stuff
    auto tp(.0);

    //-- settling time
    auto ts(.0);

    //-- overshoot
    auto overshoot(.0);
    auto peak_val(.0);

    for(auto i(0); i < _response.size(); i++){                  
                
        if(_response[i] >= 0.1*_in(0) && tr_start < 0)
            tr_start = i;
        if(_response[i] >= 0.9*_in(0) && tr_end < 0)
            tr_end = i;

        if(_response[i] > peak_val){
            peak_val = _response[i];
            tp = i;
        }
    }
    
    tr = (tr_end - tr_start)*_sampling_period;
    tp *= _sampling_period;
    if(tr_end < 0){ //-- overdamped
        overshoot = 0.;
    }else{
        overshoot = std::fabs(1. - peak_val/_in(0)) * 100.;
    } 

    //-- wayback to find when the response reach 2% of error
    for(auto i(_response.size() - 1); i >= 0; i--){
        if(std::fabs(1. - _response[i]) >= _ts_threshold){
            ts = (decltype(ts))i*_sampling_period;
            break;
        }
    }

    jacl::plot(_response, _sampling_period, std::move(_plot_name), {_plot_name});
    return std::make_tuple(tr,tp,overshoot,ts);
}

} }