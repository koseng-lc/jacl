/**
*   @author : koseng (Lintang)
*   @brief : Data of transient performance
*/

#pragma once

#include <jacl/traits.hpp>

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
template <typename _System>
static auto transient(_System _sys, double _ts_threshold=.02)
    -> typename std::enable_if_t<
        ::jacl::traits::is_siso<typename _System::state_space_t>::value,transient_data_t>{
    _sys.reset();
    constexpr auto NUM_SAMPLE_DATA(1000);
    const typename _System::input_t in(arma::fill::ones);
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
                
        if(out(0) >= 0.1 && tr_start < 0)
            tr_start = i;
        if(out(0) >= 0.9 && tr_end < 0)
            tr_end = i;

        if(out(0) > peak_val){
            peak_val = out(0);
            tp = i;
        }

        response[i] = out(0);
    }

    tr = (tr_end - tr_start)*_sys.samplingPeriod();
    tp *= _sys.samplingPeriod();
    overshoot = std::fabs(1. - peak_val) * 100.;

    //-- wayback to find when the response reach 2% of error
    for(auto i(NUM_SAMPLE_DATA - 1); i >= 0; i--){
        if(std::fabs(1. - response[i]) >= _ts_threshold){
            ts = (decltype(ts))i*_sys.samplingPeriod();
            break;
        }
    }
    
    return std::make_tuple(tr,tp,overshoot,ts);
}

} }