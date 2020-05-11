/**
*   @author : koseng (Lintang)
*   @brief : Simplified Instrument Fault Detection
*/

#pragma once

#include <jacl/diagnosis/ifd.hpp>

namespace jacl{ namespace diagnosis {

template <typename _System, std::size_t chosen_state>
class SIFD:public IFD<_System>{
public:
    SIFD(_System* _sys, std::initializer_list<double> _threshold)
        : IFD<_System>(_sys, _threshold){

    }
    ~SIFD(){}
    void init(std::initializer_list<
                 typename arma::Col<
                    typename _System::scalar_t>::template fixed<_System::n_states-1>> _poles,
              std::string&& _plot_title = "",
              std::initializer_list<std::string> _plot_name = {}) override{
        // std::cout << "SIFD POLES : " << _poles.size() << std::endl;
        // assert(("[SIFD] Number of set poles didn't matches", _poles.size() == 1));
        // static_assert(_poles.size() == 1, "[SIFD] Number of set poles didn't matches");
        std::vector<typename arma::Col<
                        typename _System::scalar_t>::template fixed<_System::n_states-1>> poles(_poles);
        IFD<_System>::template setPoleDOs<CHOSEN_STATE, decltype(this->dos_)>::set(&this->dos_, poles, typename IFD<_System>::sifd_tag());
        if(_plot_name.size()){
            this->plt_.init();
            this->plt_.setTitle(std::move(_plot_title));
            this->plt_.setPlotName(_plot_name);
        }        
    }
    auto detect(const typename _System::input_t& _in,
                const typename _System::output_t& _out)
        -> typename IFD<_System>::diag_pack_t override{        
        typename IFD<_System>::diag_pack_t res;
        typename _System::output_t y_hat;        
        IFD<_System>::template getEst<CHOSEN_STATE, decltype(this->dos_)>::get(&this->dos_, _in, _out, &y_hat);
        arma::vec delta_y( _out - y_hat );
        bool chosen_state_status(true);
        for(int i(0); i < res.size(); i++){
            if(i == CHOSEN_STATE)continue;
            std::get<0>(res[i]) = y_hat(i);
            std::get<1>(res[i]) = delta_y(i);
            if(delta_y(i) <= 1e-6)
                std::get<2>(res[i]) = false;
            else    
                std::get<2>(res[i]) = std::fabs(delta_y(i)/_out(i)) > this->threshold_[i];
            chosen_state_status &= std::get<2>(res[i]);
        }
        std::get<0>(res[CHOSEN_STATE]) = y_hat(CHOSEN_STATE);
        std::get<1>(res[CHOSEN_STATE]) = delta_y(CHOSEN_STATE);
        std::get<2>(res[CHOSEN_STATE]) = chosen_state_status;
        if(chosen_state_status){
            for(int i(0); i < res.size(); i++){
                if(i == CHOSEN_STATE)continue;
                std::get<2>(res[i]) = false;
            }
        }
        this->aux_sys_.collectSig(y_hat, _in, y_hat);
        return res;
    }
public:
    static constexpr std::size_t CHOSEN_STATE{chosen_state};
};

} } // namespace jacl::diagnosis