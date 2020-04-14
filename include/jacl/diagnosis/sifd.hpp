#pragma once

#include <jacl/diagnosis/ifd.hpp>
#include <armadillo>

namespace jacl{ namespace diagnosis {

template <typename _System, std::size_t chosen_state>
class SIFD:public IFD<_System>{
public:
    SIFD(_System* _sys, std::initializer_list<double> _threshold)
        : IFD<_System>(_sys, _threshold){

    }
    ~SIFD(){}
    auto init(std::initializer_list<arma::vec > _poles) -> void override{
        std::vector<arma::vec> poles(_poles);
        IFD<_System>::template setPoleDOs<_System::n_outputs-1, typename IFD<_System>::DOS>::set(&this->dos_, poles);
    }
    auto detect(const arma::vec& _in, const arma::vec& _out)
        -> std::array<std::pair<arma::vec, bool>, _System::n_outputs> override{        
        std::array<std::pair<arma::vec, bool>, _System::n_outputs> res;
        std::vector<arma::vec> y_hat;
        std::cout << "MASUK" << std::endl;
        IFD<_System>::template getEst<_System::n_outputs-1, typename IFD<_System>::DOS>::get(&this->dos_, _in, _out, &y_hat);  
        std::cout << "KELUAR" << std::endl;      
        arma::vec delta_y( _out - y_hat[CHOSEN_STATE] );
        bool chosen_state_status(true);
        for(int i(0); i < res.size(); i++){
            if(i == CHOSEN_STATE)continue;
            std::get<0>(res[i]) = y_hat[CHOSEN_STATE](i);
            std::get<1>(res[i]) = std::fabs(delta_y(i)) > this->threshold_[i];
            chosen_state_status &= std::get<1>(res[i]);
        }
        std::get<0>(res[CHOSEN_STATE]) = y_hat[CHOSEN_STATE](CHOSEN_STATE);
        std::get<1>(res[CHOSEN_STATE]) = chosen_state_status;
        if(chosen_state_status){
            for(int i(0); i < res.size(); i++){
                if(i == CHOSEN_STATE)continue;
                std::get<1>(res[i]) = false;
            }
        }        
        return res;
    }
public:
    static constexpr std::size_t CHOSEN_STATE{chosen_state};
};

} } // namespace jacl::diagnosis