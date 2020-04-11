#pragma once

#include <jacl/diagnosis/ifd.hpp>
#include <armadillo>

namespace jacl{ namespace diagnosis {

template <typename _System, std::size_t chosen_state>
class SIFD:public IFD<_System, 1>{
public:
    SIFD(_System* _sys, std::initializer_list<double> _threshold)
        : IFD<_System, 1>(_sys, _threshold){

    }
    ~SIFD(){}
    auto init(std::initializer_list<arma::vec > _poles) -> void override{
        this->transform(CHOSEN_STATE, typename IFD<_System, 1>::discrete_tag());        
        std::vector<arma::vec> poles(_poles);
        jacl::LinearStateSpace<_System::n_states - 1,
                                1, 1> ss(this->A22(this->An_.front()), arma::zeros(_System::n_states - 1, 1),
                                         this->A12(this->An_.front()), arma::zeros(1, 1));
        arma::mat gain;
        jacl::pole_placement::KautskyNichols(&ss, poles.front(), &gain, jacl::pole_placement::PolePlacementType::Observer);
        ss.A().print("A : ");
        ss.C().print("C : ");
        gain.print("Kautsky-Nichols gain : ");
        jacl::pole_placement::BassGura(&ss, poles.front(), &gain);
        gain.print("Bass-Gura gain : ");   
        // arma::mat gain_from_others;
        // gain_from_others << -9395.4 << arma::endr << 4567.8 << arma::endr;
        // gain = gain_from_others;
        // gain.print("Why gain : ");   
        // arma::cx_mat eigvec;
        // arma::cx_vec eigval;
        // arma::eig_gen(eigval, eigvec, ss.A());
        // eigval.print("EigVal : ");
        // if(jacl::common::observable(ss.A(), ss.C()))
        //     std::cout << "SS Observable !" << std::endl;        
        this->L_.front() = gain;
    }
    auto detect(const arma::vec& _in, const arma::vec& _out)
        -> std::array<std::pair<arma::vec, bool>, _System::n_outputs> override{        
        std::array<std::pair<arma::vec, bool>, _System::n_outputs> res;
        this->convolve(_in, _out, CHOSEN_STATE, typename IFD<_System, 1>::discrete_tag());
        arma::vec x_hat = arma::inv(this->M_.front()) * this->q_hat_.front();
        arma::vec y_hat = this->sys_->C() * x_hat + this->sys_->D() * _in;
        arma::vec delta_y = _out - y_hat;        
        bool chosen_state_status(true);
        for(int i(0); i < res.size(); i++){
            if(i == CHOSEN_STATE)continue;
            std::get<0>(res[i]) = y_hat(i);
            std::get<1>(res[i]) = std::fabs(delta_y(i)) > this->threshold_[i];
            chosen_state_status &= std::get<1>(res[i]);
        }
        std::get<0>(res[CHOSEN_STATE]) = y_hat(CHOSEN_STATE);
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