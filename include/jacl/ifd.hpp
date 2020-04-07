#pragma once

#include <jacl/system/observer.hpp>

namespace jacl{

template <typename _System>
class IFD{
public:
    IFD(_System* _sys)
        : sys_(_sys){              
    }
    ~IFD(){}
    template <typename Type>
    auto init(std::initializer_list<Type> _poles) -> void{
        arma::mat desired_form(_System::n_outputs, 1, arma::fill::eye);
        arma::mat C_t = arma::trans( sys_->C() );
        std::vector<arma::mat> M_t;
        for(int i(0); i < _System::n_outputs; i++){
            arma::mat I(_System::n_states, _System::n_states, arma::fill::eye);
            for(int j(0); j < _System::n_states; j++){
                if(C_t(j,i) != 0){
                    I = arma::shift(I, -j, 1);
                    break;
                }
            }
            M_t.push_back(I);         
        }
        arma::mat temp;
        for(int i(0); i < _System::n_outputs; i++){
            for(int j(0); j < _System::n_states; j++){
                if(C_t(j,i) != 0){
                    M_t[i](j,0) = C_t(j,i);
                }
            }
            M_.push_back(arma::trans(M_t[i]));
            temp = ss_->A() * arma::inv(M_[i]);
            An_.emplace_back(M_[i] * temp);
            Bn_.emplace_back(M_[i] * ss_->B());
        }
        // for(const auto& m:M_){
        //     m.print("M : ");
        // }
    }
    //-- A22 -> (ns-1) x (ns-1)
    //-- A12 -> 1 x (ns-1)
    //-- A21 -> (n-1) x 1
    //-- A11 -> 1 x 1
private:
    auto transform() -> void{

    }
private:
    //-- transformed A
    std::vector<arma::mat> An_;
    //-- transformed B
    std::vector<arma::mat> Bn_;
    //-- for transform coordinate
    std::vector<arma::mat> M_;
    //-- Dedicated Observer gain 
    std::vector<arma::mat> L_;
    _System* sys_;        
};

}