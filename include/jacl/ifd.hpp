#pragma once

#include <jacl/system/observer.hpp>
#include <jacl/pole_placement.hpp>

namespace jacl{

// template <typename _System>
// class IFD;

// namespace detail{
//     template <typename _System>
//     class BaseSystemClient{
//         inline static auto ss(_System* _sys) -> decltype(_sys->ss_){
//             return _sys->ss_;
//         }
//         friend class ::jacl::IFD<_System>;
//     }
// }

template <typename _System>
class IFD{
public:
    IFD(_System* _sys)
        : sys_(_sys){              
    }
    ~IFD(){}

    auto init(std::initializer_list<arma::vec > _poles) -> void{
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
            arma::mat test = desired_form.t() * M_[i];
            // test.print("TEST : ");
            temp = sys_->A() * arma::inv(M_[i]);
            An_.emplace_back(M_[i] * temp);
            Bn_.emplace_back(M_[i] * sys_->B());
        }
        std::vector<arma::vec> poles(_poles);
        for(int i(0); i < _System::n_states; i++){
            jacl::LinearStateSpace<_System::n_states - 1, 1, 1> ss(A22(An_[i]), arma::zeros(_System::n_states - 1, 1),
                                                                   A12(An_[i]), arma::zeros(1, 1));
            arma::mat gain;
            // std::cout << "TEST" << std::endl;
            ss.A().print("A : ");
            ss.C().print("C : ");
            jacl::pole_placement::KautskyNichols(&ss, poles[i], &gain, jacl::pole_placement::PolePlacementType::Observer);
            gain.print("Gain : ");
            L_.emplace_back(gain);
        }        
        //-- Debugging part
        // for(const auto& p:_poles)
        //     p.print("Poles : ");
        // for(const auto& m:M_){
        //     m.print("M : ");
        // }
        // for(const auto& a:An_){
        //     A11(a).print("A11 : ");
        //     A12(a).print("A12 : ");
        //     A21(a).print("A21 : ");
        //     A22(a).print("A22 : ");
        // }
    }

private:
    auto transform() -> void{

    }
    inline auto A11(const arma::mat& _A) -> arma::mat{
        return _A.submat(0,0,0,0);
    }
    inline auto A12(const arma::mat& _A) -> arma::mat{
        return _A.submat(0, 1, 0, _System::n_states - 1);
    }
    inline auto A21(const arma::mat& _A) -> arma::mat{
        return _A.submat(1, 0, _System::n_states - 1, 0);
    }
    inline auto A22(const arma::mat& _A) -> arma::mat{
        return _A.submat(1, 1, _System::n_states - 1, _System::n_states - 1);
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