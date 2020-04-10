#pragma once

#include <jacl/traits.hpp>
#include <jacl/system/observer.hpp>
#include <jacl/pole_placement.hpp>
#include <jacl/lti_common.hpp>

namespace jacl{ namespace diagnosis {

template <typename _System, std::size_t num_do = _System::n_outputs>
class IFD{
public:    
    IFD(_System* _sys)
        : sys_(_sys){     
        z_.fill(arma::zeros<arma::vec>(_System::n_outputs, 1));
        prev_z_.fill(arma::zeros<arma::vec>(_System::n_outputs, 1));       
    }
    ~IFD(){}

    auto init(std::initializer_list<arma::vec > _poles) -> void{
        transform();        
        std::vector<arma::vec> poles(_poles);
        for(int i(0); i < _System::n_states; i++){
            jacl::LinearStateSpace<_System::n_states - 1, 1, 1> ss(A22(An_[i]), arma::zeros(_System::n_states - 1, 1),
                                                                   A12(An_[i]), arma::zeros(1, 1));
            arma::mat gain;
            // std::cout << "TEST" << std::endl;
            // ss.A().print("A : ");
            // ss.C().print("C : ");
            // poles[i].print("Poles : ");
            jacl::pole_placement::KautskyNichols(&ss, poles[i], &gain, jacl::pole_placement::PolePlacementType::Observer);
            // gain.print("Gain : ");
            // if(jacl::common::observable(ss.A(), ss.C()))
            //     std::cout << "SS Observable !" << std::endl;
            // arma::cx_mat eigvec;
            // arma::cx_vec eigval;
            // arma::eig_gen(eigval, eigvec, ss.A() - gain*ss.C());
            // eigval.print("EigVal : ");
            L_[i] = gain;
        }
    }
    auto detect(const arma::vec& _in, const arma::vec& _out)
        -> std::array<std::pair<arma::vec, bool>, _System::n_outputs>{
        std::array<std::pair<arma::vec, bool>, _System::n_outputs> res;
        arma::vec x_hat;
        arma::vec y_hat;
        arma::vec delta_y;        
        convolve(_in, _out);
        for(int i(0); i < res.size(); i++){
            x_hat = arma::inv(M_[i]) * q_hat_[i];
            y_hat = sys_->C() * x_hat + sys_->D() * _in;
            delta_y = _out - y_hat;
            auto err(1.);
            for(int j(0); j < _System::n_outputs; j++){
                if(j != i)
                    err *= delta_y(j);
            }
            std::get<0>(res[i]) = x_hat;
            std::get<1>(res[i]) = err > threshold_;
        }
        return res;
    }
protected:
    template <typename __System = _System>
    auto convolve()
        -> typename std::enable_if<
            ::jacl::traits::is_continuous_system<__System>::value, void>::type{

    }
    template <typename __System = _System>
    auto convolve(const arma::vec& _in, const arma::vec& _out)
        -> typename std::enable_if<
            ::jacl::traits::is_discrete_system<__System>::value, void>::type{
        arma::mat term1, term2, term3;
        arma::vec yn;
        for(int i(0); i < _System::outputs; i++){
            term1 = A22(An_[i]) - L_[i]*A12(An_[i]);
            term2 = A21(An_[i]) - L_[i]*A11(An_[i]);
            term3 = B2(Bn_[i]) - L_[i]*B2(Bn_[i]);
            yn = _out(i) - sys_->D().row(i) * _in;
            z_[i] = term1*prev_z_[i] + term1*L_[i]*yn
                    + term2*yn + term3*_in;
            q_hat_[i] = arma::join_vert(yn, z_[i] + L_[i]*yn);
        }        
    }    
    auto transform() -> void{
        arma::mat desired_form(_System::n_outputs, 1, arma::fill::eye);
        arma::mat C_t = arma::trans( sys_->C() );
        std::array<arma::mat, _System::n_outputs> M_t;
        for(int i(0); i < _System::n_outputs; i++){
            arma::mat I(_System::n_states, _System::n_states, arma::fill::eye);
            for(int j(0); j < _System::n_states; j++){
                if(C_t(j,i) != 0){
                    I = arma::shift(I, -j, 1);
                    break;
                }
            }
            M_t[i] = I;  
        }
        arma::mat temp;
        for(int i(0); i < _System::n_outputs; i++){
            for(int j(0); j < _System::n_states; j++){
                if(C_t(j,i) != 0){
                    M_t[i](j,0) = C_t(j,i);
                }
            }
            M_[i] = arma::trans(M_t[i]);
            temp = sys_->A() * arma::inv(M_[i]);
            An_[i] = M_[i] * temp;
            Bn_[i] = M_[i] * sys_->B();
        }
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
    inline auto B1(const arma::mat& _B) -> arma::mat{
        return _B.head_rows(1);
    }
    inline auto B2(const arma::mat& _B) -> arma::mat{
        return _B.tail_rows(_System::n_states - 1);
    }

protected:
    //-- transformed A
    std::array<arma::mat, _System::n_outputs> An_;
    //-- transformed B
    std::array<arma::mat, _System::n_outputs> Bn_;
    //-- for transform coordinate
    std::array<arma::mat, _System::n_outputs> M_;
    //-- Dedicated Observer gain 
    std::array<arma::mat, _System::n_outputs> L_;
    
    std::array<arma::vec, _System::n_outputs> q_hat_;
    std::array<arma::vec, _System::n_outputs> z_;
    std::array<arma::vec, _System::n_outputs> prev_z_;
    std::array<double, _System::n_outputs> threshold_;
    _System* sys_;        
};

} } // namespace jacl::diagnosis