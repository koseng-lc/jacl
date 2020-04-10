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
        z_.fill(arma::zeros<arma::vec>(_System::n_outputs - 1, 1));
        prev_z_.fill(arma::zeros<arma::vec>(_System::n_outputs - 1, 1));  
        threshold_.fill(std::numeric_limits<double>::max());     
    }
    ~IFD(){}

    virtual auto init(std::initializer_list<arma::vec > _poles) -> void{
        transform(continuous_tag{});        
        std::vector<arma::vec> poles(_poles);
        arma::mat gain;
        for(int i(0); i < _System::n_states; i++){
            jacl::LinearStateSpace<_System::n_states - 1, 1, 1> ss(A22(An_[i]), arma::zeros(_System::n_states - 1, 1),
                                                                   A12(An_[i]), arma::zeros(1, 1));
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
    virtual auto detect(const arma::vec& _in, const arma::vec& _out)
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
            std::get<0>(res[i]) = y_hat;
            std::get<1>(res[i]) = err > threshold_[i];
        }
        return res;
    }
protected:
    struct continuous_tag{};
    struct discrete_tag{};
    template <typename __System = _System>
    auto convolve()
        -> typename std::enable_if<
            ::jacl::traits::is_continuous_system<__System>::value, void>::type{

    }
    template <typename __System = _System>
    auto convolve(const arma::vec& _in, const arma::vec& _out)
        -> typename std::enable_if<
            ::jacl::traits::is_discrete_system<__System>::value, void>::type{
        for(int i(0); i < _System::n_outputs; i++)
            convolve(_in, _out, i, continuous_tag{});
    }
    template <typename __System = _System>
    auto convolve(const arma::vec& _in, const arma::vec& _out, std::size_t _idx, continuous_tag)
        -> typename std::enable_if<
            ::jacl::traits::is_discrete_system<__System>::value, void>::type{
        if(_idx >= 0 && _idx < _System::n_outputs){
            arma::mat term1, term2, term3;
            arma::vec yn;
            term1 = A22(An_[_idx]) - L_[_idx]*A12(An_[_idx]);
            term2 = A21(An_[_idx]) - L_[_idx]*A11(An_[_idx]);
            term3 = B2(Bn_[_idx]) - L_[_idx]*B1(Bn_[_idx]);
            yn = _out(_idx) - sys_->D().row(_idx) * _in;
            z_[_idx] = term1*prev_z_[_idx] + term1*L_[_idx]*yn
                    + term2*yn + term3*_in;
            prev_z_[_idx] = z_[_idx];
            q_hat_[_idx] = arma::join_vert(yn, z_[_idx] + L_[_idx]*yn);
        }            
    }
    template <typename __System = _System>
    auto convolve(const arma::vec& _in, const arma::vec& _out, std::size_t _idx, discrete_tag)
        -> typename std::enable_if<
            ::jacl::traits::is_discrete_system<__System>::value, void>::type{
        if(_idx >= 0 && _idx < _System::n_outputs){
            arma::mat term1, term2, term3;
            arma::vec yn;
            term1 = A22(An_.front()) - L_.front()*A12(An_.front());
            term2 = A21(An_.front()) - L_.front()*A11(An_.front());
            term3 = B2(Bn_.front()) - L_.front()*B1(Bn_.front());
            yn = _out(_idx) - sys_->D().row(_idx) * _in;
            z_.front() = term1*prev_z_.front() + term1*L_.front()*yn
                    + term2*yn + term3*_in;
            prev_z_.front() = z_.front();
            q_hat_.front() = arma::join_vert(yn, z_.front() + L_.front()*yn);
        }            
    } 
    auto transform(continuous_tag) -> void{
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
    auto transform(std::size_t _idx, discrete_tag) -> void{
        arma::mat desired_form(_System::n_outputs, 1, arma::fill::eye);
        arma::mat C_t = arma::trans( sys_->C() );
        std::array<arma::mat, num_do> M_t;
        arma::mat I(_System::n_states, _System::n_states, arma::fill::eye);
        for(int j(0); j < _System::n_states; j++){
            if(C_t(j,_idx) != 0){
                I = arma::shift(I, -j, 1);
                break;
            }
        }
        M_t.front() = I;  
        arma::mat temp;
        for(int j(0); j < _System::n_states; j++){
            if(C_t(j,_idx) != 0){
                M_t.front()(j,0) = C_t(j,_idx);
            }
        }
        M_.front() = arma::trans(M_t.front());
        temp = sys_->A() * arma::inv(M_.front());
        An_.front() = M_.front() * temp;
        Bn_.front() = M_.front() * sys_->B();
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
    std::array<arma::mat, num_do> An_;
    //-- transformed B
    std::array<arma::mat, num_do> Bn_;
    //-- for transform coordinate
    std::array<arma::mat, num_do> M_;
    //-- Dedicated Observer gain 
    std::array<arma::mat, _System::n_outputs> L_;
    
    std::array<arma::vec, num_do> q_hat_;
    std::array<arma::vec, num_do> z_;
    std::array<arma::vec, num_do> prev_z_;
    std::array<double, _System::n_outputs> threshold_;
    _System* sys_;        
};

} } // namespace jacl::diagnosis