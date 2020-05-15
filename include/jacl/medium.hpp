/**
*   @author : koseng (Lintang)
*   @brief : Prevent circular dependencies
*/

#pragma once

namespace jacl{

namespace state_space::detail{
    template<typename _StateSpace>
    class NonLinearStateSpaceClient{
    protected:
        inline static auto setSig(_StateSpace* _ss, const arma::vec& _sig) -> void{
            _ss->sig_ = arma::conv_to<decltype(_ss->sig_)>::from(_sig);
        }
        inline static auto dstate(_StateSpace* _ss) -> arma::vec {
            arma::vec xdot(_StateSpace::n_states, 1, arma::fill::zeros);
            int i(0);
            for(const auto& fn:_ss->state_fn_){
                xdot(i) = fn(*_ss);
                ++i;
            }
            return xdot;
        }
        inline static auto output(_StateSpace* _ss) -> arma::vec {
            arma::vec out(_StateSpace::n_outputs, 1, arma::fill::zeros);
            int i(0);
            for(const auto& fn:_ss->output_fn_){
                out(i) = fn(*_ss);
                ++i;
            }
            return out;
        }
    };
}

namespace system::detail{
    template <typename _System>
    class BaseSystemClient{
    protected:
        static auto ss(_System* _sys) -> typename _System::state_space_t*{
            if(_sys)
                return _sys->ss_;
            return nullptr;
        }
        template <typename __System, typename _StateSpace>
        static auto setSubject(__System* _sys, _StateSpace* _ss) -> void{
            _sys->setSubject(_ss);
        }        
    };
}

}