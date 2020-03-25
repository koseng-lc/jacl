/**
*   @author : koseng (Lintang)
*   @brief : Observer Simulator
*/

#pragma once

#include "observer.h"
#include "sim.h"

namespace jacl{

template <class _StateSpace>
class ObserverSim final: public Sim<_StateSpace>{
public:
    ObserverSim(_StateSpace* _ss, const arma::mat& _K);
    ~ObserverSim();

    auto setInput(const arma::mat& _in) -> void;

    auto updateVariables() -> void{
        observer_.updateVariables();
    }

    auto setGain(const arma::mat& _K) -> void{
        observer_.setGain(_K);
    }

protected:
    auto signalCalc() -> arma::mat override;
    auto getSig() -> arma::mat override;

private:
    Observer<_StateSpace> observer_;
};

template <class _StateSpace>
ObserverSim<_StateSpace>::ObserverSim(_StateSpace* _ss, const arma::mat& _K)
    : Sim<_StateSpace >(_StateSpace::n_states + _StateSpace::n_outputs) // State and ouputs only
    , observer_(_ss, _K){

}

template <class _StateSpace>
ObserverSim<_StateSpace>::~ObserverSim(){

}

template <class _StateSpace>
auto ObserverSim<_StateSpace>::setInput(const arma::mat& _in) -> void{
    boost::mutex::scoped_lock lk(this->sig_mtx_);
    observer_.setInput(_in);
}

template <class _StateSpace>
auto ObserverSim<_StateSpace>::signalCalc() -> arma::mat{
    boost::mutex::scoped_lock lk(this->sig_mtx_);
    return observer_.calcStateEst();
}

template <class _StateSpace>
auto ObserverSim<_StateSpace>::getSig() -> arma::mat{
    boost::mutex::scoped_lock lk(this->sig_mtx_);
    arma::mat dummy;
    return dummy;
}


}
