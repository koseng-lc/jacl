/**
*   @author : koseng (Lintang)
*   @brief : Observer Simulator
*/

#pragma once

#include "observer.h"
#include "sim.h"

namespace jacl{

template <class _StateSpace>
class ObserverSim: public Sim<_StateSpace>{
public:
    ObserverSim(_StateSpace* _ss, const arma::mat& _K);
    ~ObserverSim();

    void setInput(const arma::mat& _in);

    void updateVariables(){
        observer_.updateVariables();
    }

    void setGain(const arma::mat& _K){
        observer_.setGain(_K);
    }

protected:
    arma::mat signalCalc();

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
void ObserverSim<_StateSpace>::setInput(const arma::mat& _in){

    observer_.setInput(_in);
}

template <class _StateSpace>
arma::mat ObserverSim<_StateSpace>::signalCalc(){
    return observer_.calcStateEst();
}

}
