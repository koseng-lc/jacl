/**
*   @author : koseng (Lintang)
*   @brief : Observer Simulator
*/

#pragma once

#include "observer.h"
#include "sim.h"

namespace jacl{

template <class SSpace>
class ObserverSim: public Sim<SSpace>{
public:
    ObserverSim(SSpace* _ss, const arma::mat& _K);
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
    Observer<SSpace> observer_;

};

template <class SSpace>
ObserverSim<SSpace>::ObserverSim(SSpace* _ss, const arma::mat& _K)
    : Sim<SSpace >(_ss->A().n_rows + _ss->C().n_rows) // State and ouputs only
    , observer_(_ss, _K){

}

template <class SSpace>
ObserverSim<SSpace>::~ObserverSim(){

}

template <class SSpace>
void ObserverSim<SSpace>::setInput(const arma::mat& _in){

    observer_.setInput(_in);
}

template <class SSpace>
arma::mat ObserverSim<SSpace>::signalCalc(){
    return observer_.calcStateEst();
}

}
