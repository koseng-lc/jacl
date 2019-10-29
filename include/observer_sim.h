/**
*   @author : koseng (Lintang)
*   @brief : Observer Simulator
*/

#pragma once

#include "observer.h"
#include "sim.h"

namespace JACL{

template <class SSpace>
class ObserverSim: public Sim<SSpace>{
public:
    ObserverSim(SSpace* _ss, const Mat& _K);
    ~ObserverSim();

    void setInput(const Mat& _in);

    void updateVariables(){
        observer_.updateVariables();
    }

    void setGain(const Mat& _K){
        observer_.setGain(_K);
    }

protected:
    Mat signalCalc();

private:
    Observer<SSpace> observer_;

};

template <class SSpace>
ObserverSim<SSpace>::ObserverSim(SSpace* _ss, const Mat& _K)
    : Sim<SSpace >(_ss->A().n_rows + _ss->C().n_rows) // State and ouputs only
    , observer_(_ss, _K){

}

template <class SSpace>
ObserverSim<SSpace>::~ObserverSim(){

}

template <class SSpace>
void ObserverSim<SSpace>::setInput(const Mat& _in){

    observer_.setInput(_in);
}

template <class SSpace>
Mat ObserverSim<SSpace>::signalCalc(){
    return observer_.calcStateEst();
}

}
