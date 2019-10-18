/**
*   @author : koseng (Lintang)
*   @brief : JACL Physical Parameter abstract
*/

#pragma once

namespace JACL{

class PhysicalParameter{
public:
    PhysicalParameter(double _nominal)
        : nominal(_nominal)
        , perturbed(nominal){}

    ~PhysicalParameter(){}

    double nominal;
    double perturbed;
};

}
