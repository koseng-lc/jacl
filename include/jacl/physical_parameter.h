/**
*   @author : koseng (Lintang)
*   @brief : jacl Physical Parameter abstract
*/

#pragma once

namespace jacl{

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
