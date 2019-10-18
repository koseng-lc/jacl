/**
*   @author : koseng (Lintang)
*   @brief : JACL Physical Parameter abstract
*/

#pragma once

namespace JACL{

class PhysicalParameter{
public:
    PhysicalParameter()
        : nominal(.0)
        , perturbed(.0){}

    ~PhysicalParameter(){}

    double nominal;
    double perturbed;
};

}
