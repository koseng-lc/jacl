/**
*   @author : koseng (Lintang)
*   @brief : jacl physical parameter abstract
*/

#pragma once

namespace jacl{

template <typename Scalar=double>
class PhysicalParameter{
public:
    using scalar_t = Scalar;

    PhysicalParameter(scalar_t _nominal)
        : nominal(_nominal)
        , perturbed(nominal){}

    ~PhysicalParameter(){}

    scalar_t nominal;
    scalar_t perturbed;
};

}
