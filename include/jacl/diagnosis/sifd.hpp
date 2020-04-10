#pragma once

#include <jacl/diagnosis/ifd.hpp>

namespace jacl{ namespace diagnosis {

template <typename _System, std::size_t dedicated_state>
class SIFD:public IFD<_System>{
public:
    SIFD(_System* _sys)
        : IFD<_System>(_sys){

    }
    ~SIFD(){}
};

} } // namespace jacl::diagnosis