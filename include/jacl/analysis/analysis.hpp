#pragma once

namespace jacl{ namespace analysis{

template <typename _StateSpace>
auto ssv(const _StateSpace& _ss){
    auto ubound(.0);
    auto lbound(.0);
    return std::make_pair(ubound, lbound);
}

} } // namespace::analysis