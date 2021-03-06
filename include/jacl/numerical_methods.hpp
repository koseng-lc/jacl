/**
*   @author : koseng (Lintang)
*   @brief : jacl numerical methods
*/

#pragma once

#include <numeric>
#include <cassert>

#define NUMERICAL_METHODS_VERBOSE

namespace jacl::numerical_methods{

//-- mf is a method that have two arguments(current estimate and that target function) and return boolean
//-- true -> less than
//-- false -> greater than
template <typename Scalar, typename TargetFunction, typename MetricFunction>
auto bisection(const TargetFunction& target,
               const MetricFunction& mf,
               Scalar ubound,
               Scalar lbound,
               Scalar tol = std::numeric_limits<Scalar>::epsilon()){
    Scalar est;
    // assert(!mval(lbound, target) && mval(ubound, target) && ubound >= lbound && "Bounding Error !");
    std::size_t num_iter(0);
    while(!((ubound - lbound)/lbound < tol)){
        est = (ubound + lbound) * .5;
        if(mf(est, target))
            ubound = est;
        else
            lbound = est;
        ++num_iter;
    }
    #ifdef NUMERICAL_METHODS_VERBOSE    
    std::cout << "[numerical_methods] Iter : " << num_iter << std::endl;
    std::cout << "[numerical_methods] Est : " << est << std::endl;
    #endif
    return est;
}

}
