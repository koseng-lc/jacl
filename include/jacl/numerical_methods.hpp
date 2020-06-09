#pragma once

#include <numeric>
#include <cassert>

// #define NUMERICAL_METHODS_VERBOSE

namespace jacl{

namespace numerical_methods{

//-- mf is a method that have two arguments(current estimate and that target function) and return boolean
//-- true -> less than
//-- false -> greater than
template <typename Scalar, typename TargetFunction, typename MetricFunction>
static auto bisection(const TargetFunction& target,
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

        #ifdef NUMERICAL_METHODS_VERBOSE
        ++num_iter;
        std::cout << "Iter : " << num_iter << std::endl;
        std::cout << "Est : " << est << std::endl;
        #endif
    }

    return est;
}

}

}
