#pragma once

#include <numeric>
#include <cassert>

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
    assert(!mval(lbound, target) && mval(ubound, target) && ubound >= lbound && "Bounding Error !");
    while((ubound - lbound)/lbound < tol){
        est = (ubound + lbound) * .5;
        if(mf(est, target))
            ubound = est;
        else
            lbound = est;
    }

    return est;
}

}

}
