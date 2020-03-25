#pragma once

#include <cassert>
#include <numeric>

namespace jacl{

namespace numerical_methods{

//-- mf is a method that have two arguments(current estimate and that target function) and return boolean
//-- true -> less than
//-- false -> greater than
template <typename TargetFunction, typename MetricFunction>
static auto bisection(const TargetFunction& target, const MetricFunction& mf, double ubound=1., double lbound=.0, double tol = std::numeric_limits<double>::epsilon()) -> double{
    auto est(.0);
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
