/**
*   @author : koseng (Lintang)
*   @brief : Some LTI Properties
*/
#pragma once

#include <armadillo>

#include <cassert>

namespace jacl{

namespace common{

static bool controlable(const arma::mat& _A, const arma::mat& _B){
    //-- must be square matrix
    assert(_A.n_rows == _A.n_cols);
    //-- must be compatible with A
    assert(_B.n_rows == _A.n_rows);

    int system_order = _A.n_cols;

    arma::mat ctrb(_B);
    for(int i(1); i < system_order; i++)
        ctrb.join_rows(_A * ctrb.col(i-1));

    return arma::rank(ctrb) == system_order;

}

static bool observable(const arma::mat& _A, const arma::mat& _C){
    //-- must be square matrix
    assert(_A.n_rows == _A.n_cols);
    //-- must be compatible with A
    assert(_C.n_cols == _A.n_cols);

    int system_order = _A.n_cols;

    arma::mat obsv(_C);
    for(int i(0); i < system_order; i++)
        obsv.join_cols(obsv.row(i) * _A);

    return arma::rank(obsv) == system_order;
}

}

}
