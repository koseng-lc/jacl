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
    int input_size = _B.n_cols;

    arma::mat ctrb(_B);
    for(int i(1); i < system_order; i++){
        int first = (i-1) * input_size;
        int last = i * input_size - 1;
        ctrb = join_rows(ctrb, _A * ctrb.cols(first, last));
    }

    return arma::rank(ctrb) == system_order;

}

static bool stabilizable(const arma::mat& _A, const arma::mat& _B){
    //-- must be square matrix
    assert(_A.n_rows == _A.n_cols);
    //-- must be compatible with A
    assert(_B.n_rows == _A.n_rows);

    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eig_gen(eigval, eigvec, _A);
    eigval.print("EigVal : ");
    eigvec.print("EigVec: ");

    arma::vec eigval_re = arma::real(eigval);

    arma::cx_mat temp;
    arma::cx_mat eye(arma::size(_A), arma::fill::eye);
    arma::cx_mat cx_B( arma::size(_B) );
    cx_B.set_real(_B);
    bool ok(true);
    for(int i(0); i < eigval_re.n_rows; i++){
        if(eigval_re(i) > .0){
            temp = _A - eigval(i)*eye;
            if(arma::rank( arma::join_horiz(temp, cx_B) ) < _A.n_rows){
                ~ok;
                break;
            }
        }
    }

    return ok;
}

static bool observable(const arma::mat& _A, const arma::mat& _C){
    //-- must be square matrix
    assert(_A.n_rows == _A.n_cols);
    //-- must be compatible with A
    assert(_C.n_cols == _A.n_cols);

    int system_order = _A.n_cols;
    int output_size = _C.n_rows;

    arma::mat obsv(_C);
    for(int i(1); i < system_order; i++){
        int first = (i-1) * output_size;
        int last = i * output_size - 1;
        obsv = join_cols(obsv, obsv.rows(first, last) * _A);
    }

    return arma::rank(obsv) == system_order;
}

static bool detectability(const arma::mat& _A, const arma::mat& _C){

}

}

}
