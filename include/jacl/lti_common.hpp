/**
*   @author : koseng (Lintang)
*   @brief : Some LTI Properties
*/
#pragma once

#include <cassert>

#include <armadillo>

namespace jacl{ namespace common{

using StateSpacePack = std::tuple<arma::mat, arma::mat, arma::mat, arma::mat >;

enum class PBHTestType{
    Column,
    Row
};

static auto controllable(const arma::mat& _A, const arma::mat& _B) -> bool{
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

static auto stabilizable(const arma::mat& _A, const arma::mat& _B) -> bool{
    //-- must be square matrix
    assert(_A.n_rows == _A.n_cols);
    //-- must be compatible with A
    assert(_B.n_rows == _A.n_rows);

    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eig_gen(eigval, eigvec, _A);
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

static auto hasUncontrollableModeInImAxis(const arma::mat& _A, const arma::mat& _B) -> bool{
    //-- must be square matrix
    assert(_A.n_rows == _A.n_cols);
    //-- must be compatible with A
    assert(_B.n_rows == _A.n_rows);

    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eig_gen(eigval, eigvec, _A);
    arma::vec eigval_re = arma::real(eigval);

    arma::cx_mat temp;
    arma::cx_mat eye(arma::size(_A), arma::fill::eye);
    arma::cx_mat cx_B( arma::size(_B) );
    cx_B.set_real(_B);
    bool ok(false);
    for(int i(0); i < eigval_re.n_rows; i++){
        if(eigval_re(i) == .0){
            temp = _A - eigval(i)*eye;
            if(arma::rank( arma::join_horiz(temp, cx_B) ) < _A.n_rows){
                ~ok;
                break;
            }
        }
    }

    return ok;
}

static auto observable(const arma::mat& _A, const arma::mat& _C, arma::mat* _O = nullptr) -> bool{
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
    if(_O != nullptr)
        *_O = obsv;
    return arma::rank(obsv) == system_order;
}

static auto detectability(const arma::mat& _A, const arma::mat& _C) -> bool{
    //-- must be square matrix
    assert(_A.n_rows == _A.n_cols);
    //-- must be compatible with A
    assert(_C.n_cols == _A.n_cols);

    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eig_gen(eigval, eigvec, _A);
    arma::vec eigval_re = arma::real(eigval);

    arma::cx_mat temp;
    arma::cx_mat eye(arma::size(_A), arma::fill::eye);
    arma::cx_mat cx_C( arma::size(_C) );
    cx_C.set_real(_C);
    bool ok(true);
    for(int i(0); i < eigval_re.n_rows; i++){
        if(eigval_re(i) > .0){
            temp = _A - eigval(i)*eye;
            if(arma::rank( arma::join_vert(temp, cx_C) ) < _A.n_cols){
                ~ok;
                break;
            }
        }
    }
    return ok;
}

static auto hasUnobservableModeInImAxis(const arma::mat& _A, const arma::mat& _C) -> bool{
    //-- must be square matrix
    assert(_A.n_rows == _A.n_cols);
    //-- must be compatible with A
    assert(_C.n_cols == _A.n_cols);

    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eig_gen(eigval, eigvec, _A);
    arma::vec eigval_re = arma::real(eigval);

    arma::cx_mat temp;
    arma::cx_mat eye(arma::size(_A), arma::fill::eye);
    arma::cx_mat cx_C( arma::size(_C) );
    cx_C.set_real(_C);
    bool ok(false);
    for(int i(0); i < eigval_re.n_rows; i++){
        if(eigval_re(i) == .0){
            temp = _A - eigval(i)*eye;
            if(arma::rank( arma::join_vert(temp, cx_C) ) < _A.n_cols){
                ~ok;
                break;
            }
        }
    }
    return ok;
}

template <typename _StateSpace>
static auto discretize(const _StateSpace& _ss, double _sampling_time) -> StateSpacePack{
    arma::mat aug = arma::join_cols(arma::join_rows(_ss.A(), _ss.B()),
                                    arma::zeros(_StateSpace::n_inputs, _StateSpace::n_states + _StateSpace::n_inputs));
    arma::mat expmAB = arma::expmat(aug * _sampling_time);
    arma::mat Ad = expmAB.submat(0, 0, _StateSpace::n_states - 1, _StateSpace::n_states - 1);
    arma::mat Bd = expmAB.submat(0, _StateSpace::n_states,
                                 _StateSpace::n_states-1, (_StateSpace::n_states + _StateSpace::n_inputs) - 1);
    return std::make_tuple(Ad,Bd, _ss.C(), _ss.D());
}

} }
