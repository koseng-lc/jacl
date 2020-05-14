/**
*   @author : koseng (Lintang)
*   @brief : Some LTI Properties
*/

#pragma once

#include <cassert>

#include <armadillo>

namespace jacl{ namespace lti_common{

namespace detail{
    template <typename StateMatrix, typename InputMatrix>
    static auto controllable(const StateMatrix& _A, const InputMatrix& _B){
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

    template <typename StateMatrix, typename InputMatrix>
    static auto stabilizable(const StateMatrix& _A, const InputMatrix& _B){
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

    template <typename StateMatrix, typename OutputMatrix>
    static auto observable(const StateMatrix& _A, const OutputMatrix& _C, arma::mat* _O = nullptr){
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

    template <typename StateMatrix, typename OutputMatrix>
    static auto detectability(const StateMatrix& _A, const OutputMatrix& _C){
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
}

using StateSpacePack = std::tuple<arma::mat, arma::mat, arma::mat, arma::mat >;

template <typename StateMatrix, typename InputMatrix>
static auto controllable(const StateMatrix& _A, const InputMatrix& _B){
    return detail::controllable(_A, _B);
}

template <typename _StateSpace>
static auto controllable(const _StateSpace& _ss){
    return detail::controllable(_ss.A(), _ss.B());
}

template <typename StateMatrix, typename InputMatrix>
static auto stabilizable(const StateMatrix& _A, const InputMatrix& _B){
    return detail::stabilizable(_A, _B);
}

template <typename _StateSpace>
static auto stabilizable(const _StateSpace& _ss){
    return detail::stabilizable(_ss.A(), _ss.B());
}

template <typename StateMatrix, typename InputMatrix>
static auto hasUncontrollableModeInImAxis(const StateMatrix& _A, const InputMatrix& _B){
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

template <typename StateMatrix, typename InputMatrix>
static auto hasUncontrollableModeInUnitCircle(const StateMatrix& _A, const InputMatrix& _B){
    //-- must be square matrix
    assert(_A.n_rows == _A.n_cols);
    //-- must be compatible with A
    assert(_B.n_rows == _A.n_rows);

    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eig_gen(eigval, eigvec, _A);

    arma::cx_mat temp;
    arma::cx_mat eye(arma::size(_A), arma::fill::eye);
    arma::cx_mat cx_B( arma::size(_B) );
    cx_B.set_real(_B);
    bool ok(false);
    for(int i(0); i < eigval.n_rows; i++){
        if(std::abs(eigval(i)) == 1.){
            temp = _A - eigval(i)*eye;
            if(arma::rank( arma::join_horiz(temp, cx_B) ) < _A.n_rows){
                ~ok;
                break;
            }
        }
    }

    return ok; 
}

template <typename StateMatrix, typename OutputMatrix>
static auto observable(const StateMatrix& _A, const OutputMatrix& _C, arma::mat* _O = nullptr){
    return detail::observable(_A, _C, _O);
}

template <typename _StateSpace>
static auto observable(const _StateSpace& _ss, arma::mat* _O = nullptr){
    return detail::observable(_ss.A(), _ss.C(), _O);
}

template <typename StateMatrix, typename OutputMatrix>
static auto detectability(const StateMatrix& _A, const OutputMatrix& _C){
    return detail::detectability(_A, _C);
}

template <typename _StateSpace>
static auto detectability(const _StateSpace& _ss){
    return detail::detectability(_ss.A(), _ss.C());
}

template <typename StateMatrix, typename OutputMatrix>
static auto hasUnobservableModeInImAxis(const StateMatrix& _A, const OutputMatrix& _C){
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

template <typename StateMatrix, typename OutputMatrix>
static auto hasUnobservableModeInUnitCircle(const StateMatrix& _A, const OutputMatrix& _C){
    //-- must be square matrix
    assert(_A.n_rows == _A.n_cols);
    //-- must be compatible with A
    assert(_C.n_cols == _A.n_cols);

    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eig_gen(eigval, eigvec, _A);

    arma::cx_mat temp;
    arma::cx_mat eye(arma::size(_A), arma::fill::eye);
    arma::cx_mat cx_C( arma::size(_C) );
    cx_C.set_real(_C);
    bool ok(false);
    for(int i(0); i < eigval.n_rows; i++){
        if(std::abs(eigval(i)) == .0){
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
static auto poles(const _StateSpace& _ss){
    arma::cx_vec p;
    arma::cx_mat eigvec;
    arma::eig_gen(p, eigvec, _ss.A());
    return p;
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
