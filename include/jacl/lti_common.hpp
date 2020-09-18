/**
*   @author : koseng (Lintang)
*   @brief : Some LTI Properties
*/

#pragma once

#include <cassert>
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>

#include <jacl/traits.hpp>
#include <jacl/numerical_methods.hpp>

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
    static auto stabilizable(const StateMatrix& _A, const InputMatrix& _B, bool _continuous = true){
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
            if(_continuous){
                if(std::real(eigval(i)) > .0){
                    temp = _A - eigval(i)*eye;
                    if(arma::rank( arma::join_horiz(temp, cx_B) ) < _A.n_rows){
                        ~ok;
                        break;
                    }
                }
            }else{
                if(std::abs(eigval(i)) - 1. > 1e-4){
                    temp = _A - eigval(i)*eye;
                    if(arma::rank( arma::join_horiz(temp, cx_B) ) < _A.n_rows){
                        ~ok;
                        break;
                    }
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
    static auto detectability(const StateMatrix& _A, const OutputMatrix& _C, bool _continuous = true){
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
            if(_continuous){
                if(std::real(eigval(i)) > .0){
                    temp = _A - eigval(i)*eye;
                    if(arma::rank( arma::join_vert(temp, cx_C) ) < _A.n_cols){
                        ~ok;
                        break;
                    }
                }
            }else{
                if(std::abs(eigval(i)) - 1. > 1e-4){
                    temp = _A - eigval(i)*eye;
                    if(arma::rank( arma::join_vert(temp, cx_C) ) < _A.n_cols){
                        ~ok;
                        break;
                    }
                }
            }
        }
        return ok;
    }
    template <typename StateMatrix>
    static auto poles(const StateMatrix& _A){
        arma::cx_vec p;
        arma::cx_mat eigvec;
        arma::eig_gen(p, eigvec, _A);
        return p;
    }

    template <typename StateMatrix>
    static auto isStable(const StateMatrix& _A, bool _continuous){
        arma::cx_vec p( poles(_A) );
        bool ok(true);
        if(_continuous){
            for(const auto& _p:p){
                if(std::real(_p) > .0){
                    ok = false;
                    break;
                }
            }
        }else{
            for(const auto& _p:p){              
                if(std::abs(_p) - 1. > .1){
                    ok = false;
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
static auto stabilizable(const StateMatrix& _A, const InputMatrix& _B, bool _continuous = true){
    return detail::stabilizable(_A, _B, _continuous);
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
static auto detectability(const StateMatrix& _A, const OutputMatrix& _C, bool _continuous = true){
    return detail::detectability(_A, _C, _continuous);
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
        if(std::abs(eigval(i)) == 1.){
            temp = _A - eigval(i)*eye;
            if(arma::rank( arma::join_vert(temp, cx_C) ) < _A.n_cols){
                ~ok;
                break;
            }
        }
    }
    return ok;
}

template <typename Scalar, std::size_t n_states>
static auto poles(const typename arma::Mat<Scalar>::template fixed<n_states,n_states>& _A){
    return detail::poles(_A);
}

template <typename _StateSpace>
static auto poles(const _StateSpace& _ss)
    -> typename std::enable_if_t<
        jacl::traits::is_state_space<_StateSpace>::value, decltype(detail::poles(_ss.A()))>{
    return detail::poles(_ss.A());
}

// template <typename Scalar, std::size_t n_states>
// static auto isStable(const typename arma::Mat<Scalar>::template fixed<n_states,n_states>& _A, bool _continuous){
//     return detail::isStable(_A, _continuous);
// }

template <typename _StateMatrix>
static auto isStable(const _StateMatrix& _A, bool _continuous){
    return detail::isStable(_A, _continuous);
}

template <typename _StateSpace>
static auto discretize(const _StateSpace& _ss, double _sampling_time) -> StateSpacePack{
    arma::mat aug = arma::join_cols(arma::join_rows(_ss.A(), _ss.B()),
                                    arma::zeros(_StateSpace::n_inputs, _StateSpace::n_states + _StateSpace::n_inputs));
    arma::mat expmAB = arma::expmat(aug*_sampling_time);
    arma::mat Ad = expmAB.submat(0, 0, _StateSpace::n_states - 1, _StateSpace::n_states - 1);
    arma::mat Bd = expmAB.submat(0, _StateSpace::n_states,
                                 _StateSpace::n_states-1, (_StateSpace::n_states + _StateSpace::n_inputs) - 1);
    return std::make_tuple(Ad,Bd, _ss.C(), _ss.D());
}

template <typename _StateSpace>
bool isInfNormLessThan(double _gam, const _StateSpace& _ss){
    arma::mat temp;
    arma::mat B = (1./_gam)*_ss.B();
    arma::mat D = (1./_gam)*_ss.D();

    arma::mat R = arma::eye(_StateSpace::n_inputs,_StateSpace::n_inputs) - arma::trans(D)*D;
    arma::mat R_inv = arma::inv(R);
    arma::mat T = arma::eye(_StateSpace::n_outputs,_StateSpace::n_outputs) - D*arma::trans(D);
    arma::mat T_inv = arma::inv(T);
    temp = B*R_inv*arma::trans(D)*_ss.C();
    arma::mat E = _ss.A() + temp;
    arma::mat E_t_inv = arma::inv(arma::trans(E));
    arma::mat G = -B*R_inv*arma::trans(B);
    arma::mat Q = arma::trans(_ss.C())*T_inv*_ss.C();

    arma::mat S11 = E + G*E_t_inv*Q;
    arma::mat S12 = -G*E_t_inv;
    arma::mat S21 = -E_t_inv*Q;
    arma::mat S22 = E_t_inv;

    //-- Symplectic matrix
    arma::mat S = arma::join_cols(
        arma::join_rows(S11, S12),
        arma::join_rows(S21, S22)
    );
    arma::mat I(arma::size(_ss.A()), arma::fill::eye);
    temp = arma::inv(I - _ss.A());
    if(arma::norm(_ss.C()*temp*B+D,2) >= 1)
        return false;  
    arma::cx_vec p( detail::poles(S) );
    bool ok(true);
    for(const auto& _p:p){
        if(std::fabs(std::abs(_p) - 1.0) < 1e-6){
            ok = false;
            break;
        }
    }
    return ok;
}

template <typename _StateSpace>
auto approxInfNorm(const _StateSpace& _ss, double _ubound, double _lbound){
    // auto metric = std::function<bool(double,_StateSpace)>(isInfNormLessThan<_StateSpace>);
    return numerical_methods::bisection(_ss, isInfNormLessThan<_StateSpace>, _ubound, _lbound);
}

} }
