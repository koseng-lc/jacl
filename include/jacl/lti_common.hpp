/**
*   @author : koseng (Lintang)
*   @brief : Some LTI Properties
*/

#pragma once

#include <cassert>

#define ARMA_DONT_USE_WRAPPER
#include <armadillo>

#include <jacl/traits.hpp>
#include <jacl/linear_algebra.hpp>
#include <jacl/numerical_methods.hpp>

namespace jacl::lti_common{

namespace detail{

    template <std::size_t idx, std::size_t order, std::size_t input_size, typename C, typename A>
    constexpr auto controllable_iter(C* c, const A& a){

        if constexpr(idx < order){
            auto first = (idx-1) * input_size;
            auto last = idx * input_size - 1;
            c->cols(first + input_size, last + input_size) = a * c->cols(first, last);
            controllable_iter<idx+1, order, input_size>(c, a);
        }
    }

    template <typename StateMatrix, typename InputMatrix, typename C=arma::mat>
    constexpr auto controllable(const StateMatrix& _A, const InputMatrix& _B, C* c = nullptr){

        //-- must be square matrix
        static_assert(StateMatrix::n_rows == StateMatrix::n_cols, "StateMatrix must be square!");
        //-- must be compatible with A
        static_assert(InputMatrix::n_rows == StateMatrix::n_rows, "InputMatrix and StateMatrix are not compatible!");

        constexpr auto system_order = StateMatrix::n_rows;
        constexpr auto input_size = InputMatrix::n_cols;

        typename arma::Mat<double>::template fixed<system_order, system_order*input_size> ctrb;
        ctrb.cols(0, input_size-1) = _B;

        controllable_iter<1, system_order, input_size>(&ctrb, _A);

        if(c != nullptr)
            *c = ctrb;

        return arma::rank(ctrb) == system_order;
    }

    template <typename A, typename B>
    constexpr auto stabilizable_c_iter(A, B){

    }

    template <typename A, typename B>
    constexpr auto stabilizable_d_iter(A, B){

    }

    template <bool continuous, typename StateMatrix, typename InputMatrix>
    auto stabilizable(const StateMatrix& _A, const InputMatrix& _B){

        //-- must be square matrix
        static_assert(StateMatrix::n_rows == StateMatrix::n_cols, "StateMatrix must be square!");
        //-- must be compatible with A
        static_assert(InputMatrix::n_rows == StateMatrix::n_rows, "InputMatrix and StateMatrix are not compatible!");

        arma::cx_vec eigval;
        arma::cx_mat eigvec;
        arma::eig_gen(eigval, eigvec, _A);

        arma::cx_mat temp;
        arma::cx_mat eye(arma::size(_A), arma::fill::eye);
        arma::cx_mat cx_B( arma::size(_B) );
        cx_B.set_real(_B);
        auto ok(true);
        if constexpr(continuous){
            for(int i(0); i < eigval.n_rows; i++){
                if(std::real(eigval(i)) > .0){
                    temp = _A - eigval(i)*eye;
                    if(arma::rank( arma::join_horiz(temp, cx_B) ) < _A.n_rows){
                        ~ok;
                        break;
                    }
                }
            }
        }else{
            for(int i(0); i < eigval.n_rows; i++){
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

    template <std::size_t idx, std::size_t order, std::size_t output_size, typename O, typename A>
    constexpr auto observable_iter(O* o, const A& a){

        if constexpr(idx < order){
            auto first = (idx-1) * output_size;
            auto last = idx * output_size - 1;
            o->rows(first + output_size, last + output_size) = o->rows(first, last) * a;
            observable_iter<idx+1, order, output_size>(o, a);
        }
    }

    template <typename StateMatrix, typename OutputMatrix, typename O=arma::mat>
    auto observable(const StateMatrix& _A, const OutputMatrix& _C, O* o = nullptr){

        //-- must be square matrix
        static_assert(StateMatrix::n_rows == StateMatrix::n_cols, "StateMatrix must be square!");
        //-- must be compatible with A
        static_assert(OutputMatrix::n_cols == StateMatrix::n_cols, "OutputMatrix and StateMatrix are not compatible!");

        constexpr auto system_order = StateMatrix::n_cols;
        constexpr auto output_size = OutputMatrix::n_rows;

        typename arma::Mat<double>::template fixed<system_order*output_size, system_order> obsv;
        obsv.rows(0,output_size-1) = _C;

        observable_iter<1,system_order,output_size>(&obsv, _A);

        if(o != nullptr)
            *o = obsv;

        return arma::rank(obsv) == system_order;
    }

    template <bool continous, typename StateMatrix, typename OutputMatrix>
    auto detectability(const StateMatrix& _A, const OutputMatrix& _C){

        //-- must be square matrix
        static_assert(StateMatrix::n_rows == StateMatrix::n_cols, "StateMatrix must be square!");
        //-- must be compatible with A
        static_assert(OutputMatrix::n_cols == StateMatrix::n_cols, "OutputMatrix and StateMatrix are not compatible!");

        arma::cx_vec eigval;
        arma::cx_mat eigvec;
        arma::eig_gen(eigval, eigvec, _A);

        arma::cx_mat temp;
        arma::cx_mat eye(arma::size(_A), arma::fill::eye);
        arma::cx_mat cx_C( arma::size(_C) );
        cx_C.set_real(_C);
        auto ok(true);
        if constexpr(continous){
            for(int i(0); i < eigval.n_rows; i++){
                if(std::real(eigval(i)) > .0){
                    temp = _A - eigval(i)*eye;
                    if(arma::rank( arma::join_vert(temp, cx_C) ) < _A.n_cols){
                        ~ok;
                        break;
                    }
                }
            }
        }else{
            for(int i(0); i < eigval.n_rows; i++){
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
    auto poles(const StateMatrix& _A){
        arma::cx_vec p;
        arma::cx_mat eigvec;
        arma::eig_gen(p, eigvec, _A);
        return p;
    }

    template <bool continuous, typename StateMatrix>
    auto isStable(const StateMatrix& _A){
        arma::cx_vec p( poles(_A) );
        // p.print("Poles : ");
        bool ok(true);
        if constexpr(continuous){
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
constexpr auto controllable(const StateMatrix& _A, const InputMatrix& _B){
    return detail::controllable(_A, _B);
}

template <typename _StateSpace>
constexpr auto controllable(const _StateSpace& _ss){
    return detail::controllable(_ss.A(), _ss.B());
}

template <bool continuous, typename StateMatrix, typename InputMatrix>
auto stabilizable(const StateMatrix& _A, const InputMatrix& _B){
    return detail::stabilizable<continuous>(_A, _B);
}

template <bool continuous, typename _StateSpace>
auto stabilizable(const _StateSpace& _ss){
    return detail::stabilizable<continuous>(_ss.A(), _ss.B());
}

template <typename StateMatrix, typename InputMatrix>
auto hasUncontrollableModeInImAxis(const StateMatrix& _A, const InputMatrix& _B){
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
auto hasUncontrollableModeInUnitCircle(const StateMatrix& _A, const InputMatrix& _B){
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
auto observable(const StateMatrix& _A, const OutputMatrix& _C, arma::mat* _O = nullptr){
    return detail::observable(_A, _C, _O);
}

template <typename _StateSpace>
auto observable(const _StateSpace& _ss, arma::mat* _O = nullptr){
    return detail::observable(_ss.A(), _ss.C(), _O);
}

template <bool continuous, typename StateMatrix, typename OutputMatrix>
auto detectability(const StateMatrix& _A, const OutputMatrix& _C){
    return detail::detectability<continuous>(_A, _C);
}

template <bool continuous, typename _StateSpace>
auto detectability(const _StateSpace& _ss){
    return detail::detectability<continuous>(_ss.A(), _ss.C());
}

template <typename StateMatrix, typename OutputMatrix>
auto hasUnobservableModeInImAxis(const StateMatrix& _A, const OutputMatrix& _C){
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
auto hasUnobservableModeInUnitCircle(const StateMatrix& _A, const OutputMatrix& _C){
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
auto poles(const typename arma::Mat<Scalar>::template fixed<n_states,n_states>& _A){
    return detail::poles(_A);
}

template <typename _StateSpace>
auto poles(const _StateSpace& _ss)
    -> typename std::enable_if_t<
        ::jacl::traits::is_state_space_v<_StateSpace>, decltype(detail::poles(_ss.A()))>{
    return detail::poles(_ss.A());
}

// template <typename Scalar, std::size_t n_states>
// static auto isStable(const typename arma::Mat<Scalar>::template fixed<n_states,n_states>& _A, bool _continuous){
//     return detail::isStable(_A, _continuous);
// }

template <bool continuous, typename _StateMatrix>
auto isStable(const _StateMatrix& _A){
    return detail::isStable<continuous>(_A);
}

template <typename _StateSpace>
auto discretize(const _StateSpace& _ss, double _sampling_time) -> StateSpacePack{
    arma::mat aug = arma::join_cols(arma::join_rows(_ss.A(), _ss.B()),
                                    arma::zeros(_StateSpace::n_inputs, _StateSpace::n_states + _StateSpace::n_inputs));
    arma::mat expmAB = arma::expmat(aug*_sampling_time);
    arma::mat Ad = expmAB.submat(0, 0, _StateSpace::n_states - 1, _StateSpace::n_states - 1);
    arma::mat Bd = expmAB.submat(0, _StateSpace::n_states,
                                 _StateSpace::n_states-1, (_StateSpace::n_states + _StateSpace::n_inputs) - 1);
    return std::make_tuple(Ad,Bd, _ss.C(), _ss.D());
}

template <typename _System>
bool isInfNormLessThan(double _gam, _System _sys){
    constexpr auto continuous = ::jacl::traits::is_continuous_system_v<_System>;
    auto ok(true);
    if constexpr(continuous){
        arma::mat temp;
        temp = arma::trans(_sys.D())*_sys.D();
        arma::mat R = (_gam*_gam)*arma::eye(_System::n_inputs, _System::n_inputs) - temp;
        arma::mat R_inv = arma::inv(R);

        temp = _sys.B()*R_inv*arma::trans(_sys.D());
        arma::mat H11 = _sys.A() + temp*_sys.C();
        arma::mat H12 = _sys.B()*R_inv*arma::trans(_sys.B());
        temp = arma::eye(_System::n_outputs, _System::n_outputs) + _sys.D()*R_inv*arma::trans(_sys.D());
        arma::mat H21 = -arma::trans(_sys.C())*temp*_sys.C();
        arma::mat H22 = -arma::trans(H11);

        arma::mat H = arma::join_cols(
            arma::join_rows(H11, H12),
            arma::join_rows(H21, H22)
        );
        arma::cx_vec p( detail::poles(H) );

        if(linear_algebra::largestSV(_sys.D()) < _gam){
            for(const auto& _p:p){
                // std::cout << "C : " << _p << std::endl;
                if(std::abs(std::real(_p)) < 1e-6){
                    ok = false;
                    break;
                }
            }
        }else{
            ok = false;
        }
    }else{
        arma::mat temp;
        arma::mat B = (1./_gam)*_sys.B();
        arma::mat D = (1./_gam)*_sys.D();

        arma::mat R = arma::eye(_System::n_inputs,_System::n_inputs) - arma::trans(D)*D;
        arma::mat R_inv = arma::inv(R);
        arma::mat T = arma::eye(_System::n_outputs,_System::n_outputs) - D*arma::trans(D);
        arma::mat T_inv = arma::inv(T);
        temp = B*R_inv*arma::trans(D)*_sys.C();
        arma::mat E = _sys.A() + temp;
        arma::mat E_t_inv = arma::inv(arma::trans(E));
        arma::mat G = -B*R_inv*arma::trans(B);
        arma::mat Q = arma::trans(_sys.C())*T_inv*_sys.C();

        arma::mat S11 = E + G*E_t_inv*Q;
        arma::mat S12 = -G*E_t_inv;
        arma::mat S21 = -E_t_inv*Q;
        arma::mat S22 = E_t_inv;

        //-- Symplectic matrix
        arma::mat S = arma::join_cols(
            arma::join_rows(S11, S12),
            arma::join_rows(S21, S22)
        );
        arma::mat I(arma::size(_sys.A()), arma::fill::eye);
        temp = arma::inv(I - _sys.A());

        if(arma::norm(_sys.C()*temp*B+D,2) < 1){
            arma::cx_vec p( detail::poles(S) );
            for(const auto& _p:p){
                // std::cout << "D : "<< _p << std::endl;
                if(std::fabs(std::abs(_p) - 1.0) < 1e-6){
                    ok = false;
                    break;
                }
            }
        }else{
            ok = false;
        }       
    }    
    
    return ok;
}

template <typename _System>
auto approxInfNorm(_System _sys){    
    return numerical_methods::bisection(_sys, isInfNormLessThan<_System> ,100., .1);
}

}
