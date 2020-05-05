/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of discrete Algebraic Riccati Equations solver
*/

#pragma once

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <numpy/ndarrayobject.h>

#include <jacl/state_space/linear.hpp>

// #define DARE_VERBOSE

namespace jacl{

namespace{
    namespace py = boost::python;
    namespace np = boost::python::numpy;
}

template <class _StateSpace>
class DARE{
public:
    DARE(_StateSpace* _ss, const arma::mat& _R, const arma::mat& _Q);
    DARE(_StateSpace* _ss);
    ~DARE();

    DARE(const DARE&) = delete;
    DARE& operator=(const DARE&) = delete;

    auto setR(const arma::mat& _R, bool update_hamiltonian = false){
        assert(arma::size(_R) == arma::size(ss_->A()));
        R_ = _R;
        if(update_hamiltonian)
            genSympleticMatrix();
    }

    auto setQ(const arma::mat& _Q, bool update_hamiltonian = false){
        assert(arma::size(_Q) == arma::size(ss_->A()));
        Q_ = _Q;
        if(update_hamiltonian)
            genSympleticMatrix();
    }

    auto setSympleticMatrix(const arma::mat& _S){
        S_ = _S;
    }

    auto solve();

private:
    auto genSympleticMatrix();

    //-- change to pointer arg for TZ
    auto auxSchur(const arma::mat& _S, std::tuple<arma::cx_mat, arma::cx_mat>* _TZ) -> int{
        ::jacl::py_stuff::AcquireGIL lk;
        try{
            py::object sys = py::import("sys");
            sys.attr("path").attr("append")("../python");

            py::object aux_schur = py::import("aux_schur");
            Py_intptr_t S_shape[2] = {(int)S_.n_rows, (int)S_.n_cols};
            PyObject* np_H = PyArray_SimpleNewFromData(2, S_shape, NPY_FLOAT64, reinterpret_cast<void*>(S_.memptr()));
            py::handle<> S_handle(np_H);
            py::object S_object(S_handle);

            py::object abc = aux_schur.attr("aux_schur")(S_object, "iuc");
            py::stl_input_iterator<py::object> begin(abc),end;
            std::vector<py::object> l1(begin, end);
            arma::mat dum[4];
            for(int i(0); i < l1.size(); i++){
                arma::mat temp(arma::size(_S));
                std::vector<py::object> layer2(py::stl_input_iterator<py::object>(l1[i]), py::stl_input_iterator<py::object>());
                for(int j(0); j < layer2.size(); j++){
                    std::vector<double> layer3(py::stl_input_iterator<double>(layer2[j]), py::stl_input_iterator<double>());
                    for(int k(0); k < layer3.size(); k++)
                        temp(j,k) = layer3[k];
                }
                dum[i] = temp;
            }
            arma::cx_mat T(arma::size(_S)),Z(arma::size(_S));
            T.set_real(dum[0]);
            T.set_imag(dum[1]);
            Z.set_real(dum[2]);
            Z.set_imag(dum[3]);
            *_TZ = std::make_tuple(T,Z);
        }catch(py::error_already_set){
            PyErr_Print();
            return 1;
        }

        return 0;
    }

private:
    _StateSpace* ss_;
    //-- Sympletic matrix
    arma::mat S_;
    //-- Symmetric matrix
    arma::mat R_;
    arma::mat Q_;
    //-- SKew-Symmetric matrix
    arma::mat J_;
    //-- Solution
    arma::mat X_;
};

template <class _StateSpace>
DARE<_StateSpace>::DARE(_StateSpace* _ss, const arma::mat& _R, const arma::mat& _Q)
    : ss_(_ss)
    , S_(_StateSpace::n_states * 2, _StateSpace::n_states * 2)
    , R_(_StateSpace::n_states, _StateSpace::n_states)
    , Q_(_StateSpace::n_states, _StateSpace::n_states){

    setR(_R);
    setQ(_Q, true);
}

template <class _StateSpace>
DARE<_StateSpace>::DARE(_StateSpace* _ss)
    : ss_(_ss)
    , S_(_StateSpace::n_states * 2, _StateSpace::n_states * 2)
    , R_(_StateSpace::n_states, _StateSpace::n_states)
    , Q_(_StateSpace::n_states, _StateSpace::n_states){

}

template <class _StateSpace>
DARE<_StateSpace>::~DARE(){

}

template <class _StateSpace>
auto DARE<_StateSpace>::genSympleticMatrix(){
    /*
    Not code yet
    */
   S_ = arma::join_rows(
       arma::join_cols(ss_->A(), R_),
       arma::join_cols(-Q_, -arma::trans(ss_->A()))
   );
}

template <class _StateSpace>
auto DARE<_StateSpace>::solve(){
    std::tuple<arma::cx_mat, arma::cx_mat> TZ;
    int ret = auxSchur(S_, &TZ);
    arma::cx_mat T = std::get<0>(TZ);
    arma::cx_mat Z = std::get<1>(TZ);
    #ifdef DARE_VERBOSE
    T.print("[DARE] T : ");
    Z.print("[DARE] Z : ");
    #endif

    arma::cx_mat ISS, T1, T2;
    ISS = Z.head_cols(Z.n_cols >> 1);
    #ifdef DARE_VERBOSE
    ISS.print("[DARE] Invariant spectral subspace : ");
    #endif

    T1 = ISS.head_rows(ISS.n_rows >> 1);
    T2 = ISS.tail_rows(ISS.n_rows >> 1);
    #ifdef DARE_VERBOSE
    T1.print("[DARE] T1 : ");
    T2.print("[DARE] T2 : ");
    #endif

    if(arma::cond(T1) > 1e6){
        //-- regularize the ill-conditioned matrix
        T1 = T1 + .001*arma::eye(arma::size(T1));
        std::cout << "[DARE] Regularization triggered" << std::endl;
    }
    arma::cx_mat T1_inv = arma::solve(T1, arma::eye<arma::cx_mat>(arma::size(T1)),
        arma::solve_opts::refine);
    arma::cx_mat solution = T2 * T1_inv;

    return solution;

}

} // namespace jacl
