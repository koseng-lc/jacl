/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of continous algebraic Ricatti equations solver
*/

#pragma once

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <numpy/ndarrayobject.h>

#include <jacl/state_space/linear.hpp>

// #define ARE_VERBOSE

namespace jacl{

namespace{
    namespace py = boost::python;
    namespace np = boost::python::numpy;
}

template <class _StateSpace>
class ARE{
public:
    ARE(_StateSpace* _ss, const arma::mat& _R, const arma::mat& _Q);
    ARE(_StateSpace* _ss);
    ~ARE();

    ARE(const ARE&) = delete;
    ARE& operator=(const ARE&) = delete;

    auto setR(const arma::mat& _R, bool update_hamiltonian = false){
        assert(arma::size(_R) == arma::size(ss_->A()));
        R_ = _R;
        if(update_hamiltonian)
            genHamiltonianMatrix();
    }

    auto setQ(const arma::mat& _Q, bool update_hamiltonian = false){
        assert(arma::size(_Q) == arma::size(ss_->A()));
        Q_ = _Q;
        if(update_hamiltonian)
            genHamiltonianMatrix();
    }

    auto setHamiltonianMatrix(const arma::mat& _H){
        H_ = _H;
    }

    auto solve() -> arma::cx_mat;

private:
    auto genHamiltonianMatrix();

    auto auxSchur(const arma::mat& _H, std::tuple<arma::cx_mat, arma::cx_mat>* _TZ){
        ::jacl::py_stuff::AcquireGIL lk;
        try{
            py::object sys = py::import("sys");
            sys.attr("path").attr("append")("../python");

            py::object aux_schur = py::import("aux_schur");
            Py_intptr_t H_shape[2] = {(int)H_.n_rows, (int)H_.n_cols};
            PyObject* np_H = PyArray_SimpleNewFromData(2, H_shape, NPY_FLOAT64, reinterpret_cast<void*>(H_.memptr()));
            py::handle<> H_handle(np_H);
            py::object H_object(H_handle);

            py::object abc = aux_schur.attr("aux_schur")(H_object, "lhp");
            py::stl_input_iterator<py::object> begin(abc),end;
            std::vector<py::object> l1(begin, end);
            arma::mat dum[4];
            for(int i(0); i < l1.size(); i++){
                arma::mat temp(arma::size(_H));
                std::vector<py::object> layer2(py::stl_input_iterator<py::object>(l1[i]), py::stl_input_iterator<py::object>());
                for(int j(0); j < layer2.size(); j++){
                    std::vector<double> layer3(py::stl_input_iterator<double>(layer2[j]), py::stl_input_iterator<double>());
                    for(int k(0); k < layer3.size(); k++)
                        temp(j,k) = layer3[k];
                }
                dum[i] = temp;
            }
            arma::cx_mat T(arma::size(_H)),Z(arma::size(_H));
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
    //-- Hamiltonian matrix
    arma::mat H_;
    //-- Symmetric matrix
    arma::mat R_;
    arma::mat Q_;
    //-- SKew-Symmetric matrix
    arma::mat J_;
    //-- Solution
    arma::mat X_;
};

template <class _StateSpace>
ARE<_StateSpace>::ARE(_StateSpace* _ss, const arma::mat& _R, const arma::mat& _Q)
    : ss_(_ss)
    , H_(_StateSpace::n_states * 2, _StateSpace::n_states * 2)
    , R_(_StateSpace::n_states, _StateSpace::n_states)
    , Q_(_StateSpace::n_states, _StateSpace::n_states){

    setR(_R);
    setQ(_Q, true);
}

template <class _StateSpace>
ARE<_StateSpace>::ARE(_StateSpace* _ss)
    : ss_(_ss)
    , H_(_StateSpace::n_states * 2, _StateSpace::n_states * 2)
    , R_(_StateSpace::n_states, _StateSpace::n_states)
    , Q_(_StateSpace::n_states, _StateSpace::n_states){

}

template <class _StateSpace>
ARE<_StateSpace>::~ARE(){

}

template <class _StateSpace>
auto ARE<_StateSpace>::genHamiltonianMatrix(){
    H_ = arma::join_cols(
        arma::join_rows(ss_->A(),R_),
        arma::join_rows(-Q_,-arma::trans(ss_->A()))
    );
}

template <class _StateSpace>
auto ARE<_StateSpace>::solve() -> arma::cx_mat{
    std::tuple<arma::cx_mat, arma::cx_mat> TZ;
    int ret = auxSchur(H_, &TZ);
    arma::cx_mat T = std::get<0>(TZ);
    arma::cx_mat Z = std::get<1>(TZ);
    #ifdef ARE_VERBOSE
    T.print("[ARE] T : ");
    Z.print("[ARE] Z : ");
    #endif

    arma::cx_mat ISS, X1, X2;
    ISS = Z.head_cols(Z.n_cols >> 1);
    #ifdef ARE_VERBOSE
    ISS.print("[ARE] Invariant spectral subspace : ");
    #endif

    X1 = ISS.head_rows(ISS.n_rows >> 1);
    X2 = ISS.tail_rows(ISS.n_rows >> 1);
    #ifdef ARE_VERBOSE
    X1.print("[ARE] X1 : ");
    X2.print("[ARE] X2 : ");
    #endif
    
    if(arma::cond(X1) > 1e6){
        //-- regularize the ill-conditioned matrix
        X1 = X1 + .001*arma::eye(arma::size(X1));
        std::cout << "[ARE] Regularization triggered" << std::endl;
    }
    arma::cx_mat X1_inv = arma::solve(X1, arma::eye<arma::cx_mat>(arma::size(X1)),
        arma::solve_opts::refine);
    arma::cx_mat solution = X2 * X1_inv;

    return solution;

}

} // namespace jacl
