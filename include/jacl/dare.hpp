/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of discrete Algebraic Riccati Equations solver
*/
#pragma once

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <numpy/ndarrayobject.h>

#include <jacl/linear_state_space.hpp>

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

    auto setQ(const arma::mat& _Q, bool update_hamiltonian = false) -> void{
        assert(arma::size(_Q) == arma::size(ss_->A()));
        Q_ = _Q;
        if(update_hamiltonian)
            genSympleticMatrix();
    }

    auto setSympleticMatrix(const arma::mat& _H){
        H_ = _H;
    }

    auto solve();

private:
    auto genSympleticMatrix();

    //-- create utils lib
    inline auto toCx(const arma::mat& _in) const -> arma::cx_mat{
        return arma::cx_mat(_in, arma::zeros<arma::mat>( arma::size(_in) ));
    }

    //-- change to pointer arg for TZ
    auto auxSchur(const arma::mat& _H, std::tuple<arma::cx_mat, arma::cx_mat>& TZ) -> int{
        ::jacl::py_stuff::AcquireGIL lk;
        try{
            py::object sys = py::import("sys");
            sys.attr("path").attr("append")("../python");

            py::object aux_schur = py::import("aux_schur");
            Py_intptr_t H_shape[2] = {(int)H_.n_rows, (int)H_.n_cols};
            PyObject* np_H = PyArray_SimpleNewFromData(2, H_shape, NPY_FLOAT64, reinterpret_cast<void*>(H_.memptr()));
            py::handle<> H_handle(np_H);
            py::object H_object(H_handle);

            py::object abc = aux_schur.attr("aux_schur")(H_object, "iuc");
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
            TZ = std::make_tuple(T,Z);
        }catch(py::error_already_set){
            PyErr_Print();
            return 1;
        }

        return 0;
    }

private:
    _StateSpace* ss_;
    //-- Sympletic matrix
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
DARE<_StateSpace>::DARE(_StateSpace* _ss, const arma::mat& _R, const arma::mat& _Q)
    : ss_(_ss)
    , H_(_StateSpace::n_states * 2, _StateSpace::n_states * 2)
    , R_(_StateSpace::n_states, _StateSpace::n_states)
    , Q_(_StateSpace::n_states, _StateSpace::n_states){

    setR(_R);
    setQ(_Q, true);
}

template <class _StateSpace>
DARE<_StateSpace>::DARE(_StateSpace* _ss)
    : ss_(_ss)
    , H_(_StateSpace::n_states * 2, _StateSpace::n_states * 2)
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
   H_ = arma::join_rows(
       arma::join_cols(ss_->A(), R_),
       arma::join_cols(-Q_, -arma::trans(ss_->A()))
   );
}

template <class _StateSpace>
auto DARE<_StateSpace>::solve(){

    std::tuple<arma::cx_mat, arma::cx_mat> TZ;
    int ret = auxSchur(H_, TZ);
    arma::cx_mat T = std::get<0>(TZ);
    arma::cx_mat Z = std::get<1>(TZ);

    arma::cx_mat ISS, X1, X2;
//    T.print("T : ");
//    Z.print("Z : ");
    ISS = Z.head_cols(Z.n_cols >> 1);
//    ISS.print("ISS : ");
    X1 = ISS.head_rows(ISS.n_rows >> 1);
    X2 = ISS.tail_rows(ISS.n_rows >> 1);
    
    arma::cx_mat solution = X2 * arma::inv(X1);

    return solution;

}

} // namespace jacl
