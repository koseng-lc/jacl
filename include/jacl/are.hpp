/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of Algebraic Ricatti Solver
*/
#pragma once

#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/bindings/lapack/computational.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <numpy/ndarrayobject.h>

#include <jacl/linear_state_space.hpp>

namespace jacl{

namespace{
    using real_t = double;
    using cx_t = std::complex<real_t>;
    namespace ublas = boost::numeric::ublas;
    namespace bindings = boost::numeric::bindings;
    using v_t = ublas::vector<real_t>;
    using cxv_t = ublas::vector<cx_t>;
    using cxm_t = ublas::matrix<cx_t>;

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

    auto setQ(const arma::mat& _Q, bool update_hamiltonian = false) -> void{
        assert(arma::size(_Q) == arma::size(ss_->A()));
        Q_ = _Q;
        if(update_hamiltonian)
            genHamiltonianMatrix();
    }

    auto setHamiltonianMatrix(const arma::mat& _H) -> void{
        H_ = _H;
    }

    auto solve() -> arma::cx_mat;

private:
    auto genHamiltonianMatrix() -> void;

    //-- create auxiliary lib !!!!!!!!!!!!!
    auto toUBLASMat(const arma::cx_mat& _in) -> cxm_t{
        cxm_t res(_in.n_rows, _in.n_cols);
        for(int i(0); i < _in.n_rows; i++){
            for(int j(0); j < _in.n_cols;j++){
                res(i,j) = _in(i,j);
            }
        }
        return res;
    }

    auto toARMAMat(const cxm_t& _in) -> arma::cx_mat{
        arma::cx_mat res(_in.size1(), _in.size2());
        for(int i(0); i < _in.size1(); i++){
            for(int j(0); j < _in.size2();j++){
                res(i,j) = _in(i,j);
            }
        }
        return res;
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
            TZ = std::make_tuple(T,Z);
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
auto ARE<_StateSpace>::genHamiltonianMatrix() -> void{

    H_.submat(                    0,                     0,     _StateSpace::n_states - 1,     _StateSpace::n_states - 1) = ss_->A();
    H_.submat(                    0, _StateSpace::n_states,     _StateSpace::n_states - 1, _StateSpace::n_states * 2 - 1) = R_;
    H_.submat(_StateSpace::n_states,                     0, _StateSpace::n_states * 2 - 1,     _StateSpace::n_states - 1) = -Q_;
    H_.submat(_StateSpace::n_states, _StateSpace::n_states, _StateSpace::n_states * 2 - 1, _StateSpace::n_states * 2 - 1) = -1*arma::trans(ss_->A());
}

template <class _StateSpace>
auto ARE<_StateSpace>::solve() -> arma::cx_mat{

    //-- Using QR-Algorithm
    /*arma::mat T, U;
    LinearAlgebra::QRAlgorithm(H_, &T, &U);

    arma::mat A, T, U;
    A << 1 << 0 << 0 << arma::endr
      << 0 << cos(M_PI/3) << -sin(M_PI/3) << arma::endr
      << 0 << sin(M_PI/3) << cos(M_PI/3) << arma::endr;

    LinearAlgebra::QRAlgorithm(A, &T, &U);

    T.print("T(ARE) : ");
    U.print("U(ARE) : ");*/
    //--

    /*arma::mat H;
    H << -3 << 2 << 0 << 0 << arma::endr
      << -2 << 1 << 0 << -1 << arma::endr
      << -0 << -0 << 3 << 2 << arma::endr
      << -0 << -0 << -2 << -1 << arma::endr;
    eig_gen(eigval, eigvec, H);*/

    //-- using Lapack Binding from Boost
    /*ublas::vector<unsigned int> sel( t.size1() );
    for(int i(0); i < sel.size(); i++){
        sel(i) = .0;
        if(t(i,i).real() > .0)
            sel(i) = 1.;
    }
    cxv_t w( t.size1() );
    fortran_int_t m;
    double s;
    double sep;
    int info = bindings::lapack::trsen('N','V',sel,t,q,w,m,s,sep);*/

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

//    X1.print("X1 : ");
//    X2.print("X2 : ");
    
    arma::cx_mat solution = X2 * arma::inv(X1);
    //-- Check, are the solution is fulfull the equation
    /*int n = H_.n_rows >> 1;
    arma::cx_mat A = toCx(H_.submat(0,0,n-1,n-1));
    arma::cx_mat R = toCx(H_.submat(0,n,n-1,(n*2)-1));
    arma::cx_mat Q = toCx(-1*H_.submat(n,0,(n*2)-1,n-1));
    arma::cx_mat temp1 = arma::trans(A)*solution;
    arma::cx_mat temp2 = solution*A;
    arma::cx_mat temp3 = solution*R*solution;
    arma::cx_mat temp4 = temp1 + temp2 + temp3;
    arma::cx_mat check = temp4 + Q;
    check.print("CHECK : ");*/

//    arma::cx_vec eigval;
//    arma::cx_mat eigvec;
//    arma::eig_gen(eigval, eigvec, H_);

//    eigval.print("EigenValue : ");
//    eigvec.print("EigenVector : ");

//    arma::mat eigval_re = arma::real(eigval);

//    eigval_re.print("EigVal Re : ");

    //-- invariant spectral subspace
//    arma::cx_mat ISS;

//    for(int i(0); i < eigval.n_rows; i++){
//        if(eigval_re(i, 0) < .0){
//            ISS = join_rows(ISS,eigvec.col(i));
//        }
//    }

//    ISS.print("ISS : ");

//    arma::cx_mat X1, X2;
//    X1 = ISS.head_rows(ISS.n_rows * .5);
//    X2 = ISS.tail_rows(ISS.n_rows * .5);

//    X1.print("X1 : ");
//    X2.print("X2 : ");

//    arma::cx_mat inv_X1( arma::inv( X1 ) );
//    arma::cx_mat solution( X2 * inv_X1 );
//    solution.print("X : ");

    return solution;

}

} // namespace jacl
