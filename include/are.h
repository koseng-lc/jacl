#pragma once

#include "state_space.h"
#include "linear_algebra.h"

namespace jacl{

template <class _StateSpace>
class ARE{
public:
    ARE(_StateSpace* _ss, const arma::mat& _R, const arma::mat& _Q);
    ARE(_StateSpace* _ss);
    ~ARE();

    auto setR(const arma::mat& _R, bool update_hamiltonian = false) -> void{
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
    H_.submat(_StateSpace::n_states, _StateSpace::n_states, _StateSpace::n_states * 2 - 1, _StateSpace::n_states * 2 - 1) = -ss_->A_.t();
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

    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eig_gen(eigval, eigvec, H_);
//    eigval.print("EigenValue : ");
//    eigvec.print("EigenVector : ");

    arma::mat eigval_re = arma::real(eigval);

//    eigval_re.print("EigVal Re : ");

    //-- invariant spectral subspace
    arma::cx_mat ISS;

    for(int i(0); i < eigval.n_rows; i++){
        if(eigval_re(i, 0) < .0){
            ISS = join_rows(ISS,eigvec.col(i));
        }
    }

//    ISS.print("ISS : ");

    arma::cx_mat X1, X2;
    X1 = ISS.head_rows(ISS.n_rows * .5);
    X2 = ISS.tail_rows(ISS.n_rows * .5);

//    X1.print("X1 : ");
//    X2.print("X2 : ");

    arma::cx_mat inv_X1( arma::inv( X1 ) );
    arma::cx_mat solution( X2 * inv_X1 );
//    solution.print("X : ");

    return solution;

}

}
