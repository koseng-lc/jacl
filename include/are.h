#pragma once

#include "state_space.h"
#include "linear_algebra.h"

namespace JACL{

namespace{
    namespace linalg = linear_algebra;
}

template <class SSpace>
class ARE{
public:
    ARE(SSpace* _ss, const arma::mat& _R, const arma::mat& _Q);
    ARE(SSpace* _ss);
    ~ARE();

    void setR(const arma::mat& _R, bool update_hamiltonian = false){
        assert(arma::size(_R) == arma::size(ss_->A()));
        R_ = _R;
        if(update_hamiltonian)
            genHamiltonianMatrix();
    }

    void setQ(const arma::mat& _Q, bool update_hamiltonian = false){
        assert(arma::size(_Q) == arma::size(ss_->A()));
        Q_ = _Q;
        if(update_hamiltonian)
            genHamiltonianMatrix();
    }

    void solve();

private:
    SSpace* ss_;

    //-- Hamiltonian arma::matrix
    arma::mat H_;

    //-- Symmetric arma::matrix
    arma::mat R_;
    arma::mat Q_;

    //-- SKew-Symmetric arma::matrix
    arma::mat J_;

    //-- Solution
    arma::mat X_;

    void genHamiltonianMatrix();

};

template <class SSpace>
ARE<SSpace>::ARE(SSpace* _ss, const arma::mat& _R, const arma::mat& _Q)
    : ss_(_ss)
    , H_(ss_->numStates() * 2, ss_->numStates() * 2)
    , R_(ss_->numStates(), ss_->numStates())
    , Q_(ss_->numStates(), ss_->numStates()){

    setR(_R);
    setQ(_Q, true);
}

template <class SSpace>
ARE<SSpace>::ARE(SSpace* _ss)
    : ss_(_ss)
    , H_(ss_->numStates() * 2, ss_->numStates() * 2)
    , R_(ss_->numStates(), ss_->numStates())
    , Q_(ss_->numStates(), ss_->numStates()){

}

template <class SSpace>
ARE<SSpace>::~ARE(){

}

template <class SSpace>
void ARE<SSpace>::genHamiltonianMatrix(){

    H_.submat(               0,                0,      ss_->numStates(),     ss_->numStates()) = ss_->A();
    H_.submat(               0, ss_->numStates(),      ss_->numStates(), ss_->numStates() * 2) = R_;
    H_.submat(ss_->numStates(),                0, ss_->numStates() * 2,      ss_->numStates()) = -Q_;
    H_.submat(ss_->numStates(), ss_->numStates(), ss_->numStates() * 2,  ss_->numStates() * 2) = -ss_->A_.t();
}

template <class SSpace>
void ARE<SSpace>::solve(){

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
    eig_gen(eigval, eigvec, H_);
//    eigval.print("EigenValue : ");
//    eigvec.print("EigenVector : ");

    arma::mat eigval_re = arma::real(eigval);

    eigval_re.print("Eval Re : ");

    //-- invariant spectral subspace
    arma::cx_mat ISS;

    for(int i(0); i < eigval.n_rows; i++){
        if(eigval_re(i, 0) < .0){
            ISS = join_rows(ISS,eigvec.col(i));
        }
    }

    ISS.print("ISS : ");

    arma::cx_mat X1, X2;
    X1 = ISS.submat(             0, 0, ISS.n_rows * 2, ISS.n_cols);
    X2 = ISS.submat(ISS.n_rows * 2, 0,     ISS.n_rows, ISS.n_cols);

    arma::mat inv_X1(arma::inv(X1));
    arma::mat solution(X2*inv_X1);
    solution.print("X : ");

}

}
