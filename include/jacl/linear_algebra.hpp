/**
*   @author : koseng (Lintang)
*   @brief : Some linear algebra tools
*/

#pragma once

#include <vector>

#define ARMA_DONT_USE_WRAPPER
#include <armadillo>

namespace jacl{ namespace linear_algebra{

// template <typename Type>
// static inline auto toCx(const arma::Mat<Type>& _in)
//     -> decltype(arma::Mat<std::complex<Type>>(_in, arma::zeros<arma::Mat<Type>>( arma::size(_in) ))){
//     return arma::Mat<std::complex<Type>>(_in, arma::zeros<arma::Mat<Type>>( arma::size(_in) ));
// }

static inline auto toCx(const arma::mat _in){
    return arma::cx_mat(_in, arma::zeros<arma::mat>( arma::size(_in) ));
}

template <typename Type>
static inline auto toReal(const arma::Mat<std::complex<Type>>& _in)
    -> decltype(arma::real(_in)){
    return arma::real(_in);
}

static void GramSchmidtProjection(int idx, const arma::mat& col, arma::mat* u, std::vector<arma::mat>* e){
    if(idx > 0){
        auto norm = arma::dot(col, (*e)[idx-1]);
        *u -= (norm * (*e)[idx-1]);
        GramSchmidtProjection(idx - 1, col, u, e);
    }
}

static auto QRDecomp(const arma::mat& in, arma::mat* Q, arma::mat* R){

    *Q = arma::mat(in.n_rows, in.n_cols, arma::fill::zeros);
    *R = arma::mat(arma::size(in), arma::fill::zeros);
    std::vector<arma::mat> e;

//    int n_rank = arma::rank(in);
    arma::mat u;

    //-- Gram-Schmidt Process
    auto norm(.0);
    arma::mat normalized_u;
    for(int i(0); i < in.n_cols; i++){
        u = in.col(i);
        GramSchmidtProjection(i, in.col(i), &u, &e);
        norm = arma::norm(u, 2);
        normalized_u = (norm == .0) ? u : u / norm;
        e.push_back(normalized_u);
        Q->col(i) = e[i];
    }

    for(int r(0); r < in.n_cols; r++){ //-- limit row iteration
        for(int c(r); c < in.n_cols; c++){
            (*R)(r, c) = arma::dot(in.col(c), Q->col(r));
        }
    }
}

static auto QRAlgorithm(const arma::mat& in, arma::mat* T, arma::mat* U, int num_iter = 20){
    arma::mat updated_in(in);
    *U = arma::mat(arma::size(in), arma::fill::eye);
    arma::mat Q, R;
    for(int i(0); i < num_iter; i++){
//        updated_in.print("UIN : ");
        QRDecomp(updated_in, &Q, &R);
//        Q.print("Q : ");
//        R.print("R : ");
        updated_in = R * Q;
        (*U) = (*U) * Q;
    }

    *T = updated_in;
}

template <typename Type>
static auto isPosSemiDefinite(const arma::Mat<Type> in){
    arma::Mat<Type> hessian_in = .5*(in + arma::trans(in));
    arma::Mat<Type> R;
    return arma::chol(R, hessian_in);
}

template <typename Type>
static auto spectralRadius(const arma::Mat<Type>& in){
    arma::Mat<Type> eigvec;
    arma::Col<Type> eigval;
    arma::eig_gen(eigval, eigvec, in);
    arma::vec eigval_re = arma::real(eigval);
    arma::vec eigval_im = arma::imag(eigval);
    auto max_mag_eigval(.0);
    auto mag_eigval(.0);
    for(int i(0); i < eigval.n_rows; i++){
        mag_eigval = std::sqrt(eigval_re(i)*eigval_re(i) + eigval_im(i)*eigval_im(i));
        if(mag_eigval > max_mag_eigval)
            max_mag_eigval = mag_eigval;
    }
    return max_mag_eigval;
}

    // template <typename _StateSpace>
    // auto isInfNormLessThan(double _gam, const _StateSpace& _ss){
    //     arma::mat C_t = arma::trans( _ss.C() );
    //     arma::mat D_t = arma::trans( _ss.D() );
    //     arma::mat D_tD = D_t * _ss.D();
    //     arma::mat R = (_gam*_gam)*arma::eye(_StateSpace::n_inputs,
    //                                         _StateSpace::n_inputs)
    //                     - D_tD;
    //     arma::mat R_inv = arma::inv(R);

    //     //-- Hamiltonian matrix
    //     arma::mat H(_StateSpace::n_states << 1,
    //                 _StateSpace::n_states << 1);

    //     arma::mat temp1, temp2;

    //     //-- block 1,1
    //     temp1 = R_inv*D_t*_ss.C();
    //     temp2 = _ss.A() + temp1;
    //     H.submat(                         0,                          0,
    //              _StateSpace::n_states - 1, _StateSpace::n_states - 1) = temp2;
    //     //-- block 2,2
    //     H.submat(_StateSpace::n_states, _StateSpace::n_states, _StateSpace::n_states * 2 - 1, _StateSpace::n_states * 2 - 1) = -arma::trans(temp2);
    //     //-- block 1,2
    //     H.submat(                         0,         _StateSpace::n_states,
    //              _StateSpace::n_states - 1, _StateSpace::n_states * 2 - 1) = _ss.B()*R_inv*arma::trans(_ss.B());
    //     //-- block 2,1
    //     temp1 = _ss.D()*R_inv*D_t;
    //     temp2 = arma::eye(_StateSpace::n_outputs,
    //                       _StateSpace::n_outputs)
    //              - temp1;
    //     H.submat(_StateSpace::n_states, 0, _StateSpace::n_states * 2 - 1, _StateSpace::n_states - 1) = -C_t*temp2*_ss.C();

    //     arma::cx_mat eigvec;
    //     arma::cx_vec eigval;
    //     arma::eig_gen(eigval, eigvec, H);
    //     arma::mat eigval_re = arma::real(eigval);
    //     auto ok(true);
    //     for(int i(0); i < eigval.n_rows; i++){
    //         if(std::fabs(eigval_re(i, 0)) < std::numeric_limits<double>::epsilon()){
    //             ~ok;
    //             break;
    //         }
    //     }
    //     return ok;
    // }

// template <typename _StateSpace>
    // auto isInfNormLessThan(double _gam, const _StateSpace& _ss){
    //     arma::mat C_t = arma::trans( _ss.C() );
    //     arma::mat D_t = arma::trans( _ss.D() );
    //     arma::mat D_tD = D_t * _ss.D();
    //     arma::mat R = (_gam*_gam)*arma::eye(_StateSpace::n_inputs,
    //                                         _StateSpace::n_inputs)
    //                     - D_tD;
    //     arma::mat R_inv = arma::inv(R);        

    //     arma::mat temp1, temp2, temp3, temp4;

    //     //-- block 1,1
    //     temp1 = R_inv*D_t*_ss.C();
    //     temp2 = _ss.A() + temp1;
    //     temp3 = _ss.D()*R_inv*D_t;
    //     temp4 = arma::eye(_StateSpace::n_outputs,
    //                       _StateSpace::n_outputs)
    //              - temp3;
    //     //-- Sympletic matrix
    //     arma::mat H = arma::join_cols(
    //         arma::join_rows(temp2, _ss.B()*R_inv*arma::trans(_ss.B())),
    //         arma::join_rows(-C_t*temp4*_ss.C(), -arma::trans(temp2))
    //     );

    //     arma::cx_mat eigvec;
    //     arma::cx_vec eigval;
    //     arma::eig_gen(eigval, eigvec, H);
    //     arma::mat eigval_re = arma::real(eigval);
    //     auto ok(true);
    //     for(int i(0); i < eigval.n_rows; i++){
    //         if(std::fabs(eigval_re(i, 0)) < std::numeric_limits<double>::epsilon()){
    //             ~ok;
    //             break;
    //         }
    //     }
    //     return ok;
    // }

} } // namespace::linear_algebra
