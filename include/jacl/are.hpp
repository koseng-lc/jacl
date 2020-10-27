/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of continous algebraic Ricatti equations solver
*/

#pragma once

#include <jacl/state_space/linear.hpp>

// #define ARE_VERBOSE

namespace jacl::are{

namespace detail{

    template <typename H>
    static auto solve(H h){
        std::tuple<arma::cx_mat, arma::cx_mat> tz;
        int ret = ::jacl::linear_algebra::auxSchur(h, &tz, "lhp");
        arma::cx_mat t = std::get<0>(tz);
        arma::cx_mat z = std::get<1>(tz);
        #ifdef ARE_VERBOSE
        t.print("[ARE] T : ");
        z.print("[ARE] Z : ");
        #endif

        arma::cx_mat iss, x1, x2;
        iss = z.head_cols(z.n_cols >> 1);
        #ifdef ARE_VERBOSE
        iss.print("[ARE] Invariant spectral subspace : ");
        #endif

        x1 = iss.head_rows(iss.n_rows >> 1);
        x2 = iss.tail_rows(iss.n_rows >> 1);
        #ifdef ARE_VERBOSE
        x1.print("[ARE] X1 : ");
        x2.print("[ARE] X2 : ");
        #endif
        
        if(arma::cond(x1) > 1e6){
            //-- regularize the ill-conditioned matrix
            x1 = x1 + .001*arma::eye(arma::size(x1));
            std::cout << "[ARE] Regularization triggered" << std::endl;
        }
        arma::cx_mat x1_inv = arma::solve(x1, arma::eye<arma::cx_mat>(arma::size(x1)),
            arma::solve_opts::refine);
        arma::cx_mat solution = x2 * x1_inv;

        return solution;
    }
}

template <typename H>
static auto solve(H h){
    return detail::solve(h);
}

template <typename A, typename R, typename Q>
static auto solve(A a, R r, Q q){
    return detail::solve(
        arma::join_rows(
            arma::join_cols(a,r),
            arma::join_cols(-q,-arma::trans(a))
        )
    );
}

} // namespace jacl::are
