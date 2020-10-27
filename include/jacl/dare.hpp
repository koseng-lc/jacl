/**
*   @author : koseng (Lintang)
*   @brief : Simple implementation of discrete algebraic Riccati equations solver
*/

#pragma once

#include <jacl/state_space/linear.hpp>

// #define DARE_VERBOSE

namespace jacl::dare{

namespace detail{

    template <typename S>
    static auto solve(S s){
        std::tuple<arma::cx_mat, arma::cx_mat> tz;
        int ret = ::jacl::linear_algebra::auxSchur(s, &tz, "iuc");
        arma::cx_mat t = std::get<0>(tz);
        arma::cx_mat z = std::get<1>(tz);
        #ifdef DARE_VERBOSE
        t.print("[DARE] T : ");
        z.print("[DARE] Z : ");
        #endif

        arma::cx_mat iss, t1, t2;
        iss = z.head_cols(z.n_cols >> 1);
        #ifdef DARE_VERBOSE
        iss.print("[DARE] Invariant spectral subspace : ");
        #endif

        t1 = iss.head_rows(iss.n_rows >> 1);
        t2 = iss.tail_rows(iss.n_rows >> 1);
        #ifdef DARE_VERBOSE
        t1.print("[DARE] T1 : ");
        t2.print("[DARE] T2 : ");
        #endif

        if(arma::cond(t1) > 1e6){
            //-- regularize the ill-conditioned matrix
            t1 = t1 + .001*arma::eye(arma::size(t1));
            std::cout << "[DARE] Regularization triggered" << std::endl;
        }
        arma::cx_mat t1_inv = arma::solve(t1, arma::eye<arma::cx_mat>(arma::size(t1)),
            arma::solve_opts::refine + arma::solve_opts::equilibrate + arma::solve_opts::allow_ugly);
        arma::cx_mat solution = t2 * t1_inv;

        return solution;
    }
}

template <typename S>
static auto solve(S s){
    return detail::solve(s);
}

template <typename A, typename G, typename Q>
static auto solve(A a, G g, Q q){
    auto a_t_inv = arma::inv( arma::trans(a) );
    return detail::solve(
        arma::join_rows(
            arma::join_cols(a + g*a_t_inv*q,-g*a_t_inv),
            arma::join_cols(-a_t_inv*q,a_t_inv)
        )
    );
}

} // namespace jacl::dare
