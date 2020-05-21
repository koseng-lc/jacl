/**
*   @author : koseng (Lintang)
*   @brief : jacl random generator wrapper
*/

#pragma once

#include <random>
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>

namespace jacl{

namespace detail{
    template <typename T>
    struct is_prob_dist:std::false_type{};

    template <typename Scalar>
    struct is_prob_dist<std::normal_distribution<Scalar>>:std::true_type{};
}

template<typename Scalar, std::size_t dim, typename Distribution>
class Random{
public:
    using type_t = arma::vec::fixed<dim>;
    using cov_t = arma::mat::fixed<dim, dim>;

public:
    explicit Random(const type_t& _param1, const cov_t _param2)
        : rand_gen_(std::random_device{}())
        , dist_(0, 1)
        , mean_(_param1)
        , cov_(_param2){}

    ~Random(){}

    inline auto operator()() -> type_t{
        genNormal();
        transform();

        return result_;
    }

private:
    std::mt19937 rand_gen_;
    Distribution dist_;
    type_t result_;
    type_t mean_;
    cov_t cov_;

    auto genNormal() -> type_t{
        for(auto& r:result_)
            r = dist_(rand_gen_);
    }

    auto transform() -> type_t{
        arma::cx_vec::fixed<dim> eigval;
        arma::cx_mat::fixed<dim,dim> eigvec;
        arma::eig_gen(eigval, eigvec, cov_);
        arma::cx_vec::fixed<dim> temp1; temp1.set_real(result_);
        arma::cx_vec::fixed<dim> temp2 = eigvec * temp1;
        result_ = arma::real(arma::sqrtmat(arma::diagmat(eigval)) * temp2);
        result_ = result_ + mean_;
    }    

};

template<typename Scalar, typename Distribution>
class Random<Scalar,1,Distribution>{
public:
    using type_t = Scalar;
    using cov_t = Scalar;

public:
    explicit Random(Scalar _param1)
        : rand_gen_(std::random_device{}())
        , dist_(_param1){}

    explicit Random(Scalar _param1, Scalar _param2)
        : rand_gen_(std::random_device{}())
        , dist_(_param1, _param2){}

    Random(){}
    ~Random(){}    

    inline auto operator()() -> Scalar{
        return dist_(rand_gen_);
    }

private:
    std::mt19937 rand_gen_;
    Distribution dist_;

};

}