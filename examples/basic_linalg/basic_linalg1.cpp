/**
*   @author : koseng (Lintang)
*   @brief : Basic linear algebra 1
*/

#include <jacl/jacl.hpp>

int main(int argc, char** argv){
    using a_t = jacl::linear_algebra::Matrix<3,3>;
    using b_t = jacl::linear_algebra::Matrix<3,4>;
    constexpr a_t a1(1.,2.,3.,
                    4.,5.,6.,
                    7.,8.,9.);
    constexpr a_t a2(2.,3.,4.,
                    5.,6.,7.,
                    8.,9.,10.);
    constexpr b_t b(1.,2.,3.,4.,
                    5.,6.,7.,8.,
                    9.,10.,11.,12.);
    constexpr b_t mul_res(a1*b);
    constexpr a_t add_res(a1+a2);
    constexpr a_t sub_res(a2-a1);

    jacl::linear_algebra::show(a1, "A1 : ");
    jacl::linear_algebra::show(mul_res, "Mult. Res : ");
    jacl::linear_algebra::show(add_res, "Add. Res : ");
    jacl::linear_algebra::show(sub_res, "Sub. Res : ");

    return 0.;
}