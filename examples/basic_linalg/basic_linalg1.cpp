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

    using c_t = jacl::linear_algebra::Matrix<3,1>;
    using d_t = jacl::linear_algebra::Matrix<1,4>;
    constexpr c_t pos2_col(jacl::linear_algebra::getCol<2>(mul_res));
    constexpr d_t pos1_row(jacl::linear_algebra::getRow<1>(mul_res));

    constexpr auto split_res = jacl::linear_algebra::split<0,1,1,3>(mul_res);

    using e_t = b_t::transpose_t;
    constexpr e_t transpose_res(jacl::linear_algebra::transpose(mul_res));

    jacl::linear_algebra::show(a1, "A1 : ");
    jacl::linear_algebra::show(mul_res, "Mult. Res : ");
    jacl::linear_algebra::show(add_res, "Add. Res : ");
    jacl::linear_algebra::show(sub_res, "Sub. Res : ");
    jacl::linear_algebra::show(pos2_col, "Pos.2 Col : ");
    jacl::linear_algebra::show(pos1_row, "Pos.1 Row : ");
    jacl::linear_algebra::show(split_res, "Split Res : ");
    jacl::linear_algebra::show(transpose_res, "Transpose of Mult. Res : ");

    return 0.;
}