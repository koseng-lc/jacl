/**
*   @author : koseng (Lintang)
*   @brief : basic linear algebra 1
*/

#include <jacl/jacl.hpp>

int main(int argc, char** argv){
    using a_t = jacl::linear_algebra::Matrix<3,3>;
    using b_t = jacl::linear_algebra::Matrix<3,4>;
    constexpr a_t a(1.,2.,3.,4.,5.,6.,7.,8.,9.);
    constexpr b_t b(1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.);
    constexpr b_t res(a*b);
    if constexpr(res(0)){
        
    }
    for(int i(0); i < b_t::n_elems; i++){
        std::cout << res(i) << "\t" << std::endl;
    }
    for(int i(0); i < a_t::n_elems; i++){
        std::cout << a(i) << "\t" << std::endl;
    }  
    return 0.;
}