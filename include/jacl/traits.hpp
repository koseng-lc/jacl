/**
*   @author : koseng (Lintang)
*   @brief : jacl type traits
*/

#pragma once

#include <type_traits>

#include <jacl/system/discrete.hpp>
#include <jacl/system/continuous.hpp>
#include <jacl/state_space/linear.hpp>
#include <jacl/state_space/nonlinear.hpp>

namespace jacl::traits{

    template <typename T>
    struct is_state_space:std::false_type{};

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no>
    struct is_state_space<::jacl::state_space::Linear<Scalar,ns,ni,no> >:std::true_type{};    

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no,
                class PhysicalParam, class  ...Rest>
    struct is_state_space<::jacl::state_space::Linear<Scalar,ns,ni,no,PhysicalParam,Rest...> >:std::true_type{};

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no>
    struct is_state_space<::jacl::state_space::NonLinear<Scalar,ns,ni,no> >:std::true_type{};

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no,
                class PhysicalParam, class  ...Rest>
    struct is_state_space<::jacl::state_space::NonLinear<Scalar,ns,ni,no,PhysicalParam,Rest...> >:std::true_type{};

    //-- helper variable
    template <typename T>
    constexpr bool is_state_space_v = is_state_space<T>::value;

    //--
    template <typename T>
    struct is_linear_state_space:std::false_type{};

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no>
    struct is_linear_state_space<::jacl::state_space::Linear<Scalar,ns,ni,no> >:std::true_type{};

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no,
                class PhysicalParam, class  ...Rest>
    struct is_linear_state_space<::jacl::state_space::Linear<Scalar,ns,ni,no,PhysicalParam,Rest...> >:std::true_type{};

    template <typename T>
    constexpr bool is_linear_state_space_v = is_linear_state_space<T>::value;

    //--
    template <typename T>
    struct is_nonlinear_state_space:std::false_type{};

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no>
    struct is_nonlinear_state_space<::jacl::state_space::NonLinear<Scalar,ns,ni,no> >:std::true_type{};

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no,
                class PhysicalParam, class  ...Rest>
    struct is_nonlinear_state_space<::jacl::state_space::NonLinear<Scalar,ns,ni,no,PhysicalParam,Rest...> >:std::true_type{};

    template <typename T>
    constexpr bool is_nonlinear_state_space_v = is_nonlinear_state_space<T>::value;

    //--
    template <typename T>
    struct is_system:std::is_base_of<::jacl::system::BaseSystem<typename T::state_space_t>, T>{};

    template <typename T>
    constexpr bool is_system_v = is_system<T>::value;

    //--
    template <typename T>
    struct is_continuous_system:std::false_type{};

    template <typename _StateSpace>
    struct is_continuous_system<::jacl::system::Continuous<_StateSpace> >:std::true_type{};

    template <typename T>
    constexpr bool is_continuous_system_v = is_continuous_system<T>::value;

    //--
    template <typename T>
    struct is_discrete_system:std::false_type{};

    template <typename _StateSpace>
    struct is_discrete_system<::jacl::system::Discrete<_StateSpace> >:std::true_type{};

    template <typename T>
    constexpr bool is_discrete_system_v = is_discrete_system<T>::value;

    //--
    template <typename T>
    struct is_siso:std::false_type{};

    template <typename Scalar, std::size_t ns>
    struct is_siso<::jacl::state_space::Linear<Scalar,ns,1,1>>:std::true_type{};        

    template <typename Scalar, std::size_t ns, class PhysicalParam, class  ...Rest>
    struct is_siso<::jacl::state_space::Linear<Scalar,ns,1,1,PhysicalParam,Rest...>>:std::true_type{};

    template <typename Scalar, std::size_t ns>
    struct is_siso<::jacl::state_space::NonLinear<Scalar,ns,1,1>>:std::true_type{};

    template <typename Scalar, std::size_t ns, class PhysicalParam, class  ...Rest>
    struct is_siso<::jacl::state_space::NonLinear<Scalar,ns,1,1,PhysicalParam,Rest...>>:std::true_type{};

    template <typename T>
    constexpr bool is_siso_v = is_siso<T>::value;

    //--
    template <typename T>
    struct is_complex:std::false_type{};

    template <typename T>
    struct is_complex<std::complex<T>>:std::true_type{};

    template <typename T>
    constexpr bool is_complex_v = is_complex<T>::value;

    //-- convert to underlying type
    //-- the use of function instead of struct to prevent explicitly write the T
    template <typename T>
    constexpr auto toUType(T t) -> typename std::underlying_type<T>::type{
        return static_cast<typename std::underlying_type<T>::type>(t);
    }
}
