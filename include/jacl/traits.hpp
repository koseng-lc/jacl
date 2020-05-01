#pragma once

#include <type_traits>

#include <jacl/linear_state_space.hpp>
#include <jacl/nonlinear_state_space.hpp>
#include <jacl/system/continuous_system.hpp>
#include <jacl/system/discrete_system.hpp>

namespace jacl{ namespace traits{
    template <typename T>
    struct is_state_space:std::false_type{};

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no>
    struct is_state_space<LinearStateSpace<Scalar,ns,ni,no> >:std::true_type{};        

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no,
                class PhysicalParam, class  ...Rest>
    struct is_state_space<LinearStateSpace<Scalar,ns,ni,no,PhysicalParam,Rest...> >:std::true_type{};

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no>
    struct is_state_space<NonLinearStateSpace<Scalar,ns,ni,no> >:std::true_type{};

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no,
                class PhysicalParam, class  ...Rest>
    struct is_state_space<NonLinearStateSpace<Scalar,ns,ni,no,PhysicalParam,Rest...> >:std::true_type{};

    //--
    template <typename T>
    struct is_linear_state_space:std::false_type{};

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no>
    struct is_linear_state_space<LinearStateSpace<Scalar,ns,ni,no> >:std::true_type{};

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no,
                class PhysicalParam, class  ...Rest>
    struct is_linear_state_space<LinearStateSpace<Scalar,ns,ni,no,PhysicalParam,Rest...> >:std::true_type{};

    //--
    template <typename T>
    struct is_nonlinear_state_space:std::false_type{};

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no>
    struct is_nonlinear_state_space<NonLinearStateSpace<Scalar,ns,ni,no> >:std::true_type{};

    template <typename Scalar, std::size_t ns, std::size_t ni, std::size_t no,
                class PhysicalParam, class  ...Rest>
    struct is_nonlinear_state_space<NonLinearStateSpace<Scalar,ns,ni,no,PhysicalParam,Rest...> >:std::true_type{};

    //--
    template <typename T>
    struct is_continuous_system:std::false_type{};

    template <typename _StateSpace>
    struct is_continuous_system<::jacl::system::ContinuousSystem<_StateSpace> >:std::true_type{};

    //--
    template <typename T>
    struct is_discrete_system:std::false_type{};

    template <typename _StateSpace>
    struct is_discrete_system<::jacl::system::DiscreteSystem<_StateSpace> >:std::true_type{};

    //--
    template <typename T>
    struct is_siso:std::false_type{};

    template <typename Scalar, std::size_t ns>
    struct is_siso<LinearStateSpace<Scalar,ns,1,1>>:std::true_type{};        

    template <typename Scalar, std::size_t ns, class PhysicalParam, class  ...Rest>
    struct is_siso<LinearStateSpace<Scalar,ns,1,1,PhysicalParam,Rest...>>:std::true_type{};

    template <typename Scalar, std::size_t ns>
    struct is_siso<NonLinearStateSpace<Scalar,ns,1,1>>:std::true_type{};

    template <typename Scalar, std::size_t ns, class PhysicalParam, class  ...Rest>
    struct is_siso<NonLinearStateSpace<Scalar,ns,1,1,PhysicalParam,Rest...>>:std::true_type{};

    //-- convert to underlying type
    //-- the use of function instead of struct to prevent explicitly write the T
    template <typename T>
    constexpr auto toUType(T t) -> typename std::underlying_type<T>::type{
        return static_cast<typename std::underlying_type<T>::type>(t);
    }
} }
