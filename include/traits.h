#pragma once

#include <type_traits>

#include "state_space.h"

namespace jacl{
    namespace traits{
        template <typename T>
        struct is_state_space:std::false_type{};

        template <std::size_t ns, std::size_t ni, std::size_t no>
        struct is_state_space<StateSpace<ns,ni,no> >:std::true_type{};

        //-- convert to underlying type
        //-- the use of function instead of struct to prevent explicitly write the T
        template <typename T>
        constexpr auto toUType(T t) -> typename std::underlying_type<T>::type{
            return static_cast<typename std::underlying_type<T>::type>(t);
        }
    }
}
