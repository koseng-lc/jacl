#pragma once

#include <type_traits>

#include "state_space.h"

namespace jacl{
    namespace traits{
        template <typename T>
        struct is_state_space:std::false_type{};

        template <int ns, int ni, int no>
        struct is_state_space<StateSpace<ns,ni,no> >:std::true_type{};        
    }
}
