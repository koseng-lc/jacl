/**
*   @author : koseng (Lintang)
*   @brief : Singleton design pattern
*/

#pragma once

#include <boost/thread.hpp>

namespace jacl::pattern{

template <typename T>
class Singleton{
public:
    static auto getInstance() -> T&{
        static T instance;     
        return instance;
    }
    Singleton(const Singleton&) = delete;
    Singleton& operator=(const Singleton&) = delete;
protected:
    Singleton() = default;
    virtual ~Singleton(){}
};

} // namespace::pattern