#pragma once

#include <vector>

namespace jacl{ namespace pattern{

class Observer;

class Subject{
public:
    Subject(){}
    ~Subject(){}
    auto attach(Observer* _c) -> void{        
        observer_.push_back(_c);
    }
    inline auto notify() -> void;

private:
    std::vector<Observer*> observer_;
};

class Observer{
public:
    Observer(Subject* _s)
        : s_(_s){
        if(s_)
            s_->attach(this);
    }
    ~Observer(){}
    virtual auto update() -> void = 0;    
protected:
    Subject* s_;
};

auto Subject::notify() -> void{
    if(observer_.size()){
        std::cout << "Observer size : " << observer_.size() << std::endl;
        for(const auto& o:observer_)
            o->update();
    }        
}

} } // namespace jacl::pattern