#include <iostream>

#include <boost/python.hpp>

void greet(){
    std::cout << "TEST !!!" << std::endl;
}

BOOST_PYTHON_MODULE(hello_ext){

    boost::python::def("greet", &greet);
}
