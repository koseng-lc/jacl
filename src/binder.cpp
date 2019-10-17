#include <boost/python.hpp>

using boost::python;

int main(int argc, char** argv){

    python::object sim_module = import("sim");

    return 0;
}
