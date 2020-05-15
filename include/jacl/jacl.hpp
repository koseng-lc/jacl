/**
*   @author : koseng (Lintang)
*   @brief : jacl
*   @details : jacl stands for Just Another Control Library
*/

#pragma once

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <jacl/physical_parameter.hpp>
#include <jacl/defs.hpp>
#include <jacl/transfer_matrix.hpp>
#include <jacl/state_space/linear.hpp>
#include <jacl/state_space/nonlinear.hpp>
#include <jacl/system/base_system.hpp>
#include <jacl/system/continuous.hpp>
#include <jacl/system/discrete.hpp>
#include <jacl/system/observer.hpp>
#include <jacl/system/continuous_observer.hpp>
#include <jacl/system/discrete_observer.hpp>
#include <jacl/py_stuff.hpp>
#include <jacl/diagnosis/ifd.hpp>
#include <jacl/diagnosis/sifd.hpp>
#include <jacl/pole_placement.hpp>
#include <jacl/analysis/transient.hpp>
#include <jacl/synthesis/h_inf.hpp>
#include <jacl/synthesis/dh_inf.hpp>
#include <jacl/plotter.hpp>
#include <jacl/parser.hpp>
#include <jacl/traits.hpp>
#include <jacl/random.hpp>