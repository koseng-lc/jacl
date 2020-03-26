#pragma once

#include <fstream>

#include <boost/lexical_cast.hpp>

#include "state_space.h"
#include "traits.h"

namespace jacl{

namespace parser{

namespace detail{
    static auto serialize(const arma::mat& _m) -> std::string{
        std::string res;
        for(int i(0); i < _m.n_rows * _m.n_cols; i++){
            res += std::to_string(_m(i)) + " ";
        }
        return res;
    }

    static auto deserialize(const std::string& _str, int _rows, int _cols) -> arma::mat{
        arma::mat res(_rows, _cols, arma::fill::zeros);
        std::string elem;
        int j(0);
        for(int i(0); i < _str.size(); i++){
            if(_str.at(i) == ' '){
                res(j) = boost::lexical_cast<double>(elem);
                elem.clear();
                j++;
                continue;
            }
            elem += _str.at(i);
        }
        return res;
    }
}

template <typename _StateSpace>
static auto saveStateSpace(const _StateSpace& _ss, std::string _file_path)
    -> typename std::enable_if<traits::is_state_space<_StateSpace>::value, int>::type{
    std::ofstream ss_file(_file_path);
    if(ss_file.is_open()){
        ss_file << detail::serialize(_ss.A()) << "\n";
        ss_file << detail::serialize(_ss.B()) << "\n";
        ss_file << detail::serialize(_ss.C()) << "\n";
        ss_file << detail::serialize(_ss.D());
        ss_file.close();
    }
    return 0;
}

template <typename _StateSpace>
static auto readStateSpace(_StateSpace* _ss, std::string _file_path)
    -> typename std::enable_if<traits::is_state_space<_StateSpace>::value, int>::type{
    std::ifstream ss_file(_file_path);
    if(ss_file.is_open()){
        std::string A,B,C,D;
        std::getline(ss_file, A);
        _ss->setA(detail::deserialize(A, _StateSpace::n_states, _StateSpace::n_states));
        std::getline(ss_file, B);
        _ss->setB(detail::deserialize(B, _StateSpace::n_states, _StateSpace::n_inputs));
        std::getline(ss_file, C);
        _ss->setC(detail::deserialize(C, _StateSpace::n_outputs, _StateSpace::n_states));
        std::getline(ss_file, D);
        _ss->setD(detail::deserialize(D, _StateSpace::n_outputs, _StateSpace::n_inputs));
        ss_file.close();
    }
    return 0;
}

}

}
