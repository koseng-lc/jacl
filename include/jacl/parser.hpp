#pragma once

#include <fstream>

#include <boost/lexical_cast.hpp>

#include "state_space.hpp"
#include "traits.hpp"

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

    static auto deserialize(const std::string& _str, std::size_t _rows, std::size_t _cols) -> arma::mat{
        arma::mat res(_rows, _cols, arma::fill::zeros);
        std::string elem;
        int j(0);
        for(int i(0); i < _str.size(); i++){
            if(_str.at(i) == ' '){                
                res(j) = std::atof(elem.c_str());//boost::lexical_cast<double>(elem);
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
        ss_file << "state_space\n";
        ss_file << std::to_string(_StateSpace::n_states) << "\n";
        ss_file << std::to_string(_StateSpace::n_inputs) << "\n";
        ss_file << std::to_string(_StateSpace::n_outputs) << "\n";
        ss_file << detail::serialize(_ss.A()) << "\n";
        ss_file << detail::serialize(_ss.B()) << "\n";
        ss_file << detail::serialize(_ss.C()) << "\n";
        ss_file << detail::serialize(_ss.D());
        ss_file.close();
    }
    return 0;
}

template <typename Type>
static auto saveGain(const arma::Mat<Type>& _G, std::string _file_path)
    -> int{
    std::ofstream gain_file(_file_path);
    if(gain_file.is_open()){
        gain_file << "gain\n";
        gain_file << std::to_string(_G.n_rows) << "\n";
        gain_file << std::to_string(_G.n_cols) << "\n";
        gain_file << detail::serialize(_G);
        gain_file.close();
    }
    return 0;
}

template <typename _StateSpace>
static auto readStateSpace(_StateSpace* _ss, std::string _file_path)
    -> typename std::enable_if<traits::is_state_space<_StateSpace>::value, int>::type{
    std::ifstream ss_file(_file_path);
    if(ss_file.is_open()){
        std::string type;
        std::getline(ss_file, type);
        if(type != "state_space")
            throw "Invalid file";
        std::string ns, ni, no;
        std::getline(ss_file, ns);
        std::getline(ss_file, ni);
        std::getline(ss_file, no);
        if(boost::lexical_cast<std::size_t>(ns) != _StateSpace::n_states
                && boost::lexical_cast<std::size_t>(ni) != _StateSpace::n_inputs
                && boost::lexical_cast<std::size_t>(no) != _StateSpace::n_outputs)
            throw "Invalid state space size";
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

template <typename Type>
static auto readGain(arma::Mat<Type>* _G, std::string _file_path)
    -> int{
    std::ifstream gain_file(_file_path);
    if(gain_file.is_open()){
        std::string type;
        std::getline(gain_file, type);
        if(type != "gain")
            throw "Invalid file";
        std::string r,c;
        std::getline(gain_file, r);
        std::getline(gain_file, c);
        std::size_t nr( boost::lexical_cast<std::size_t>(r) )
                ,nc( boost::lexical_cast<std::size_t>(c) );
        std::string G;
        std::getline(gain_file, G);
        *_G = detail::deserialize(G, nr, nc);
        gain_file.close();
    }
    return 0;
}

}

}
