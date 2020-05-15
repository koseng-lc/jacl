/**
*   @author : koseng (Lintang)
*   @brief : jacl transfer matrix module
*/

#pragma once

#include <functional>
#include <complex>
#include <vector>

namespace jacl{

template <typename T> // temporary template, to obligate One Definition Rule
class TransferMatrix{
public:
    using S = std::complex<double >;
    using TransferFunction = std::function<S(S) >;
    using TMElem = std::vector<TransferFunction >;

public:
    TransferMatrix(int _rows, int _cols);
    ~TransferMatrix(){}
    void setElements(std::initializer_list<TransferFunction > _elem);
    TMElem getElements();    

private:
    TMElem elem_;

    int rows_;
    int cols_;

};

template <typename T>
TransferMatrix<T>::TransferMatrix(int _rows, int _cols)
    : rows_(_rows)
    , cols_(_cols){

}

template <typename T>
void TransferMatrix<T>::setElements(std::initializer_list<TransferFunction> _elem){
    elem_.insert(elem_.begin(), _elem.begin(), _elem.end());
}

template <typename T>
typename TransferMatrix<T>::TMElem TransferMatrix<T>::getElements(){
    return elem_;
}

}
