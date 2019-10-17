/**
*   @author : koseng (Lintang)
*   @brief : JACL Transfer Matrix module
*/

#pragma once

#include <functional>
#include <complex>
#include <vector>

namespace JACL{

class TransferMatrix{

public:
    typedef std::complex<double > S;
    typedef std::function<S(S) > TransferFunction;
    typedef std::vector<TransferFunction > TMElem;

public:
    void setElements(std::initializer_list<TransferFunction > _elem);
    TMElem getElements();
    TransferMatrix(int _rows, int _cols);

private:
    TMElem elem_;

    int rows_;
    int cols_;

};

}
