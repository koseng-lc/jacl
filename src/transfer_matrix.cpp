#include "transfer_matrix.h"

using namespace JACL;

TransferMatrix::TransferMatrix(int _rows, int _cols)
    : rows_(_rows)
    , cols_(_cols){

}

void TransferMatrix::setElements(std::initializer_list<TransferFunction> _elem){
    elem_.insert(elem_.begin(), _elem.begin(), _elem.end());
}

TransferMatrix::TMElem TransferMatrix::getElements(){
    return elem_;
}
