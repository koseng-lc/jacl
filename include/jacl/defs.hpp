/**
*   @author : koseng (Lintang)
*   @brief : JACL Macros
*/

#pragma once

#define JSS_VAR 0
#define JSS_TYPE 0

#define JACL_TF(x) [](JACL::TransferMatrix::S s)->JACL::TransferMatrix::S{return x;}
#define JACL_CONST_TF(x) [](JACL::TransferMatrix::S s)->JACL::TransferMatrix::S{(void)s; return x;}

#define JE(p1,p2) [](decltype(p1) _p)->double{return p2;}
#define JC(p1,p2) [](decltype(p1) _p)->double{(void)_p; return p2;}