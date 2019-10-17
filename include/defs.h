/**
*   @author : koseng (Lintang)
*   @brief : JACL Macros
*/

#pragma once

#define JACL_TF(x) [](JACL::TransferMatrix::S s)->JACL::TransferMatrix::S{return x;}
#define JACL_CONST_TF(x) [](JACL::TransferMatrix::S s)->JACL::TransferMatrix::S{(void)s; return x;}

#define JACL_SS(x) [](SS ss)->double{return x;}
#define JACL_CONST_SS(x) [](SS ss)->double{(void)ss; return x;}
