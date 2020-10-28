/**
*   @author : koseng (Lintang)
*   @brief : Such a python embed handler
*/

#pragma once

#include <mutex>

#include <Python.h>
#include <numpy/ndarrayobject.h>

namespace jacl::py_stuff{

//-- more safe with RAII style
//-- it is required because the thread being used exclusive from Python API
class AcquireGIL{
public:
    AcquireGIL():state_(PyGILState_Ensure()){}
    ~AcquireGIL(){PyGILState_Release(state_);}

private:
    PyGILState_STATE state_;
};

//-- this is one is temporary
class PyEmbedHandler{
public:
    PyEmbedHandler(){
        // std::call_once(init_,[&](){
            Py_Initialize();
            PyEval_InitThreads();
            state_ = PyEval_SaveThread();
            int ret = initImportArray();
        // });        
    }   
    ~PyEmbedHandler(){
        // std::call_once(end_,[&](){
            PyEval_RestoreThread(state_);
            Py_Finalize();
        // });         
    }
    auto initImportArray() -> int{
        AcquireGIL lk;
        import_array1(1);
        return 1;
    }
private:    
    PyThreadState* state_;
    // std::once_flag init_;
    // std::once_flag end_;
};

} // namespace jacl::parser