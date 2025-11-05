#include "../include/decorators.h"
#include "core.h"
namespace QP_NNLS {

    QPNNLS::QPNNLS():
        isInitialized(false),
        core(new Core())
    {}
    QPNNLS::~QPNNLS()
    {
        delete core;
    }
    void QPNNLS::Init(const Configuration& config) {
        core->Set(config);
        isInitialized = true;
    }
    void QPNNLS::SetCallback(Callback* callback) {
        if (callback) {
            core->SetCallback(callback);
        }
    }
    const Output& QPNNLS::GetOutput() {
        return (output = core->GetOutput());
    }   
    bool QPNNLS::SetProblem(const Input& problem) {
        if (!isInitialized) {
            return false;
        }
        return core->InitProblem(problem);
    }
    void QPNNLS::Solve() {
        core->Solve();
    }
    unsigned char QPNNLS::GetInitStatus() {
        InitStageStatus coreInitStatus = core->GetInitStatus();
        return static_cast<unsigned char>(coreInitStatus);
    }


}

