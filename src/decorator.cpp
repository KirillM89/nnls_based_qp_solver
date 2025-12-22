#include "decorator.h"
#include "core.h"
namespace QP_NNLS {

    QPNNLS::QPNNLS():
        core(new Core())
    {}
    QPNNLS::~QPNNLS()
    {
        delete core;
    }
    void QPNNLS::SetCallback(Callback* callback) {
        if (callback) {
            core->SetCallback(callback);
        }
    }
    const Output& QPNNLS::GetOutput() {
        return (output = core->GetOutput());
    }   
    void QPNNLS::Solve(const Input& problem, const Configuration& config = Configuration()) {
        core->Set(config);
        ErrorCode err = core->InitProblem(problem);
        core->Solve();
    }



}

