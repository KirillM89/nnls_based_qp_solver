#ifndef QP_NNLS_DECORATORS_H
#define QP_NNLS_DECORATORS_H
#include <memory>
#include "types.h"
#include "callback.h"
namespace QP_NNLS {
    class Core;
    class QPNNLS {
    public:
        QPNNLS();
         ~QPNNLS();
        QPNNLS(const QPNNLS& other) = delete;
        QPNNLS(QPNNLS&& other) = delete;
        QPNNLS& operator=(const QPNNLS& other) = delete;
        QPNNLS& operator=(QPNNLS&& other) = delete;
        void SetCallback(Callback* callback = nullptr);
        void Solve(const Input& problem, const Configuration& config = Configuration());
        const Output& GetOutput();
    protected:
        bool isInitialized;
        Core* core;
        Output output;
    };
}

#endif // DECORATORS_H
