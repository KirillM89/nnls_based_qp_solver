#ifndef QP_NNLS_CALLBACK_H
#define QP_NNLS_CALLBACK_H
#include <set>
#include <string>
#include "types.h"
namespace QP_NNLS {
    struct InitializationData {
        unsg_t nVariables;
        unsg_t nConstraints;
        unsg_t nEqConstraints;
        fp_t scaleDB;
        std::string tChol;
        std::string tInv;
        std::string tM;
        std::string initStatus;
        std::vector<fp_t>* s;
        std::vector<fp_t>* b;
        std::vector<fp_t>* c;
        std::vector<fp_t>* Chol;
        std::vector<fp_t>* CholInv;
        matrix_t* M;
    };
    struct IterationData {
       bool singular;
       unsg_t newIndex;
       unsg_t iteration;
       fp_t gamma;
       fp_t dualTol;
       fp_t rsNorm;
       std::vector<fp_t>* dual;
       std::vector<fp_t>* primal;
       std::vector<fp_t>* zp;
       std::set<unsg_t>* activeSet;
    };
    struct FinalData {
        unsigned char primalStatus;
        unsigned char dualStatus;
        unsg_t nIterations;
        fp_t cost;
        std::vector<fp_t>* x;
        std::vector<fp_t>* lambda;
        std::vector<fp_t>* lambdaUp;
        std::vector<fp_t>* lambdaLw;
    };

    class Callback {
    public:
        Callback() = default;
        virtual ~Callback() = default;
        virtual void Init() {
            return;
        }
        virtual void SetLogLevel(unsg_t logLevel) {
            return;
        }
        virtual void ProcessData(int stage) {
            return;
        };
        InitializationData initData;
        IterationData iterData;
        FinalData finalData;
    };

    class Logger;
    class Callback1 : public Callback {
    public:
        Callback1(const std::string& filePath);
        virtual ~Callback1() override;
        void ProcessData(int stage) override;
        void SetLogLevel(unsg_t logLevel) override {
            this->logLevel = logLevel;
        }
        void Init() override;
    private:
        const std::string& filePath;
        Logger* logger;
        unsg_t logLevel = 0u;
    };
}

#endif // CALLBACK_H
