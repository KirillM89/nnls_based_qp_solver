#ifndef TIMER_H
#define TIMER_H
#include <chrono>
namespace TIMER {
using ticks_t = unsigned long long; //microseconds
enum class TIMER_TYPES {
    WALL_CLOCK = 0,
    CPU,
};

class iTimer {
public:
    virtual ~iTimer() = default;
    virtual void Start() = 0;
    virtual void Stop() = 0;
    virtual void Reset() = 0;
    virtual ticks_t TimeFromFirstStart() = 0;
    virtual ticks_t TimeFromLastStart() = 0;
protected:
    iTimer() = default;
};

class wcTimer : public iTimer {
public:
    wcTimer();
    virtual ~wcTimer() = default;
    void Start() override;
    void Stop() override;
    void Reset() override;
    ticks_t TimeFromFirstStart() override { return ticksFromFirstStart;}
    ticks_t TimeFromLastStart() override { return ticksFromLastStart;}
protected:
    ticks_t ticksFromLastStart = 0;
    ticks_t ticksFromFirstStart = 0;
    std::chrono::time_point<std::chrono::steady_clock> tStart{};
    bool active = false;
};
}
#endif // TIMER_H
