#include "timer.h"
namespace TIMER {
wcTimer::wcTimer() {
    active = false;
}
void wcTimer::Start() {
    tStart = std::chrono::steady_clock::now();
    if (!active) {
        Reset();
        active = true;
    }
}

void wcTimer::Stop() {
    if (!active) {
        return;
    }
    ticksFromLastStart = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now()
                                                                               - tStart).count();
    ticksFromFirstStart += ticksFromLastStart;
}

void wcTimer::Reset() {
    ticksFromFirstStart = 0;
    ticksFromLastStart = 0;
    active = false;
}

} //TIMER
