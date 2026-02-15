#include "CycleTimer.h"
#include <iostream>

template <typename Func, typename... Args>
double timeItNS(const std::string& label, const bool output, Func&& func, Args&&... args) {
    double start = CycleTimer::currentSeconds();

    func(std::forward<Args>(args)...);

    double end = CycleTimer::currentSeconds();

    double ns = (end - start) * 1e9;
    
    if (output){
        std::cout << label << " took "
                << ns << " ns ("
                << ns / 1e6 << " ms)"
                << std::endl;
    }

    return ns;
}