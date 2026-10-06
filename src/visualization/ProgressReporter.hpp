#ifndef AURORA_PROGRESSREPORTER_HPP
#define AURORA_PROGRESSREPORTER_HPP

#include <cstddef>
#include <mutex>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>

class ProgressReporter {
private:
    std::string stage;
    std::size_t total;
    std::size_t completed{};
    std::size_t lastReportedPercent{};
    std::ostream& output;
    std::mutex mutex;

public:
    ProgressReporter(std::string stage, std::size_t total, std::ostream& output)
        : stage(std::move(stage)), total(total), output(output) {
        if (total == 0) {
            throw std::invalid_argument("ProgressReporter: total must be positive");
        }
        output << this->stage << ": 0% (0/" << total << ")\n" << std::flush;
    }

    void advance() {
        std::lock_guard lock(mutex);
        if (completed >= total) {
            throw std::logic_error("ProgressReporter: completed exceeds total");
        }
        ++completed;
        const std::size_t percent = completed * 100 / total;
        if (percent == 100 || percent / 5 > lastReportedPercent / 5) {
            output << stage << ": " << percent << "% (" << completed
                   << "/" << total << ")\n" << std::flush;
            lastReportedPercent = percent;
        }
    }
};

#endif // AURORA_PROGRESSREPORTER_HPP
