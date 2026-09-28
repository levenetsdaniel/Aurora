#ifndef AURORA_PARALLELFOR_HPP
#define AURORA_PARALLELFOR_HPP

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <stdexcept>
#include <thread>
#include <vector>

template <typename Action>
void parallelFor(std::size_t count, std::size_t threads, Action action) {
    if (threads == 0) {
        throw std::invalid_argument("parallelFor: thread count must be positive");
    }
    if (count == 0) return;
    threads = std::min(count, threads);
    std::atomic<std::size_t> next{0};
    std::vector<std::exception_ptr> failures(count);
    std::vector<std::jthread> workers;
    workers.reserve(threads);
    for (std::size_t worker = 0; worker < threads; ++worker) {
        workers.emplace_back([&] {
            while (true) {
                const std::size_t index = next.fetch_add(1, std::memory_order_relaxed);
                if (index >= count) break;
                try {
                    action(index);
                } catch (...) {
                    failures[index] = std::current_exception();
                }
            }
        });
    }
    for (std::jthread& worker : workers) worker.join();
    for (const auto& failure : failures) {
        if (failure) std::rethrow_exception(failure);
    }
}

#endif // AURORA_PARALLELFOR_HPP
