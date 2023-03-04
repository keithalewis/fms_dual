// fms_timer.h - time function execution
#pragma once
#include <chrono>

namespace fms {

	template<typename Duration = std::chrono::milliseconds, class F, class... Xs>
	inline typename Duration::rep timer(F&& f, Xs&&... xs)
	{
		const auto beg = std::chrono::high_resolution_clock::now();
		std::forward<F>(f)(std::forward<Xs>(xs)...);
		const auto end = std::chrono::high_resolution_clock::now();

		return std::chrono::duration_cast<Duration>(end - beg).count();
	}

} // namespace fms