/*************************************************************************
*  2020      Janek Kozicki                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

/*
 * This simple class measures time since last .restart() call which is also executed in constructor.
 *
 * The function bool check(howOften) returns if the howOften period has passed since last restart() and calls restart() if true.
 *
 * The howOften argument uses C++14 / C++20 time units defined in namespace std::chrono_literals;
 * So the arguments can be any of the units defined in https://en.cppreference.com/w/cpp/header/chrono#Literals
 *
 * Typical usage is to use namespace std::chrono_literals; Then write directly: 10s or 200ms.
 * For example usage see file py/_log.cpp function testTimedLevels();
 *
 * See also:
 * - https://en.cppreference.com/w/cpp/header/chrono#Literals , note: by default LOG_TIMED_* in lib/base/LoggingUtils.hpp only accepts seconds and milliseconds (we can change that if necessary)
 * - https://en.cppreference.com/w/cpp/thread/sleep_for
 * - https://en.cppreference.com/w/cpp/chrono/duration
 */

#pragma once
#include <chrono>
#include <type_traits>

namespace yade {
struct Timer {
	std::chrono::time_point<std::chrono::high_resolution_clock> start;
	bool                                                        first { true };

	// the Timer starts in constructor
	// firstAlwaysTrue == true means that the first call of check() always returns true; even if the howOften is much shorter than the time since construction of Timer.
	// this flag ↓ is useful when we want to see the first log message
	Timer(bool firstAlwaysTrue)
	{
		first = firstAlwaysTrue;
		restart();
	};

	// restart stopwatch to current time.
	inline void restart(std::chrono::time_point<std::chrono::high_resolution_clock> now = std::chrono::high_resolution_clock::now()) { start = now; }

	// check if howOften has passed. If yes, then restart() and return true. Otherwise return false. Example call:
	//    using namespace std::chrono_literals;
	//    if(t.check(10ms)) // 10 milliseconds.
	// The 10ms is a compile-time rational constant from namespace std::chrono_literals, so the compiler optimizes it into a single integer comparison.
	template <typename A, typename B> inline bool check(const std::chrono::duration<A, B>& howOften)
	{
		auto now = std::chrono::high_resolution_clock::now();
		if (first or ((now - start) > howOften)) {
			first = false;
			restart(now);
			return true;
		} else
			return false;
	}
};
}
