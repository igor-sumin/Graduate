#pragma once

#include <iostream>
#include <chrono>
#include <string>

using namespace std;
using namespace std::chrono;

// Позволяет объявить уникальный идентификатор в пределах заданного cpp-файла
#define UNIQ_ID(line) UNIQ_ID_IMPL(line)
#define UNIQ_ID_IMPL(line) loc_var_##line
#define LOG_DURATION(message) LogDuration UNIQ_ID(__LINE__){message};

class LogDuration {
private:
	// Момент начала замера времени выполнения программы
	steady_clock::time_point start;
	string message;

public:
	// Добавление explicit - защита от неявного преобразования string в LogDuration
	explicit LogDuration(const string& msg = "") : start(steady_clock::now()), message(msg + ": ")
	{
	}

	~LogDuration() {
		auto finish = steady_clock::now();
		cerr << message << duration_cast<milliseconds>(finish - start).count() << " ms ("
			<< duration_cast<seconds>(finish - start).count() << " s)" << endl;
	}

};