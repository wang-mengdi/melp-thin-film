//////////////////////////////////////////////////////////////////////////
// Timer cpp
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Timer.h"

void TimerRecord::Add(double cur_time)
{
	total += cur_time;
	cur = cur_time;
	num++;
}

std::pair<double, double> TimerRecord::Profile(void)
{
	return std::pair<double, double>(cur, total / num);
}


void Timer::Reset()
{
	start=std::chrono::system_clock::now();
}

double Timer::Elapse(const real unit)
{
	std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
	std::chrono::duration<double, std::ratio<1> > elapse = end - start;
	return elapse.count() / unit;
}

double Timer::Total_Elapsed(const real unit)
{
	std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
	std::chrono::duration<double, std::ratio<1> > elapse = end - total_start;//in seconds
	return elapse.count() / unit;
}

double Timer::Elapse_And_Reset()
{
	double elapse=Elapse();
	Reset();
	return elapse;
}

void Timer::Elapse_And_Output(const std::string& message)
{
	double elapse=Elapse();
	Info("{} cost {:.2f}ms", message, elapse);
}

void Timer::Elapse_And_Output_And_Reset(const std::string& message)
{
	Elapse_And_Output(message);
	Reset();
}


void Timer::Begin_Loop(void)
{
	loop_start = std::chrono::system_clock::now();
	Reset();
}

void Timer::Record(const std::string& name)
{
	int k = -1;
	auto iter = name_index.find(name);
	if (iter == name_index.end()) {
		k = records.size();
		name_index[name] = k;
		records.push_back(TimerRecord(name));//push a zero record
	}
	else {
		k = iter->second;
	}
	records[k].Add(Elapse_And_Reset() / 1000.0);//Elapse_And_Reset() returns in ms
}

void Timer::Output_Profile(std::ostream& out)
{
	std::string str = "#     Time record cur(avg) in seconds: ";
	for (int i = 0; i < records.size(); i++) {
		std::pair<double, double> rec = records[i].Profile();
		str += fmt::format("{}: {:.3f}({:.3f}), ", records[i].name, rec.first, rec.second);
	}
	auto loop_end = std::chrono::system_clock::now();
	std::chrono::duration<double, std::ratio<1> > loop_elapse = loop_end - loop_start;//in seconds
	double loop_time_sec = loop_elapse.count() / PhysicalUnits::s;
	str += fmt::format("total: {:.3f}", loop_time_sec);
	out << str << "\n";
}