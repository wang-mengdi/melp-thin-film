//////////////////////////////////////////////////////////////////////////
// Timer
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Timer_h__
#define __Timer_h__
#include <chrono>
#include <ctime>
#include <string>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include "Constants.h"

//everything happens in seconds here
//assume your main loop is: A-B-C. 
//then you'll have TimerRecord("A"), TimerRecord("B"), TimerRecord("C").
class TimerRecord {
public:
	double total;//total time
	double cur;//time of the current loop
	int num;//number of loops
	std::string name;
	TimerRecord(std::string _name) :total(0), cur(0), num(0),name(_name) {}
	void Add(double cur_time);//add a loop
	std::pair<double, double> Profile(void);//first:current time. second: avg time
};

class Timer {
public:
	Timer() {
		name_index.clear();
		records.clear();
		total_start = std::chrono::system_clock::now();
		Reset();
	}

	void Reset();
	double Elapse(const real unit = PhysicalUnits::ms);
	double Total_Elapsed(const real unit = PhysicalUnits::ms);
	double Elapse_And_Reset();
	void Elapse_And_Output(const std::string& message);
	void Elapse_And_Output_And_Reset(const std::string& message);

	//Now suppose your main loop is A->B->C. You call as follow:
	//Begin_Loop() at the beginning
	//Record("A");Record("B");Record("C"); after corresponding parts are done
	//Output_Profile(std::cout); at the end of the loop (an ofstream is also OK)
	void Begin_Loop(void);
	void Record(const std::string& name);
	void Output_Profile(std::ostream& out = std::cout);

protected:
	std::chrono::time_point<std::chrono::system_clock> total_start;
	std::chrono::time_point<std::chrono::system_clock> start;
	std::chrono::time_point<std::chrono::system_clock> loop_start;
	std::map<std::string, int> name_index;
	std::vector<TimerRecord> records;
};
#endif
