//////////////////////////////////////////////////////////////////////////
// Json adapter
// NOTE: DO NOT INCLUDE THIS IN ANY CUDA FILE
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "Common.h"
#include <nlohmann/json.hpp>
#include <fmt/color.h>
#include <fmt/format.h>

using json = nlohmann::json;

//nlohmann::json adapter for eigen vector
namespace Eigen {
	//NOTE: we may need to write some specific functions -- Mengdi
	template<class T, int d>
	void to_json(json& j, const Vector<T, d>& v) {
		j = std::vector<T>(v.data(), v.data() + d);
	}
	template<class T, int d>
	void from_json(const json& j, Vector<T, d>& v) {
		Assert(!j.is_object() && j.size() == d, "Illigal json when trying to parse a {}-dimension Vector: \n{}\n", d, j.dump(4));
		for (int i = 0; i < j.size(); i++) j.at(i).get_to(v[i]);
	}

	template<class T>
	void to_json(json& j, const Quaternion<T>& q) {
		Vector<T, 4> v = q.coeffs();
		to_json(j, v);
	}
	template<class T>
	void from_json(const json& j, Quaternion<T>& q) {
		Vector<T, 4> v;
		from_json(j, v);
		q = Quaternion<T>(v.data());
	}

	template<class T>
	void to_json(json& j, const Matrix<T, 3, 3>& mat) {
		Vector<T, 3 * 3> vec(mat.data());
		to_json(j, vec);
	}
	template<class T>
	void from_json(const json& j, Matrix<T, 3, 3>& mat) {
		Vector<T, 3* 3> vec;
		from_json(j, vec);
		mat = Matrix<T, 3, 3>(vec.data());
	}
}

namespace Json {
	template<class T>
	T Value(json& j, const std::string key, const T default_value) {
		if (j.contains(key)) {
			T value = j.at(key);
			fmt::print(fg(fmt::color::green), "#     [=] Parse key ");
			fmt::print("{}", key);
			fmt::print(fg(fmt::color::green), " from json: ");
			fmt::print("{}\n", value);
			return value;
		}
		else {
			j[key] = default_value;
			fmt::print(fg(fmt::color::yellow), "#     [+] Can't parse key ");
			fmt::print("{}", key);
			fmt::print(fg(fmt::color::yellow), " in json, set to default value: ");
			fmt::print("{}\n", default_value);
			return default_value;
		}
	}

	template<class T>
	void Set_Non_Override(json& j, const std::string key, const T value) {
		if (!j.contains(key)) {
			j[key] = value;
		}
	}
}