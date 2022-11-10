#pragma once
#include <math.hpp>
#include <integrand.hpp>
#include <integral.hpp>
#include <spheroid.hpp>
#include <sqlitedb.hpp>
#include <csvdb.hpp>
#include <boost/multiprecision/float128.hpp>

namespace database {
	const std::string sqlite_name = "../data/database/integrals.db";
}

template <typename T>
T get_An_spheroidal_integral(float radius, float height, unsigned n, bool save = true) {

	const std::string shape_name = "prolate_spheroid";
	// Check if integral is already calculated.
	{
		Data<T> res = SQLite::get_from_database<T>(database::sqlite_name, shape_name, radius, height, __func__, n, 0);
		if (res.found == true) return res.result;
	}

	// Check if integral is in CSV database.
    {
        std::stringstream ss;
		ss << "../data/database/spheroidal/integrals_r" << radius << "h" << height << ".csv";
		Data<T> res = CSVDB::get_from_database<T>(ss.str(), shape_name, radius, height, __func__, n, 0);
        // if (res.found == true) SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, n, 0, res.result);
        if (res.found == true) return res.result;
    }



	// Check if it was calculated in the deprecated textfiles.
	// Please delete this block as soon as possible.
	/*{
		std::stringstream ss;
		ss << "../data/database/relic/an_spheroidal_h" << height;
		std::string name = ss.str();
		
		using SmallFloat = boost::multiprecision::float128;
		std::ifstream file(name, std::ios::in);
		std::string content = "";

		// Load lines.
		unsigned line = 0;
		while (getline(file, content) and line != n/2) ++line;
		file.close();

		// The file contains the number. Load it, save, and return.
		if (line == n/2 and !content.empty()) {
			SmallFloat result = SmallFloat(content);

			if (save) {
				SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, n, 0, result);
			}

			return T(content);
		}
	}*/


	// Calculate integral.
	T tradius = static_cast<T>(radius);
	T theight = static_cast<T>(height);
	T result = integrate<T>(-1.0, 1.0, An_integrand<T>(tradius, theight, n));

	// Save result on database.
	if (save) {
        std::stringstream ss;
		ss << "../data/database/spheroidal/integrals_r" << radius << "h" << height << ".csv";
		CSVDB::save_at_database(ss.str(), shape_name, radius, height, __func__, n, 0, result);
		SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, n, 0, result);
	}

	// Return.
	return result;
}

template <typename T>
T get_Gl_spheroidal_integral(float radius, float height, unsigned l, bool save = true) {
	// Check if integral is already calculated.
	const std::string shape_name = "prolate_spheroid";
	{
		Data<T> res = SQLite::get_from_database<T>(database::sqlite_name, shape_name, radius, height, __func__, l, 0);
		if (res.found == true) return res.result;
	}

	// Check if integral is in relic CSV database.
    {
        std::stringstream ss;
		ss << "../data/database/spheroidal/integrals_r" << radius << "h" << height << ".csv";
		Data<T> res = CSVDB::get_from_database<T>(ss.str(), shape_name, radius, height, __func__, l, 0);
        // if (res.found == true) SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, l, 0, res.result);
        if (res.found == true) return res.result;
    }

	// Calculate integral.
	T tradius = static_cast<T>(radius);
	T theight = static_cast<T>(height);
	T result = integrate<T>(-1.0, 1.0, Gl_spheroidal_integrand<T>(tradius, theight, l));

	// Save result on database.
	if (save) {
        std::stringstream ss;
		ss << "../data/database/spheroidal/integrals_r" << radius << "h" << height << ".csv";
		CSVDB::save_at_database(ss.str(), shape_name, radius, height, __func__, l, 0, result);
		SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, l, 0, result);
	}

	// Return.
	return result;
}


template <typename T>
T get_Iij_spheroidal_integral(float radius, float height, unsigned i, unsigned j, bool save = true) {
	// Check if integral is already calculated.
	const std::string shape_name = "prolate_spheroid";
	{
		Data<T> res = SQLite::get_from_database<T>(database::sqlite_name, shape_name, radius, height, __func__, i, j);
		if (res.found == true) return res.result;
	}

	// Check if integral is in relic CSV database.
    {
        std::stringstream ss;
		ss << "../data/database/spheroidal/integrals_r" << radius << "h" << height << ".csv";
		Data<T> res = CSVDB::get_from_database<T>(ss.str(), shape_name, radius, height, __func__, i, j);
        // if (res.found == true) SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, i, j, res.result);
        if (res.found == true) return res.result;
    }

	// Calculate integral.
	T tradius = static_cast<T>(radius);
	T theight = static_cast<T>(height);
	T result = integrate<T>(-1.0, 1.0, Iij_spheroidal_integrand<T>(tradius, theight, i, j));

	// Save result on database.
	if (save) {
        std::stringstream ss;
		ss << "../data/database/spheroidal/integrals_r" << radius << "h" << height << ".csv";
		CSVDB::save_at_database(ss.str(), shape_name, radius, height, __func__, i, j, result);
		SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, i, j, result);
	}

	// Return.
	return result;
}


template <typename T>
T get_Gl_spheroidal_summation(float radius, float height, unsigned l, bool save = true) {
	// Check if integral is already calculated.
	const std::string shape_name = "prolate_spheroid";
	{
		Data<T> res = SQLite::get_from_database<T>(database::sqlite_name, shape_name, radius, height, __func__, l, 0);
		if (res.found == true) return res.result;
	}

	// Check if integral is in relic CSV database.
    {
        std::stringstream ss;
		ss << "../data/database/spheroidal/integrals_r" << radius << "h" << height << ".csv";
		Data<T> res = CSVDB::get_from_database<T>(ss.str(), shape_name, radius, height, __func__, l, 0);
        // if (res.found == true) SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, l, 0, res.result);
        if (res.found == true) return res.result;
    }

	// Get the An values.
	using SmallFloat = boost::multiprecision::float128;
	unsigned num = static_cast<unsigned>(std::ceil(height / radius));
	unsigned total = l + 1;
	unsigned an_total = 4*total*num*num + 30;

	std::vector<T> an;
	for (unsigned n = 0; n < an_total; ++n) {
		an.push_back(
			static_cast<T>(
				get_An_spheroidal_integral<SmallFloat>(radius, height, n)
			)
		);
	}

	// Calculate summation.
	T result = Gl_spheroidal_summation_exact<T>(radius, height, l, an);

	// Save result on database.
	if (save) {
        std::stringstream ss;
		ss << "../data/database/spheroidal/integrals_r" << radius << "h" << height << ".csv";
		CSVDB::save_at_database(ss.str(), shape_name, radius, height, __func__, l, 0, result);
		SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, l, 0, result);
	}

	// Return.
	return result;
}


template <typename T>
T get_Iij_spheroidal_summation(float radius, float height, unsigned i, unsigned j, bool save = true) {
	// Check if integral is already calculated.
	const std::string shape_name = "prolate_spheroid";
	{
		Data<T> res = SQLite::get_from_database<T>(database::sqlite_name, shape_name, radius, height, __func__, i, j);
		if (res.found == true) return res.result;
	}

	// Check if integral is in relic CSV database.
    {
        std::stringstream ss;
		ss << "../data/database/spheroidal/integrals_r" << radius << "h" << height << ".csv";
		Data<T> res = CSVDB::get_from_database<T>(ss.str(), shape_name, radius, height, __func__, i, j);
        // if (res.found == true) SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, i, j, res.result);
        if (res.found == true) return res.result;
    }

	// Get the An values.
	using SmallFloat = boost::multiprecision::float128;
	unsigned num = static_cast<unsigned>(std::ceil(height / radius));
	unsigned total = i + j + 1;
	unsigned an_total = 4*total*num*num + 30;

	std::vector<T> an;
	for (unsigned n = 0; n < an_total; ++n) {
		an.push_back(
			static_cast<T>(
				get_An_spheroidal_integral<SmallFloat>(radius, height, n)
			)
		);
	}

	// Calculate integral.
	T result = Iij_spheroidal_summation_exact<T>(radius, height, i, j, an);

	// Save result on database.
	if (save) {
        std::stringstream ss;
		ss << "../data/database/spheroidal/integrals_r" << radius << "h" << height << ".csv";
		CSVDB::save_at_database(ss.str(), shape_name, radius, height, __func__, i, j, result);
		SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, i, j, result);
	}

	// Return.
	return result;
}

template <typename T>
T get_Gl_suspended_hemispherical_integral(float radius, float height, unsigned l, bool save = true) {
	// Check if integral is already calculated.
	const std::string shape_name = "suspended hemisphere";
	{
		Data<T> res = SQLite::get_from_database<T>(database::sqlite_name, shape_name, radius, height, __func__, l, 0);
		if (res.found == true) return res.result;
	}

	// Check if integral is in relic CSV database.
    {
        std::stringstream ss;
		ss << "../data/database/suspended/integrals_r" << radius << "h" << height << ".csv";
		Data<T> res = CSVDB::get_from_database<T>(ss.str(), shape_name, radius, height, __func__, l, 0);
        // if (res.found == true) SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, l, 0, res.result);
        if (res.found == true) return res.result;
    }

	// Calculate integral.
	T tradius = static_cast<T>(radius);
	T theight = static_cast<T>(height);
	T tlength = theight - tradius;
	T x_limit = tlength / math::sqrt(tlength * tlength + tradius * tradius);
	T upper_x = integrate<T>(x_limit, 1.0, Gl_suspended_hemispherical_integrand<T>(tradius, tlength, l));
	T lower_x = integrate<T>(-1.0, -x_limit, Gl_suspended_hemispherical_integrand<T>(tradius, -tlength, l));
	T result = upper_x + lower_x;

	// Save result on database.
	if (save) {
        std::stringstream ss;
		ss << "../data/database/suspended/integrals_r" << radius << "h" << height << ".csv";
		CSVDB::save_at_database(ss.str(), shape_name, radius, height, __func__, l, 0, result);
		SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, l, 0, result);
	}

	// Return.
	return result;
}


template <typename T>
T get_Iij_suspended_hemispherical_integral(float radius, float height, unsigned i, unsigned j, bool save = true) {
	// Check if integral is already calculated.
	const std::string shape_name = "suspended hemisphere";
	{
		Data<T> res = SQLite::get_from_database<T>(database::sqlite_name, shape_name, radius, height, __func__, i, j);
		if (res.found == true) return res.result;
	}


	// Check if integral is in relic CSV database.
    {
        std::stringstream ss;
		ss << "../data/database/suspended/integrals_r" << radius << "h" << height << ".csv";
		Data<T> res = CSVDB::get_from_database<T>(ss.str(), shape_name, radius, height, __func__, i, j);
        // if (res.found == true) SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, l, 0, res.result);
        if (res.found == true) return res.result;
    }

	// Calculate integral.
	T tradius = static_cast<T>(radius);
	T theight = static_cast<T>(height);
	T tlength = theight - tradius;
	T x_limit = tlength / math::sqrt(tlength * tlength + tradius * tradius);
	T upper_x = integrate<T>(x_limit, 1.0, Iij_suspended_hemispherical_integrand<T>(tradius, tlength, i, j));
	T lower_x = integrate<T>(-1.0, -x_limit, Iij_suspended_hemispherical_integrand<T>(tradius, -tlength, i, j));
	T result = upper_x + lower_x;

	// Save result on database.
	if (save) {
        std::stringstream ss;
		ss << "../data/database/suspended/integrals_r" << radius << "h" << height << ".csv";
		CSVDB::save_at_database(ss.str(), shape_name, radius, height, __func__, i, j, result);
		SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, i, j, result);
	}

	// Return.
	return result;
}


template <typename T>
T get_Gl_cylindrical_integral(float radius, float height, unsigned l, bool save = true) {
	// Check if integral is already calculated.
	const std::string shape_name = "cylinder";
	{
		Data<T> res = SQLite::get_from_database<T>(database::sqlite_name, shape_name, radius, height, __func__, l, 0);
		if (res.found == true) return res.result;
	}

	// Check if integral is in relic CSV database.
    {
        std::stringstream ss;
		ss << "../data/database/cylinder/integrals_r" << radius << "h" << height << ".csv";
		Data<T> res = CSVDB::get_from_database<T>(ss.str(), shape_name, radius, height, __func__, l, 0);
        // if (res.found == true) SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, l, 0, res.result);
        if (res.found == true) return res.result;
    }

	// Calculate integral.
	T tradius = static_cast<T>(radius);
	T theight = static_cast<T>(height);
	T tlength = theight - tradius;
	T x_limit = tlength / math::sqrt(tlength * tlength + tradius * tradius);
	T result = integrate<T>(-x_limit, x_limit, Gl_cylindrical_integrand<T>(tradius, tlength, l));

	// Save result on database.
	if (save) {
        std::stringstream ss;
		ss << "../data/database/cylinder/integrals_r" << radius << "h" << height << ".csv";
		CSVDB::save_at_database(ss.str(), shape_name, radius, height, __func__, l, 0, result);
		SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, l, 0, result);
	}

	// Return.
	return result;
}

template <typename T>
T get_Iij_cylindrical_integral(float radius, float height, unsigned i, unsigned j, bool save = true) {
	// Check if integral is already calculated.
	const std::string shape_name = "cylinder";
	{
		Data<T> res = SQLite::get_from_database<T>(database::sqlite_name, shape_name, radius, height, __func__, i, j);
		if (res.found == true) return res.result;
	}


	// Check if integral is in relic CSV database.
    {
        std::stringstream ss;
		ss << "../data/database/cylinder/integrals_r" << radius << "h" << height << ".csv";
		Data<T> res = CSVDB::get_from_database<T>(ss.str(), shape_name, radius, height, __func__, i, j);
        // if (res.found == true) SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, l, 0, res.result);
        if (res.found == true) return res.result;
    }


	// Calculate integral.
	T tradius = static_cast<T>(radius);
	T theight = static_cast<T>(height);
	T tlength = theight - tradius;
	T x_limit = tlength / math::sqrt(tlength * tlength + tradius * tradius);
	T result = integrate<T>(-x_limit, x_limit, Iij_cylindrical_integrand<T>(tradius, tlength, i, j));

	// Save result on database.
	if (save) {
        std::stringstream ss;
		ss << "../data/database/cylinder/integrals_r" << radius << "h" << height << ".csv";
		CSVDB::save_at_database(ss.str(), shape_name, radius, height, __func__, i, j, result);
		SQLite::save_at_database(database::sqlite_name, shape_name, radius, height, __func__, i, j, result);
	}

	// Return.
	return result;
}

template <typename T>
T get_Gl_hcp_integral(float radius, float height, unsigned l, bool save = true) {
	T cylinder = get_Gl_cylindrical_integral<T>(radius, height, l, save);
	T susp = get_Gl_suspended_hemispherical_integral<T>(radius, height, l, save);
	return cylinder + susp;
}

template <typename T>
T get_Iij_hcp_integral(float radius, float height, unsigned i, unsigned j, bool save = true) {
	T cylinder = get_Iij_cylindrical_integral<T>(radius, height, i, j, save);
	T susp = get_Iij_suspended_hemispherical_integral<T>(radius, height, i, j, save);
	return cylinder + susp;
}

