#include <integral.hpp>
#include <integrand.hpp>
#include <polynomial.hpp>
#include <calcs.hpp>
#include <leastsquares.hpp>
#include <cte.hpp>

#include <csvdb.hpp>
#include <sqlitedb.hpp>

#include <functional>
#include <iostream>
#include <vector>

#include <boost/multiprecision/float128.hpp> 
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <Eigen/Dense>

/*
POSTGRESQL:
https://www.digitalocean.com/community/tutorials/sqlite-vs-mysql-vs-postgresql-a-comparison-of-relational-database-management-systems
https://github.com/malnvenshorn/OctoPrint-FilamentManager/wiki/Setup-PostgreSQL-on-Arch-Linux
https://wiki.archlinux.org/index.php/PostgreSQL
http://pqxx.org/development/libpqxx/
https://www.tutorialspoint.com/postgresql/postgresql_c_cpp.htm

MYSQL:
https://dev.mysql.com/doc/connector-cpp/1.1/en/connector-cpp-examples-complete-example-2.html
*/

/*
BEFORE PROMOTING TO SPHEROID.HPP:
	- Identify and correct overflow errors.
	- Differentiate between calculate Gl and get Gl. (or I, or A).
		- Calc: Raw calcualtions.
		- Get: Search on DB. If doens't exist: Calculate, and save on DB, and return.
		- Perhaps Calc and Get should be in different .hpps.:w
	- Create area_element_moment function calculating *in T* (really important) (do not calc at float128 and convert). Also, do calc and get functions.
*/

const std::string sqlite_name = "../data/database/integrals.db";
const std::string shape_name = "prolate_spheroid";

template <typename T>
T Gl_spheroidal_summation_approx(float r, float h, unsigned l, bool save = true) {
	if (l % 2 == 0) return T(0);

	{
		Data<T> res = SQLite::get_from_database<T>(sqlite_name, shape_name, r, h, __func__, l, 0);
		if (res.found == true) return res.result;
	}

	unsigned p = (l-1) / 2;
	T radius = T(r);
	T height = T(h);

	T pt = static_cast<T>(p);
	T pi = T(math::cte::pi_string);
	T cte = 2 * pi / boost::multiprecision::pow(radius, 2*p-1);
	T ee2 = 1.0 - radius * radius / (height * height);
	math::Polynomial<T> p_legendre = math::Polynomial<T>::legendre(2*p+1);
	T result = 0.0;
	T bin = 1.0;

	for (unsigned n = 0; n < 10; ++n) {
		T mul = 1;
		T nt = static_cast<T>(n);

		if (n >= 1) bin *= (pt - 0.5 - nt+1.0) / nt;
		if (n % 2 == 0) mul = 1;
		else mul = -1; 

		for (unsigned k = 0; k <= p; ++k) {
			T kt = static_cast<T>(k);
			T main = p_legendre[2*k+1] * bin * boost::multiprecision::pow(ee2, n);
			T moment = 1.0 / (nt+kt+1.0+0.5) + ee2 * ee2 / 2.0 / (nt+kt+1.0+1.5);

			std::cout << "(" << pt-0.5 << ";" << nt << ") = " << bin << std::endl;
			

			result += mul * main * moment;
		}
	}

	if (save) {
		SQLite::save_at_database(sqlite_name, shape_name, r, h, __func__, l, 0, cte * result);
	}

	return cte * result;
}

template <typename T>
T Gl_spheroidal_summation_exact(float r, float h, unsigned l, const std::vector<T>& An = {}, bool save = true) {
	if (l % 2 == 0) return T(0);

	{
		Data<T> res = SQLite::get_from_database<T>(sqlite_name, shape_name, r, h, __func__, l, 0);
		if (res.found == true) return res.result;
	}

	unsigned p = (l-1) / 2;
	T radius = T(r);
	T height = T(h);

	T pt = static_cast<T>(p);
	T pi = T(math::cte::pi_string);
	T cte = 2 * pi / boost::multiprecision::pow(radius, 2*p-1);
	T ee2 = 1.0 - radius * radius / (height * height);
	math::Polynomial<T> p_legendre = math::Polynomial<T>::legendre(2*p+1);
	T result = 0.0;
	T partial_sum = 0.0;
	T bin = 1.0;
	T error = 0.0;

	unsigned n = 0;
	do {
		T mul = 1;
		T nt = static_cast<T>(n);

		if (n % 2 == 0) mul = 1;
		else mul = -1; 

		for (unsigned k = 0; k <= p; ++k) {
			T main = p_legendre[2*k+1] * bin;
			T pe = boost::multiprecision::pow(ee2, n);

			T moment;
			if (An.empty()) moment = integrate<T>(-1.0, 1.0, An_integrand<T>(radius, height, 2*(k+n+1)));
			else if (An.size() <= 2*(k+n+1)) moment = integrate<T>(-1.0, 1.0, An_integrand<T>(radius, height, 2*(k+n+1)));
			else moment = An[2*(k+n+1)];

		/*	
			std::cout << "l=" << l << ", n=" << n << ", k=" << k << std::endl;
			std::cout << "(" << pt-0.5 << ";" << nt << ") = " << bin << std::endl;
			std::cout << "mul: " << mul << std::endl;
			std::cout << "bin: " << bin << std::endl;
			std::cout << "main: " << main << std::endl;
			std::cout << "moment: " << moment << std::endl;
			std::cout << "pe: " << pe << std::endl;
			std::cout << "p result: " << mul * main * moment * pe << std::endl;
			std::cout << std::endl;
		*/

			partial_sum += mul * main * (moment * pe);
		}

		bin *= (pt - 0.5 - nt) / (nt + 1.0);
		result += partial_sum;
		error = boost::multiprecision::abs(partial_sum / result);
		partial_sum = 0.0;
		++n;
	} while (2*p > n+1  or  error > 1e-5);

	// std::cout << "n=" << n << "|" << std::flush;
	if (save) {
		// CSVDB::save_at_database(csv_name, shape_name, r, h, __func__, l, 0, cte * result);
		SQLite::save_at_database(sqlite_name, shape_name, r, h, __func__, l, 0, cte * result);
	}

	return cte * result;
}

template <typename T>
T Iij_spheroidal_summation_approx(float r, float h, unsigned i, unsigned j, bool save = true) {
	if ((i+j) % 2 == 1) return T(0);

	{
		Data<T> res = SQLite::get_from_database<T>(sqlite_name, shape_name, r, h, __func__, i, j);
		if (res.found == true) return res.result;
	}

	/*{
		Data<T> res = CSVDB::get_from_database<T>(csv_name, shape_name, r, h, __func__, i, j);
		if (res.found == true) SQLite::save_at_database(sqlite_name, shape_name, r, h, __func__, i, j, res.result);
		if (res.found == true) return res.result;
	}*/

	T radius = static_cast<T>(r);
	T height = static_cast<T>(h);

	T ij = static_cast<T>((i+j) / 2);
	T pi = T(math::cte::pi_string);
	T cte = 2.0 * pi / boost::multiprecision::pow(radius, i+j);
	T ee2 = 1.0 - radius * radius / (height * height);
	math::Polynomial<T> i_legendre = math::Polynomial<T>::legendre(i);
	math::Polynomial<T> j_legendre = math::Polynomial<T>::legendre(j);
	T result = 0.0;
	T bin = 1.0;

	for (unsigned n = 0; n <= (i+j) / 2; ++n) {
		T mul = 1;
		T nt = static_cast<T>(n);

		if (n % 2 == 0) mul = 1;
		else mul = -1; 

		for (unsigned u = 0; u <= i; ++u) {
			for (unsigned v = 0; v <= j; ++v) {
				if ((u + v) % 2 == 1) continue;

				T uv = static_cast<T>((u+v) / 2);
				T main = i_legendre[u] * j_legendre[v] * bin * boost::multiprecision::pow(ee2, n);
				T moment = 1.0 / (nt+uv+0.5) + ee2 * ee2 / 2.0 / (nt+uv+1.5);

				std::cout << "(" << ij << ";" << nt << ") = " << bin << std::endl;

				result += mul * main * moment;
			}
		}

		bin *= (ij - nt) / (nt + 1.0);
	}

	if (save) {
		// CSVDB::save_at_database(csv_name, shape_name, r, h, __func__, i, j, cte * result);
		SQLite::save_at_database(sqlite_name, shape_name, r, h, __func__, i, j, cte * result);
	}

	return cte * result;
}

template <typename T>
T Iij_spheroidal_summation_exact(float r, float h, unsigned i, unsigned j, const std::vector<T>& An = {}, bool save = true) {
	if ((i+j) % 2 == 1) return T(0);

	{
		Data<T> res = SQLite::get_from_database<T>(sqlite_name, shape_name, r, h, __func__, i, j);
		if (res.found == true) return res.result;
	}

	/*{
		Data<T> res = CSVDB::get_from_database<T>(csv_name, shape_name, r, h, __func__, i, j);
		if (res.found == true) SQLite::save_at_database(sqlite_name, shape_name, r, h, __func__, i, j, res.result);
		if (res.found == true) return res.result;
	}*/

	T radius = static_cast<T>(r);
	T height = static_cast<T>(h);


	T ij = static_cast<T>((i+j) / 2);
	T pi = T(math::cte::pi_string);
	T cte = 2 * pi / boost::multiprecision::pow(radius, i+j);
	T ee2 = 1.0 - radius * radius / (height * height);
	math::Polynomial<T> i_legendre = math::Polynomial<T>::legendre(i);
	math::Polynomial<T> j_legendre = math::Polynomial<T>::legendre(j);
	T result = 0.0;
	T bin = 1.0;

	for (unsigned n = 0; n <= (i+j) / 2; ++n) {
		T mul = 1;
		T nt = static_cast<T>(n);

		if (n >= 1) bin *= (ij - nt+1.0) / nt;
		if (n % 2 == 0) mul = 1;
		else mul = -1; 

		for (unsigned u = 0; u <= i; ++u) {
			for (unsigned v = 0; v <= j; ++v) {
				T main = i_legendre[u] * j_legendre[v] * bin * boost::multiprecision::pow(ee2, n);

				T moment;
				if (An.empty()) moment = integrate<T>(-1, 1, An_integrand<T>(radius, height, u+v+2*n));
				else if (An.size() <= u+v+2*n) moment = integrate<T>(-1, 1, An_integrand<T>(radius, height, u+v+2*n));
				else moment = An[u+v+2*n];

				// std::cout << "(" << ij << ";" << nt << ") = " << bin << std::endl;

				result += mul * main * moment;
			}
		}
	}

	if (save) {
		// CSVDB::save_at_database(csv_name, shape_name, r, h, __func__, i, j, cte * result);
		SQLite::save_at_database(sqlite_name, shape_name, r, h, __func__, i, j, cte * result);
	}

	return cte * result;
}

/*
61) -7.12592
62) 0
63) -10.2899
64) 0
65) -14.6473
66) 0
67) -20.5532
68) 0
69) -28.6212
70) 0
71) -38.462
72) 0
73) -57.2676
74) 0
75) -45.1968
76) 0
77) -250.788
78) 0
79) 806.525
80) 0
81) -5928.18
82) 0
83) 37925.5
84) 0
85) -271662
86) 0
87) 2.04369e+06
88) 0
89) -1.59791e+07
90) 0
91) 1.28107e+08
92) 0
93) -1.04646e+09
94) 0
95) 8.6402e+09
96) 0
97) -7.12399e+10
98) 0
99) 5.78138e+11


59) n=30|2.06789e-05
60) 0
61) n=31|-0.000170966
62) 0
63) n=31|0.00140328
64) 0
65) n=31|-0.0116627
66) 0
67) n=31|0.10178
68) 0
69) n=31|-0.968016
70) 0
71) n=31|10.1402
72) 0
73) n=32|-113.846
74) 0
75) n=33|1308.52
76) 0
77) n=34|-14766.9
78) 0
79) n=35|158754
80) 0
81) n=36|-1.59118e+06
82) 0
83) n=37|1.45932e+07
84) 0
85) n=38|-1.19654e+08
86) 0
87) n=39|8.40992e+08
88) 0
89) n=40|-4.53517e+09
90) 0
91) n=41|9.94845e+09
92) 0
93) n=42|1.69708e+11
94) 0
95) n=43|-3.2921e+12
96) 0
97) n=44|3.84219e+13
98) 0
99) n=45|-3.64689e+14
*/

template <typename T>
void populate_an_data(const std::string& name, const T& radius, const T& height, unsigned group = 20) {
	using SmallFloat = boost::multiprecision::float128;
	std::ifstream file(name, std::ios::in);
	std::vector<std::string> lines;
	std::string content;

	// Load lines.
	while (getline(file, content)) lines.push_back(content);
	file.close();

	// Calculate 20 more.
	unsigned next = 2 * static_cast<unsigned>(lines.size());
	for (unsigned i = 0; i < group; ++i) {
		SmallFloat sradius = static_cast<SmallFloat>(radius);
		SmallFloat sheight = static_cast<SmallFloat>(height);
		SmallFloat an = integrate<SmallFloat>(-1.0, 1.0, An_integrand(sradius, sheight, next));
		next += 2;
		
		std::stringstream ss;
		ss << std::setprecision(std::numeric_limits<SmallFloat>::max_digits10) << an;
		
		std::string san = ss.str();
		lines.push_back(san);
	}

	// Save on file.
	std::ofstream output(name, std::ios::out);
	for (const std::string& num : lines) output << num << std::endl;
	output.close();
}

template <typename T>
T get_An(const std::string& name, const T& radius, const T& height, unsigned n) {
	if (n % 2 == 1) return T(0);

	using SmallFloat = boost::multiprecision::float128;
	std::ifstream file(name, std::ios::in);
	std::string content = "";

	// Load lines.
	unsigned line = 0;
	while (getline(file, content) and line != n/2) ++line;
	file.close();

	// The file contains the number. Load it and return.
	if (line == n/2 and !content.empty()) return T(content);

	// The file doesn't contain the number. Calculate an.
	SmallFloat sradius = static_cast<SmallFloat>(radius);
	SmallFloat sheight = static_cast<SmallFloat>(height);
	SmallFloat an = integrate<SmallFloat>(-1.0, 1.0, An_integrand(sradius, sheight, n));

	// Convert to string
	std::stringstream ss;
	ss << std::setprecision(std::numeric_limits<SmallFloat>::max_digits10) << an;
	std::string san = ss.str();

	// Load lines.
	std::vector<std::string> lines;
	file.open(name, std::ios::in);
	while (getline(file, content)) lines.push_back(content);
	file.close();

	// Save computed an value at the lines.
	unsigned next = 2*lines.size();
	if (n < next) lines[n/2] = san;
	else {
		while (next < n) {
			lines.push_back("");
			next += 2;
		}
		lines.push_back(san);
	}
	
	// Output at the file.
	std::ofstream output(name, std::ios::out);
	for (const std::string& content : lines) output << content << std::endl;
	output.close();

	// Return the value.
	return static_cast<T>(an);
}

template <typename T>
void get_detail(const T& height, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& I, const Eigen::Matrix<T, Eigen::Dynamic, 1>& G, const T& real_fef = 0.0) {
	std::ofstream file("data/aareport.dat", std::ios::out);
	std::ofstream plot("data/aaplotter.m", std::ios::out);
	
	std::vector<T> orders;
	std::vector<T> fefs;
	std::vector<T> errors;
	bool has_real_fef = (real_fef > 1e-4);
	for (unsigned dim = 1; dim <= I.rows(); ++dim) {
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Ic(dim, dim);
		Eigen::Matrix<T, Eigen::Dynamic, 1> Gc(dim);
		
		for (unsigned i = 0; i < dim; ++i) {
			Gc(i) = G(i);
			for (unsigned j = 0; j < dim; ++j) Ic(i, j) = I(i, j);
		}

		Eigen::Matrix<T, Eigen::Dynamic, 1> Al = Ic.colPivHouseholderQr().solve(Gc);
		T FEF = get_FEF<T>(Al, height);

		std::cout << dim-1 << ") FEF: " << FEF;
		if (has_real_fef) std::cout << ",\terr: " << (real_fef - FEF) / real_fef;
		std::cout << std::endl;

		if ((dim-1) % 2 == 1) {
			file << dim-1 << "  " << FEF;
			orders.push_back(dim-1);
			fefs.push_back(FEF);
			
			if (has_real_fef) {
				T error = (real_fef - FEF) / real_fef;
				file << "  " << error;
				errors.push_back(error);
			}

			file << std::endl;
		}
	}
	file.close();

	plot << "orders = [" << orders[0];
	for (unsigned i = 1; i < orders.size(); ++i) plot << ", " << orders[i];
	plot << "];\n" << std::endl;
	
	plot << "fefs = [" << fefs[0];
	for (unsigned i = 1; i < fefs.size(); ++i) plot << ", " << fefs[i];
	plot << "];\n" << std::endl;

	if (has_real_fef) {
		plot << "errors = [" << errors[0];
		for (unsigned i = 1; i < errors.size(); ++i) plot << ", " << errors[i];
		plot << "];\n" << std::endl;
	}

	if (has_real_fef) {
		std::function<T(T)> logarithm = [=](const T& value){return std::log(value);};
		// std::function<T(T)> logarithm = [=](const T& value){return boost::multiprecision::log(value);};
		auto data = leastSquares(orders, errors, logarithm);

		std::cout << "log(error) = " << data.independent() << " + " << data.linear() << " order" << std::endl;
		std::cout << "error = " << std::exp(data.independent()) << " exp(" << data.linear() << " order)" << std::endl;
		std::cout << "error = " << std::exp(data.independent()) << " exp(-order/" << -1.0 / data.linear() << ")" << std::endl;
		std::cout << "Relative Error: " << std::sqrt(data.relative_quadratic_error()) << std::endl;
		std::cout << "Absolute Error: " << std::sqrt(data.absolute_quadratic_error()) << std::endl;

		plot << "error_fit = [" << std::exp(data.independent() + orders[0] * data.linear());
		for (unsigned i = 1; i < errors.size(); ++i) plot << ", " << std::exp(data.independent() + orders[i] * data.linear());
		plot << "];\n" << std::endl;

		plot << "hold on\n";
		plot << "semilogy(orders, errors, '+r')\n";
		plot << "semilogy(orders, error_fit)\n";
		plot << std::endl;
	}

	plot.close();
}

void legendre_test() {
	//using TestFloat = double;
	using TestFloat = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<8192>>;
	math::Polynomial<TestFloat> legendre = math::Polynomial<TestFloat>::legendre(501);
	std::cout << legendre << std::endl;
	std::cout << std::endl;

	TestFloat x = -1;
	do {
		std::cout << legendre.at(x) << ", " << std::flush;
		x += 0.05;
	} while (x <= 1);
	std::cout << legendre.at(x) << std::endl;
}
/*
int main() {
	using Float = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<8192>>;
	using SmallFloat = boost::multiprecision::float128;

	float radius = 1.0;
	float height = 2.0;
	unsigned max_order = 100;

	SmallFloat sradius = static_cast<SmallFloat>(radius);
	SmallFloat sheight = static_cast<SmallFloat>(height);
	std::string name = "data/an_spheroidal_h" + std::to_string(static_cast<unsigned>(height));

	std::vector<SmallFloat> an;
	unsigned num = static_cast<unsigned>(boost::multiprecision::ceil(sheight / sradius));
	unsigned total = max_order + 1;
	unsigned an_total = 4*total*num*num + 30;
	for (unsigned i = 0; i < an_total; ++i) {
		std::cout << "\r                                              ";
		std::cout << "\rIntegrating: " << i << "/" << an_total-1 << "\t" << double(i) / (an_total-1) * 100.0 << "%" << std::flush;
		if (i % 2 == 0) an.push_back(get_An(name, sradius, sheight, i));
		else an.push_back(0.0);
	}

	std::cout << std::endl;
	std::vector<Float> An;
	for (const SmallFloat& value : an) An.push_back(static_cast<Float>(value));

	std::cout << Gl_spheroidal_summation_exact<Float>(radius, height, 37, An, false) << std::endl;
	std::cout << Gl_spheroidal_summation_exact<Float>(radius, height, 39, An, false) << std::endl;
	std::cout << Gl_spheroidal_summation_exact<Float>(radius, height, 41, An, false) << std::endl;
}*/


int main() {
	using Float = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<8192>>;
	using SmallFloat = boost::multiprecision::float128;
	using EigenFloat = double;

	float radius = 1.0;
	float height = 1.1;
	unsigned max_order = 29;

	EigenFloat aspect_ratio = static_cast<EigenFloat>(height / radius);
	EigenFloat xi = std::sqrt(aspect_ratio * aspect_ratio - 1.0);
	EigenFloat real_fef = xi * xi * xi / (aspect_ratio * std::log(aspect_ratio + xi) - xi);
	std::cout << "Theoretical FEF: " << real_fef << std::endl;

	SmallFloat sradius = static_cast<SmallFloat>(radius);
	SmallFloat sheight = static_cast<SmallFloat>(height);

	SmallFloat ecc = boost::multiprecision::sqrt(1 - sradius*sradius / sheight/sheight);
	std::cout << "Eccentricity: " << ecc << std::endl;

	std::stringstream ss;
	ss << "data/an_spheroidal_h" << height;

	std::string name = ss.str();
	std::cout << name << std::endl;
	std::cout << sqlite_name << std::endl;

	std::vector<SmallFloat> an;
	unsigned num = static_cast<unsigned>(boost::multiprecision::ceil(sheight / sradius));
	unsigned total = max_order + 1;
	unsigned an_total = 4*total*num*num + 30;
	for (unsigned i = 0; i < an_total; ++i) {
		std::cout << "\r                                              ";
		std::cout << "\rIntegrating: " << i << "/" << an_total-1 << "\t" << double(i) / (an_total-1) * 100.0 << "%" << std::flush;
		if (i % 2 == 0) an.push_back(get_An(name, sradius, sheight, i));
		else an.push_back(0.0);
	}

	std::cout << std::endl;
	std::vector<Float> An;
	for (const SmallFloat& value : an) An.push_back(static_cast<Float>(value));

	Eigen::Matrix<EigenFloat, Eigen::Dynamic, 1> Gl(total);
	for (unsigned l = 0; l < total; ++l) {
		std::cout << "\r                                              ";
		std::cout << "\rG summation: " << l << "/" << total-1 << "\t" << double(l) / (total-1) * 100.0 << "%" << std::flush;
		Float current = Gl_spheroidal_summation_exact<Float>(radius, height, l, An);
		Gl(l) = static_cast<EigenFloat>(current);
	}

	std::cout << std::endl;
	Eigen::Matrix<EigenFloat, Eigen::Dynamic, Eigen::Dynamic> Iij(total, total);
	for (unsigned i = 0; i < total; ++i) {
		for (unsigned j = 0; j < total; ++j) {
			unsigned c = i*total + j;
			std::cout << "\r                                              ";
			std::cout << "\rI summation: " << c << "/" << total*total-1 << "\t" << double(c) / (total*total -1) * 100.0 << "%" << std::flush;
			Float current = Iij_spheroidal_summation_exact<Float>(radius, height, i, j, An);
			Iij(i, j) = static_cast<EigenFloat>(current);
		}
	}

	std::cout << std::endl;	
	Eigen::Matrix<EigenFloat, Eigen::Dynamic, 1> Al = Iij.colPivHouseholderQr().solve(Gl);
	std::vector<EigenFloat> contribution = get_FEF_contributions(Al, static_cast<EigenFloat>(height));
	EigenFloat FEF = get_FEF<EigenFloat>(Al, static_cast<EigenFloat>(height));

	std::cout << "[" << contribution[0];
	for (unsigned i = 1; i < contribution.size(); ++i) std::cout << ", " << contribution[i];
	std::cout << "]" << std::endl;
	std::cout << std::endl;

	std::cout << "FEF: " << FEF << std::endl;
	std::cout << std::endl;

	std::cout << Gl << std::endl;
	std::cout << std::endl;

	std::cout << Iij << std::endl;
	std::cout << std::endl;

	std::cout << Al << std::endl;
	std::cout << std::endl;

	// Convergence test.
	for (unsigned k = 0; k < total-2; ++k) {
		if (std::abs(Gl(k+2)) > std::abs(Gl(k))) {
			std::cout << "Convergence fail G: " << k << std::endl;
		}
	}

	// Convergence test.
	for (unsigned k = 0; k < total-1; ++k) {
		if (std::abs(Iij(k+1, k+1)) > std::abs(Iij(k, k))) {
			std::cout << "Convergence fail I: " << k << std::endl;
		}
	}

	// Display results.
	// get_detail<EigenFloat>(static_cast<EigenFloat>(height), Iij, Gl, 5.76156);
	get_detail<EigenFloat>(static_cast<EigenFloat>(height), Iij, Gl, real_fef);
	std::cout << "Theoretical FEF: " << real_fef << std::endl;
}
