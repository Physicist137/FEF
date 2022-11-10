#include <math.hpp>
#include <integrand.hpp>
#include <integral.hpp>
#include <spheroid.hpp>
#include <sqlitedb.hpp>
#include <csvdb.hpp>
#include <get_data.hpp>

int main() {
	unsigned n = 2;
	std::string shape_name = "prolate_spheroid";
	float radius = 1.0;
	float height = 1.1;
	std::cout << std::to_string(height) << std::endl;

	using Float128 = boost::multiprecision::float128;
	Data<Float128> res = CSVDB::get_from_database<Float128>(
		"../data/database/spheroidal/integrals_r1h1.1.csv",
		shape_name,
		radius,
		height,
		"get_An_spheroidal_integral",
		2,
		0
	);

	std::cout << res.found << std::endl;
	std::cout << res.result << std::endl;



    std::stringstream ss;
	ss << "../data/database/spheroidal/integrals_r"  << radius << "h" << height << ".csv";
	Data<Float128> r = CSVDB::get_from_database<Float128>(ss.str(), shape_name, radius, height, "get_An_spheroidal_integral", n, 0);
	std::cout << r.found << std::endl;
	std::cout << r.result << std::endl;
	
	

	std::cout << get_An_spheroidal_integral<Float128>(1.0, 1.1, 20) << std::endl;
	
}
