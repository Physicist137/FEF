#include <get_data.hpp>
#include <calcs.hpp>
#include <leastsquares.hpp>

#include <boost/multiprecision/float128.hpp> 
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <Eigen/Dense>


int main() {
	using Float = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<8192>>;
	using SmallFloat = boost::multiprecision::float128;
	using EigenFloat = double;

	float radius = 1.0;
	float height = 2.0;
	unsigned max_order = 29;
	unsigned total = max_order + 1;

	// Calculate Theoretical FEF.
	double aspect_ratio = static_cast<EigenFloat>(height / radius);
	double xi = std::sqrt(aspect_ratio * aspect_ratio - 1.0);
	double real_fef = xi * xi * xi / (aspect_ratio * std::log(aspect_ratio + xi) - xi);
	std::cout << "Aspect Ratio: " << aspect_ratio << std::endl;
	std::cout << "Theoretical FEF: " << real_fef << std::endl;

	// Display eccentricity.
	float ecc = math::sqrt(1 - radius*radius / height / height);
	std::cout << "Eccentricity: " << ecc << std::endl;

	// Calculate An sequence.
	unsigned num = static_cast<unsigned>(std::ceil(height / radius));
	unsigned an_total = 4*total*num*num + 30;
	for (unsigned i = 0; i < an_total; ++i) {
		std::cout << "\r                                              ";
		std::cout << "\rIntegrating: " << i << "/" << an_total-1 << "\t" << double(i) / (an_total-1) * 100.0 << "%" << std::flush;
		get_An_spheroidal_integral<SmallFloat>(radius, height, i);
	}

	// Calculate Gl vectors.
	std::cout << std::endl;
	Eigen::Matrix<EigenFloat, Eigen::Dynamic, 1> Gl(total);
	for (unsigned l = 0; l < total; ++l) {
		std::cout << "\r                                              ";
		std::cout << "\rG summation: " << l << "/" << total-1 << "\t" << double(l) / (total-1) * 100.0 << "%" << std::flush;
		Float current = get_Gl_spheroidal_summation<Float>(radius, height, l);
		Gl(l) = static_cast<EigenFloat>(current);
	}

	// Calculate Iij matrix.
	std::cout << std::endl;
	Eigen::Matrix<EigenFloat, Eigen::Dynamic, Eigen::Dynamic> Iij(total, total);
	for (unsigned i = 0; i < total; ++i) {
		for (unsigned j = 0; j < total; ++j) {
			unsigned c = i*total + j;
			std::cout << "\r                                              ";
			std::cout << "\rI summation: " << c << "/" << total*total-1 << "\t" << double(c) / (total*total -1) * 100.0 << "%" << std::flush;
			Float current = get_Iij_spheroidal_summation<Float>(radius, height, i, j);
			Iij(i, j) = static_cast<EigenFloat>(current);
		}
	}

	// Prepare to save all data.
	std::stringstream ss;
	ss << "../data/results/ellipsoid_r" << radius;
	ss << "h" << height;
	ss << ".dat";

	// Calculate FEF sequence and save.
	std::vector<double> orders;
	std::vector<double> fefs;
	std::vector<double> errors;
	std::ofstream file(ss.str());
	for (unsigned dim = 1; dim <= Iij.rows(); ++dim) {
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Ic(dim, dim);
		Eigen::Matrix<double, Eigen::Dynamic, 1> Gc(dim);
		
		for (unsigned i = 0; i < dim; ++i) {
			Gc(i) = static_cast<double>(Gl(i));
			for (unsigned j = 0; j < dim; ++j) Ic(i, j) = static_cast<double>(Iij(i, j));
		}

		Eigen::Matrix<double, Eigen::Dynamic, 1> Al = Ic.colPivHouseholderQr().solve(Gc);
		double FEF = get_FEF<double>(Al, height);


		if ((dim-1) % 2 == 1) {
			file << dim-1 << "  " << std::setprecision(std::numeric_limits<double>::max_digits10) << FEF;

			double error = (real_fef - FEF) / real_fef;
			orders.push_back(dim-1);
			fefs.push_back(FEF);
			errors.push_back(error);

			file << "  " << std::setprecision(std::numeric_limits<double>::max_digits10) << error;
			file << std::endl;

			std::cout << dim-1 << ") FEF: " << FEF;
			std::cout << ",\terr: " << error;
			std::cout << std::endl;
		}
	}

	file.close();


	// Least squares calculation.
	std::function<double(double)> logarithm = [=](double value){return std::log(value);};
	auto data = leastSquares(orders, errors, logarithm);


	// Show exponencial fit.
	std::cout << "log(error) = " << data.independent() << " + " << data.linear() << " order" << std::endl;
	std::cout << "error = " << std::exp(data.independent()) << " exp(" << data.linear() << " order)" << std::endl;
	std::cout << "error = " << std::exp(data.independent()) << " exp(-order/" << -1.0 / data.linear() << ")" << std::endl;
	std::cout << "Relative Error: " << std::sqrt(data.relative_quadratic_error()) << std::endl;
	std::cout << "Absolute Error: " << std::sqrt(data.absolute_quadratic_error()) << std::endl;
}
