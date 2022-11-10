#include <get_data.hpp>
#include <boost/multiprecision/float128.hpp> 
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <math.hpp>
#include <calcs.hpp>
#include <leastsquares.hpp>

#include <Eigen/Dense>
#include <fstream>

template <typename T>
void calculate_coeficients(float radius, float height, unsigned max) {
	float length = height - radius;
	float r0 = std::sqrt(radius * radius + length * length);

	std::cout << "Length: " << length << std::endl;
	std::cout << "Height: " << height << std::endl;

	for (unsigned l = 0; l < max; ++l) {
		if (l % 2 == 0) continue;
		std::cout << l << ") " << get_Gl_hcp_integral<T>(radius, height, l) << std::endl;
	}

	for (unsigned i = 0; i < max; ++i) {
		std::cout << i << "," << i << ") " << get_Iij_hcp_integral<T>(radius, height, i, i) << std::endl;
		for (unsigned j = 0; j < i; ++j) {
			if ((i+j) % 2 == 1) continue;
			std::cout << i << "," << j << ") " << get_Iij_hcp_integral<T>(radius, height, i, j) << std::endl;
			std::cout << j << "," << i << ") " << get_Iij_hcp_integral<T>(radius, height, j, i) << std::endl;
		}
	}
}

/*
MATRIX:
	n=0		n=1		n=2		n=3		n=4		...
A0
A1	
A2
A3
A4
A5
...
*/
template <typename T>
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
calculate_result_matrix(float radius, float height, unsigned order, bool fef = false, double real_fef = -1.0) {
	unsigned size;
	if (fef) size = order + 1;
	else size = order;
	if (real_fef > 0.0) size = order + 2;

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> numeric_result(size, order);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> I(order, order);
	Eigen::Matrix<T, Eigen::Dynamic, 1> G(order);

	for (unsigned l = 0; l < order; ++l) {
		if (l % 2 == 0) G(l) = 0;
		else G(l) = get_Gl_hcp_integral<T>(radius, height, l);
		std::cout << l << " \r" << std::flush;
	}

	for (unsigned i = 0; i < order; ++i) {
		for (unsigned j = 0; j < order; ++j) {
			if ((i+j) % 2 == 1) I(i, j) = T(0.0);
			else I(i, j) = get_Iij_hcp_integral<T>(radius, height, i, j);
			std::cout << "(" << i << ", " << j << ")   \r" << std::flush;
		}
	}

	for (unsigned dim = 1; dim <= I.rows(); ++dim) {
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Ic(dim, dim);
		Eigen::Matrix<double, Eigen::Dynamic, 1> Gc(dim);
		
		for (unsigned i = 0; i < dim; ++i) {
			Gc(i) = static_cast<double>(G(i));
			for (unsigned j = 0; j < dim; ++j) Ic(i, j) = static_cast<double>(I(i, j));
		}

		Eigen::Matrix<double, Eigen::Dynamic, 1> Al = Ic.colPivHouseholderQr().solve(Gc);
	
		double FEF = get_FEF<double>(Al, height);
		for (unsigned n = 0; n < dim; ++n) numeric_result(n, dim-1) = Al(n);
		for (unsigned n = dim; n < order; ++n) numeric_result(n, dim-1) = 0;
		if (fef == true) numeric_result(order, dim-1) = FEF;
		if (real_fef > 0.0) numeric_result(order+1, dim-1) = std::abs(real_fef - FEF) / real_fef;
	}

	return numeric_result;
}

// To modify this section using the calculate_result_matrix function.
template <typename T>
void calculate_fef(float radius, float height, unsigned order, double real_fef = -1.0) {
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> I(order, order);
	Eigen::Matrix<T, Eigen::Dynamic, 1> G(order);

	for (unsigned l = 0; l < order; ++l) {
		if (l % 2 == 0) G(l) = 0;
		else G(l) = get_Gl_hcp_integral<T>(radius, height, l);
		std::cout << l << " \r" << std::flush;
	}

	for (unsigned i = 0; i < order; ++i) {
		for (unsigned j = 0; j < order; ++j) {
			if ((i+j) % 2 == 1) I(i, j) = T(0.0);
			else I(i, j) = get_Iij_hcp_integral<T>(radius, height, i, j);
			std::cout << "(" << i << ", " << j << ")   \r" << std::flush;
		}
	}
	std::cout << std::endl;

	std::stringstream ss;
	ss << "../data/results/hcp_r" << radius;
	ss << "h" << height;
	ss << ".dat";

	bool has_real_fef = (real_fef > 1e-4);
	std::ofstream file(ss.str());
	for (unsigned dim = 1; dim <= I.rows(); ++dim) {
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Ic(dim, dim);
		Eigen::Matrix<double, Eigen::Dynamic, 1> Gc(dim);
		
		for (unsigned i = 0; i < dim; ++i) {
			Gc(i) = static_cast<double>(G(i));
			for (unsigned j = 0; j < dim; ++j) Ic(i, j) = static_cast<double>(I(i, j));
		}

		Eigen::Matrix<double, Eigen::Dynamic, 1> Al = Ic.colPivHouseholderQr().solve(Gc);
		T FEF = get_FEF<double>(Al, height);

		std::cout << dim-1 << ") FEF: " << FEF;
		if (has_real_fef) {
			std::cout << ",\terr: " << (real_fef - FEF) / real_fef;
		}
		std::cout << std::endl;

		if ((dim-1) % 2 == 1) {
			file << dim-1 << "  " << std::setprecision(std::numeric_limits<double>::max_digits10) << FEF;

			if (has_real_fef) {
				T error = (real_fef - FEF) / real_fef;
				file << "  " << error;
			}

			file << std::endl;
		}
	}

	file.close();
}

template <typename T>
std::string python_write(const std::vector<T>& vec, const std::string& name) {
	std::stringstream result;
	result << name << " = [" << vec[0];
	for (unsigned i = 0; i < vec.size(); ++i) result << ", " << vec[i];
	result << "]";
	return result.str();
}

template <typename T>
void calculate_ab_error_coeficients(float radius, const std::vector<float>& heights, unsigned order, const std::vector<double>& mds_fefs) {
	unsigned size = heights.size();

	std::vector<double> a, b;
	for (unsigned k = 0; k < size; ++k) {
		auto matrix = calculate_result_matrix<T>(radius, heights[k], order, true, mds_fefs[k]);

		std::vector<double> orders, fefs, errors;
		for (unsigned i = 0; i < order; ++i) {
			if (i % 2 == 1) {
				orders.push_back(i);
				fefs.push_back(matrix(order, i));
				errors.push_back(matrix(order+1, i));
			}
		}
		
		// Least squares calculation.
		std::function<double(double)> logarithm = [=](double value){return std::log(value);};
		auto data = leastSquares(orders, errors, logarithm);
		a.push_back(std::exp(data.independent()));
		b.push_back(data.linear());
		std::cout << "(" << heights[k] << ", " << std::exp(data.independent()) << ", " << data.linear() << ")" << std::endl;
	}

	std::cout << std::endl;
	std::cout << python_write(heights, "r") << std::endl;
	std::cout << python_write(a, "a") << std::endl;
	std::cout << python_write(b, "b") << std::endl;
}

/*
h/r      Apex-FEF (MDS)
 2            4.20577
 3            5.28913
 4            6.30749
 5            7.28259
 6            8.22638
11           12.65419

r = [3, 4, 5, 6, 11]
a = [0.783828, 0.813784, 0.846263, 0.868566, 0.919471]
b = [-0.0388824, -0.0088113, -0.00323957, -0.00153679, -0.000160214]

*/
int main() {
	using SmallFloat = boost::multiprecision::float128;
	float radius = 1.0;
	float height = 2.0;	
	unsigned order = 20;

	std::vector<double> r = {3, 4, 5, 6, 11};
	std::vector<double> b = {0.0388824, 0.0088113,0.00323957,0.00153679,0.000160214};
	std::function<double(double)> logarithm = [=](double value){return std::log(value);};

	auto data = leastSquares(r, b, logarithm, logarithm);
	std::cout << data.independent() << std::endl;
	std::cout << data.linear() << std::endl;
	

	/*std::vector<float> heights = {2.0,     3.0,     4.0,     5.0,     6.0,     11.0};
	std::vector<double> mdsfef = {4.20577, 5.28913, 6.30749, 7.28259, 8.22638, 12.65419};
	calculate_ab_error_coeficients<SmallFloat>(radius, heights, order, mdsfef);*/
	
	/*std::vector<float> calculated_heights = {1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 11.0};
	for (float calculated_height : calculated_heights) {
		std::cout << "Height: " << calculated_height << std::endl;
		calculate_coeficients<SmallFloat>(radius, calculated_height, order);
		calculate_fef<SmallFloat>(radius, calculated_height, order);
	}*/


	// calculate_coeficients<SmallFloat>(radius, height, order);
	// calculate_fef<SmallFloat>(radius, height, order);
}	

void cylinder_test() {
	// using Float = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<8192>>;
	using SmallFloat = boost::multiprecision::float128;
	float pi = 3.1415926535;

	unsigned l = 1;
	unsigned i = 1;
	unsigned j = 1;
	float radius = 1.0;
	float length = 0.5;
	float height = radius + length;
	float r0 = std::sqrt(radius * radius + length * length);

	/*SmallFloat g = get_Gl_cylindrical_integral<SmallFloat>(radius, height, l, false);
	std::cout << g << std::endl;

	float analytic = 4.0 * 3.1415926535 * radius * (std::atanh(length / r0) - length / r0);
	std::cout << analytic << std::endl;

	float analytic_approx = 4.0 / 3.0 * 3.1415926535 * radius * std::pow(length/r0, 3);
	std::cout << analytic_approx << std::endl;*/

	// --

	/*
	SmallFloat I = get_Iij_cylindrical_integral<SmallFloat>(radius, height, i, j, false);
	std::cout << I << std::endl;

	float analytic00 = 2*pi * (pi - 2*std::acos(length / r0));
	std::cout << analytic00 << std::endl;

	float analytic11 = 2*pi / radius / radius * (pi/8.0 - 1.0/4.0 * std::acos(length / r0) + 1.0/4.0 * (radius*length/r0/r0) * (length*length - radius*radius)/r0/r0);
	std::cout << analytic11 << std::endl;*/
}
