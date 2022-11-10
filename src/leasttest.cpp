#include <iostream>
#include <leastsquares.hpp>
#include <cmath>
#include <functional>


int main() {
	std::vector<double> x = {1, 2, 3, 4};
	std::vector<double> y = {1, 4, 9, 16};
	std::vector<double> z = {2, 4, 8, 16};

	std::function<double(double)> asd = [](double a){return std::log(a);};

	auto data = leastSquares(x, z, asd);
	std::cout << "y = " << data.linear() << "x + " << data.independent() << std::endl;

	Eigen::Matrix<double, Eigen::Dynamic, 1> xx(4);
	for (unsigned i = 0; i < 4; ++i) xx(i) = x[i];

	Eigen::Matrix<double, Eigen::Dynamic, 1> yy(4);
	for (unsigned i = 0; i < 4; ++i) yy(i) = y[i];

	Eigen::Matrix<double, Eigen::Dynamic, 1> zz(4);
	for (unsigned i = 0; i < 4; ++i) zz(i) = z[i];

	data = leastSquares(xx, zz, asd);
	std::cout << "y = " << data.linear() << "x + " << data.independent() << std::endl;
}
