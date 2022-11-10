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
#include <math.hpp>

#include <Eigen/Dense>


template <typename T>
std::function<T(const T&)>
Gl_suspended_hemispherical_integrand_theta(const T& radius, const T& length, unsigned l) {
	return [=](const T& theta) {
		int p = static_cast<int>(l);

		T pi = 3.1415926535;
		T cos = math::cos(theta);
		T sin = math::sin(theta);
		T sq1 = math::sqrt(radius * radius - length * length * sin * sin);

		math::Polynomial<T> l_legendre = math::Polynomial<T>::legendre(l);
		
		T r;
		if (theta < pi/2) r = length * cos + sq1;
		else r = -length * cos + sq1;

		T factor1 = length / r * sin;
		T factor2 = 1.0 + length * cos / sq1;
		T JA = math::sqrt(1.0 + factor1 * factor1 * factor2 * factor2);

		return 2.0 * pi * math::pow(r, -p+2) * l_legendre.at(cos) * cos * sin * JA;
	};
}


template <typename T>
std::function<T(const T&)>
Gl_suspended_hemispherical_integrand_approx(const T& radius, const T& length, unsigned l) {
	if (l % 2 == 0) {
		return [=](const T& x) {
			return T(0.0);
		};
	}

	return [=](const T& x) {
		int p = static_cast<int>(l);
		int n = -p+2;

		T pi = 3.1415926535;
		math::Polynomial<T> l_legendre = math::Polynomial<T>::legendre(l);

		T nt = static_cast<T>(n);
		T term1 = math::pow(radius, n);
		T term2 = nt * x * math::pow(radius, n-1) * length;
		T term3 = ((nt * nt - 2.0) * x*x - (nt - 2.0)) * math::pow(radius, n-2) * length * length;

		//T factor1 = length * length / r / r * s;
		//T factor2 = 1.0 + length * x / sq1;
		//T JA = math::sqrt(1.0 + factor1 * factor2 * factor2);
		T rnJA = term1 + term2 + term3;
		return 4.0 * pi * l_legendre.at(x) * x * rnJA;
	};
}

template <typename T>
std::function<T(const T&)>
Iij_suspended_hemispherical_integrand_approx(const T& radius, const T& length, unsigned i, unsigned j) {
	if ((i+j) % 2 == 1) {
		return [=](const T& x) {
			return T(0.0);
		};
	}

	return [=](const T& x) {
		int p = static_cast<int>(i+j);
		int n = -p;

		T pi = 3.1415926535;
		math::Polynomial<T> i_legendre = math::Polynomial<T>::legendre(i);
		math::Polynomial<T> j_legendre = math::Polynomial<T>::legendre(j);

		T nt = static_cast<T>(n);
		T term1 = math::pow(radius, n);
		T term2 = nt * x * math::pow(radius, n-1) * length;
		T term3 = ((nt * nt - 2.0) * x*x - (nt - 2.0)) * math::pow(radius, n-2) * length * length;

		//T factor1 = length * length / r / r * s;
		//T factor2 = 1.0 + length * x / sq1;
		//T JA = math::sqrt(1.0 + factor1 * factor2 * factor2);
		T rnJA = term1 + term2 + term3;
		return 4.0 * pi * i_legendre.at(x) * j_legendre.at(x) * rnJA;
	};
}

template <typename T>
std::function<T(const T&)>
Gl_cylindrical_integrand(const T& radius, const T& length, unsigned l) {
	return [=](const T& x) {
		
	};
}


void checking_Gl();
void checking_Iij();

int main() {
	checking_Gl();
	// checking_Iij();
}

void checking_Iij() {
	double radius = 1.0;
	double length = 0.5;
	unsigned i = 1;
	unsigned j = 1;
	
	double pi = 3.1415926535;
	double x_limit = length / math::sqrt(length * length + radius * radius);

	// Generate an upper sphere.
	TriangularMesh<double> mesh;
	mesh.generate_uniform_hemisphere(1.0, 1000);
	mesh.translate({0.0, 0.0, length});
	
	// Integrate upper sphere only.
	double upper_num = surface_integral<double>(mesh, I_integrand<double>(i, j));
	double upper_x = integrate<double>(x_limit, 1.0, Iij_suspended_hemispherical_integrand<double>(radius, length, i, j));
	double full_apprx = integrate<double>(x_limit, 1.0, Iij_suspended_hemispherical_integrand_approx<double>(radius, length, i, j));

	// Integrate lower sphere.
	mesh.z_scale(-1.0);
	double lower_num = surface_integral<double>(mesh, I_integrand<double>(i, j));
	double lower_x = integrate<double>(-1.0, -x_limit, Iij_suspended_hemispherical_integrand<double>(radius, -length, i, j));

	// Display results.
	std::cout << "Upper Mesh:  " << upper_num << std::endl;	
	std::cout << "Upper Integ: " << upper_x << std::endl;

	std::cout << "Lower Mesh:  " << lower_num << std::endl;	
	std::cout << "Lower Integ: " << lower_x << std::endl;

	std::cout << "Total Mesh:  " << upper_num + lower_num << std::endl;
	std::cout << "Total Integ: " << upper_x + lower_x << std::endl;
	std::cout << "Total Apprx: " << full_apprx << std::endl;

	double i00term1 = 4*pi * (1 - x_limit);
	double i00term3 = -8*pi / 3.0 / radius / radius * length * length * (1 - math::pow(x_limit, 3));
	double i00term4 = 8*pi / radius / radius * length * length * (1 - x_limit);
	double I00 = i00term1 + i00term3 + i00term4;
	std::cout << "I00: " << I00 << std::endl;

	double i11term1 = 4.0 / 3.0 * pi / radius / radius * (1 - math::pow(x_limit, 3));
	double i11term2 = -2.0 * pi / radius / radius / radius * length * (1 - math::pow(x_limit, 4));
	double i11term3 = 8.0 / 5.0 * pi / math::pow(radius, 4) * length * length * (1 - math::pow(x_limit, 5));
	double i11term4 = 16.0 / 3.0 * pi / math::pow(radius, 4) * length * length * (1 - math::pow(x_limit, 3));
	double I11 = i11term1 + i11term2 + i11term3 + i11term4;
	std::cout << "I11: " << I11 << std::endl;
}

void checking_Gl() {
	double radius = 1.0;
	double length = 0.1;
	unsigned l = 1;
	
	double pi = 3.1415926535;
	double theta_zero = std::atan(radius / length);
	double x_limit = length / math::sqrt(length * length + radius * radius);

	// Generate an upper sphere.
	TriangularMesh<double> mesh;
	mesh.generate_uniform_hemisphere(1.0, 1000);
	mesh.translate({0.0, 0.0, length});
	
	// Integrate upper sphere only.
	double upper_num = surface_integral<double>(mesh, G_integrand<double>(l));
	double upper_theta = integrate<double>(0, theta_zero, Gl_suspended_hemispherical_integrand_theta<double>(radius, length, l));
	double upper_x = integrate<double>(x_limit, 1.0, Gl_suspended_hemispherical_integrand<double>(radius, length, l));
	double full_apprx = integrate<double>(x_limit, 1.0, Gl_suspended_hemispherical_integrand_approx<double>(radius, length, l));

	// Integrate lower sphere.
	mesh.z_scale(-1.0);
	double lower_num = surface_integral<double>(mesh, G_integrand<double>(l));
	double lower_theta = integrate<double>(pi-theta_zero, pi, Gl_suspended_hemispherical_integrand_theta<double>(radius, length, l));
	double lower_x = integrate<double>(-1.0, -x_limit, Gl_suspended_hemispherical_integrand<double>(radius, -length, l));

	// Display results.
	std::cout << "Upper Mesh:  " << upper_num << std::endl;	
	std::cout << "Upper Theta: " << upper_theta << std::endl;
	std::cout << "Upper Integ: " << upper_x << std::endl;

	std::cout << "Lower Mesh:  " << lower_num << std::endl;	
	std::cout << "Lower Theta: " << lower_theta << std::endl;
	std::cout << "Lower Integ: " << lower_x << std::endl;

	std::cout << "Total Mesh:  " << upper_num + lower_num << std::endl;
	std::cout << "Total Integ: " << upper_x + lower_x << std::endl;
	std::cout << "Total Apprx: " << full_apprx << std::endl;
	
	// G1
	double ratio = x_limit;
	double g1term1 = 4.0 / 3.0 * pi * radius * (1.0 - math::pow(x_limit, 3));
	double g1term2 = pi * length * (1.0 - math::pow(x_limit, 4));
	double g1term3 = -4.0 / 5.0 * pi / radius * length * length * (1.0 - math::pow(x_limit, 5));
	double g1term4 = 4.0 / 3.0 * pi / radius * length * length * (1.0 - math::pow(x_limit, 3));
	double g1 = g1term1 + g1term2 + g1term3 + g1term4;
	std::cout << "G1: " << g1 << std::endl;
}
