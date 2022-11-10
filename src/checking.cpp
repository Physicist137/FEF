#include <integral.hpp>
#include <integrand.hpp>
#include <polynomial.hpp>
#include <calcs.hpp>
#include <cte.hpp>

#include <boost/multiprecision/float128.hpp> 
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/gmp.hpp>

void check_area_element() {
	// using Float = double;
	using Float = boost::multiprecision::float128;	

	Float radius = 1.0;
	Float height = 2.0;
	Float ecc = boost::multiprecision::sqrt(1 - radius*radius / height / height);	

	// The elipse.
	std::function<Float(Float)> r = [=](Float t){
		Float value = ecc * boost::multiprecision::cos(t);
		return 1.0 / boost::multiprecision::sqrt(1 - value * value);
	};

	// The derivative.
	std::function<Float(Float)> drdt = [=](Float t) {
		Float rv = r(t);
		Float value = ecc / radius;
		Float cos = boost::multiprecision::cos(t);
		Float sin = boost::multiprecision::sin(t):
		return -rv*rv*rv * value*value * cos * sin;
	};

	// Numerical J
	std::function<Float(Float, Float)> Jn = [=](Float theta, Float phi){
		std::function<Float(Float, Float)> x = [=](Float theta, Float phi) {
			Float cos = boost::multiprecision::cos(phi);
			Float sin = boost::multiprecision::sin(phi):
			return r(theta) * cos * sin;
		};
		std::function<Float(Float, Float)> y = [=](Float theta, Float phi) {
			return r(theta) * std::sin(phi) * std::sin(theta);
		};
		std::function<Float(Float, Float)> z = [=](Float theta, Float phi) {
			return r(theta) * std::cos(theta);
		};

		math::Vector<Float, 3> azimuthal = {
			differentiate_first<Float>(theta, phi, x),
			differentiate_first<Float>(theta, phi, y),
			differentiate_first<Float>(theta, phi, z),
		};

		math::Vector<Float, 3> polar = {
			differentiate_second<Float>(theta, phi, x),
			differentiate_second<Float>(theta, phi, y),
			differentiate_second<Float>(theta, phi, z),
		};
		
		math::Vector<Float, 3> normal = azimuthal.cross(polar);
		return normal.length();
	};

	// Analytical J
	std::function<Float(Float, Float)> Ja = [=](Float theta, Float phi) {
		Float rv = r(theta);
		Float value = drdt(theta) / rv;
		return rv * rv * std::sin(theta) * boost::multiprecision::sqrt(1 + value * value);
	};

	std::cout << Jn(1.0, 0.8) << std::endl;
	std::cout << Ja(1.0, 0.8) << std::endl;	

}

void GL_comparer() {
	using Float = boost::multiprecision::float128;

	constexpr unsigned l = 1;

	constexpr Float pi = 3.1415926535;
	constexpr Float radius = 1.0;
	constexpr Float height = 2.0;
	constexpr Float scale = height / radius;
	constexpr Float ecc = boost::multiprecision::sqrt(1 - radius*radius / height/height);
	
	// The elipse (spherical).
	std::function<Float(Float)> r_elipse = [=](Float t){
		Float value = ecc * std::cos(t);
		return 1.0 / boost::multiprecision::sqrt(1 - value * value);
	};

	// The mesh.
	TriangularMesh<Float> mesh;
	mesh.generate_uniform_sphere(radius, 50000.0);
	mesh.z_scale(scale);

	// Computations.
	Float numerical_area = surface_integral<Float>(mesh, [=](Float r, Float cos){return 1.0;});
	Float numerical_result = surface_integral(mesh, G_integrand<Float>(l));

	Float semi1 = integrate<Float>(
		-1.0, 1.0, Gl_spheroidal_integrand<Float>(radius, height, l)
	);

	Float semi2 = integrate<Float>(
		0.0, pi, Gl_spheroidal_integrand_theta<Float>(radius, height, r_elipse, l)
	);

	Float semi3 = Gl_spheroidal_summation_exact<Float>(radius, height, l);
	Float semi4 = Gl_spheroidal_summation_approx<Float>(radius, height, l);

	std::cout << "Numerical Result: " << numerical_result << std::endl;
	std::cout << "Semi-Analytical Result: " << semi1 << std::endl;
	std::cout << "Semi-Analytical Result: " << semi2 << std::endl;
	std::cout << "Semi-Analytical Result: " << semi3 << std::endl;
	std::cout << "Semi-Analytical Result: " << semi4 << std::endl;
	std::cout << "Ratio: " << numerical_result / semi1 << std::endl;
	std::cout << "Ratio: " << numerical_result / semi2 << std::endl;
	std::cout << "Ratio: " << numerical_result / semi3 << std::endl;
	std::cout << "Ratio: " << numerical_result / semi4 << std::endl;
}

void IIJ_comparator() {
	using Float = boost::multiprecision::float128;

	constexpr unsigned i = 1;
	constexpr unsigned j = 1;

	constexpr Float pi = 3.1415926535;
	constexpr Float radius = 1.0;
	constexpr Float height = 2.0;
	constexpr Float scale = height / radius;
	constexpr Float ecc = boost::multiprecision::sqrt(1 - radius*radius / height/height);
	
	// The elipse (spherical).
	std::function<Float(Float)> r_elipse = [=](Float t){
		Float value = ecc * std::cos(t);
		return 1.0 / boost::multiprecision::sqrt(1 - value * value);
	};

	// The mesh.
	TriangularMesh<Float> mesh;
	mesh.generate_uniform_sphere(radius, 50000.0);
	mesh.z_scale(scale);

	// Computations.
	Float numerical_area = surface_integral<Float>(mesh, [=](Float r, Float cos){return 1.0;});
	Float numerical_result = surface_integral(mesh, I_integrand<Float>(i, j));

	Float semi1 = integrate<Float>(
		-1.0, 1.0, Iij_spheroidal_integrand<Float>(radius, height, i, j)
	);

	Float semi2 = integrate<Float>(
		0.0, pi, Iij_spheroidal_integrand_theta<Float>(radius, height, r_elipse, i, j)
	);

	Float semi3 = Iij_spheroidal_summation_exact<Float>(radius, height, i, j);
	Float semi4 = Iij_spheroidal_summation_approx<Float>(radius, height, i, j);

	std::cout << "Numerical Result: " << numerical_result << std::endl;
	std::cout << "Semi-Analytical Result: " << semi1 << std::endl;
	std::cout << "Semi-Analytical Result: " << semi2 << std::endl;
	std::cout << "Semi-Analytical Result: " << semi3 << std::endl;
	std::cout << "Semi-Analytical Result: " << semi4 << std::endl;
	std::cout << "Ratio: " << numerical_result / semi1 << std::endl;
	std::cout << "Ratio: " << numerical_result / semi2 << std::endl;
	std::cout << "Ratio: " << numerical_result / semi3 << std::endl;
	std::cout << "Ratio: " << numerical_result / semi4 << std::endl;
