#pragma once
#include <Eigen/Dense>
#include <polynomial.hpp>
#include <math.hpp>

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
calculate_I_matrix(const TriangularMesh<T>& surface, unsigned order) {
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(order, order);
	for (unsigned i = 0; i < order; ++i)
		for (unsigned j = 0; j < order; ++j)
			result(i, j) = surface_integral(surface, I_integrand<T>(i, j));
		
	return result;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
calculate_G_vector(const TriangularMesh<T>& surface, unsigned order) {
	Eigen::Matrix<T, Eigen::Dynamic, 1> result(order);
	for (unsigned l = 0; l < order; ++l)
		result(l, 0) = surface_integral(surface, G_integrand<T>(l));
		
	return result;
}

template <typename T>
T calculate_potential(const Eigen::Matrix<T, Eigen::Dynamic, 1>& A_vector, const T& radius, const T& cossine) {
	unsigned order = A_vector.rows();
	std::vector<math::Polynomial<T>> base = math::Polynomial<T>::legendre_base(order);

	T result = radius * cossine;
	for (unsigned l = 0; l < order; ++l) {
		int p = static_cast<int>(l);
		result += A_vector(l) * math::pow(radius, -p) * base[l].at(cossine);
	}

	return result;
}

template <typename T>
T calculate_error(const TriangularMesh<T>& surface, const Eigen::Matrix<T, Eigen::Dynamic, 1>& A_vector) {
	T result = T();
	unsigned size = surface.amount_faces();
	for (unsigned i = 0; i < size; ++i) {
		TriangularFace<T> face = surface.face(i);
		math::Vector<T, 3> central_point = face.center();
		T r = central_point.length();
		T cos = central_point.z() / r;
		T value = calculate_potential(A_vector, r, cos);
		result += value * value;
	}

	return result;
}

template <typename T>
T calculate_error(
const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& I_matrix,
const Eigen::Matrix<T, Eigen::Dynamic, 1> G_vector,
const Eigen::Matrix<T, Eigen::Dynamic, 1> A_vector) {
	T linear = 2.0 * A_vector.dot(G_vector); 
	T quadratic = (A_vector.transpose() * I_matrix * A_vector)(0, 0);
	return linear + quadratic;
}

template <typename T>
std::vector<T> get_FEF_contributions(const Eigen::Matrix<T, Eigen::Dynamic, 1>& A_vector, const T& radius) {
	unsigned order = A_vector.rows();
	std::vector<T> result;
	result.push_back(1.0);
	for (unsigned l = 0; l < order; ++l) {
		T value = (l+1) / math::pow(radius, l+2) * A_vector(l, 0);
		result.push_back(value);
	}

	return result;
}

template <typename T>
T get_FEF(const Eigen::Matrix<T, Eigen::Dynamic, 1>& A_vector, const T& radius) {
	unsigned order = A_vector.rows();
	T result = 1.0;
	for (unsigned l = 0; l < order; ++l) {
		T value = (l+1) / math::pow(radius, l+2) * A_vector(l, 0);
		result += value;
	}

	return result;
}



