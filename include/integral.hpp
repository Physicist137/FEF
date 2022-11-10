#pragma once
#include <functional>
#include <cmath>

#include <vector.hpp>
#include <mesh.hpp>

template <typename T>
T integrate(const T& inf, const T& sup, const std::function<T(const T&)>& integrand, const T& div = 1e-6) {
	// https://en.wikipedia.org/wiki/Trapezoidal_rule
	T result = T();

	// Find lower and upper integrals.
	T lower, upper, sgn;
	if (sup >= inf) {
		lower = inf;
		upper = sup;
		sgn = 1.0;
	}

	else {
		lower = sup;
		upper = inf;
		sgn = -1.0;
	}

	// Solve by trapezoidal rule.
	// Do not include the limits themselves (they could be assymptoptic).
	T dx = div * (upper - lower);
	T x = lower+dx;
	while (x < upper) {
		result += (integrand(x + dx) + integrand(x)) * dx / 2.0;
		x += dx;
	}

	return sgn * result;
}

template <typename T>
T differentiate(const T& x, const std::function<T(const T&)>& function) {
	T div = 1e-6;
	T dx = div;
	T df = (function(x + dx) - function(x - dx)) / 2.0;
	return df / dx;
}

template <typename T>
T differentiate_first(const T& a, const T& b, const std::function<T(const T&, const T&)>& function) {
	T div = 1e-6;
	T da = div;
	T df = (function(a + da, b) - function(a - da, b)) / 2.0;
	return df / da;
}

template <typename T>
T differentiate_second(const T& a, const T& b, const std::function<T(const T&, const T&)>& function) {
	T div = 1e-6;
	T db = div;
	T df = (function(a, b+db) - function(a, b-db)) / 2.0;
	return df / db;
}


template <typename T>
T surface_integral(const TriangularMesh<T>& surface, const std::function<T(const T&, const T&)>& integrand) {
	T result = T();
	unsigned size = surface.amount_faces();
	for (unsigned i = 0; i < size; ++i) {
		TriangularFace<T> face = surface.face(i);
		math::Vector<T, 3> central_point = face.center();
		T r = central_point.length();
		T cos = central_point.z() / r;
		T value = integrand(r, cos);
		result += value * face.area();
	}

	return result;
}


