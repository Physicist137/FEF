#pragma once
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

namespace math {

// Sqrt -----------------------------
template <typename T>
T sqrt(const T& number) {
	return boost::multiprecision::sqrt(number);
}

template<>
double sqrt<double>(const double& number) {
	return std::sqrt(number);
}

template<>
float sqrt<float>(const float& number) {
	return std::sqrt(number);
}

// Pow --------------------------
template <typename T>
T pow(const T& number, int p) {
	return boost::multiprecision::pow(number, p);
}

template<>
double pow<double>(const double& number, int p) {
	return std::pow(number, p);
}

template<>
float pow<float>(const float& number, int p) {
	return std::pow(number, p);
}

// log ---------------------------
template <typename t>
t log(const t& number) {
	return boost::multiprecision::log(number);
}

template<>
double log<double>(const double& number) {
	return std::log(number);
}

template<>
float log<float>(const float& number) {
	return std::log(number);
}

// Sin ---------------------------
template <typename t>
t sin(const t& number) {
	return boost::multiprecision::sin(number);
}

template<>
double sin<double>(const double& number) {
	return std::sin(number);
}

template<>
float sin<float>(const float& number) {
	return std::sin(number);
}


// Cos ---------------------------
template <typename t>
t cos(const t& number) {
	return boost::multiprecision::cos(number);
}

template<>
double cos<double>(const double& number) {
	return std::cos(number);
}

template<>
float cos<float>(const float& number) {
	return std::cos(number);
}
}
