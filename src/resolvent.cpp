#include <mesh.hpp>
#include <integral.hpp>
#include <functional>
#include <Eigen/Dense>

template <typename T>
T kernel(const math::Vector<T, 3>& x) {
	const T pi = 3.1415926535;
	return 1.0 / (4.0 * pi) / x.length();
}

template <typename T>
T kernel_directional_derivative_numerical(const math::Vector<T, 3>& x, const math::Vector<T, 3>& normal) {
	T dt = 1e-5;
	T dk = kernel(x + normal * dt) - kernel(x);
	return dk / dt;
}

template <typename T>
T kernel_directional_derivative(const math::Vector<T, 3>& x, const math::Vector<T, 3>& normal) {
	const T pi = 3.1415926535;
	T len = x.length();
	math::Vector<T, 3> gradient = -x / (4.0 * pi) / (len * len * len);
	return gradient.dot(normal);
}

template <typename T>
T resolvent_potential(const TriangularMesh<T>& surface, const math::Vector<T, 3>& point) {
	T result = T();
	unsigned size = surface.amount_faces();
	for (unsigned i = 0; i < size; ++i) {
		TriangularFace<T> face = surface.face(i);
		math::Vector<T, 3> central_point = face.center();
		T cte = 2.0;
		T value = kernel_directional_derivative(point - central_point, face.normal());
		result += cte * central_point.z() * value * face.area();
	}

	for (unsigned i = 0; i < size; ++i) {
		for (unsigned j = 0; j < size; ++j) {
			if (i == j) continue;
			TriangularFace<T> face1 = surface.face(i);
			TriangularFace<T> face2 = surface.face(j);

			math::Vector<T, 3> central_point1 = face1.center();
			math::Vector<T, 3> central_point2 = face2.center();

			T cte = 2.0 * (-2.0);

			T value1 = kernel_directional_derivative(point - central_point1, face1.normal());
			T value2 = kernel_directional_derivative(central_point1 - central_point2, face2.normal());

			result += cte * central_point2.z() * value1 * value2 * face1.area() * face2.area();
		}
	}

	for (unsigned i = 0; i < size; ++i) {
		for (unsigned j = 0; j < size; ++j) {
			if (i == j) continue;
			for (unsigned k = 0; k < size; ++k) {
				if (k == j or k == i) continue;
				TriangularFace<T> face1 = surface.face(i);
				TriangularFace<T> face2 = surface.face(j);
				TriangularFace<T> face3 = surface.face(k);

				math::Vector<T, 3> central_point1 = face1.center();
				math::Vector<T, 3> central_point2 = face2.center();
				math::Vector<T, 3> central_point3 = face3.center();

				T cte = 2.0 * (-2.0) * (-2.0);

				T value1 = kernel_directional_derivative(point - central_point1, face1.normal());
				T value2 = kernel_directional_derivative(central_point1 - central_point2, face2.normal());
				T value3 = kernel_directional_derivative(central_point2 - central_point3, face3.normal());

				T area = face1.area() * face2.area() * face3.area();
				T value = value1 * value2 * value3;

				result += cte * central_point3.z() * value * area;
			}
		}
	}

	for (unsigned i = 0; i < size; ++i) {
		for (unsigned j = 0; j < size; ++j) {
			if (i == j) continue;
			for (unsigned k = 0; k < size; ++k) {
				if (k == j or k == i) continue;
				for (unsigned l = 0; l < size; ++l) {
					TriangularFace<T> face1 = surface.face(i);
					TriangularFace<T> face2 = surface.face(j);
					TriangularFace<T> face3 = surface.face(k);
					TriangularFace<T> face4 = surface.face(l);

					math::Vector<T, 3> central_point1 = face1.center();
					math::Vector<T, 3> central_point2 = face2.center();
					math::Vector<T, 3> central_point3 = face3.center();
					math::Vector<T, 3> central_point4 = face4.center();

					T cte = 2.0 * (-2.0) * (-2.0) * (-2.0);

					T value1 = kernel_directional_derivative(point - central_point1, face1.normal());
					T value2 = kernel_directional_derivative(central_point1 - central_point2, face2.normal());
					T value3 = kernel_directional_derivative(central_point2 - central_point3, face3.normal());
					T value4 = kernel_directional_derivative(central_point3 - central_point4, face4.normal());

					T area = face1.area() * face2.area() * face3.area() * face4.area();
					T value = value1 * value2 * value3 * value4;

					result += cte * central_point3.z() * value * area;
				}
			}
		}
	}

	std::cout << "                   \r";
	return result;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> generate_kernel_matrix(const TriangularMesh<T>& mesh) {
	// Define few constants.
	T pi = 3.1415926535;
	T cte = 1.0 / (4.0 * pi);
	
	// Generate the (assymetric) kernel matrix.
	unsigned faces = mesh.amount_faces();
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> kernel(faces, faces);
	for (unsigned i = 0; i < faces; ++i) {
		for (unsigned j = 0; j < faces; ++j) {
			if (i == j) {
				kernel(i, j) = 0;
				continue;	
			}

			math::Vector<T, 3> ri = mesh.face(i).center();
			math::Vector<T, 3> rj = mesh.face(j).center();

			kernel(i, j) = cte * kernel_directional_derivative(ri - rj, mesh.face(j).normal()) * mesh.face(j).area();
		}
	}

	return kernel;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> generate_potential_vector(const TriangularMesh<T>& mesh) {
	unsigned faces = mesh.amount_faces();
	Eigen::Matrix<T, Eigen::Dynamic, 1> gvec(faces);
	for (unsigned i = 0; i < faces; ++i) {
		gvec(i) = -mesh.face(i).center().z();
	}

	return gvec;
}

template <typename T>
T potential(const TriangularMesh<T>& mesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& hvector, const math::Vector<T, 3>& point) {
	T result = T();

	unsigned faces = mesh.amount_faces();
	for (unsigned i = 0; i < faces; ++i) {
		math::Vector<T, 3> pos = mesh.face(i).center();
		result -= hvector(i) * kernel_directional_derivative(point - pos, mesh.face(i).normal()) * mesh.face(i).area();
	}

	return result;
}

template <typename T>
T potential_spherical_theoretical(const math::Vector<T, 3>& point, const T& radius = 1.0) {
	T r = point.length();
	T z = point.z();
	return -(1.0 * 0.0 - radius * radius * radius / (r * r * r)) * z;
}
/*
int main() {
	TriangularMesh<double> mesh;
	mesh.generate_uniform_sphere(1.0, 100);
	std::cout << mesh.amount_faces() << std::endl;

	auto kernel = generate_kernel_matrix<double>(mesh);
	auto gvec = generate_potential_vector<double>(mesh);
	auto identity = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Identity(mesh.amount_faces(), mesh.amount_faces());

	auto matrix = -(0.5 * identity + kernel);
	auto inverse = matrix.inverse();
	auto hvec = inverse * gvec;

	math::Vector<double, 3> point = {1.0, 1.0, 1.0};
	std::cout << potential<double>(mesh, hvec, point) << std::endl;
	std::cout << potential_spherical_theoretical<double>(point) << std::endl;
}*/

/*
int main() {
	TriangularMesh<double> mesh;
	mesh.generate_uniform_sphere(1.0, 100);

	unsigned f = 40;	
	std::cout << kernel_directional_derivative_numerical(mesh.face(f).center(), mesh.face(f).normal()) << std::endl;
	std::cout << kernel_directional_derivative(mesh.face(f).center(), mesh.face(f).normal()) << std::endl;
}
*/

int main() {
	TriangularMesh<double> mesh;
	mesh.generate_uniform_sphere(1.0, 100);
	double dx = 1e-5;
	double value1 = resolvent_potential(mesh, {0, 0, 1.01});
	double value2 = resolvent_potential(mesh, {0, 0, 1.01 + dx});
	std::cout << (value2 - value1) / dx << std::endl;
}
