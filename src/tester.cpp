#include <polynomial.hpp>
#include <vector.hpp>
#include <mesh.hpp>
#include <iostream>
#include <string>

int main() {
	
}


/*using namespace math;
int main() {
	TriangularMesh<double> mesh;
	mesh.generate_uniform_sphere(1.0, 100);
	mesh.save_as_wavefront_object("sphere.obj");
}*/


/*
int main() {
	double height = 5.0;
	double radius = 1.0;
	TriangularMesh<double> mesh = TriangularMesh<double>::generate_cylinder(radius, height, 10.0);
	std::cout << "Theoretical Area:  " << 2.0 * 3.1415926535 * radius * height << std::endl;
	std::cout << "Mesh Area:         " << mesh.area() << std::endl;
	std::cout << "Face/Area density: " << mesh.face_area_density() << std::endl;
	mesh.save_as_wavefront_object("cylinder.obj");
	
	TriangularMesh<double> cube;
	cube.load_from_wavefront_object("cube.obj");
	std::cout << "Theoretical Area:  " << 2.0 * 1.0 + 4.0 * 2.0 * 1.0 << std::endl;
	std::cout << "Mesh Area:         " << cube.area() << std::endl;
	std::cout << "Face/Area density: " << cube.face_area_density() << std::endl;

	TriangularFace<double> face = cube.face(0);
	std::cout << face[0] << std::endl;
	std::cout << face[1] << std::endl;
	std::cout << face[2] << std::endl;
	std::cout << face.area() << std::endl;

	TriangularMesh<double> square;
	unsigned id1 = square.addVertex({0.0, 0.0, 0.0});
	unsigned id2 = square.addVertex({1.0, 0.0, 0.0});
	unsigned id3 = square.addVertex({0.0, 1.0, 0.0});
	unsigned id4 = square.addVertex({1.0, 1.0, 0.0});
	
	square.addFace({id1, id2, id3});
	square.addFace({id1, id4, id3});

	std::cout << square.area() << std::endl;
}*/

/*
int main() {
	Vector3 v1 = {0.0, 0.0, 0.0};
	Vector3 v2 = {1.0, 0.0, 0.0};
	Vector3 v3 = {0.5, 0.5, 0.0};
	TriangularFace<double> face(v1, v2, v3);
	std::cout << face.area() << std::endl;

	TriangularMesh<double> mesh;
	unsigned id1 = mesh.addVertex(v1);
	unsigned id2 = mesh.addVertex(v2);
	unsigned id3 = mesh.addVertex(v3);
	mesh.addFace({id1, id2, id3});

	std::cout << mesh.face(0).area() << std::endl;
	mesh.save_as_wavefront_object("test.obj");
	mesh.load_from_wavefront_object("test.obj");

	TriangularMesh<double> cube;
	cube.load_from_wavefront_object("/home/physicist/mesh/cube.obj");
	
	for (unsigned i = 0; i < cube.amount_vertices(); ++i) std::cout << cube.vertex(i) << std::endl;
	for (unsigned i = 0; i < cube.amount_faces(); ++i) {
		std::cout << cube.face_from_vertices_id(i)[0] << " ";
		std::cout << cube.face_from_vertices_id(i)[1] << " ";
		std::cout << cube.face_from_vertices_id(i)[2] << " ";
		std::cout << "\tNormal: " << cube.face(i).area_normal() << std::endl;
	}
}*/


/*
// Tests the polynomials. Builds the legendre polynomials.
int main() {
	constexpr unsigned num = 10;
	auto base = math::Polynomial<double>::legendre_base(num);
	for (unsigned i = 0; i < num; ++i) std::cout << i << ") " << base[i] * 2.0 << std::endl;
}*/
