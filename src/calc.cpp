#include <integral.hpp>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <sstream>

/*
AR   FEFAd      FEFAn
 2    6.25494    5.76156
 3    8.75592    9.19883
 4    10.8197    13.26133
 5    12.578      17.91441
 6    14.1039    23.13227
 7    15.4476    28.89458
 8    16.6456    35.18479
 9    17.7249    41.98901
11  19.6048    57.09354
31  29.3076    306.80503
*/

double prolate_spheroid(double radius, double height, unsigned order, std::ostream& os, std::ostream& dat) {
	double scale = height / radius;
	
	TriangularMesh<double> mesh;
	mesh.generate_uniform_sphere(radius, 1000.0);
	mesh.z_scale(scale);

	auto I_matrix = calculate_I_matrix<double>(mesh, order);
	auto G_vector = calculate_G_vector<double>(mesh, order);
	Eigen::Matrix<double, Eigen::Dynamic, 1> A_vector = I_matrix.colPivHouseholderQr().solve(-G_vector);
	double error1 = calculate_error(I_matrix, G_vector, A_vector);
	double error2 = calculate_error(mesh, A_vector);
	double FEF = get_FEF(A_vector, height);

	os << "Order: " << order << std::endl;
	os << "Radius: " << radius << std::endl;
	os << "Height: " << height << std::endl;
	os << "Error 01: " << error1 << std::endl;
	os << "Error 02: " << error2 << std::endl;
	os << "FEF:      " << FEF << std::endl;
	os << std::endl;
	dat << order << "      " << FEF << std::endl;

	return FEF;
}

int main() {
	std::ofstream file("report_test_order.txt");
	std::ofstream dat("report_order.dat");

	unsigned i = 2;
	double theory = 5.76156;
	double error;
	do {
		double FEF = prolate_spheroid(1.0, 2.0, ++i, file, dat);
		error = std::abs(FEF - theory) / theory;
		std::cout << error << std::endl;
	} while (error > 0.01);

	prolate_spheroid(1.0, 2.0, i, file, dat);
	file.close();
}


/*
int main() {
	constexpr unsigned order = 9;
	std::vector<unsigned> scales = {
		1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 31
	};

	std::stringstream ss;
	ss << "order" << order << ".dat";
	std::ofstream file(ss.str());
	for (unsigned i : scales) {
		double scale = static_cast<double>(i);	

		TriangularMesh<double> mesh;
		mesh.generate_uniform_sphere(1.0, 1000.0);
		mesh.z_scale(scale);

		auto I_matrix = calculate_I_matrix<double>(mesh, order);
		auto G_vector = calculate_G_vector<double>(mesh, order);
		Eigen::Matrix<double, Eigen::Dynamic, 1> A_vector = I_matrix.colPivHouseholderQr().solve(-G_vector);
		double FEF = get_FEF(A_vector, scale);

		std::cout << "Radius: 1. Height: " << i << ". ";
		std::cout << "FEF: " << FEF << std::endl;

		std::string sep = "     ";
		file << i << sep << FEF << std::endl;
	}

	file.close();
}*/
