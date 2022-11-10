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

int main() {
	constexpr double radius = 1.0;
	constexpr double height = 2.0;
	constexpr double scale = height / radius;
	constexpr double theory = 5.76156;

	const std::string filename = "report_orderh2r3.txt";
	const std::string dat1name = "report_order1h2r3.dat";
	const std::string dat2name = "report_order2h2r3.dat";

	// Generate the mesh.
	TriangularMesh<double> mesh;
	mesh.generate_uniform_sphere(radius, 1000.0);
	mesh.z_scale(scale);

	// Prepare files to receive data.
	std::ofstream file(filename);
	std::ofstream dat1(dat1name);
	std::ofstream dat2(dat2name);
	file << "Shape: Prolate Spheroid" << std::endl;
	file << "Radius:     " << radius << std::endl;
	file << "Height:     " << height << std::endl;
	file << "Faces:      " << mesh.amount_faces() << std::endl;
	file << "Faces/Area: " << mesh.faces_per_area() << std::endl;
	file << "Theory FEF: " << theory << std::endl;
	file << "---------------------------------------------" << std::endl;
	file << std::endl;
	
	file.close();
	dat1.close();
	dat2.close();

	// Begin calcualtions.
	unsigned order = 2;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> I_matrix = calculate_I_matrix<double>(mesh, order);
	Eigen::Matrix<double, Eigen::Dynamic, 1> G_vector = calculate_G_vector<double>(mesh, order);
	Eigen::Matrix<double, Eigen::Dynamic, 1> A_vector = I_matrix.colPivHouseholderQr().solve(-G_vector);
	double error1 = calculate_error(I_matrix, G_vector, A_vector);
	double error2 = calculate_error(mesh, A_vector);
	double FEF = get_FEF(A_vector, height);
	double error = std::abs(FEF - theory) / theory;
	while (error > 0.01 * 0.01) {
		++order;

		I_matrix.conservativeResize(order, order);
		G_vector.conservativeResize(order);

		for (unsigned i = 0; i < order; ++i) {
			unsigned j = order-1;
			I_matrix(i, j) = surface_integral(mesh, I_integrand<double>(i, j));
			I_matrix(j, i) = surface_integral(mesh, I_integrand<double>(j, i));
			G_vector(j) = surface_integral(mesh, G_integrand<double>(j));
		}	
		
		A_vector = I_matrix.colPivHouseholderQr().solve(-G_vector);
		error1 = calculate_error(I_matrix, G_vector, A_vector);
		error2 = calculate_error(mesh, A_vector);
		FEF = get_FEF(A_vector, height);
		error = std::abs(FEF - theory) / theory;

		std::ofstream file(filename, std::ios_base::app);
		std::ofstream dat1(dat1name, std::ios_base::app);
		std::ofstream dat2(dat2name, std::ios_base::app);

		std::string sep = "     ";
		dat1 << order << sep << FEF << std::endl;
		dat2 << order << sep << FEF << sep << error << std::endl;

		file << "Order: " << order << std::endl;
		file << "FEF:   " << FEF << std::endl;
		file << "Theoretical Error:     " << error << std::endl;
		file << "Matrix Relative Error: " << error1 << std::endl;
		file << "Absolute Error:        " << error2 << std::endl;
		file << std::endl;

		dat1.close();
		dat2.close();
		file.close();

		std::cout << order << ") " << error << std::endl;
	}
}



