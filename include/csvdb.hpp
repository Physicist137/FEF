#pragma once
#include <string>
#include <fstream>
#include <vector>
#include <db.hpp>

namespace CSVDB {

template <typename T>
Data<T> get_from_database(const std::string& db, const std::string& shape, float radius, float height, const std::string& function, unsigned i, unsigned j) {
	std::ifstream file(db, std::ios::in);
	std::string value_type = std::string(typeid(T).name());

	// One by one tuple search.
	std::string content;
	while (std::getline(file, content)) {
		std::vector<std::string> tuple;
		std::string value;
		for (char c : content) {
			if (c != ',') value.push_back(c);
			else {
				tuple.push_back(value);
				value = "";
			}
		}
		tuple.push_back(value);

		// Display.
/*      std::cout << (tuple[0] == shape) << ") " << tuple[0] << " == " << shape << std::endl;
        std::cout << (std::stof(tuple[1]) == radius) << ") " << T(tuple[1]) << " == " << radius << std::endl;
        std::cout << (std::stof(tuple[2]) == height) << ") " << T(tuple[2]) << " == " << height << std::endl;
        std::cout << (tuple[3] == function) << ") " << tuple[3] << " == " <<  function << std::endl;
        std::cout << (tuple[4] == std::to_string(i)) << ") " << tuple[4] << " == " << std::to_string(i) << std::endl;
        std::cout << (tuple[5] == std::to_string(j)) << ") " << tuple[5] << " == " <<  std::to_string(j) << std::endl;
        std::cout << (tuple[6] == value_type) << ") " << tuple[6] << value_type << std::endl;
*/
		// Tuple has been loaded. Search.
		if (tuple[0] != shape) continue;
		if (std::stof(tuple[1]) != radius) continue;
		if (std::stof(tuple[2]) != height) continue;
		if (tuple[3] != function) continue;
		if (tuple[4] != std::to_string(i)) continue;
		if (tuple[5] != std::to_string(j)) continue;
		if (tuple[6] != value_type) continue;
		
		// Item has been found.
		file.close();
		return Data<T>(T(tuple[7]), true);
	}

	// Item was not not found.
	return Data<T>();
}

template <typename T>
bool save_at_database(const std::string& db, const std::string& shape, float radius, float height, const std::string& function, unsigned i, unsigned j, const T& result) {
	std::ofstream file(db, std::ios::app);
	if (file.is_open() == false) return false;
	
	file << shape << ",";
	file << radius << ",";
	file << height << ",";
	file << function << ",";
	file << i << ",";
	file << j << ",";
	file << typeid(result).name() << ",";
	file << std::setprecision(std::numeric_limits<T>::max_digits10) << result;
	file << std::endl;
	
	file.close();
	return true;
}


} // CSVDB Namespace.
