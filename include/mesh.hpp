#pragma once
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <map>
#include <initializer_list>
#include <cmath>
#include <vector.hpp>
#include <iostream>
#include <algorithm>

template <typename T>
class TriangularFace {
	std::array<math::Vector<T, 3>, 3> _vertices;

public:
	TriangularFace(std::initializer_list<math::Vector<T, 3>> list) : _vertices(list) {}

	TriangularFace(
		const math::Vector<T, 3>& v1,
		const math::Vector<T, 3>& v2,
		const math::Vector<T, 3>& v3
	) : _vertices({v1, v2, v3}) {}

	inline const math::Vector<T, 3>& operator[](unsigned i) const {return _vertices[i];}
	inline const math::Vector<T, 3>& vertex(unsigned i) const {return _vertices[i];}
	inline math::Vector<T, 3> edge(unsigned i, unsigned j) const {return _vertices[j] - _vertices[i];}
	inline const math::Vector<T, 3>& edge(unsigned i) const {
		switch (i) {
			case 0: return _vertices[1] - _vertices[0];
			case 1: return _vertices[2] - _vertices[0];
			case 2: return _vertices[2] - _vertices[1];
		}	
	}
	
	T area() const;
	math::Vector<T, 3> normal() const;	
	math::Vector<T, 3> area_normal() const;	
	math::Vector<T, 3> center() const;
};


template <typename T>
class TriangularMesh {
	std::vector<math::Vector<T, 3>> _vertices;
	std::vector<std::array<unsigned, 3>> _faces;

public:
	TriangularMesh() : _vertices(), _faces() {}

	// Accessor functions
	inline unsigned amount_vertices() const {return _vertices.size();}
	inline unsigned amount_faces() const {return _faces.size();}
	inline unsigned lastID() const {return _vertices.size()-1;}
	inline const math::Vector<T, 3>& vertex(unsigned id) const {return _vertices[id];}
	inline const std::array<unsigned, 3>& face_from_vertices_id(unsigned id) const {return _faces[id];}
	inline TriangularFace<T> face(unsigned id) const {
		return TriangularFace<T>(
			_vertices[_faces[id][0]],
			_vertices[_faces[id][1]],
			_vertices[_faces[id][2]]
		);
	}

	// Mesh functions.
	T area() const;
	T face_area_density() const;
	inline T faces_per_area() const {return face_area_density();}

	// Creating mesh.
	inline unsigned addVertex(const math::Vector<T, 3>& vertex) {_vertices.push_back(vertex); return lastID();}
	inline unsigned addFace(const std::array<unsigned, 3>& face) {_faces.push_back(face); return lastID();}

	// Handling the mesh.
	TriangularMesh& z_scale(const T& scale);
	TriangularMesh& translate(const math::Vector<T, 3>& point);

	// Generating functions.
	TriangularMesh& generate_uniform_cylinder(const T& radius, const T& height, const T& density);
	TriangularMesh& generate_uniform_sphere(const T& radius, const T& density);
	TriangularMesh& generate_uniform_hemisphere(const T& radius, const T& density);

	// Static functions.
	static TriangularMesh generate_cylinder(const T& radius, const T& height, const T& density);
	static TriangularMesh generate_sphere(const T& radius, const T& density);

	// Export functions.
	void save_as_wavefront_object(const std::string& file) const;

	// Import functions.
	void load_from_wavefront_object(const std::string& file);
};


template <typename T>
T TriangularFace<T>::area() const {
	return area_normal().length();
}

template <typename T>
math::Vector<T, 3> TriangularFace<T>::normal() const {
	return area_normal().unit();
}

template <typename T>
math::Vector<T, 3> TriangularFace<T>::area_normal() const {
	math::Vector<T, 3> edge01 = _vertices[1] - _vertices[0];
	math::Vector<T, 3> edge02 = _vertices[2] - _vertices[0];
	return 0.5 * edge01.cross(edge02);
}

template <typename T>
math::Vector<T, 3> TriangularFace<T>::center() const {
	return _vertices[0] / 3.0 + _vertices[1] / 3.0 + _vertices[2] / 3.0;
}


template <typename T>
T TriangularMesh<T>::area() const {
	unsigned amount = amount_faces();
	
	T result = 0;
	for (unsigned i = 0; i < amount; ++i) result += face(i).area();
	return result;
}

template <typename T>
T TriangularMesh<T>::face_area_density() const {
	return amount_faces() / area();
}


template <typename T>
TriangularMesh<T>& TriangularMesh<T>::z_scale(const T& scale) {
	unsigned sizes = _vertices.size();
	for (unsigned i = 0; i < sizes; ++i) _vertices[i][2] = _vertices[i][2] * scale;
	return *this;
}


template <typename T>
TriangularMesh<T>& TriangularMesh<T>::translate(const math::Vector<T, 3>& point) {
	unsigned sizes = _vertices.size();
	for (unsigned i = 0; i < sizes; ++i) _vertices[i] += point;
	return *this;
}

// Mesh Generating functions ---------------------------------------------------------

template <typename T>
TriangularMesh<T>& TriangularMesh<T>::generate_uniform_cylinder(const T& radius, const T& height, const T& density) {
	constexpr T pi = 3.14159265358979323;
	T length = 2.0 * pi * radius;

	// The density of vertices is 1/d^2, where d is the distance between two horizontal/vertical vertices.
	T distance =  1.0 / std::sqrt(density);

	// Calculate the amount of points going towards the base.
	T famount = std::ceil(length / distance);
	long unsigned amount = static_cast<long unsigned>(famount);

	// Calculate the amount of points at the cylinder height.
	T fheight_amount = std::ceil(height / distance);
	T offset = height / 2.0;
	long unsigned height_amount = static_cast<long unsigned>(fheight_amount);

	// Generate the planar grid mapping. It uses the fact cylinders are geometrically isomorphic to a plane.
	std::map<std::tuple<long unsigned, long unsigned>, math::Vector<T, 3>> grid;
	for (long unsigned j = 0; j < amount; ++j) {
		T fj = static_cast<T>(j);
		T arg = fj * distance / radius;

		for (long unsigned i = 0; i < height_amount; ++i) {
			T fi = static_cast<T>(i);
			math::Vector<T, 3> point = {radius * std::cos(arg), radius * std::sin(arg), fi * distance - offset};
			grid[std::make_tuple(i, j)] = point;
		}
	}

	// Clear the mesh to open space for the cylinder.
	_faces.clear();
	_vertices.clear();

	// Create the mesh and a mapping for the ID of each mesh.
	std::map<std::tuple<long unsigned, long unsigned>, unsigned> gridid;
	for (long unsigned j = 0; j < amount; ++j)
		for (long unsigned i = 0; i < height_amount; ++i)
			gridid[std::make_pair(i, j)] = this->addVertex(grid[std::make_pair(i, j)]);
	

	// Assemble the cylinder itself by connecting the planar grid into faces.
	for (long unsigned i = 0; i < height_amount-1; ++i) {
		for (long unsigned j = 0; j < amount-1; ++j) {
			this->addFace({
				gridid[std::make_pair(i, j)],
				gridid[std::make_pair(i+1, j)],
				gridid[std::make_pair(i, j+1)]
			});

			this->addFace({
				gridid[std::make_pair(i+1, j)],
				gridid[std::make_pair(i+1, j+1)],
				gridid[std::make_pair(i, j+1)]
			});
		}
	}
	
	// Now proceed with final phase, connecting the 0th and Jth planar grid columns.
	long unsigned j = amount - 1;
	long unsigned jp = 0;
	for (long unsigned i = 0; i < height_amount-1; ++i) {
		this->addFace({
			gridid[std::make_pair(i, j)],
			gridid[std::make_pair(i+1, j)],
			gridid[std::make_pair(i, jp)]
		});

		this->addFace({
			gridid[std::make_pair(i+1, j)],
			gridid[std::make_pair(i+1, jp)],
			gridid[std::make_pair(i, jp)]
		});
	}

	// With all faces added, return the mesh.	
	return *this;
}

template <typename T>
TriangularMesh<T>& TriangularMesh<T>::generate_uniform_sphere(const T& radius, const T& density) {
	// TO ASK ON STACK OVERFLOW! WHY TETRAHEDRON WAS UGLY.. WHILE ICOASHEDRON IT WORKED?
	// Create regular tetrahedron (pyramid with triangular base) inside sphere of said radius.
	// https://en.wikipedia.org/wiki/Tetrahedron#Formulas_for_a_regular_tetrahedron
	/*math::Vector<T, 3> v1 = {+std::sqrt(8.0 / 9.0), 0, -1.0 / 3.0};
	math::Vector<T, 3> v2 = {-std::sqrt(2.0 / 9.0), +std::sqrt(2.0 / 3.0), -1.0 / 3.0};
	math::Vector<T, 3> v3 = {-std::sqrt(2.0 / 9.0), -std::sqrt(2.0 / 3.0), -1.0 / 3.0};
	math::Vector<T, 3> v4 = {0, 0, 1};

	unsigned id1 = this->addVertex(radius * v1);
	unsigned id2 = this->addVertex(radius * v2);
	unsigned id3 = this->addVertex(radius * v3);
	unsigned id4 = this->addVertex(radius * v4);

	this->addFace({id1, id2, id3});	
	this->addFace({id1, id2, id4});	
	this->addFace({id1, id3, id4});	
	this->addFace({id2, id3, id4});*/
	// --------------------------------------------------------------

	// Clear current mesh.
	_faces.clear();
	_vertices.clear();

	// Create a regualr icosahedron.
	// https://schneide.blog/2016/07/15/generating-an-icosphere-in-c/
	const T X=.525731112119133606 * radius;
	const T Z=.850650808352039932 * radius;
	const T N=0;

	_vertices = {
		{-X,N,Z}, {X,N,Z}, {-X,N,-Z}, {X,N,-Z},
		{N,Z,X}, {N,Z,-X}, {N,-Z,X}, {N,-Z,-X},
		{Z,X,N}, {-Z,X, N}, {Z,-X,N}, {-Z,-X, N}
	};

	_faces = {
		{0,4,1},{0,9,4},{9,5,4},{4,5,8},{4,8,1},
		{8,10,1},{8,3,10},{5,3,8},{5,2,3},{2,7,3},
		{7,10,3},{7,6,10},{7,11,6},{11,0,6},{0,1,6},
		{6,1,10},{9,0,11},{9,11,2},{9,2,5},{7,2,11}
	};

	// At each iteration, divide each face into four faces.
	// Once done, push the new created vertices towards a sphere.
	T current_density = face_area_density();
	while (current_density < density) {
		std::vector<std::array<unsigned, 3>> new_faces;
	
		unsigned amount = this->amount_faces();
		for (unsigned i = 0; i < amount; ++i) {
			TriangularFace<T> current_face = this->face(i);
			std::array<unsigned, 3> vts = this->face_from_vertices_id(i);

			math::Vector<T, 3> v1 = current_face.vertex(0) + current_face.edge(0, 1) / 2.0;
			math::Vector<T, 3> v2 = current_face.vertex(0) + current_face.edge(0, 2) / 2.0;
			math::Vector<T, 3> v3 = current_face.vertex(1) + current_face.edge(1, 2) / 2.0;

			unsigned id1 = this->addVertex(radius * v1.unit());
			unsigned id2 = this->addVertex(radius * v2.unit());
			unsigned id3 = this->addVertex(radius * v3.unit());

			new_faces.push_back({vts[0], id2, id1});
			new_faces.push_back({id1, id3, vts[1]});
			new_faces.push_back({id3, id2, vts[2]});
			new_faces.push_back({id1, id2, id3});
		}

		_faces.clear();
		_faces.swap(new_faces);
		current_density = face_area_density();
	}

	// Return the mesh.
	return *this;
}


template <typename T>
TriangularMesh<T>& TriangularMesh<T>::generate_uniform_hemisphere(const T& radius, const T& density) {
	// Generate a sphere.
	this->generate_uniform_sphere(radius, density);

	// The plan: Cut every vertex below the z=0 plane. Delete faces belonging to such vertices.
	// Flaws: It creates some teeths. To fix that.
	
	// Count vertices.
	unsigned vertices = _vertices.size();

	// Mark all those to exclude and those not to.
	std::vector<unsigned> to_exclude;
	std::vector<unsigned> not_exclude;
	for (unsigned i = 0; i < vertices; ++i) {
		if (_vertices[i].z() < 0.0) to_exclude.push_back(i);
		else not_exclude.push_back(i);
	}

	// Map current IDs with new IDs.	
	std::map<unsigned, unsigned> newids;
	unsigned size = not_exclude.size();
	for (unsigned i = 0; i < size; ++i) newids[not_exclude[i]] = i;

	// Erase vertices.	
	unsigned counter = 0;
	for (auto it = _vertices.begin(); it != _vertices.end(); /*no condition*/ ) {

		// Check if current vertice *it is to be excluded. If yes, exclude it and update iterator.
		if (std::find(to_exclude.begin(), to_exclude.end(), counter) != to_exclude.end()) {
			it = _vertices.erase(it);
		}

		// If not, just update the iterator.
		else ++it;

		// Update the counter
		++counter;
	}


	// Erase faces. Same thing as with vertices.
	counter = 0;
	for (auto it = _faces.begin(); it != _faces.end(); /* no condition */ ) {
		if (std::find(to_exclude.begin(), to_exclude.end(), (*it)[0]) != to_exclude.end()) {
			it = _faces.erase(it);
		}

		else if (std::find(to_exclude.begin(), to_exclude.end(), (*it)[1]) != to_exclude.end()) {
			it = _faces.erase(it);
		}

		else if (std::find(to_exclude.begin(), to_exclude.end(), (*it)[2]) != to_exclude.end()) {
			it = _faces.erase(it);
		}

		else ++it;

		++counter;
	}

	// Update the faces vector with the new faces.
	for (unsigned i = 0; i < _faces.size(); ++i) {
		_faces[i][0] = newids[_faces[i][0]];
		_faces[i][1] = newids[_faces[i][1]];
		_faces[i][2] = newids[_faces[i][2]];
	}

	// Return the mesh.
	return *this;
}

// Static functions ------------------------------------------------------------------
template <typename T>
TriangularMesh<T> TriangularMesh<T>::generate_sphere(const T& radius, const T& density) {
	TriangularMesh<T> mesh;
	mesh.generate_uniform_sphere(radius, density);
	return mesh;
}

template <typename T>
TriangularMesh<T> TriangularMesh<T>::generate_cylinder(const T& radius, const T& height, const T& density) {
	TriangularMesh<T> mesh;
	mesh.generate_uniform_cylinder(radius, height, density);
	return mesh;
}


// Importing/Export functions --------------------------------------------------------
template <typename T>
void TriangularMesh<T>::save_as_wavefront_object(const std::string& filename) const {
	std::ofstream file(filename);

	// file << "# Vertices of the mesh" << std::endl;
	for (const math::Vector<T, 3>& vertex : _vertices) {
		file << "v ";
		file << vertex[0] << " ";
		file << vertex[1] << " ";
		file << vertex[2] << std::endl;
	}

	file << std::endl;
	//file << "# Faces of the mesh" << std::endl;
	for (const std::array<unsigned, 3>& face : _faces) {
		file << "f ";
		file << face[0]+1 << " ";
		file << face[1]+1 << " ";
		file << face[2]+1 << std::endl;
	}

	file.close();
}

template <typename T>
void TriangularMesh<T>::load_from_wavefront_object(const std::string& filename) {
	std::ifstream file(filename);

	_vertices.clear();
	_faces.clear();
	
	struct internal {
		static std::vector<std::string> split(const std::string& str, char separator) {
			std::vector<std::string> result;
			std::string temp;
			for (char c : str) {
				if (c == separator and !temp.empty()) {
					result.push_back(temp);
					temp.clear();
				} else temp.push_back(c);
			}
	
			result.push_back(temp);	
			return result;
		}
	};

	for (std::string line; std::getline(file, line); ) {
		auto split = internal::split(line, ' ');
		if (split.empty()) continue;
		if (split[0] == "#") continue;

		if (split[0] == "v") {
			math::Vector<T, 3> vertex = {
				static_cast<T>(std::stod(split[1])),
				static_cast<T>(std::stod(split[2])),
				static_cast<T>(std::stod(split[3]))
			};

			_vertices.push_back(vertex);
		}
		
		if (split[0] == "f") {
			std::array<unsigned, 3> face = {
				static_cast<unsigned>(std::stoi(split[1]))-1,
				static_cast<unsigned>(std::stoi(split[2]))-1,
				static_cast<unsigned>(std::stoi(split[3]))-1
			};
		
			_faces.push_back(face);
		}
	}
}
