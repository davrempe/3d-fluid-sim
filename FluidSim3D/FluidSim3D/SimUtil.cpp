#include "SimUtil.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <sstream>

namespace SimUtil {

	
	template<typename T>
	T** initMat2D(int m, int n) {
		T** mat = new T*[m];
		for (int i = 0; i < m; i++) {
			mat[i] = new T[n];
		}

		return mat;
	}

	template<typename T>
	void deleteMat2D(int m, int n, T ** mat) {
		for (int i = 0; i < m; i++) {
			delete[] mat[i];
		}
		delete[] mat;
	}

	template <typename T>
	void printMat2D(int m, int n, T** mat) {
		for (int i = 0; i < m; i++) {
			std::cout << "[";
			for (int j = 0; j < n - 1; j++) {
				std::cout << mat[i][j] << ", ";
			}
			std::cout << mat[i][n - 1] << "]" << std::endl;
		}
	}

	template<typename T>
	T** initGrid2D(int x, int y) {
		return initMat2D<T>(x, y);
	}

	template<typename T>
	void deleteGrid2D(int x, int y, T ** grid) {
		deleteMat2D<T>(x, y, grid);
	}

	template <typename T>
	void printGrid2D(int x, int y, T** grid) {
		std::cout << "====================================================================================\n";
		for (int i = y - 1; i >= 0; i--) {
			std::cout << "[";
			for (int j = 0; j < x - 1; j++) {
				std::cout << grid[j][i] << ", ";
			}
			std::cout << grid[x - 1][i] << "]" << std::endl;
		}
		std::cout << "====================================================================================\n";
	}

	template<typename T>
	T*** initMat3D(int m, int n, int l) {
		T*** mat = new T**[m];
		for (int i = 0; i < m; i++) {
			mat[i] = new T*[n];
			for (int j = 0; j < n; j++) {
				mat[i][j] = new T[l];
			}
		}

		return mat;
	}

	template<typename T>
	void deleteMat3D(int m, int n, int l, T *** mat) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				delete[] mat[i][j];
			}
			delete[] mat[i];
		}
		delete[] mat;
	}

	template<typename T>
	T*** initGrid3D(int x, int y, int z) {
		return initMat3D<T>(x, y, z);
	}

	template<typename T>
	void deleteGrid3D(int x, int y, int z, T *** grid) {
		deleteMat3D<T>(x, y, z, grid);
	}

	void readInGeom2D(int x, int y, std::string geomFileName, Mat2Di grid) {
		// open the geometry file
		std::ifstream geomFile(geomFileName);
		if (geomFile.is_open()) {
			std::string lineStr;
			// parse file based on given dimensions, will error if file does not match these
			// fills grid so that [0][0] is at bottom left corner of simulation
			for (int i = y - 1; i >= 0; i--) {
				std::getline(geomFile, lineStr);
				for (int j = 0; j < x; j++) {
					switch (lineStr[j]) {
					case 'f':
						grid[j][i] = SimUtil::FLUID;
						break;
					case 's':
						grid[j][i] = SimUtil::SOLID;
						break;
					case 'a':
						grid[j][i] = SimUtil::AIR;
						break;
					}
				}
			}
			geomFile.close();
		}
	}

	void readInGeom3D(int x, int y, int z, float dx, std::string geomFileName, SimUtil::Mat3Di grid) {
		// all vertices
		std::vector<Vec3> verts;
		// tris that make up fluid, these are vec3 of indices in verts list
		std::vector<Vec3> fluidTris;
		// tris that make up solid
		std::vector<Vec3> solidTris;
		
		// open the geometry file
		std::ifstream geomFile(geomFileName);
		if (geomFile.is_open()) {
			std::string lineStr;
			bool inFluidGroup = false;
			bool inSolidGroup = false;
			while (std::getline(geomFile, lineStr)) {
				std::vector<std::string> tokens;
				strSplit(lineStr, ' ', tokens);
				std::string lineType = "";
				if (tokens.size() > 0) {
					lineType = tokens[0];
				}
				
				// we loop through and collect vertices until we hit a group, which defines faces
				// then break these faces into tris and depending on group name store in corresponding tri list
				if (lineType == "v") {
					// it's a vertex definition
					// so we're no longer in a group, set both to false just in case
					inFluidGroup = false;
					inSolidGroup = false;
					// now take in vertex
					Vec3 newVert(std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3]));
					verts.push_back(newVert);
				} else if (lineType == "g") {
					// we're starting a group
					// figure out if it's the solid or fluid
					std::string groupName = tokens[1];
					if (groupName == "solid") {
						inSolidGroup = true;
					} else if (groupName == "fluid") {
						inFluidGroup = true;
					}
				} else if (lineType == "f") {
					// it's a face definition
					// go through each "v/vt" to get vertex indices
					std::vector<int> vertInds(tokens.size() - 1);
					for (int i = 1; i < tokens.size(); i++) {
						std::vector<std::string> v_vt;
						strSplit(tokens[i], '/', v_vt);
						// first one is vertex index
						vertInds[i-1] = std::stoi(v_vt[0]);
					}

					// create 2 tris for face, indices are in order [SW, SE, NE, NW]
					// must subtract 1 from each index b/c file is 1-indexed and we are 0-indexed
					Vec3 tri1(vertInds[0] - 1, vertInds[1] -1, vertInds[2] - 1);
					Vec3 tri2(vertInds[0] - 1, vertInds[2] - 1, vertInds[3] - 1);

					if (inSolidGroup) {
						solidTris.push_back(tri1);
						solidTris.push_back(tri2);
					} else if (inFluidGroup) {
						fluidTris.push_back(tri1);
						fluidTris.push_back(tri2);
					}
				} else {
					// it's something we don't care about, skip to next line
				}

			}
			geomFile.close();
		}

		// now build level sets for fluid and solid


		// now label grid based on level sets
	}


	Vec2 getGridCellPosition(float i, float j, float dx) {
		return Vec2(i*dx + 0.5f*dx, j*dx + 0.5f*dx);
	}

	int* getGridCellIndex(Vec2 pos, float dx) {
		int index[2] = { (int)(pos.x / dx), (int)(pos.y / dx) };
		return index;
	}

	Vec2 add(Vec2 vec1, Vec2 vec2) {
		return Vec2(vec1.x + vec2.x, vec1.y + vec2.y);
	}

	Vec3 add(Vec3 vec1, Vec3 vec2) {
		return Vec3(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z);
	}

	Vec2 sub(Vec2 vec1, Vec2 vec2) {
		return Vec2(vec1.x - vec2.x, vec1.y - vec2.y);
	}
	
	Vec3 sub(Vec3 vec1, Vec3 vec2) {
		return Vec3(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z);
	}

	Vec2 scale(Vec2 vec1, float scalar) {
		return Vec2(scalar * vec1.x, scalar * vec1.y);
	}

	Vec3 scale(Vec3 vec1, float scalar) {
		return Vec3(scalar * vec1.x, scalar * vec1.y, scalar * vec1.z);
	}

	float norm(Vec2 vec) {
		return sqrt(vec.x * vec.x + vec.y * vec.y);
	}

	float norm(Vec3 vec) {
		return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
	}

	template <typename T>
	double dot(T** grid1, T** grid2, int x, int y) {
		double dotProd = 0.0;
		for (int i = 0; i < x; i++) {
			for (int j = 0; j < y; j++) {
				dotProd += grid1[i][j] * grid2[i][j];
			}
		}

		return dotProd;
	}

	template <typename T>
	T max(T** grid1, int x, int y) {
		T maxVal = std::numeric_limits<T>::lowest();
		for (int i = 0; i < x; i++) {
			for (int j = 0; j < y; j++) {
				if (grid1[i][j] > maxVal) {
					maxVal = grid1[i][j];
				}
			}
		}

		return maxVal;
	}

	void strSplit(const std::string &s, char delim, std::vector<std::string> &elems) {
		std::stringstream ss;
		ss.str(s);
		std::string item;
		while (std::getline(ss, item, delim)) {
			elems.push_back(item);
		}
	}

}

// explicit instantiation of template functions for compilation
template int** SimUtil::initMat2D<int>(int, int);
template float** SimUtil::initMat2D<float>(int, int);
template double** SimUtil::initMat2D<double>(int, int);
template void SimUtil::deleteMat2D<int>(int, int, int**);
template void SimUtil::deleteMat2D<float>(int, int, float**);
template void SimUtil::deleteMat2D<double>(int, int, double**);
template void SimUtil::printMat2D<int>(int, int, int**);
template void SimUtil::printMat2D<float>(int, int, float**);
template int** SimUtil::initGrid2D<int>(int, int);
template float** SimUtil::initGrid2D<float>(int, int);
template double** SimUtil::initGrid2D<double>(int, int);
template void SimUtil::deleteGrid2D<int>(int, int, int**);
template void SimUtil::deleteGrid2D<float>(int, int, float**);
template void SimUtil::deleteGrid2D<double>(int, int, double**);
template void SimUtil::printGrid2D<int>(int, int, int**);
template void SimUtil::printGrid2D<float>(int, int, float**);
template void SimUtil::printGrid2D<double>(int, int, double**);
template int*** SimUtil::initMat3D<int>(int, int, int);
template float*** SimUtil::initMat3D<float>(int, int, int);
template double*** SimUtil::initMat3D<double>(int, int, int);
template void SimUtil::deleteMat3D<int>(int, int, int, int***);
template void SimUtil::deleteMat3D<float>(int, int, int, float***);
template void SimUtil::deleteMat3D<double>(int, int, int, double***);
template int*** SimUtil::initGrid3D<int>(int, int, int);
template float*** SimUtil::initGrid3D<float>(int, int, int);
template double*** SimUtil::initGrid3D<double>(int, int, int);
template void SimUtil::deleteGrid3D<int>(int, int, int, int***);
template void SimUtil::deleteGrid3D<float>(int, int, int, float***);
template void SimUtil::deleteGrid3D<double>(int, int, int, double***);

template double SimUtil::dot(int**, int**, int, int);
template double SimUtil::dot(float**, float**, int, int);
template double SimUtil::dot(double**, double**, int, int);

template int SimUtil::max(int**, int, int);
template float SimUtil::max(float**, int, int);
template double SimUtil::max(double**, int, int);