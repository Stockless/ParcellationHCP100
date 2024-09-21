#pragma once

#include <iostream>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

const float dotProduct(const float a[], float *&b);
float* crossProduct(const float a[], const float b[]);
bool ray_triangle_intersection(const float ray_near[], const float ray_dir[], const float Points[][3], float &t);

float** multiple_vertices(float **&triangle);
float*** multiple_triangles(float ***&triangles, uint16_t &len, const uint8_t polys[][3]);
float** triangle_interpolation(float ***&triangles, const uint8_t &N);

const bool getMeshAndFiberEndIntersection(float *&fiberP0, float *&fiberP1, const uint16_t &nPoints, const uint8_t &nPtsLine, const uint8_t &N, const uint8_t &npbp,
	float **&index, const float &step, bool ***&cubeNotEmpty, const std::vector<std::vector<std::vector<std::vector<uint32_t>>>> &centroidIndex,
	const std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> &almacen, float **&vertex, uint32_t **&polygons, uint32_t &Ind, float *&ptInt);

const bool getMeshAndFiberIntersection(float **&fiber, const uint16_t &nPoints, const uint8_t &nPtsLine, const uint8_t &N, const uint8_t &npbp, float **&index,
	const float &step, bool ***&cubeNotEmpty, const std::vector<std::vector<std::vector<std::vector<uint32_t>>>> &centroidIndex,
	const std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> &almacen, float **&vertex, uint32_t **&polygons,
	uint32_t &InInd, uint32_t &FnInd, float *&InPtInt, float *&FnPtInt);

void meshAndBundlesIntersection(float **&vertex, const uint32_t &n_vertex, uint32_t **&polygons, const uint32_t &n_polygons,
	const uint16_t &nBundles, uint32_t *&nFibers, uint16_t **&nPoints, float ****&Points, const uint8_t& nPtsLine,
	std::vector<std::vector<uint32_t>> &InTri, std::vector<std::vector<uint32_t>> &FnTri, std::vector<std::vector<std::vector<float>>> &InPoints,
	std::vector<std::vector<std::vector<float>>> &FnPoints, std::vector<std::vector<uint32_t>> &fib_index);