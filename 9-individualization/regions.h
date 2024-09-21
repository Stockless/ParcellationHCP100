#pragma once

#include <set>
#include <vector>
#include <iostream>
#include <algorithm>

std::vector<uint32_t> get_edges(const std::vector<uint32_t> &Tri, uint32_t **&polygons);

void Thresholding(const uint8_t &th, std::vector<std::vector<uint32_t>> &Tri);
void Dilation(const uint8_t &dil, std::vector<std::vector<uint32_t>> &Tri, uint32_t **&polygons, uint32_t &len_polygons);
void Erosion(const uint8_t &ero, std::vector<std::vector<uint32_t>> &Tri, uint32_t **&polygons);