#pragma once

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include "dirent.h"
#include <sys/stat.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

std::string GetCurrentWorkingDir( void );

void read_mesh_vtk(const std::string &filename, float **&vertex, uint32_t **&polygons, uint32_t &len_vertex, uint32_t &len_polygons);
void read_mesh(const std::string &filename, float **&vertex, uint32_t **&polygons, uint32_t &len_vertex, uint32_t &len_polygons);

std::vector<std::string> listDir(const std::string &path);

void read_parcels(const std::string &path, std::vector<std::vector<uint32_t>> &Tri);

void write_parcels(const std::string &path, const std::string &pathBundles,
	const std::vector<std::vector<uint32_t>> &Tri);

