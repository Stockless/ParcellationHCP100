#include <iostream>
#include "alignment.h"
#include <cmath>

// ======== obtiene el directorio actual ========
std::string GetCurrentWorkingDir( void ) {
	char buff[FILENAME_MAX];
	return GetCurrentDir( buff, FILENAME_MAX );
}

// función que lee el mallado (vértices y polígonos de los triángulos) en vtk
void read_mesh_vtk(const std::string &filename, float **&vertex, uint32_t **&polygons, uint32_t &len_vertex, uint32_t &len_polygons) {

	std::vector<std::vector<float>> vertex_vector;
	std::vector<std::vector<int>> polygons_vector;
	int dim = 3;
	len_vertex = 0; 
	len_polygons = 0;
	
	std::ifstream file(filename);
	if(file.is_open()) {

		std::string id;
		std::string s1,s2,s3;
		std::string line;
		while(getline(file, line)){  
		    std::stringstream iss1(line);
		    iss1>>id;
		    if (id=="v"){
		    	std::vector<float> vertex_coord(dim);
			    iss1>>vertex_coord[0];
			    iss1>>vertex_coord[1];
			    iss1>>vertex_coord[2];
			    len_vertex++;
			    vertex_vector.push_back(vertex_coord);
			    //std::cout<<"id "<<id<<" x: "<<vertex_coord[0]<<" y:"<<vertex_coord[1]<<" z: "<<vertex_coord[2]<<std::endl;
			    id = "n";
			} else if (id=="f"){
				std::vector<int> poly_vertex(dim);
			    iss1>>s1;
  				poly_vertex[0] = std::stoi(s1)-1;     
			    iss1>>s2;
  				poly_vertex[1] = std::stoi(s2)-1;     
			    iss1>>s3;
  				poly_vertex[2] = std::stoi(s3)-1;     
			    len_polygons++;
			    //std::cout<<"id "<<id<<" v1: "<<poly_vertex[0]<<" v2:"<<poly_vertex[1]<<" v3: "<<poly_vertex[2]<<std::endl;
			    polygons_vector.push_back(poly_vertex);
			    id = "n";
			}

		}
		// std::cout<<"Número de vértices"<<len_vertex<<std::endl;
		// std::cout<<"Número de triángulos"<<len_polygons<<std::endl;

	}

	file.close();

	vertex = new float*[len_vertex]; //Inicializa el vector de vértices
	for(uint32_t i=0; i<len_vertex; i++) {
		vertex[i] = new float[dim]; 	//Inicializa cada vértice con la dimensión... (x,y,z)
		vertex[i][0] = vertex_vector[i][0]; // inicializa cada vértice con coordenadas (x,y,z)
		vertex[i][1] = vertex_vector[i][1];
		vertex[i][2] = vertex_vector[i][2];
	}
	polygons = new uint32_t*[len_polygons];	
	for(uint32_t i=0; i<len_polygons; i++) {
		polygons[i] = new uint32_t[dim]; 	//Inicializa cada vértice con la dimensión... (x,y,z)
		polygons[i][0] = polygons_vector[i][0]; // inicializa cada vértice con coordenadas (x,y,z)
		polygons[i][1] = polygons_vector[i][1];
		polygons[i][2] = polygons_vector[i][2];
	}

}

// función que lee el mallado (vértices y polígonos de los triángulos) de un archivo *.mesh
void read_mesh(const std::string &filename, float **&vertex, uint32_t **&polygons, uint32_t &len_vertex, uint32_t &len_polygons) {

	uint32_t option = 2;

	if(option == 1){
		std::ifstream file(filename, std::ios::in | std::ios::binary);

		if(file.is_open()) {

			file.seekg(17);

			uint32_t dim;

			file.read ((char *) & dim, sizeof(uint32_t)); // dimensiones

			file.seekg(29);

			file.read ((char *) & len_vertex, sizeof(uint32_t)); // largo de vértices

			vertex = new float*[len_vertex];

			for(uint32_t i=0; i<len_vertex; i++) {
				vertex[i] = new float[3];

				for(uint32_t j=0; j<3; j++)
					file.read ((char *) & vertex[i][j], sizeof(float)); // lee cada punto (x,y,z)
			}

			file.seekg(4*3*len_vertex + 41);

			file.read ((char*)&len_polygons, sizeof(uint32_t)); // largo de polígonos
			polygons = new uint32_t*[len_polygons];

			for(uint32_t i = 0; i < len_polygons; i++) {
				polygons[i] = new uint32_t[3];

				for(uint32_t j = 0; j < 3; j++) 
					file.read ((char*)&polygons[i][j], sizeof(uint32_t)); // lee cada índice del triángulo (a,b,c)
			}
		}

		file.close();
	}

	else if(option == 2){

		std::vector<std::vector<float>> vertex_vector;
		std::vector<std::vector<int>> polygons_vector;
		int dim = 3;
		len_vertex = 0; 
		len_polygons = 0;
		
		std::ifstream file(filename);
		if(file.is_open()) {

			std::string id;
			std::string s1,s2,s3;
			std::string line;
			while(getline(file, line)){  
			    std::stringstream iss1(line);
			    iss1>>id;
			    if (id=="v"){
			    	std::vector<float> vertex_coord(dim);
				    iss1>>vertex_coord[0];
				    iss1>>vertex_coord[1];
				    iss1>>vertex_coord[2];
				    len_vertex++;
				    vertex_vector.push_back(vertex_coord);
				    id = "n";
				} else if (id=="f"){
					std::vector<int> poly_vertex(dim);
				    iss1>>s1; 
	  				poly_vertex[0] = std::stoi(s1)-1;
				    iss1>>s2;
	  				poly_vertex[1] = std::stoi(s2)-1;     
				    iss1>>s3;
	  				poly_vertex[2] = std::stoi(s3)-1;     
				    len_polygons++;
				    polygons_vector.push_back(poly_vertex);
				    id = "n";
				}

			}
		}

		file.close();

		vertex = new float*[len_vertex]; 
		for(uint32_t i=0; i<len_vertex; i++) {
			vertex[i] = new float[dim]; 
			vertex[i][0] = vertex_vector[i][0]; 
			vertex[i][1] = vertex_vector[i][1];
			vertex[i][2] = vertex_vector[i][2];
		}
		polygons = new uint32_t*[len_polygons];	
		for(uint32_t i=0; i<len_polygons; i++) {
			polygons[i] = new uint32_t[dim]; 	
			polygons[i][0] = polygons_vector[i][0]; 
			polygons[i][1] = polygons_vector[i][1];
			polygons[i][2] = polygons_vector[i][2];
		}
	}
}


// ========= función que retorna una lista de string con los nombres de los archivos de una ruta determinada ==========
std::vector<std::string> listDir(const std::string &path) {

	std::vector<std::string> Dir;

	DIR *dir;
	struct dirent *ent;
	if((dir = opendir(path.c_str())) != NULL) {
	  /* get all the file names and directories within directory */
	  while((ent = readdir (dir)) != NULL) {
	    const std::string file_dir = ent->d_name;

	    if(*file_dir.rbegin() == 'a')
	    	Dir.emplace_back(path + file_dir);
	  }
	  closedir (dir);
	}
	else {
	  /* could not open directory */
	  perror("Error");
	  exit( EXIT_FAILURE );
	}

	std::sort(Dir.begin(), Dir.end());
	return Dir;
}

// ================ lee la carpeta "intersection" con los datos de intersección de cada archivo ================
void read_parcels(const std::string &path, std::vector<std::vector<uint32_t>> &Tri) {

	std::vector<std::string> intersectionDir;

	DIR *dir;
	struct dirent *ent;
	if((dir = opendir(path.c_str())) != NULL) {
	  /* get all the file names and directories within directory */
	  while((ent = readdir (dir)) != NULL) {
	    const std::string intersection_dir = ent->d_name;

	    if(intersection_dir[0] != '.')
	    	intersectionDir.emplace_back(path + intersection_dir);
	  }
	  closedir (dir);
	}
	else {
	  /* could not open directory */
	  perror("Error");
	  exit( EXIT_FAILURE );
	}


	std::sort(intersectionDir.begin(), intersectionDir.end());
	const uint16_t nIntersections = intersectionDir.size();

	Tri.resize(nIntersections);

	for(uint16_t i=0; i<nIntersections; i++) {

		std::string filename = intersectionDir[i];
		std::ifstream file(filename, std::ios::in | std::ios::binary);

		if(file.is_open()) {
			uint32_t total_triangles;
			file.read ((char*)&total_triangles, sizeof(uint32_t));

			Tri[i].resize(total_triangles);

			for(uint32_t j=0; j<total_triangles; j++)
				file.read ((char*)&Tri[i][j], sizeof(uint32_t));
		}
		file.close();
	}
}

// ======= genera la carpeta "regions" con los datos de la generación de regiones corticales ===========
void write_parcels(const std::string &path, const std::string &pathBundles,
	const std::vector<std::vector<uint32_t>> &Tri) {
	
	std::vector<std::string> bundlesDir;

	DIR *dirInt;
	struct dirent *ent;
	if((dirInt = opendir(pathBundles.c_str())) != NULL) {
	  /* get all the file names and directories within directory */
	  while((ent = readdir (dirInt)) != NULL) {
	    const std::string bundle_dir = ent->d_name;
	    if(*bundle_dir.rbegin() == 'a')
	    	bundlesDir.emplace_back(bundle_dir.substr(0,bundle_dir.size() - 12));
	  }
	  closedir (dirInt);
	}
	else {
	  /* could not open directory */
	  perror("Error");
	  exit( EXIT_FAILURE );
	}

	std::sort(bundlesDir.begin(), bundlesDir.end());

	DIR *dir;
	if((dir = opendir(path.c_str())) == NULL) { // Checks if a directory path exists

		const int dir_err = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if(dir_err == -1) {
			perror("Error creating directory!");
			exit( EXIT_FAILURE );
		}
	}
	closedir (dir);

	for(uint16_t i=0; i<Tri.size(); i++) {

		// std::cout << path + bundlesDir[i] << std::endl;
		std::ofstream file(path + bundlesDir[i] + ".parcelsdata", std::ios::out | std::ios::binary);


		if(file.is_open()) {
			const uint32_t len_Tri = Tri[i].size();
			file.write(reinterpret_cast<const char*>( &len_Tri ), sizeof( uint32_t )); // número de triágulos iniciales

			for(uint32_t j=0; j<len_Tri; j++) // escribe el índice cada triángulo inicial{
				file.write(reinterpret_cast<const char*>( &Tri[i][j] ), sizeof( uint32_t ));
		}
		file.close();
	}
}

