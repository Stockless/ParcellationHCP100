#include "intersection.h"
#include "alignment.h"
#include "regions.h"
#include <iostream>

int main(int argc, char const *argv[])
{

	// ====================== Se lee el mallado cortical de ambos hemisferios =========================
	float **Lvertex, **Rvertex;
	uint32_t **Lpolygons, **Rpolygons;
	uint32_t n_Lvertex, n_Lpolygons;
	uint32_t n_Rvertex, n_Rpolygons;

	read_mesh_vtk(argv[1], Lvertex, Lpolygons, n_Lvertex, n_Lpolygons);
	read_mesh_vtk(argv[2], Rvertex, Rpolygons, n_Rvertex, n_Rpolygons);

	//const uint8_t th = atoi(argv[4]); // parámetro de umbralización
	const uint8_t ero = atoi(argv[3]); // parámetro de erosión
	const uint8_t dil = atoi(argv[4]); // parámetro de dilatación
	std::vector<std::vector<uint32_t>> LTri, RTri;

	// ========== lee la carpeta "intersection" ==============

	read_parcels(argv[5], LTri);
	read_parcels(argv[6], RTri);
	// Umbralización

	//Thresholding(th, LTri);
	//Thresholding(th, RTri);

	// Erosión

	Erosion(ero, LTri, Lpolygons);
	Erosion(ero, RTri, Rpolygons);

	// Dilatación

	Dilation(dil, LTri, Lpolygons, n_Lpolygons);	
	Dilation(dil, RTri, Rpolygons, n_Rpolygons);

	// ============ escribe la carpeta "regions" =============

	//std::cout << LTri.size() << std::endl;
	write_parcels(argv[7], argv[5], LTri);
	write_parcels(argv[8], argv[6], RTri);

	// Se libera el espacio reservado en memoria

	//Delete(Lvertex, Lpolygons, n_Lvertex, n_Lpolygons);
	//Delete(Rvertex, Rpolygons, n_Rvertex, n_Rpolygons);

	return 0;
}
