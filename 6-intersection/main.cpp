#include "intersection.h"
#include "Tools.h"
#include <cmath>

int main(int argc, char const *argv[])
{
	float **Lvertex, **Rvertex;
	uint32_t **Lpolygons, **Rpolygons;
	uint32_t n_Lvertex, n_Lpolygons;
	uint32_t n_Rvertex, n_Rpolygons;

	// ========================== Se lee el mallado =======================================
	std::cout <<argv[1]<< std::endl;
	read_mesh(argv[1], Lvertex, Lpolygons, n_Lvertex, n_Lpolygons);
	read_mesh(argv[2], Rvertex, Rpolygons, n_Rvertex, n_Rpolygons);

	// ========================== Se leen las fibras ======================================
	uint32_t nLBundles, nRBundles;
	uint32_t *nLFibers, *nRFibers;
	uint32_t **nLPoints, **nRPoints;
	float ****LPoints, ****RPoints;
	read_bundles(argv[3], nLBundles, nLFibers, nLPoints, LPoints);
	std::cout<<"Readed Left Bundles"<<std::endl;
	std::cout<<"nLBundles: "<<nLBundles<<std::endl;
	read_bundles(argv[4], nRBundles, nRFibers, nRPoints, RPoints);
	std::cout<<"Readed Right Bundles"<<std::endl;
	std::cout<<"nRBundles: "<<nRBundles<<std::endl;

	// =========================== Intersección ===========================================
	const uint32_t nPtsLine = 2;

	std::vector<std::vector<uint32_t>> Lfib_index, Rfib_index;
	std::vector<std::vector<uint32_t>> LInTri, LFnTri, RInTri, RFnTri;
	std::vector<std::vector<std::vector<float>>> LInPoints, LFnPoints, RInPoints, RFnPoints;

	meshAndBundlesIntersection(Lvertex, n_Lvertex, Lpolygons, n_Lpolygons, nLBundles, nLFibers,
							   nLPoints, LPoints, nPtsLine, LInTri, LFnTri, LInPoints, LFnPoints, Lfib_index);
	std::cout<<"Intersection Left"<<std::endl;
	meshAndBundlesIntersection(Rvertex, n_Rvertex, Rpolygons, n_Rpolygons, nRBundles, nRFibers,
							   nRPoints, RPoints, nPtsLine, RInTri, RFnTri, RInPoints, RFnPoints, Rfib_index);

	write_intersection(argv[5], argv[3], LInTri, LFnTri, LInPoints, LFnPoints, Lfib_index);
	write_intersection(argv[6], argv[4], RInTri, RFnTri, RInPoints, RFnPoints, Rfib_index);
	std::cout<<"Wrote Intersection"<<std::endl;
	Delete(Lvertex, Lpolygons, n_Lvertex, n_Lpolygons);
	Delete(Rvertex, Rpolygons, n_Rvertex, n_Rpolygons);

	return 0;
}