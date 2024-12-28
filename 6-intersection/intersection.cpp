#include "intersection.h"
#include <cmath>

// ============ Producto Punto =============
const float dotProduct(const float a[], float *&b) {
	float c = 0;

	#pragma omp simd reduction(+:c)
	for(uint32_t i=0; i<3; i++)
		c += a[i]*b[i];

	return c;
}

// ============ Producto Cruz ==============
float* crossProduct(const float a[], const float b[]) {
	float *c = new float[3];
	uint32_t i=1, j=2;

	#pragma omp simd
	for(uint32_t k=0; k<3; k++) {
		c[k] = a[i]*b[j] - a[j]*b[i];

		i = (i + 1) % 3;
		j = (j + 1) % 3;
	}

	return c;
}

// =========== Función que calcula la intersección entre un rayo con un triángulo ====================
bool ray_triangle_intersection(const float ray_near[], const float ray_dir[], const float Points[][3], float &t) {

	const float eps = 0.000001;
	float edges[2][3];

	#pragma omp simd collapse(2)
	for(uint32_t i=0; i<2; i++) {
		for(uint32_t j=0; j<3; j++) {
			edges[i][j] = Points[i+1][j] - Points[0][j];
		}
	}

	float *pvec = crossProduct(ray_dir, edges[1]);
	const float det = dotProduct(edges[0], pvec);

	if(fabs(det) < eps) {
		delete[] pvec;
		return false;
	}

  	const float inv_det = 1. / det;

  	float tvec[3];
  	#pragma omp simd
  	for(uint32_t i=0; i<3; i++)
  		tvec[i] = ray_near[i] - Points[0][i];

  	const float u = dotProduct(tvec, pvec) * inv_det;
  	delete[] pvec;
  	
  	if((u < 0.) || (u > 1.))
  		return false;

  	float *qvec = crossProduct(tvec, edges[0]);
  	const float v = dotProduct(ray_dir, qvec) * inv_det;

  	if((v < 0.) || (u + v > 1.)) {
  		delete[] qvec;
  		return false;
  	}

  	t = dotProduct(edges[1], qvec) * inv_det;
  	delete[] qvec;

  	if(t < eps)
  		return false;

	return true;
}

float** multiple_vertices(float **&triangle) {

	float **pt = new float*[6];
	for(uint32_t i=0; i<6; i++)
		pt[i] = new float[3];

	for(uint32_t i=0; i<3; i++) {
		pt[0][i] = triangle[0][i];
		pt[2][i] = triangle[1][i];
		pt[4][i] = triangle[2][i];

		pt[1][i] = (triangle[0][i] + triangle[1][i]) / 2.0;
		pt[3][i] = (triangle[1][i] + triangle[2][i]) / 2.0;
		pt[5][i] = (triangle[0][i] + triangle[2][i]) / 2.0;
	}

    delete[] triangle;
    return pt;
}

float*** multiple_triangles(float ***&triangles, uint32_t &len, const uint32_t polys[][3]) {
	float ***new_triangles = new float**[len*4];

	for(uint32_t i=0; i<len; i++) {
		float **tri = multiple_vertices(triangles[i]);

		for(uint32_t j=0; j<4; j++) {
			new_triangles[i*4 + j] = new float*[3];

			for(uint32_t k=0; k<3; k++) {
				new_triangles[i*4 + j][k] = new float[3];

				for(uint32_t l=0; l<3; l++)
					new_triangles[i*4 + j][k][l] = tri[polys[j][k]][l];
			}
		}

		for(uint32_t j=0; j<6; j++)
			delete[] tri[j];
		delete[] tri;
	}
	len = len * 4;
	delete[] triangles;
	return new_triangles;
}

float** triangle_interpolation(float ***&triangles, const uint32_t &N) {
	float **centroid = new float*[int(pow(4,N))];
	uint32_t len = 1;
	uint32_t polys[4][3] = {{0,1,5},{1,2,3},{5,3,4},{1,3,5}};

	for(uint32_t i=0; i<N; i++)
		triangles = multiple_triangles(triangles, len, polys);

	for(uint32_t i=0; i<len; i++) {
		centroid[i] = new float[3];

		for(uint32_t j=0; j<3; j++) {
			float sum = 0;

			#pragma omp simd reduction(+:sum)
			for(uint32_t k=0; k<3; k++)
				sum += triangles[i][k][j];

			centroid[i][j] = sum / 3.0;
		}
		for(uint32_t j=0; j<3; j++)
			delete[] triangles[i][j];
		delete[] triangles[i];
	}

	delete[] triangles;
	return centroid;
}

const bool getMeshAndFiberEndIntersection(float *&fiberP0, float *&fiberP1, const uint32_t &nPoints, const uint32_t &nPtsLine, const uint32_t &N, const uint32_t &npbp,
	float **&index, const float &step, bool ***&cubeNotEmpty, const std::vector<std::vector<std::vector<std::vector<uint32_t>>>> &centroidIndex,
	const std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> &almacen, float **&vertex, uint32_t **&polygons, uint32_t &Ind, float *&ptInt) {
	
    float dd[3];

    #pragma omp simd
    for(uint32_t i=0; i<3; i++)
    	dd[i] = (fiberP0[i] - fiberP1[i]) / float(npbp);

    std::vector<std::vector<uint32_t>> indexes;

    for(uint32_t i=0; i <= nPtsLine*npbp + npbp; i++) {

    	uint32_t I[3];

    	#pragma omp simd
    	for(uint32_t j=0; j<3; j++)
    		I[j] = ( (fiberP1[j] + i * dd[j]) - index[j][0] ) / step;

        #pragma omp simd collapse(3)
        for(int32_t a=-1; a<2; a++) {
        	for(int32_t b=-1; b<2; b++) {
        		for(int32_t c=-1; c<2; c++) {

        			if(cubeNotEmpty[ I[0]+a ][ I[1]+b ][ I[2]+c ]) {
        				std::vector<uint32_t> INDEX(3);

        				int32_t abc[3] = {a, b, c};

        				for(uint32_t k=0; k<3; k++)
        					INDEX[k] = I[k] + abc[k];

        				indexes.emplace_back(INDEX);
        			}
        		}
        	}
        }
    }

    if(indexes.empty())
    	return false;

    std::sort(indexes.begin(), indexes.end());
    indexes.erase( std::unique( indexes.begin(), indexes.end() ), indexes.end() );
    std::vector<std::vector<double>> listDist;

    for(const std::vector<uint32_t>& I : indexes) {
    	for(uint32_t u=0; u<centroidIndex[I[0]][I[1]][I[2]].size(); u++) {
    		double cen[3], dist = 0;

    		#pragma omp simd reduction(+:dist)
    		for(uint32_t i=0; i<3; i++) {
    			cen[i] = almacen[I[0]][I[1]][I[2]][u][i];
    			dist += fabs(cen[i] - fiberP0[i]);
    		}

    		const uint32_t& c_index = centroidIndex[I[0]][I[1]][I[2]][u];
    		listDist.emplace_back((std::vector<double>){dist, (double)c_index});
    	}
    }
    std::sort(listDist.begin(), listDist.end());
    std::vector<uint32_t> listIndex;
    ptInt = new float[3];
    
	for(const std::vector<double>& ind : listDist) {
		if(std::find(listIndex.begin(), listIndex.end(), (uint32_t)ind[1]) == listIndex.end())
			listIndex.emplace_back((uint32_t)ind[1]);
		else
			continue;

		float ray_dir[3], ray_near[3];

		#pragma omp simd
		for(uint32_t i=0; i<3; i++) {
			ray_dir[i] = fiberP1[i] - fiberP0[i];
			ray_near[i] = fiberP1[i];
		}
        
        uint32_t *Triangle = polygons[(uint32_t)ind[1]];
        float Pts[3][3];

        #pragma omp simd collapse(2)
        for(uint32_t i=0; i<3; i++) {
        	for(uint32_t j=0; j<3; j++) {
        		Pts[i][j] = vertex[Triangle[i]][j];
        	}
        }



		float t;

		// ========== Verifica la intersección entre el rayo y el triángulo =============
        if(ray_triangle_intersection(ray_near, ray_dir, Pts, t)) {

        	//std::cout << "hola" << std::endl;
        	#pragma omp simd
        	for(uint32_t i=0; i<3; i++)
        		ptInt[i] = fiberP1[i] + (fiberP1[i] - fiberP0[i])*t;

        	Ind = (uint32_t)ind[1];
        	return true;
        }
        else {
        	// ============= Ídem, pero con el rayo apuntando hacia el sentido contrario =============

        	float ray_invert[3];
        	#pragma omp simd
        	for(uint32_t i=0; i<3; i++)
        		ray_invert[i] = -ray_dir[i];

            if(ray_triangle_intersection(ray_near, ray_invert, Pts, t)) {
            	//std::cout << "hola" << std::endl;
            	#pragma omp simd
	        	for(uint32_t i=0; i<3; i++)
	        		ptInt[i] = fiberP1[i] - (fiberP1[i] - fiberP0[i])*t;

	        	Ind = (uint32_t)ind[1];
	        	return true;
            }
        }
	}
	delete[] ptInt;
	return false;
}

const bool getMeshAndFiberIntersection(float **&fiber, const uint32_t &nPoints, const uint32_t &nPtsLine, const uint32_t &N, const uint32_t &npbp, float **&index,
	const float &step, bool ***&cubeNotEmpty, const std::vector<std::vector<std::vector<std::vector<uint32_t>>>> &centroidIndex,
	const std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> &almacen, float **&vertex, uint32_t **&polygons,
	uint32_t &InInd, uint32_t &FnInd, float *&InPtInt, float *&FnPtInt) {

	bool findInt;

	findInt = getMeshAndFiberEndIntersection(fiber[0], fiber[1], nPoints, nPtsLine, N, npbp, index, step, cubeNotEmpty,
								   			 centroidIndex, almacen, vertex, polygons, InInd, InPtInt);

	//std::cout << "got mesh" << std::endl;

	if(!findInt) return false;


	findInt = getMeshAndFiberEndIntersection(fiber[nPoints-1], fiber[nPoints-2], nPoints, nPtsLine, N, npbp, index, step, cubeNotEmpty,
								   			 centroidIndex, almacen, vertex, polygons, FnInd, FnPtInt);

//	std::cout << findInt << std::endl;

	if(!findInt) {
		delete[] InPtInt;
		return false;
	}

	return true;
}

void meshAndBundlesIntersection(float **&vertex, const uint32_t &n_vertex, uint32_t **&polygons, const uint32_t &n_polygons,
	const uint32_t &nBundles, uint32_t *&nFibers, uint32_t **&nPoints, float ****&Points, const uint32_t& nPtsLine,
	std::vector<std::vector<uint32_t>> &InTri, std::vector<std::vector<uint32_t>> &FnTri, std::vector<std::vector<std::vector<float>>> &InPoints,
	std::vector<std::vector<std::vector<float>>> &FnPoints, std::vector<std::vector<uint32_t>> &fib_index) {

	uint32_t N = 1;

	float mdbv = 0; // maximum distance between vertices
	//#pragma omp parallel for schedule(dynamic)
	for(uint32_t i=0; i<n_polygons; i++) {
		float **pts = new float*[3];
		pts[0] = vertex[polygons[i][0]];
		pts[1] = vertex[polygons[i][1]];
		pts[2] = vertex[polygons[i][2]];

		float dists[3] = {0,0,0};
		if(n_polygons-i<12){
			std::cout << i << std::endl;
			std::cout << vertex[polygons[i-1][1]] << std::endl;
			std::cout << vertex[polygons[i][1]] << std::endl;
			std::cout << pts[2][0] << std::endl;
		}
		for(uint32_t k=0; k<3; k++) dists[0] += pow(pts[0][k] - pts[1][k], 2);
		for(uint32_t k=0; k<3; k++) dists[1] += pow(pts[1][k] - pts[2][k], 2);
		for(uint32_t k=0; k<3; k++) dists[2] += pow(pts[2][k] - pts[0][k], 2);
		for(uint32_t k=0; k<3; k++) dists[k] = sqrt(dists[k]);

		delete[] pts;

		const float newMax = *std::max_element(std::begin(dists), std::end(dists));
		//#pragma omp critical
		if(newMax > mdbv)
			mdbv = newMax;
	}

	const float step = mdbv / pow(2, N + 1);

	float mdbp = 0; // maximum distance between points
	for(uint32_t i=0; i<nBundles; i++) {
		for(uint32_t j=0; j<nFibers[i]; j++) {

			float dist = 0;

			for(uint32_t k=0; k<3; k++)
				dist += pow(Points[i][j][0][k] - Points[i][j][1][k], 2);

			dist = sqrt(dist);
			if(dist > mdbp)
				mdbp = dist;

			dist = 0;
			for(uint32_t k=0; k<3; k++)
				dist += pow(Points[i][j][0][k] - Points[i][j][1][k], 2);

			dist = sqrt(dist);
			if(dist > mdbp)
				mdbp = dist;
		}
	}

	uint32_t npbp = mdbp / step; // number of points between points

	std::vector<float> vx(n_vertex);
	std::vector<float> vy(n_vertex);
	std::vector<float> vz(n_vertex);
	// ===== Find outermost vertex ===========
	#pragma omp parallel for schedule(dynamic)
	for(uint32_t i=0; i<n_vertex; i++) {
		vx[i] = vertex[i][0];
		vy[i] = vertex[i][1];
		vz[i] = vertex[i][2];
	}
	
	
	std::sort(vx.begin(), vx.end());
	std::sort(vy.begin(), vy.end());
	std::sort(vz.begin(), vz.end());

	const float minx = *vx.begin()  - (nPtsLine + 1) * mdbp - 4 * step;		   // coordenada x minima
	const float maxx = *vx.rbegin() + (nPtsLine + 1) * mdbp + 4 * step;		   // coordenada x maxima
	const float miny = *vy.begin()  - (nPtsLine + 1) * mdbp - 4 * step;		   // coordenada y minima
	const float maxy = *vy.rbegin() + (nPtsLine + 1) * mdbp + 4 * step;		   // coordenada y maxima
	const float minz = *vz.begin()  - (nPtsLine + 1) * mdbp - 4 * step;		   // coordenada z minima
	const float maxz = *vz.rbegin() + (nPtsLine + 1) * mdbp + 4 * step;		   // coordenada z maxima

	//std::cout << "[" << minx << "," << maxx << "],[" << miny << "," << maxy << "],[" << minz << "," << maxz << "]." << std::endl;

	std::vector<std::vector<float**>> Bundles;
	std::vector<std::vector<uint32_t>> new_nBundles;

	for(uint32_t i=0; i<nBundles; i++) {
		std::vector<float**> new_Points;
		std::vector<uint32_t> new_nPoints;

		for(uint32_t j=0; j<nFibers[i]; j++) {

			const bool exi = ((*vx.begin() - mdbp - 2*step) <= Points[i][j][0][0]) && (Points[i][j][0][0] <= (*vx.rbegin() + mdbp + 2*step));
			const bool eyi = ((*vy.begin() - mdbp - 2*step) <= Points[i][j][0][1]) && (Points[i][j][0][1] <= (*vy.rbegin() + mdbp + 2*step));
			const bool ezi = ((*vz.begin() - mdbp - 2*step) <= Points[i][j][0][2]) && (Points[i][j][0][2] <= (*vz.rbegin() + mdbp + 2*step));

			const bool exf = ((*vx.begin() - mdbp - 2*step) <= Points[i][j][nPoints[i][j]-1][0]) && (Points[i][j][nPoints[i][j]-1][0] <= (*vx.rbegin() + mdbp + 2*step));
			const bool eyf = ((*vy.begin() - mdbp - 2*step) <= Points[i][j][nPoints[i][j]-1][1]) && (Points[i][j][nPoints[i][j]-1][1] <= (*vy.rbegin() + mdbp + 2*step));
			const bool ezf = ((*vz.begin() - mdbp - 2*step) <= Points[i][j][nPoints[i][j]-1][2]) && (Points[i][j][nPoints[i][j]-1][2] <= (*vz.rbegin() + mdbp + 2*step));

			if(exi && eyi && ezi && exf && eyf && ezf) {
				new_Points.emplace_back(Points[i][j]);
				new_nPoints.emplace_back(nPoints[i][j]);
			}
		}
		Bundles.emplace_back(new_Points);
		new_nBundles.emplace_back(new_nPoints);
	}

	for(uint32_t i=0; i<nBundles; i++) {
		nFibers[i] = Bundles[i].size();
		delete[] Points[i];
		delete[] nPoints[i];

		Points[i] = new float**[nFibers[i]];
		nPoints[i] = new uint32_t[nFibers[i]];

		for(uint32_t j=0; j<nFibers[i]; j++) {
			Points[i][j] = Bundles[i][j];
			nPoints[i][j] = new_nBundles[i][j];
		}
	}

	// ================ Obtiene la cantidad de intervalos por eje ==================
	const float mins[3] = {minx, miny, minz};
	const float maxs[3] = {maxx, maxy, maxz};
	uint32_t counts[3] = {0, 0, 0};

	for(uint32_t i=0; i<3; i++) {
		float ini = mins[i];
		while(ini < maxs[i]) {
			counts[i]++;
			ini += step;
		}
	}

	// ====== Generacion de intervalos (coordenadas de los vertices de cada cubo) ===================
	float** index = new float*[3];

	for(uint32_t i=0; i<3; i++) {
		index[i] = new float[counts[i] + 1];

		#pragma omp simd
		for(uint32_t j=0; j<counts[i] + 1; j++) {
			index[i][j] = mins[i] + j * step;
		}
	}

	std::vector<std::vector<float>> centroids(n_polygons*pow(4,N), std::vector<float>(3));

	#pragma omp parallel for schedule(dynamic)
	for(uint32_t i=0; i<n_polygons; i++) {
		float ***triangles = new float**[1];
		triangles[0] = new float*[3];

		for(uint32_t j=0; j<3; j++)
			triangles[0][j] = vertex[polygons[i][j]];

		float **centroid = triangle_interpolation(triangles, N);

		for(uint32_t j=0; j<pow(4,N); j++) {
			for(uint32_t k=0; k<3; k++)
				centroids[i*pow(4,N) + j][k] = centroid[j][k];
			delete[] centroid[j];
		}
		delete[] centroid;
	}

	std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> almacen;
	std::vector<std::vector<std::vector<std::vector<uint32_t>>>> centroidIndex;
	bool ***cubeNotEmpty = new bool**[counts[0]];

	almacen.resize(counts[0]);
	centroidIndex.resize(counts[0]);

	for(uint32_t ix = 0; ix < counts[0]; ix++) {

		almacen[ix].resize(counts[1]);
		centroidIndex[ix].resize(counts[1]);
		cubeNotEmpty[ix] = new bool*[counts[1]];

		for(uint32_t iy = 0; iy < counts[1]; iy++) {

			almacen[ix][iy].resize(counts[2]);
			centroidIndex[ix][iy].resize(counts[2]);
			cubeNotEmpty[ix][iy] = new bool[counts[2]];

			for(uint32_t iz = 0; iz < counts[2]; iz++) {
				cubeNotEmpty[ix][iy][iz] = false;
			}
		}
	}

	// =========================== En cada cubo almacena una lista de centroides y su indice ======================

	for(uint32_t i=0; i<n_polygons*pow(4,N); i++) {

		uint32_t I[3];

		#pragma omp simd
		for(uint32_t j=0; j<3; j++)
			I[j] = (centroids[i][j] - index[j][0]) / step;

		almacen[I[0]][I[1]][I[2]].emplace_back(centroids[i]);
		centroidIndex[I[0]][I[1]][I[2]].emplace_back(uint32_t(i / pow(4,N)));
		cubeNotEmpty[I[0]][I[1]][I[2]] = true;
	}

	// ===============================================================================================================

	InTri.resize(nBundles);
	FnTri.resize(nBundles);
	InPoints.resize(nBundles);
	FnPoints.resize(nBundles);
	fib_index.resize(nBundles);
	for(uint32_t i=0; i<nBundles; i++) {

		std::cout << "Bundle: " << i+1 << "/" << nBundles;
		std::cout << ", Num. Fibers: " << nFibers[i] << std::endl;

		std::vector<uint32_t> listFibInd;
		std::vector<std::vector<uint32_t>> listTri;
		std::vector<std::vector<std::vector<float>>> listPtInt;
		#pragma omp parallel for schedule(dynamic)
		for(uint32_t j=0; j<nFibers[i]; j++) {
			if(j > 100000) continue;
			bool findInt; // find intersection
			uint32_t InT, FnT; // initial and final triangle
			float *InPtInt, *FnPtInt; // Initial and Final point intersection
			//std::cout<<i<<" "<<j<<std::endl;
			findInt = getMeshAndFiberIntersection(Points[i][j], nPoints[i][j], nPtsLine, N, npbp, index, step, cubeNotEmpty,
				 						    centroidIndex, almacen, vertex, polygons, InT, FnT, InPtInt, FnPtInt);
			//std::cout<<i<<" "<<j<<" here"<<std::endl;

			if(!findInt) {
				continue;
			}
			

			//#pragma omp critical

			#pragma omp critical
			listFibInd.emplace_back(j);
			#pragma omp critical
			listTri.emplace_back((std::vector<uint32_t>){InT, FnT});
			#pragma omp critical
			listPtInt.emplace_back((std::vector<std::vector<float>>){{InPtInt[0], InPtInt[1], InPtInt[2]}, {FnPtInt[0], FnPtInt[1], FnPtInt[2]}});

			delete[] InPtInt;
			delete[] FnPtInt;
		}

		for(uint32_t j=0; j<listTri.size(); j++) {
			InTri[i].emplace_back(listTri[j][0]);
			FnTri[i].emplace_back(listTri[j][1]);
			InPoints[i].emplace_back(listPtInt[j][0]);
			FnPoints[i].emplace_back(listPtInt[j][1]);
		}
		fib_index[i] = listFibInd;
	}

	for(uint32_t i=0; i<3; i++)
		delete[] index[i];
	delete[] index;

	for(uint32_t i=0; i<counts[0]; i++) {
		for(uint32_t j=0; j<counts[1]; j++) {
			delete[] cubeNotEmpty[i][j];
		}
		delete[] cubeNotEmpty[i];
	}
	delete[] cubeNotEmpty;

	for(uint32_t i=0; i<nBundles; i++) {
		for(uint32_t j=0; j<nFibers[i]; j++) {
			for(uint32_t k=0; k<nPoints[i][j]; k++)
				delete[] Points[i][j][k];
			delete[] Points[i][j];
		}
		delete[] Points[i];
		delete[] nPoints[i];
	}
	
	delete[] Points;
	delete[] nPoints;
	delete[] nFibers;
}
