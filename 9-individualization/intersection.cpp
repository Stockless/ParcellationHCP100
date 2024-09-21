#include "intersection.h"
#include <cmath>

// ============ Producto Punto =============
const float dotProduct(const float a[], float *&b) {
	float c = 0;

	#pragma omp simd reduction(+:c)
	for(uint16_t i=0; i<3; i++)
		c += a[i]*b[i];

	return c;
}

// ============ Producto Cruz ==============
float* crossProduct(const float a[], const float b[]) {
	float *c = new float[3];
	uint16_t i=1, j=2;

	#pragma omp simd
	for(uint16_t k=0; k<3; k++) {
		c[k] = a[i]*b[j] - a[j]*b[i];

		i = (i + 1) % 3;
		j = (j + 1) % 3;
	}

	return c;
}

// =========== Función que calcula la intersección entre un rayo con un triángulo ====================
// Para más información, visitar: https://cadxfem.org/inf/Fast%20MinimumStorage%20RayTriangle%20Intersection.pdf
bool ray_triangle_intersection(const float ray_near[], const float ray_dir[], const float Points[][3], float &t) {

	const float eps = 0.000001;
	float edges[2][3];

	#pragma omp simd collapse(2)
	for(uint16_t i=0; i<2; i++) {
		for(uint16_t j=0; j<3; j++) {
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
  	for(uint16_t i=0; i<3; i++)
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

// ========= Generación de 6 vértices a partir de 1 triángulo =============
float** multiple_vertices(float **&triangle) {

	float **pt = new float*[6];
	for(uint8_t i=0; i<6; i++)
		pt[i] = new float[3];

	for(uint8_t i=0; i<3; i++) {
		// ======== los 3 vértices principales se mantienen =======
		pt[0][i] = triangle[0][i];
		pt[2][i] = triangle[1][i];
		pt[4][i] = triangle[2][i];

		// ======== se calculan los 3 nuevos vértices como los puntos intermedios entre los 3 vértices anteriores =======
		pt[1][i] = (triangle[0][i] + triangle[1][i]) / 2.0;
		pt[3][i] = (triangle[1][i] + triangle[2][i]) / 2.0;
		pt[5][i] = (triangle[0][i] + triangle[2][i]) / 2.0;
	}

    delete[] triangle;
    return pt;
}

// ========= Se subdivide el triángulo "len" veces ==================
float*** multiple_triangles(float ***&triangles, uint16_t &len, const uint8_t polys[][3]) {
	float ***new_triangles = new float**[len*4];

	for(uint16_t i=0; i<len; i++) { // se obtienen los nuevos vértices del triángulo
		float **tri = multiple_vertices(triangles[i]);

		// se definen los nuevos triángulos con sus aristas y vértices
		for(uint8_t j=0; j<4; j++) {
			new_triangles[i*4 + j] = new float*[3];

			for(uint8_t k=0; k<3; k++) {
				new_triangles[i*4 + j][k] = new float[3];

				for(uint8_t l=0; l<3; l++)
					new_triangles[i*4 + j][k][l] = tri[polys[j][k]][l];
			}
		}

		for(uint8_t j=0; j<6; j++)
			delete[] tri[j];
		delete[] tri;
	}
	len = len * 4;
	delete[] triangles;
	return new_triangles;
}

// ======== Se subdivide 1 triángulo N veces y retorna sus centroides ============
float** triangle_interpolation(float ***&triangles, const uint8_t &N) {
	float **centroid = new float*[int(pow(4,N))];
	uint16_t len = 1; // factor que incrementa en múltiplos de 4
	uint8_t polys[4][3] = {{0,1,5},{1,2,3},{5,3,4},{1,3,5}}; // índice de los vértices de los 4 nuevos triángulos

	for(uint8_t i=0; i<N; i++) // se subdivide el triángulo N veces de manera recursiva
		triangles = multiple_triangles(triangles, len, polys);

	// se calcula el centroide de cada subtriángulo
	for(uint16_t i=0; i<len; i++) {
		centroid[i] = new float[3];

		for(uint8_t j=0; j<3; j++) {
			float sum = 0;

			#pragma omp simd reduction(+:sum)
			for(uint8_t k=0; k<3; k++)
				sum += triangles[i][k][j];

			centroid[i][j] = sum / 3.0;
		}
		for(uint8_t j=0; j<3; j++)
			delete[] triangles[i][j];
		delete[] triangles[i];
	}

	delete[] triangles;
	return centroid;
}

// =========== Proyecta los extremos de la fibra y calcula la intersección con el mallado ==============
const bool getMeshAndFiberEndIntersection(float *&fiberP0, float *&fiberP1, const uint16_t &nPoints, const uint8_t &nPtsLine, const uint8_t &N, const uint8_t &npbp,
	float **&index, const float &step, bool ***&cubeNotEmpty, const std::vector<std::vector<std::vector<std::vector<uint32_t>>>> &centroidIndex,
	const std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> &almacen, float **&vertex, uint32_t **&polygons, uint32_t &Ind, float *&ptInt) {
	
    float dd[3]; // delta entre los 2 últimos puntos del extremo de la fibra

    #pragma omp simd
    for(uint8_t i=0; i<3; i++)
    	dd[i] = (fiberP0[i] - fiberP1[i]) / float(npbp);

    std::vector<std::vector<uint16_t>> indexes; // índice del cubo que contiene al posible triángulo que intersecta con la fibra

    for(uint16_t i=0; i <= nPtsLine*npbp + npbp; i++) { // el extremo de la fibra es proyectado

    	uint16_t I[3]; // índice del cubo que contiene al punto de proyección

    	#pragma omp simd
    	for(uint16_t j=0; j<3; j++)
    		I[j] = ( (fiberP1[j] + i * dd[j]) - index[j][0] ) / step;

    	// ======== se recorre el cubo junto a sus vecinos ==========
        #pragma omp simd collapse(3)
        for(int8_t a=-1; a<2; a++) {
        	for(int8_t b=-1; b<2; b++) {
        		for(int8_t c=-1; c<2; c++) {

        			// si el cubo no está vacío, entonces se guarda el índice de ese cubo en "indexes"
        			if(cubeNotEmpty[ I[0]+a ][ I[1]+b ][ I[2]+c ]) {
        				std::vector<uint16_t> INDEX(3);

        				int8_t abc[3] = {a, b, c};

        				for(uint8_t k=0; k<3; k++)
        					INDEX[k] = I[k] + abc[k];

        				indexes.emplace_back(INDEX);
        			}
        		}
        	}
        }
    }

    if(indexes.empty()) // si la lista está vacía, no habrá intersección con la fibra
    	return false;

    std::sort(indexes.begin(), indexes.end());
    indexes.erase( std::unique( indexes.begin(), indexes.end() ), indexes.end() );
    std::vector<std::vector<double>> listDist; // lista de distancias entre el extremo de la fibra y los centroides más cercanos

    // =========== se recorren los índices de estos cubos ==========
    for(const std::vector<uint16_t>& I : indexes) {
    	for(uint16_t u=0; u<centroidIndex[I[0]][I[1]][I[2]].size(); u++) {
    		double cen[3], dist = 0;

    		// se calcula la distancia entre cada centroide dentro del cubo con el extremo de la fibra
    		#pragma omp simd reduction(+:dist)
    		for(uint8_t i=0; i<3; i++) {
    			cen[i] = almacen[I[0]][I[1]][I[2]][u][i];
    			dist += fabs(cen[i] - fiberP0[i]);
    		}

    		const uint32_t& c_index = centroidIndex[I[0]][I[1]][I[2]][u];
    		listDist.emplace_back((std::vector<double>){dist, (double)c_index}); // se guarda la distancia junto con su índice
    	}
    }
    std::sort(listDist.begin(), listDist.end()); // se ordenan las distancias de menor a mayor
    std::vector<uint32_t> listIndex; // lista de los índices de los centroides sin repetición
    ptInt = new float[3]; // punto de intersección
    
    // ======= se evalúa la intersección desde los centroides más cercanos a los más lejanos =========
	for(const std::vector<double>& ind : listDist) {
		// si ya se evaluó el índice de un centroide, entonces continúa con el siguiente
		if(std::find(listIndex.begin(), listIndex.end(), (uint32_t)ind[1]) == listIndex.end())
			listIndex.emplace_back((uint32_t)ind[1]);
		else
			continue;

		float ray_dir[3], ray_near[3];

		#pragma omp simd
		for(uint8_t i=0; i<3; i++) {
			ray_dir[i] = fiberP1[i] - fiberP0[i]; // sentido del rayo
			ray_near[i] = fiberP1[i]; // punto de donde sale el rayo
		}
        
        uint32_t *Triangle = polygons[(uint32_t)ind[1]]; // triángulo que contiene el índice del centroide
        float Pts[3][3]; // vértices de dicho triángulo

        #pragma omp simd collapse(2)
        for(uint8_t i=0; i<3; i++) {
        	for(uint8_t j=0; j<3; j++) {
        		Pts[i][j] = vertex[Triangle[i]][j];
        	}
        }

		float t; // distancia desde el rayo de origen al punto de intersección en el triángulo

		// ========== Verifica la intersección entre el rayo y el triángulo =============
        if(ray_triangle_intersection(ray_near, ray_dir, Pts, t)) {

        	#pragma omp simd
        	for(uint8_t i=0; i<3; i++)
        		ptInt[i] = fiberP1[i] + (fiberP1[i] - fiberP0[i])*t;

        	Ind = (uint32_t)ind[1];
        	return true;
        }
        else {
        	// ============= Ídem, pero con el rayo apuntando hacia el sentido contrario =============

        	float ray_invert[3];
        	#pragma omp simd
        	for(uint8_t i=0; i<3; i++)
        		ray_invert[i] = -ray_dir[i];

            if(ray_triangle_intersection(ray_near, ray_invert, Pts, t)) {

            	#pragma omp simd
	        	for(uint8_t i=0; i<3; i++)
	        		ptInt[i] = fiberP1[i] - (fiberP1[i] - fiberP0[i])*t;

	        	Ind = (uint32_t)ind[1];
	        	return true;
            }
        }
	}

	delete[] ptInt;
	return false; // retorna falso en caso que no haya encontrado ninguna intersección
}

// ======== Obtiene los datos de intersección entre 1 fibra y el mallado cortical ==========
const bool getMeshAndFiberIntersection(float **&fiber, const uint16_t &nPoints, const uint8_t &nPtsLine, const uint8_t &N, const uint8_t &npbp, float **&index,
	const float &step, bool ***&cubeNotEmpty, const std::vector<std::vector<std::vector<std::vector<uint32_t>>>> &centroidIndex,
	const std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> &almacen, float **&vertex, uint32_t **&polygons,
	uint32_t &InInd, uint32_t &FnInd, float *&InPtInt, float *&FnPtInt) {

	bool findInt; // verifica si encuentra o no intersección

	// Primer extremo de la fibra
	findInt = getMeshAndFiberEndIntersection(fiber[0], fiber[1], nPoints, nPtsLine, N, npbp, index, step, cubeNotEmpty,
								   			 centroidIndex, almacen, vertex, polygons, InInd, InPtInt);

	if(!findInt) return false;

	// Segundo extremo de la fibra
	findInt = getMeshAndFiberEndIntersection(fiber[nPoints-1], fiber[nPoints-2], nPoints, nPtsLine, N, npbp, index, step, cubeNotEmpty,
								   			 centroidIndex, almacen, vertex, polygons, FnInd, FnPtInt);

	if(!findInt) {
		delete[] InPtInt;
		return false;
	}

	return true;
}

// ============= Calcula la intersección entre un fascículo y el mallado ===================
void meshAndBundlesIntersection(float **&vertex, const uint32_t &n_vertex, uint32_t **&polygons, const uint32_t &n_polygons,
	const uint16_t &nBundles, uint32_t *&nFibers, uint16_t **&nPoints, float ****&Points, const uint8_t& nPtsLine,
	std::vector<std::vector<uint32_t>> &InTri, std::vector<std::vector<uint32_t>> &FnTri, std::vector<std::vector<std::vector<float>>> &InPoints,
	std::vector<std::vector<std::vector<float>>> &FnPoints, std::vector<std::vector<uint32_t>> &fib_index) {

	uint8_t N = 1; // cantidad de subdivisiones de cada triángulo del mallado

	float mdbv = 0; // maximum distance between vertices (distancia máxima entre los vertices)

	#pragma omp parallel for schedule(dynamic)
	for(uint32_t i=0; i<n_polygons; i++) {
		float **pts = new float*[3];
		pts[0] = vertex[polygons[i][0]];
		pts[1] = vertex[polygons[i][1]];
		pts[2] = vertex[polygons[i][2]];

		float dists[3] = {0,0,0};
		for(uint8_t k=0; k<3; k++) dists[0] += pow(pts[0][k] - pts[1][k], 2);
		for(uint8_t k=0; k<3; k++) dists[1] += pow(pts[1][k] - pts[2][k], 2);
		for(uint8_t k=0; k<3; k++) dists[2] += pow(pts[2][k] - pts[0][k], 2);
		for(uint8_t k=0; k<3; k++) dists[k] = sqrt(dists[k]);

		delete[] pts;

		const float newMax = *std::max_element(std::begin(dists), std::end(dists));
		#pragma omp critical
		if(newMax > mdbv)
			mdbv = newMax;
	}

	const float step = mdbv / pow(2, N + 1); // tamaño de los cubos

	float mdbp = 0; // maximum distance between points (distancia máxima entre los puntos de las fibras)

	for(uint16_t i=0; i<nBundles; i++) {
		for(uint32_t j=0; j<nFibers[i]; j++) {

			float dist = 0; // distancia entre los 2 puntos extremos de una fibra

			for(uint8_t k=0; k<3; k++)
				dist += pow(Points[i][j][0][k] - Points[i][j][1][k], 2);

			dist = sqrt(dist);
			if(dist > mdbp)
				mdbp = dist;

			dist = 0; // distancia entre los 2 puntos del otro extremo de la fibra
			for(uint8_t k=0; k<3; k++)
				dist += pow(Points[i][j][nPoints[i][j]-1][k] - Points[i][j][nPoints[i][j]-2][k], 2);

			dist = sqrt(dist);
			if(dist > mdbp)
				mdbp = dist;
		}
	}

	uint8_t npbp = mdbp / step; // number of points between points (número de puntos entre los puntos de proyección)

	std::vector<float> vx(n_vertex);
	std::vector<float> vy(n_vertex);
	std::vector<float> vz(n_vertex);
	// ===== Encuentra los vértices más externos ===========
	#pragma omp parallel for schedule(dynamic)
	for(uint32_t i=0; i<n_vertex; i++) {
		vx[i] = vertex[i][0];
		vy[i] = vertex[i][1];
		vz[i] = vertex[i][2];
	}
	
	std::sort(vx.begin(), vx.end());
	std::sort(vy.begin(), vy.end());
	std::sort(vz.begin(), vz.end());

	const float minx = *vx.begin()  - (nPtsLine + 1) * mdbp - 4 * step;		   // coordenada x mínima
	const float maxx = *vx.rbegin() + (nPtsLine + 1) * mdbp + 4 * step;		   // coordenada x máxima
	const float miny = *vy.begin()  - (nPtsLine + 1) * mdbp - 4 * step;		   // coordenada y mínima
	const float maxy = *vy.rbegin() + (nPtsLine + 1) * mdbp + 4 * step;		   // coordenada y máxima
	const float minz = *vz.begin()  - (nPtsLine + 1) * mdbp - 4 * step;		   // coordenada z mínima
	const float maxz = *vz.rbegin() + (nPtsLine + 1) * mdbp + 4 * step;		   // coordenada z máxima

	std::vector<std::vector<float**>> Bundles;
	std::vector<std::vector<uint16_t>> new_nBundles;

	// ========= Considera solo aquellas fibras que estén dentro de los márgenes establecidos ========
	for(uint16_t i=0; i<nBundles; i++) {
		std::vector<float**> new_Points;
		std::vector<uint16_t> new_nPoints;

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

	// ======== Actualiza los otros parámetros como el número de fibras por fascículo, etc ========
	for(uint16_t i=0; i<nBundles; i++) {
		nFibers[i] = Bundles[i].size();
		delete[] Points[i];
		delete[] nPoints[i];

		Points[i] = new float**[nFibers[i]];
		nPoints[i] = new uint16_t[nFibers[i]];

		for(uint32_t j=0; j<nFibers[i]; j++) {
			Points[i][j] = Bundles[i][j];
			nPoints[i][j] = new_nBundles[i][j];
		}
	}

	// ================ Obtiene la cantidad de intervalos por eje ==================
	const float mins[3] = {minx, miny, minz};
	const float maxs[3] = {maxx, maxy, maxz};
	uint16_t counts[3] = {0, 0, 0};

	for(uint16_t i=0; i<3; i++) {
		float ini = mins[i];
		while(ini < maxs[i]) {
			counts[i]++;
			ini += step;
		}
	}

	// ====== Generación de intervalos (coordenadas de los vértices de cada cubo) ===================
	float** index = new float*[3];

	for(uint8_t i=0; i<3; i++) {
		index[i] = new float[counts[i] + 1];

		#pragma omp simd
		for(uint16_t j=0; j<counts[i] + 1; j++) {
			index[i][j] = mins[i] + j * step;
		}
	}

	std::vector<std::vector<float>> centroids(n_polygons*pow(4,N), std::vector<float>(3)); // centroide de los triángulos

	// ======== Subdivide los triángulos del mallado N veces ==========
	#pragma omp parallel for schedule(dynamic)
	for(uint32_t i=0; i<n_polygons; i++) {
		float ***triangles = new float**[1];
		triangles[0] = new float*[3];

		for(uint8_t j=0; j<3; j++)
			triangles[0][j] = vertex[polygons[i][j]];

		float **centroid = triangle_interpolation(triangles, N); // obtiene los centroides de los triángulos

		for(uint16_t j=0; j<pow(4,N); j++) {
			for(uint8_t k=0; k<3; k++)
				centroids[i*pow(4,N) + j][k] = centroid[j][k];
			delete[] centroid[j];
		}
		delete[] centroid;
	}

	// =========== Se definen las variables que almacenarán los centroides en cada cubo ===========

	std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> almacen;
	std::vector<std::vector<std::vector<std::vector<uint32_t>>>> centroidIndex;
	bool ***cubeNotEmpty = new bool**[counts[0]];

	almacen.resize(counts[0]);
	centroidIndex.resize(counts[0]);

	for(uint16_t ix = 0; ix < counts[0]; ix++) {

		almacen[ix].resize(counts[1]);
		centroidIndex[ix].resize(counts[1]);
		cubeNotEmpty[ix] = new bool*[counts[1]];

		for(uint16_t iy = 0; iy < counts[1]; iy++) {

			almacen[ix][iy].resize(counts[2]);
			centroidIndex[ix][iy].resize(counts[2]);
			cubeNotEmpty[ix][iy] = new bool[counts[2]];

			for(uint16_t iz = 0; iz < counts[2]; iz++) {
				cubeNotEmpty[ix][iy][iz] = false;
			}
		}
	}

	// =========== En cada cubo se almacena una lista de centroides, su índice y se encuentra o no vacío =============

	for(uint32_t i=0; i<n_polygons*pow(4,N); i++) {

		uint32_t I[3];

		#pragma omp simd
		for(uint16_t j=0; j<3; j++)
			I[j] = (centroids[i][j] - index[j][0]) / step;

		almacen[I[0]][I[1]][I[2]].emplace_back(centroids[i]);
		centroidIndex[I[0]][I[1]][I[2]].emplace_back(uint32_t(i / pow(4,N)));
		cubeNotEmpty[I[0]][I[1]][I[2]] = true;
	}

	// ====================== Calcula la intersección con cada fascículo =============================

	InTri.resize(nBundles);
	FnTri.resize(nBundles);
	InPoints.resize(nBundles);
	FnPoints.resize(nBundles);
	fib_index.resize(nBundles);

	for(uint16_t i=0; i<nBundles; i++) {

		std::cout << "Bundle: " << i << "/" << nBundles;
		std::cout << ", Num. Fibers: " << nFibers[i] << std::endl;

		std::vector<uint32_t> listFibInd; // lista de los índices de las fibras
		std::vector<std::vector<uint32_t>> listTri; // lista de los índices de los triángulos
		std::vector<std::vector<std::vector<float>>> listPtInt; // lista con los puntos de intersección (inicial y final)

		#pragma omp parallel for schedule(dynamic)
		for(uint32_t j=0; j<nFibers[i]; j++) {
			bool findInt; // find intersection (encuentra intersección)
			uint32_t InT, FnT; // initial and final triangle (triángulo inicial y final)
			float *InPtInt, *FnPtInt; // Initial and Final point intersection (punto de intersección inicial y final)

			findInt = getMeshAndFiberIntersection(Points[i][j], nPoints[i][j], nPtsLine, N, npbp, index, step, cubeNotEmpty,
				 						    centroidIndex, almacen, vertex, polygons, InT, FnT, InPtInt, FnPtInt);

			// si no existe intersección, continúa con la siguiente fibra
			if(!findInt)
				continue;

			#pragma omp critical
			listFibInd.emplace_back(j);
			#pragma omp critical
			listTri.emplace_back((std::vector<uint32_t>){InT, FnT});
			#pragma omp critical
			listPtInt.emplace_back((std::vector<std::vector<float>>){{InPtInt[0], InPtInt[1], InPtInt[2]}, {FnPtInt[0], FnPtInt[1], FnPtInt[2]}});

			delete[] InPtInt;
			delete[] FnPtInt;
		}

		// ======== Almacena los datos de intersección ordenadamente ========
		for(uint32_t j=0; j<listTri.size(); j++) {
			InTri[i].emplace_back(listTri[j][0]);
			FnTri[i].emplace_back(listTri[j][1]);
			InPoints[i].emplace_back(listPtInt[j][0]);
			FnPoints[i].emplace_back(listPtInt[j][1]);
		}
		fib_index[i] = listFibInd;
	}

	// ======== Se libera el espacio reservado en memoria =========

	for(uint8_t i=0; i<3; i++)
		delete[] index[i];
	delete[] index;

	for(uint16_t i=0; i<counts[0]; i++) {
		for(uint16_t j=0; j<counts[1]; j++) {
			delete[] cubeNotEmpty[i][j];
		}
		delete[] cubeNotEmpty[i];
	}
	delete[] cubeNotEmpty;

	for(uint16_t i=0; i<nBundles; i++) {
		for(uint32_t j=0; j<nFibers[i]; j++) {
			for(uint16_t k=0; k<nPoints[i][j]; k++)
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
