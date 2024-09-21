#include "regions.h"

// =========== Calcula los vértices que contienen el borde de la región =============
std::vector<uint32_t> get_edges(const std::vector<uint32_t> &Tri, uint32_t **&polygons) {

	std::multiset<std::set<uint32_t>> edges;
	std::vector<uint32_t> newEdges;

	// ========= obtiene las 3 aristas de un triángulo =========
	for(const uint32_t &i : Tri) {
		const std::set<uint32_t> a = {polygons[i][0], polygons[i][1]};
		const std::set<uint32_t> b = {polygons[i][0], polygons[i][2]};
		const std::set<uint32_t> c = {polygons[i][1], polygons[i][2]};

		edges.emplace(a);
		edges.emplace(b);
		edges.emplace(c);
	}

	// ========= se guardan los vértices de aquellas aristas que nunca se repiten =======
	for(const std::set<uint32_t> &edge : edges) {
		if(edges.count(edge) == 1) {
			newEdges.emplace_back(*edge.begin());
			newEdges.emplace_back(*edge.rbegin());
		}
	}

	// ====== elimina los vértices repetidos =======
	std::sort(newEdges.begin(), newEdges.end());
	newEdges.erase( std::unique( newEdges.begin(), newEdges.end() ), newEdges.end() );
	return newEdges;
}

// ========== Umbraliza dejando solo aquellos triángulos intersectados que se repitan más de 1 vez ========
void Thresholding(const uint8_t &th, std::vector<std::vector<uint32_t>> &Tri) {

	std::vector<std::vector<uint32_t>> newTri;
	newTri.resize(Tri.size());

	for(uint16_t i=0; i<Tri.size(); i++) {

		std::multiset<uint32_t> mymultiset(Tri[i].begin(), Tri[i].end()); // índices de los triángulos de la región
		
		// se calculan los índices de todos los triángulos de la región sin repetición
		std::sort(Tri[i].begin(), Tri[i].end());
		Tri[i].erase( std::unique( Tri[i].begin(), Tri[i].end() ), Tri[i].end() );

		// se evalúa que cada índice se repita "th" veces
		for(uint32_t j=0; j<Tri[i].size(); j++) {

			if(mymultiset.count(Tri[i][j]) > th)
				newTri[i].emplace_back(Tri[i][j]);
		}
	}

	Tri = newTri;
}

// ============ Se realiza una dilatación morfológica con los triángulos de la región ============
void Dilation(const uint8_t &dil, std::vector<std::vector<uint32_t>> &Tri, uint32_t **&polygons, uint32_t &len_polygons) {

	// ===== se dilata la región "dil" veces ========
	for(uint8_t i=0; i<dil; i++) {

		std::vector<std::vector<uint32_t>> newTri;
		newTri.resize(Tri.size()); // conjunto de regiones de un hemisferio

		// ======== se aplica dilatación para cada grupo de triángulo =======
		for(uint16_t j=0; j<Tri.size(); j++) {

			const std::vector<uint32_t> edges = get_edges(Tri[j], polygons); // se obtienen los vértices del borde de la región

			// ====== encuentra mediante búsqueda binaria a los triángulos que coincidan con los vértices de los bordes ======
			#pragma omp parallel for schedule(dynamic)
			for(uint32_t k=0; k<len_polygons; k++) {

				for(uint8_t l=0; l<3; l++) {

					if(std::binary_search(edges.begin(), edges.end(), polygons[k][l])) {
						#pragma omp critical
						newTri[j].emplace_back(k); // almacena dichos triángulos
						break;
					}
				}
			}

			// ======== se fusionan y eliminan los triángulos repetidos al dilatar ==========
			std::vector<uint32_t> temp;
			temp.reserve(Tri[j].size() + newTri[j].size());
			temp.insert(temp.end(), Tri[j].begin(), Tri[j].end());
			temp.insert(temp.end(), newTri[j].begin(), newTri[j].end());

			std::sort(temp.begin(), temp.end());
			temp.erase( std::unique( temp.begin(), temp.end() ), temp.end() );

			newTri[j] = temp;
		}

		Tri = newTri;
	}
}

// ============ Se realiza una erosión morfológica con los triángulos de la región ============
void Erosion(const uint8_t &ero, std::vector<std::vector<uint32_t>> &Tri, uint32_t **&polygons) {

	// ========= se erosiona la región "ero" veces ===========
	for(uint8_t i=0; i<ero; i++) {
	
		std::vector<std::vector<uint32_t>> newTri;
		newTri.resize(Tri.size()); // conjunto de regiones de un hemisferio

		// ======== se aplica erosión para cada grupo de triángulo =======
		for(uint16_t j=0; j<Tri.size(); j++) {

			const std::vector<uint32_t> edges = get_edges(Tri[j], polygons); // se obtienen los vértices del borde de la región

			// encuentra mediante búsqueda binaria a los triángulos que coincidan con los vértices de los bordes
			for(const uint32_t &ind : Tri[j]) {

				for(uint8_t k=0; k<3; k++) {

					if(std::binary_search(edges.begin(), edges.end(), polygons[ind][k]))
						break;

					if(k == 2) newTri[j].emplace_back(ind); // los almacena solamente si no pertenecen al borde
				}
			}		
		}

		Tri = newTri;
	}
}