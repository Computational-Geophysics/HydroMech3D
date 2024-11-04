//
// Created by Federico Ciardo on 28.07.21. All rights reserved.
//

// Inclusion from the project
#include "Mesh.h"

#include "Utils.cpp"

namespace EQSim {

EQSim::ElementData Mesh::getElementData(il::int_t ne) const {
  // Get node coordinates of element ne
  // CoorEltNe -> 4x3 array 2D with columns corresponding to x, y, z of the node
  il::Array2D<double> CoorEltNe{4, 3, 0.};
  CoorEltNe(0, 0) = coordinates_(nodes_connectivity_(ne, 0), 0);
  CoorEltNe(0, 1) = coordinates_(nodes_connectivity_(ne, 0), 1);
  CoorEltNe(0, 2) = coordinates_(nodes_connectivity_(ne, 0), 2);

  CoorEltNe(1, 0) = coordinates_(nodes_connectivity_(ne, 1), 0);
  CoorEltNe(1, 1) = coordinates_(nodes_connectivity_(ne, 1), 1);
  CoorEltNe(1, 2) = coordinates_(nodes_connectivity_(ne, 1), 2);

  CoorEltNe(2, 0) = coordinates_(nodes_connectivity_(ne, 2), 0);
  CoorEltNe(2, 1) = coordinates_(nodes_connectivity_(ne, 2), 1);
  CoorEltNe(2, 2) = coordinates_(nodes_connectivity_(ne, 2), 2);

  CoorEltNe(3, 0) = coordinates_(nodes_connectivity_(ne, 3), 0);
  CoorEltNe(3, 1) = coordinates_(nodes_connectivity_(ne, 3), 1);
  CoorEltNe(3, 2) = coordinates_(nodes_connectivity_(ne, 3), 2);

  // Get centroid coordinate of element ne
  il::Array<double> CentroidNe{3, 0.};
  for (int I = 0; I < CentroidNe.size(); ++I) {
    CentroidNe[I] = centroids_(ne, I);
  }

  EQSim::ElementData ElementData(CoorEltNe, CentroidNe);
  return ElementData;
}

il::Array2D<il::int_t> Mesh::getNeighbourElements_NonUniformMesh() const {
  // Get neighbour elements, i.e. for each element find its surrounding elements
  // Max number of elements surrounding one element is 9!
  // Not very efficient, but it can be used for non-uniform meshes
  il::Array2D<il::int_t> neigh_elts{
      nodes_connectivity_.size(0), 9,
      -1};  // The maximum element patch is composed of nine elements!

  // Aux int
  il::int_t K1 = 0;

  // Get for each element, its surrounding elements indexes
  // Loop over all the elements
  //#pragma omp parallel
  //  {
  //#pragma omp for
  for (il::int_t I = 0; I < nodes_connectivity_.size(0); ++I) {
    // Auxiliary arrays
    il::Array2D<il::int_t> aux{};
    il::Array<il::int_t> aux2{};
    il::Array<il::int_t> aux3{};
    il::int_t K;

    // Loop over all the nodes of a given element, i.e. in total they are 4
    //#pragma omp parallel
    //      {
    //#pragma omp for reduction(+ : K1)
    for (il::int_t J = 0; J < nodes_connectivity_.size(1); ++J) {
      aux = EQSim::position_2d_array(nodes_connectivity_,
                                     nodes_connectivity_(I, J));

      for (K = 0; K < aux.size(0); ++K) {
        aux2.Resize(K1 + aux.size(0));
        aux2[K1 + K] = aux(K, 0);
      }
      K1 += K;
    }
    //}

    // Make sure to delete duplicates and sort the result
    aux3 = EQSim::deleteDuplicates(aux2);
    std::sort(aux3.begin(), aux3.end());

    // Fill the neigh_elts 2D array
    // Note that all the negative entries must be discarded
    for (il::int_t j = 0; j < aux3.size(); ++j) {
      neigh_elts(I, j) = aux3[j];
    }

    // Restart from zero
    K1 = 0;
  }
  //}

  return neigh_elts;
}

il::Array2D<il::int_t> Mesh::getNeighbourElements_UniformMesh(
    EQSim::Mesh &Mesh) {
  // Get neighbour elements, i.e. for each element find its surrounding elements
  // Max number of elements surrounding one element is 9!
  // Very efficient, but it can be used for uniform meshes only!!

  il::Array<il::int_t> boundary_elts{};
  il::Array<il::int_t> aux{};
  il::int_t K;
  il::int_t kk = 0;
  il::Array2D<il::int_t> res{centroids_.size(0), 9, -1};

  // Get element data of zero element -> this is why it works only for uniform
  // meshes such that hx = hy
  ElementData elm_data_0 = Mesh.getElementData(0);

  double hx = (2. * elm_data_0.getA());  // This is the element size on x
                                         // direction which is equal to hy

  // Loop over all the elements
  for (il::int_t I = 0; I < centroids_.size(0); ++I) {
    aux.Resize(0);
    K = 0.;

      // Loop over all the elements
      for (il::int_t J = 0; J < centroids_.size(0); ++J) {
        // Check Euclidean distance between every element to select neighbour
        // elements
        if ((EQSim::euclideanDistance(centroids_(I, 0), centroids_(I, 1),
                                      centroids_(I, 2), centroids_(J, 0),
                                      centroids_(J, 1),
                                      centroids_(J, 2)) <= 10e-4) ||
            (EQSim::euclideanDistance(centroids_(I, 0), centroids_(I, 1),
                                      centroids_(I, 2), centroids_(J, 0),
                                      centroids_(J, 1), centroids_(J, 2)) -
                 hx <=
             10e-4) ||
            (EQSim::euclideanDistance(centroids_(I, 0), centroids_(I, 1),
                                      centroids_(I, 2), centroids_(J, 0),
                                      centroids_(J, 1), centroids_(J, 2)) -
                 (sqrt(2.) * hx) <=
             10e-4)) {
          aux.Resize(K + 1);
          aux[K] = J;
          K = K + 1;
        }
      }

    std::sort(aux.begin(), aux.end());
    for (il::int_t J1 = 0; J1 < aux.size(); ++J1) {
      res(I, J1) = aux[J1];
    }

    if (aux.size() == 4 || aux.size() == 6) {
      boundary_elts.Resize(kk + 1);
      boundary_elts[kk] = I;
      ++kk;
    }
  }

  boundary_elements_ = boundary_elts;

  return res;
}

il::Array2D<double> Mesh::getInnerCentroids(
    il::Array2D<il::int_t> &neigh_elts) {
  // Get the coordinates of the inner centroids, i.e. the centroids of all the
  // elements that are not on the boundary of the mesh
  il::Array2D<double> InnerCentroids{};

  il::int_t Nelts = neigh_elts.size(0);
  il::int_t MaxNumEltsSurrOneElt = neigh_elts.size(1);

  // Auxiliary integers
  il::int_t counter;
  il::int_t n = 0;

  // Get indexes of all the inner elements / centroids
  il::Array<il::int_t> IndexesInnerCentroids{};

  // For a given element I, count how many elements surround it
  // If the number is equal to 9, then fill aux array with an
  // integer value equal to 1
  for (il::int_t I = 0; I < Nelts; ++I) {
    counter = 0;
    for (il::int_t J = 0; J < MaxNumEltsSurrOneElt; ++J) {
      if (neigh_elts(I, J) >= 0) {
        counter += 1;
      }
    }

    if (counter == 9) {
      InnerCentroids.Resize(n + 1, 3);
      InnerCentroids(n, 0) = centroids_(I, 0);
      InnerCentroids(n, 1) = centroids_(I, 1);
      InnerCentroids(n, 2) = centroids_(I, 2);
      IndexesInnerCentroids.Resize(n + 1);
      IndexesInnerCentroids[n] = I;
      n += 1;
    }
  }

  indexes_inner_centroids_ = IndexesInnerCentroids;

  return InnerCentroids;
}

}  // namespace EQSim
