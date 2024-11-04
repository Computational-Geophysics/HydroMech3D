//
// Created by Federico Ciardo on 28.07.21. All rights reserved.
//

// Inclusion from InsideLibrary
#include <il/Array.h>
#include <il/Array2D.h>

// Import from the project
#include "ElementData.h"

// Import from standard library
#include <algorithm>

#ifndef INC_3DEQSIM_SRC_MESH_H
#define INC_3DEQSIM_SRC_MESH_H

namespace EQSim {

class Mesh {
 private:
  il::Array2D<double> coordinates_{};
  il::Array2D<il::int_t> nodes_connectivity_{};
  il::Array2D<il::int_t> edges_connectivity_{};
  il::int_t interpolation_order_ = 0;
  il::Array2D<double> centroids_{};
  il::Array<il::int_t> indexes_inner_centroids_;
  il::Array<il::int_t> boundary_elements_;

 public:
  // Constructors
  Mesh() = default;

  Mesh(const il::Array2D<double> &Coordinates,
       const il::Array2D<il::int_t> &Connectivity,
       const il::int_t &interpolationOrder) {
    IL_EXPECT_FAST(Coordinates.size(0) > 1 && Coordinates.size(1) == 3);
    IL_EXPECT_FAST(interpolationOrder == 0 || interpolationOrder == 4);

    coordinates_ = Coordinates;
    nodes_connectivity_ = Connectivity;
    interpolation_order_ = interpolationOrder;

    il::int_t Nelts = nodes_connectivity_.size(0);
    il::Array2D<double> Centroids{Nelts, 3, 0.};

    for (il::int_t I = 0; I < Nelts; ++I) {
      Centroids(I, 0) = (Coordinates(Connectivity(I, 0), 0) +
                         Coordinates(Connectivity(I, 1), 0) +
                         Coordinates(Connectivity(I, 2), 0) +
                         Coordinates(Connectivity(I, 3), 0)) /
                        4.;
      Centroids(I, 1) = (Coordinates(Connectivity(I, 0), 1) +
                         Coordinates(Connectivity(I, 1), 1) +
                         Coordinates(Connectivity(I, 2), 1) +
                         Coordinates(Connectivity(I, 3), 1)) /
                        4.;
      Centroids(I, 2) = (Coordinates(Connectivity(I, 0), 2) +
                         Coordinates(Connectivity(I, 1), 2) +
                         Coordinates(Connectivity(I, 2), 2) +
                         Coordinates(Connectivity(I, 3), 2)) /
                        4.;
    }
    centroids_ = Centroids;

    il::Array2D<il::int_t> EdgeConnectivity{Nelts, 8, 0};
      for (il::int_t I = 0; I < Nelts; ++I) {
        EdgeConnectivity(I, 0) = Connectivity(I, 0);
        EdgeConnectivity(I, 1) = Connectivity(I, 1);
        EdgeConnectivity(I, 2) = Connectivity(I, 1);
        EdgeConnectivity(I, 3) = Connectivity(I, 2);
        EdgeConnectivity(I, 4) = Connectivity(I, 2);
        EdgeConnectivity(I, 5) = Connectivity(I, 3);
        EdgeConnectivity(I, 6) = Connectivity(I, 3);
        EdgeConnectivity(I, 7) = Connectivity(I, 0);
      }
    edges_connectivity_ = EdgeConnectivity;
  }

  // Getter methods
  il::Array2D<double> getCoordinates() const { return coordinates_; };
  double getCoordinates(il::int_t k, il::int_t i) const {
    return coordinates_(k, i);
  };

  il::Array2D<il::int_t> getNodesConnettivity() const {
    return nodes_connectivity_;
  };
  il::int_t getNodesConnettivity(il::int_t k, il::int_t i) const {
    return nodes_connectivity_(k, i);
  };

  il::Array2D<il::int_t> getEdgesConnettivity() const {
    return edges_connectivity_;
  };
  il::int_t getEdgesConnettivity(il::int_t k, il::int_t i) const {
    return edges_connectivity_(k, i);
  };

  il::Array2D<double> getCentroids() const { return centroids_; };
  double getCentroids(il::int_t k, il::int_t i) const {
    return centroids_(k, i);
  };

  il::Array2D<double> getInnerCentroids(il::Array2D<il::int_t> &neigh_elts);

  il::Array<il::int_t> getIndexesInnerCentroids() const {
    return indexes_inner_centroids_;
  };

  il::int_t getIndexesInnerCentroids(il::int_t centroid_i) const {
    return indexes_inner_centroids_[centroid_i];
  };

  il::int_t getInterpolationOrder() const { return interpolation_order_; };

  il::int_t getNumberOfElts() const { return nodes_connectivity_.size(0); };

  il::int_t getNumberOfNodes() const { return coordinates_.size(0); };

  il::int_t getNumberOfDofs() const { return 3 * nodes_connectivity_.size(0); };

  EQSim::ElementData getElementData(il::int_t ne) const;

  il::Array2D<il::int_t> getNeighbourElements_NonUniformMesh() const;
  il::Array2D<il::int_t> getNeighbourElements_UniformMesh(
      EQSim::Mesh &Mesh);

  il::Array<il::int_t> getBoundaryElements() const {
    return boundary_elements_;
  };

  // Setter methods
  void setCentroids(il::Array2D<double> &centroids) { centroids_ = centroids; };
};

}  // namespace EQSim
#endif  // INC_3DEQSIM_SRC_MESH_H
