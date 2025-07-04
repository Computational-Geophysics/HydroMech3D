//
// Created by Federico Ciardo on 11.08.21. All rights reserved.
//

// Inclusion from IL library
#include <armadillo>

// Inclusion from the project
#include "ElementData.h"
#include "Mesh.h"

#ifndef INC_3DEQSIM_SRC_UTILS_H
#define INC_3DEQSIM_SRC_UTILS_H

namespace EQSim {

arma::mat RotationMatrix3D(arma::vec &normal_vector,
                                            double &theta);

double euclideanDistance(double x1, double x2, double x3, double y1,
                                double y2, double y3);

template <typename T>
inline arma::Mat<T> position_2d_array(const arma::Mat<T> &arr2D, T seek);

template <class T>
inline arma::Col<T> deleteDuplicates(const arma::Col<T> &arr);

arma::vec CalculatePressureLaplacianUniformMeshOnly(
    arma::imat &neigh_elts, EQSim::Mesh &Mesh,
    const arma::vec &my_vector, bool use_neumann_bc = false);

template <typename T>
inline arma::Mat<T> position_2d_array(const arma::Mat<T> &arr2D, T seek) {
  arma::uword arr2D_size1 = arr2D.n_cols;
  arma::uword arr2D_size0 = arr2D.n_rows;
  arma::Mat<T> M(arr2D_size1 * arr2D_size0, 2, arma::fill::value(-1));
  arma::uword k = 0;

  for (arma::uword i = 0; i < arr2D_size0; ++i) {
    for (arma::uword j = 0; j < arr2D_size1; ++j) {
      if (arr2D(i, j) == seek) {
        M(k, 0) = i;
        M(k, 1) = j;
        k = k + 1;
      }
    }
  }

  arma::imat outp(k, 2, arma::fill::zeros);

  for (arma::uword l = 0; l < k; ++l) {
    for (arma::uword k2 = 0; k2 < 2; ++k2) {
      outp(l, k2) = M(l, k2);
    }
  }

  return outp;
}

template <class T>
inline arma::Col<T> deleteDuplicates(const arma::Col<T> &arr) {
  arma::Col<T> res;

  for (arma::uword i = 0; i < arr.n_elem; ++i) {
    bool already_there = false;
    for (arma::uword j = 0; j < res.n_elem; ++j) {
      if (arr[i] == res[j]) {
        already_there = true;
      }
    }
    if (!already_there) {
      res.insert_rows(res.n_elem, 1);
      res(res.n_elem - 1) = arr[i];
    }
  }

  return res;
}
}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_UTILS_H
