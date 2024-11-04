//
// Created by Federico Ciardo on 11.08.21. All rights reserved.
//

#include "Utils.h"

namespace EQSim {

inline il::Array2D<double> RotationMatrix3D(il::Array<double> &normal_vector,
                                            double &theta) {
  il::Array2D<double> R{3, 3, 0.};

  double n1, n2, n3, n1square, n2square, n3square;
  n1 = normal_vector[0];
  n2 = normal_vector[1];
  n3 = normal_vector[2];
  n1square = n1 * n1;
  n2square = n2 * n2;
  n3square = n3 * n3;

  R(0, 0) = cos(theta) + (n1square * (1 - cos(theta)));
  R(0, 1) = (n1 * n2 * (1 - cos(theta))) - (n3 * sin(theta));
  R(0, 2) = (n1 * n3 * (1 - cos(theta))) + (n2 * sin(theta));

  R(1, 0) = (n1 * n2 * (1 - cos(theta))) + n3 * sin(theta);
  R(1, 1) = cos(theta) + (n2square * (1 - cos(theta)));
  R(1, 2) = (n2 * n3 * (1 - cos(theta))) - (n1 * sin(theta));

  R(2, 0) = (n1 * n3 * (1 - cos(theta))) - (n2 * sin(theta));
  R(2, 1) = (n2 * n3 * (1 - cos(theta))) + (n1 * sin(theta));
  R(2, 2) = cos(theta) + (n3square * (1 - cos(theta)));

  return R;
}

/////////
// This function calculates the 3D euclidean distance between two points;
// Point 1=(x1 = x1, x2 = y1, x3 = z1) -- Point 2=(y1 = x2, y2 = y2, y3 = z2)
inline double euclideanDistance(double x1, double x2, double x3, double y1,
                                double y2, double y3) {
  double res = 0.;

  res = sqrt(pow(x1 - y1, 2) + pow(x2 - y2, 2) + pow(x3 - y3, 2));

  return res;
}

/////////
template <typename T>
inline il::Array2D<T> position_2d_array(const il::Array2D<T> &arr2D, T seek) {
  il::int_t arr2D_size1 = arr2D.size(1);
  il::int_t arr2D_size0 = arr2D.size(0);
  il::Array2D<T> M{arr2D_size1 * arr2D_size0, 2, -1};
  il::int_t k = 0;

  for (il::int_t i = 0; i < arr2D_size0; ++i) {
    for (il::int_t j = 0; j < arr2D_size1; ++j) {
      if (arr2D(i, j) == seek) {
        M(k, 0) = i;
        M(k, 1) = j;
        k = k + 1;
      }
    }
  }

  il::Array2D<T> outp{k, 2, 0};

  for (il::int_t l = 0; l < k; ++l) {
    for (il::int_t k2 = 0; k2 < 2; ++k2) {
      outp(l, k2) = M(l, k2);
    }
  }

  return outp;
}

/////////
template <class T>
inline il::Array<T> deleteDuplicates(const il::Array<T> &arr) {
  il::Array<T> res{};

  for (il::int_t i = 0; i < arr.size(); ++i) {
    bool already_there = false;
    for (il::int_t j = 0; j < res.size(); ++j) {
      if (arr[i] == res[j]) {
        already_there = true;
      }
    }
    if (!already_there) {
      res.Append(arr[i]);
    }
  }

  return res;
}

////////
inline il::Array<double> CalculatePressureLaplacianUniformMeshOnly(
    il::Array2D<il::int_t> &neigh_elts, EQSim::Mesh &Mesh,
    const il::Array<double> &my_vector) {
  il::Array<double> res{my_vector.size(), 0.};
  il::Array<double> resX{my_vector.size(), 0.};
  il::Array<double> resY{my_vector.size(), 0.};
  il::Array<il::int_t> neighElts{};
  il::int_t t;
  EQSim::ElementData elemI;
  il::int_t N_eltsX, N_eltsY;

  for (il::int_t I = 0; I < neigh_elts.size(0); ++I) {
    elemI = Mesh.getElementData(I);

    // Find the real neighbour elements, i.e. the ones different than -1.
    t = 0;
    for (il::int_t J = 0; J < neigh_elts.size(1); ++J) {
      if (neigh_elts(I, J) != -1) {
        neighElts.Resize(t + 1);
        neighElts[t] = neigh_elts(I, J);
        ++t;
      }
    }

    N_eltsX = 0;
    N_eltsY = 0;
    // Loop over the real neighbour elements
    for (il::int_t J = 0; J < neighElts.size(); ++J) {
      if ((Mesh.getCentroids(neighElts[J], 1) == Mesh.getCentroids(I, 1)) &&
          (neighElts[J] != I)) {
        resX[I] += my_vector[neighElts[J]];
        ++N_eltsX;
      }

      if ((Mesh.getCentroids(neighElts[J], 0) == Mesh.getCentroids(I, 0)) &&
          (neighElts[J] != I)) {
        resY[I] += my_vector[neighElts[J]];
        ++N_eltsY;
      }
    }

    if (N_eltsX == 2) {
      resX[I] -= 2. * my_vector[I];
      resX[I] = (resX[I] / (pow((2. * elemI.getA()), 2)));
    } else {
      resX[I] -= my_vector[I];
      resX[I] = (resX[I] / (pow((2. * elemI.getA()), 2)));
    }

    if (N_eltsY == 2) {
      resY[I] -= 2. * my_vector[I];
      resY[I] = (resY[I] / (pow((2. * elemI.getB()), 2)));
    } else {
      resY[I] -= my_vector[I];
      resY[I] = (resY[I] / (pow((2. * elemI.getB()), 2)));
    }

    res[I] = resX[I] + resY[I];
  }

  return res;
}

}  // namespace EQSim