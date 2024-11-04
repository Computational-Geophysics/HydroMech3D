//
// Created by Federico Ciardo on 11.08.21. All rights reserved.
//

// Inclusion from IL library
#include <il/Array.h>
#include <il/Array2D.h>

// Inclusion from the project
#include "ElementData.h"
#include "Mesh.h"

#ifndef INC_3DEQSIM_SRC_UTILS_H
#define INC_3DEQSIM_SRC_UTILS_H

namespace EQSim {

inline il::Array2D<double> RotationMatrix3D(il::Array<double> &normal_vector,
                                            double &theta);

inline double euclideanDistance(double x1, double x2, double x3, double y1,
                                double y2, double y3);

template <typename T>
inline il::Array2D<T> position_2d_array(const il::Array2D<T> &arr2D, T seek);

template <class T>
inline il::Array<T> deleteDuplicates(const il::Array<T> &arr);

il::Array<double> CalculatePressureLaplacianUniformMeshOnly(
    il::Array2D<il::int_t> &neigh_elts, EQSim::Mesh &Mesh,
    const il::Array<double> &my_vector);

}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_UTILS_H
