//
// Created by Federico Ciardo on 11.08.21. All rights reserved.
//

// Inclusion from IL library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/linearAlgebra.h>

// Inclusion from the project
#include "SolidMatrixProperties.h"
#include "src/StressKernelsP0/StressKernelsDxP0.h"

#ifndef INC_3DEQSIM_SRC_FULLSPACEELASTICITY_H
#define INC_3DEQSIM_SRC_FULLSPACEELASTICITY_H

namespace EQSim {

il::Array2D<double> TractionsDueToDDsOnSingleEltP0(
    double &a, double &b, il::Array<double> &coor_centroid_elt,
    il::Array<double> &shear1_vector, il::Array<double> &shear2_vector,
    il::Array<double> &normal_vector, SolidMatrixProperties Matrix_Prop, il::io_t);

}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_FULLSPACEELASTICITY_H
