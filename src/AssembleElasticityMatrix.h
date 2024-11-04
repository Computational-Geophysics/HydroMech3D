//
// Created by Federico Ciardo on 11.10.22. All rights reserved.
//

#include <il/il/linearAlgebra/dense/blas/dot.h>
#include <il/Array2D.h>
#include <src/Mesh.h>
#include "FullSpaceElasticity.h"

#ifndef INC_3DEQSIM_SRC_ASSEMBLEELASTICITYMATRIX_H
#define INC_3DEQSIM_SRC_ASSEMBLEELASTICITYMATRIX_H

namespace EQSim {

il::Array2D<double> AssembleElastMat(Mesh &mesh, EQSim::SolidMatrixProperties &matrix_prop);

il::Array2D<double> RotationMatrix3D(il::Array<double> &normal_vector,
                                     double &theta);

}
#endif  // INC_3DEQSIM_SRC_ASSEMBLEELASTICITYMATRIX_H
