//
// Created by Federico Ciardo on 11.10.22. All rights reserved.
//
#include <armadillo>
#include <omp.h>
#include <src/Mesh.h>
#include "FullSpaceElasticity.h"
#include "Utils.h"
#ifndef INC_3DEQSIM_SRC_ASSEMBLEELASTICITYMATRIX_H
#define INC_3DEQSIM_SRC_ASSEMBLEELASTICITYMATRIX_H

namespace EQSim {

arma::mat AssembleElastMat(Mesh &mesh, EQSim::SolidMatrixProperties &matrix_prop);


}
#endif  // INC_3DEQSIM_SRC_ASSEMBLEELASTICITYMATRIX_H
