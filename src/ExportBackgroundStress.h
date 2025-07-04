#ifndef EXPORT_BACKGROUND_STRESS_H
#define EXPORT_BACKGROUND_STRESS_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "io/bigwham_io.h"                // BigWham IO interface (H-matrix builder)

// Forward declare BigWham IO class to avoid heavy includes in header
//namespace bigwham { class BigWhamIO; }

namespace EQSim {

class ExportBackgroundStress {
public:
    // Constructor takes fault mesh and material properties. 
    // snapshotIndex selects which preset snapshot (0 for top plane by default).
    ExportBackgroundStress(const Mesh &faultMesh, double E, double nu);
    ~ExportBackgroundStress();
    // Compute stress from current slip and export to CSV file.
    void ComputeAndExport(const std::vector<double> &slip, const std::string &filename);
    arma::vec multiply(const arma::vec& x) const;

private:
    // Grid dimensions for observation points (e.g. NX x NY plane at fixed k index)
    int n_openMP_threads_;

    std::string kernel_name_;

   
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace EQSim

#endif
