#ifndef EXPORT_BACKGROUND_STRESS_DIRECT_H
#define EXPORT_BACKGROUND_STRESS_DIRECT_H

#include <vector>
#include <string>
#include "Mesh.h"

namespace EQSim {

class ExportBackgroundStressDirect {
public:
    ExportBackgroundStressDirect(const Mesh& faultMesh, double E, double nu);
    ~ExportBackgroundStressDirect();

    void ComputeAndExport(const std::vector<double>& slip, const std::string& filename);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace EQSim

#endif // EXPORT_BACKGROUND_STRESS_DIRECT_H