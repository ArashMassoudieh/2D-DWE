// =============================================================================
// vti_writer.hpp -- VTK ImageData (.vti) output, written directly as XML.
//
// WHY NOT LINK VTK FOR THIS
//   A structured grid maps exactly onto vtkImageData, and the .vti container is
//   plain XML with ASCII arrays. Writing it by hand keeps the solver free of any
//   VTK link dependency, so the core builds and runs anywhere, and the files
//   open natively in ParaView and in any VTK reader.
//
//   The interactive VTK viewer (CONFIG += with_vtk in the .pro) links the real
//   library and reads these same files. Keeping the two apart means a broken VTK
//   install can never stop the solver from running.
//
// ORIENTATION
//   Row 0 of the arrays is the NORTH row (grid.hpp). VTK ImageData indexes +y
//   upward, so rows are emitted bottom-up and the origin is the SOUTH-WEST
//   corner. Get this wrong and every result is mirrored.
// =============================================================================
#pragma once

#include <cstdio>
#include <string>
#include <vector>

#include "grid.hpp"

namespace dwe {

inline bool write_vti(const std::string& path,
                      const Grid& g,
                      const std::vector<double>& depth,
                      const std::vector<double>& bed,
                      const std::vector<double>& wse)
{
    std::FILE* fp = std::fopen(path.c_str(), "w");
    if (!fp) return false;

    const double oy = g.y0 - g.ny * g.dx;   // south edge

    std::fprintf(fp,
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
        "  <ImageData WholeExtent=\"0 %d 0 %d 0 0\" Origin=\"%.6f %.6f 0\" "
        "Spacing=\"%.6f %.6f 1\">\n"
        "    <Piece Extent=\"0 %d 0 %d 0 0\">\n"
        "      <CellData Scalars=\"depth\">\n",
        g.nx, g.ny, g.x0, oy, g.dx, g.dx, g.nx, g.ny);

    auto emit = [&](const char* name, const std::vector<double>& a) {
        std::fprintf(fp, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n",
                     name);
        for (int j = g.ny - 1; j >= 0; --j) {          // bottom-up for VTK
            for (int i = 0; i < g.nx; ++i)
                std::fprintf(fp, "%.4f ", float(a[g.k(i, j)]));
            std::fputc('\n', fp);
        }
        std::fprintf(fp, "        </DataArray>\n");
    };

    emit("depth", depth);
    emit("bed",   bed);
    emit("wse",   wse);

    std::fprintf(fp,
        "      </CellData>\n"
        "    </Piece>\n"
        "  </ImageData>\n"
        "</VTKFile>\n");

    std::fclose(fp);
    return true;
}

} // namespace dwe
