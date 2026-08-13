# =============================================================================
# 2D-DWE.pro -- qmake project (Qt Creator)
#
# Build with qmake6 on this machine: plain `qmake` is a broken Qt5 symlink.
#   In Qt Creator: Projects -> Build -> make sure the kit's qmake is qmake6.
#   From a shell:  qmake6 2D-DWE.pro && make -j
#
# The solver core has NO third-party dependency. GDAL and VTK are opt-in so a
# broken install of either can never stop the solver building or running:
#
#   qmake6 CONFIG+=with_gdal          terrain / roughness raster I/O
#   qmake6 CONFIG+=with_vtk           interactive viewer (links real VTK)
#
# Result files are written as .vti (VTK ImageData XML) by a hand-rolled writer
# regardless, so ParaView works out of the box with neither flag set.
# =============================================================================

TEMPLATE = app
TARGET   = 2D-DWE
CONFIG  += console c++17
CONFIG  -= app_bundle
CONFIG  -= qt                       # no Qt runtime -- this is a solver

DESTDIR     = $$PWD/bin
OBJECTS_DIR = $$PWD/build/obj
MOC_DIR     = $$PWD/build/moc

# --- warnings and optimisation ----------------------------------------------
QMAKE_CXXFLAGS         += -Wall -Wextra -Wpedantic
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3 -march=native -funroll-loops

# --- OpenMP ------------------------------------------------------------------
# The face and cell loops are the GPU kernels in waiting; keep them threaded.
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS   += -fopenmp
LIBS           += -lgomp

INCLUDEPATH += $$PWD/src

HEADERS += \
    src/params.hpp \
    src/grid.hpp \
    src/fields.hpp \
    src/boundary.hpp \
    src/massbalance.hpp \
    src/scheme_inertial.hpp \
    src/vti_writer.hpp

SOURCES += \
    src/scheme_inertial.cpp \
    src/main.cpp

# --- optional: GDAL ----------------------------------------------------------
with_gdal {
    DEFINES     += DWE_WITH_GDAL
    INCLUDEPATH += /usr/include/gdal
    LIBS        += -lgdal
    HEADERS     += src/io_gdal.hpp
    SOURCES     += src/io_gdal.cpp
    message("GDAL enabled")
}

# --- optional: VTK (interactive viewer) --------------------------------------
# VTK 9.4 is a CMake-first project; under qmake the modules must be named
# explicitly. Only the handful needed to render an ImageData scene are linked.
with_vtk {
    DEFINES     += DWE_WITH_VTK
    VTK_DIR      = /usr/local
    VTK_VER      = 9.4
    INCLUDEPATH += $$VTK_DIR/include/vtk-$$VTK_VER
    LIBS        += -L$$VTK_DIR/lib \
                   -lvtkCommonCore-$$VTK_VER \
                   -lvtkCommonDataModel-$$VTK_VER \
                   -lvtkCommonExecutionModel-$$VTK_VER \
                   -lvtkCommonTransforms-$$VTK_VER \
                   -lvtkCommonMath-$$VTK_VER \
                   -lvtkFiltersCore-$$VTK_VER \
                   -lvtkFiltersGeometry-$$VTK_VER \
                   -lvtkIOXML-$$VTK_VER \
                   -lvtkIOCore-$$VTK_VER \
                   -lvtkRenderingCore-$$VTK_VER \
                   -lvtkRenderingOpenGL2-$$VTK_VER \
                   -lvtkInteractionStyle-$$VTK_VER \
                   -lvtksys-$$VTK_VER
    HEADERS     += src/viewer_vtk.hpp
    SOURCES     += src/viewer_vtk.cpp
    message("VTK enabled from $$VTK_DIR (vtk-$$VTK_VER)")
}

OTHER_FILES += DESIGN.md README.md
