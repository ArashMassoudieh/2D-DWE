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
#
# Results are written as .vti (VTK ImageData XML) by a hand-rolled writer and
# opened in ParaView. Nothing links VTK.
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
    src/vti_writer.hpp \
    src/case_config.hpp

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

# --- visualisation -----------------------------------------------------------
# No VTK link. Results are written as .vti (VTK ImageData XML) by src/vti_writer.hpp
# and opened in ParaView, which is what the .vti container is for. Linking VTK would
# buy an in-process viewer nobody asked for and a fragile module list with it.

OTHER_FILES += DESIGN.md README.md
