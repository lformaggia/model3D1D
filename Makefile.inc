# Project Path (this file path)
PROJECT_DIR = /home/pacs_student/model3D1D

# Path to BGLgeom Library
export mkBGLgeomHome = /media/More_SPACE/Qua/Pacs_BGLgeom_Ilaria_Mattia

# Path to Boost Graph Library
export mkBGLInc = /media/More_SPACE/boost_1_61_0/boost

# Path to GetFEM main folder
export mkGetFEMHome = /opt/getfem




## EDIT THESE ONLY IF YOU DO NOT USE MODULE SYSTEM OF MOX CLUSTER ##
# Path to VTK includes
export mkVtkInc   = /u/geo2/sw/Toolchains/gcc+system/4.9/Packages/libs/vtk/5/include/vtk-5.10
# Path to VTK libraries
export mkVtkLib   = /u/geo2/sw/Toolchains/gcc+system/4.9/Packages/libs/vtk/5/lib/vtk-5.10
# Path to Eigen includes
export mkEigenInc = /u/geo2/sw/Packages/libs/eigen/3/include/eigen3


### NOW DO NOT EDIT ANYTHING ###


# Path to GetFEM library
GETFEM_LIB_DIR     = $(mkGetFEMHome)/lib
# Path to GetFEM includes
GETFEM_INC_DIR     = $(mkGetFEMHome)/include
# libBGLgeom include directory
BGLGEOM_INC_DIR    = $(mkBGLgeomHome)/libBGLgeom/include
# libBGLgeom library directory
BGLGEOM_LIB_DIR    = $(mkBGLgeomHome)/libBGLgeom/lib
# Path to BGL includes
BGL_INC_DIR        = $(mkBGLInc)
# Path to BGL Netdiff includes
BGL_NETW_INC       = $(mkBGLgeomHome)/network_diffusion/include
# Path to Eigen include
EIGEN_INC          = $(mkEigenInc)
# Path to VTK includes
VTK_INC_DIR      = $(mkVtkInc)
# Path to VTK library
VTK_LIB_DIR      = $(mkVtkLib)

# Parsing vtk libraries. In the modules are present also libraries for Python. We filter them out
VTK_LIBS1 = $(wildcard $(VTK_LIB_DIR)/libvtk*.so.5.10.1) #collecting all libraries
VTK_LIBS2 = $(notdir $(VTK_LIBS1))
VTK_LIBS3 = $(basename $(basename $(basename $(basename $(VTK_LIBS2))))) #Discarding .so.5.10.1
VTK_LIBS4 = $(foreach lib, $(VTK_LIBS3), $(if $(findstring Python, $(lib)),,$(lib)))
VTK_LIBS  = $(VTK_LIBS4:lib%=%)



# Path to Curved M3D1D include
C_MANWORKS_INC_DIR = $(PROJECT_DIR)/CurvedFormulation/include
# Path to Curved M3D1D library
C_MANWORKS_LIB_DIR = $(PROJECT_DIR)/CurvedFormulation/lib
# Path to Linear M3D1D includes
MANWORKS_INC_DIR   = $(PROJECT_DIR)/MixedFormulation/include 
# Path to Linear M3D1D library
MANWORKS_LIB_DIR   = $(PROJECT_DIR)/MixedFormulation/lib
# Path to Graph Generator include
GRAPH_INC_DIR      = $(PROJECT_DIR)/GraphGenerator/include
# Path to Graph Generator library
GRAPH_LIB_DIR      = $(PROJECT_DIR)/GraphGenerator/lib



# Specific Graph Generator Include list
GRAPH_INCLUDE = -I$(GRAPH_INC_DIR) -I$(BGLGEOM_INC_DIR) -I$(EIGEN_INC) -I$(GETFEM_INC_DIR) -I$(C_MANWORKS_INC_DIR) -I$(MANWORKS_INC_DIR) -I$(BGL_NETW_INC) -I$(VTK_INC_DIR)
# Specific Graph Generator Library list (for test case)
GRAPH_LIBRARY = -lBGLgeom -lc_problem3d1d -lproblem3d1d -lgetfem $(foreach lib, $(VTK_LIBS), -l$(lib))
# Specific Graph Generator Library directory
GRAPH_LIBRARY_DIR = -L$(BGLGEOM_LIB_DIR) -L$(MANWORKS_LIB_DIR) -L$(GETFEM_LIB_DIR) -L$(C_MANWORKS_LIB_DIR) -L$(VTK_LIB_DIR)
# Specific Linear M3D1D Include list
MANWORKS_INCLUDE = -I$(MANWORKS_INC_DIR) -I$(GETFEM_INC_DIR) 
# Specific Linear M3D1D Library list
MANWORKS_LIBRARY = -L$(GETFEM_LIB_DIR) -lgetfem
# Specific Curved M3D1D Include list
C_MANWORKS_INCLUDE = -I$(MANWORKS_INC_DIR) -I$(GETFEM_INC_DIR) -I$(C_MANWORKS_INC_DIR)
# Specific Curved M3D1D Library list
C_MANWORKS_LIBRARY =   -lproblem3d1d -lgetfem
# Specific Curved M3D1D Library directory
C_MANWORK_LIBRARY_DIR = -L$(GETFEM_LIB_DIR) -L$(MANWORKS_LIB_DIR)

# Compiler and standard for the compiler
CXX = g++
STDFLAG = -std=c++11
WARNS = -Wall

ifeq ($(DEBUG),yes)
  OPTFLAG =-g
else
  OPTFLAG =-O3 -march=native -DNDEBUG
endif
