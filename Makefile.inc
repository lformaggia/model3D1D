# Project Path (this file path)
PROJECT_DIR = /home/pacs_student/model3D1D

# Path to BGLgeom Library
export mkBGLgeomHome = /media/More_SPACE/BGLgeom

# Path to Boost Graph Library
export mkBGLInc = /media/More_SPACE/boost_1_61_0/boost

# Path to GetFEM main folder
export mkGetFEMHome = /opt/getfem

# Path to Eigen includes
export mkEigenInc = /u/geo2/sw/Packages/libs/eigen/3/include/eigen3



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
GRAPH_INCLUDE = -I$(GRAPH_INC_DIR) -I$(BGLGEOM_INC_DIR) -I$(EIGEN_INC) -I$(GETFEM_INC_DIR) -I$(C_MANWORKS_INC_DIR) -I$(MANWORKS_INC_DIR) -I$(BGL_NETW_INC)
# Specific Graph Generator Library list (for test case)
GRAPH_LIBRARY = -lBGLgeom -lc_problem3d1d -lproblem3d1d -lgetfem
# Specific Graph Generator Library directory
GRAPH_LIBRARY_DIR = -L$(BGLGEOM_LIB_DIR) -L$(MANWORKS_LIB_DIR) -L$(GETFEM_LIB_DIR) -L$(C_MANWORKS_LIB_DIR)
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
