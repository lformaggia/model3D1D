# Project Path (this file path)
PROJECT_DIR = /home/pacs_student/model3D1D

# Installation path
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
# libBGLgeom is stored here (we consider the installation path)
BGLGEOM_DIR        = $(mkBGLgeomHome)/libBGLgeom
# Example on network diffusion is stored here
NETW_DIR           = $(mkBGLgeomHome)/network_diffusion
# Path to BGL includes
BGL_INC_DIR        = $(mkBGLInc)
# Path to BGL Netdiff includes
BGL_NETW_INC       = $(mkBGLgeomHome)/include
# Path to Eigen include
EIGEN_INC          = $(mkEigenInc)



C_MANWORKS_INC_DIR = $(PROJECT_DIR)/CurvedFormulation/include
# Path to Curved M3D1D library
C_MANWORKS_LIB_DIR = $(PROJECT_DIR)/CurvedFormulation/lib
# Path to Linear M3D1D includes
MANWORKS_INC_DIR   = $(PROJECT_DIR)/MixedFormulation/include 
# Path to Linear M3D1D library
MANWORKS_LIB_DIR   = $(PROJECT_DIR)/MixedFormulation/lib
# Path to Graph Generator include
GRAPH_INC_DIR      = $(PROJECT_DIR)/GraphInterface/include
# Path to Graph Generator library
GRAPH_LIB_DIR      = $(PROJECT_DIR)/GraphInterface/lib



# Specific Graph Generator Include list
GRAPH_INCLUDE = -I$(GRAPH_INC_DIR) -I$(NETW_DIR)/include -I$(BGL_NETW_INC) -I$(BGLGEOM_DIR)/include -I$(EIGEN_INC) -I$(GETFEM_INC_DIR) -I$(C_MANWORKS_INC_DIR) -I$(MANWORKS_INC_DIR)
# Specific Graph Generator Library list (for test case)
GRAPH_LIBRARY = -L$(GETFEM_LIB_DIR) -lgetfem -L$(MANWORKS_LIB_DIR) -lproblem3d1d -L$(C_MANWORKS_LIB_DIR) -lc_problem3d1d
# Specific Linear M3D1D Include list
MANWORKS_INCLUDE = -I$(MANWORKS_INC_DIR) -I$(GETFEM_INC_DIR)
# Specific Linear M3D1D Library list
MANWORKS_LIBRARY = -L$(GETFEM_LIB_DIR) -lgetfem
# Specific Curved M3D1D Include list
C_MANWORKS_INCLUDE = -I$(MANWORKS_INC_DIR) -I$(GETFEM_INC_DIR) -I$(C_MANWORKS_INC_DIR)
# Specific Curved M3D1D Library list
C_MANWORKS_LIBRARY =  -lgetfem -lproblem3d1d
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
