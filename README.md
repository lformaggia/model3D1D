# model3D1D
3D1D code
# Mixed Finite Element Methods for Coupled 3D/1D Fluid Curved Problems
#### *Politecnico di Milano* (ITALY)

**Author** : Giorgio Raimondi

**Mailto** : <giorgio3.raimondi@mail.polimi.it>

**Date**   : May 2017

-------------------------------------------------------
## How to install and run the program
-------------------------------------------------------
## THE PACKAGE

- 'MixedFormulation/': Directory for the library 'libproblem3d1d'

- 'CurvedFormulation/': Directory for the library 'libc_problem3d1d'

- 'GraphGenerator/': Directory for the library 'libgraphgenerator'

All this directories have the same structure:

- `doc/`     : Code documentation

- `include/` : General include files

- `lib/`     : Main library (to be generated)

- `src/`     : Example sources

- `Doxyfile` : Instruction to build the code documentation

- `Makefile` : Instruction to install the project (see INSTALLATION)



## INSTALLATION
### Prerequisites

#You need the open source finite element library "GetFEM++".

See <http://download.gna.org/getfem/html/homepage>

Version >= 4.2 is preferible


##If you want to install also the graphgenerator library you need to install before:

#"BGLgeom" library.

See <https://github.com/MOX-PolitecnicoMilano/Pacs_BGLgeom_Ilaria_Mattia>

# Boost Graph Library
Version >= 1.61.0 required.

You can download it from here: <https://sourceforge.net/projects/boost/files/>.

You don't need to compile anything, since the BGL is a header only library. 

##### VTK library

Version >= 5.10 required

For more information see: <http://www.vtk.org/download/>

Alternatively, at Mox cluster just load the correspondent module:
```
module load vtk/5

# Eigen

Version >= 3 required.

For more information see: <http://eigen.tuxfamily.org/index.php?title=Main_Page>

You do not have to compile anything, since the Eigen is a header only library.


#Before the installation must also modify the library path in Makefile.inc
Examples:
  # Project Path (this file path)
   PROJECT_DIR = /home/pacs_student/model3D1D

  # Path to BGLgeom
   export mkBGLgeomHome = /media/More_SPACE/BGLgeom

  # Path to Boost Graph Library
   export mkBGLInc = /media/More_SPACE/boost_1_61_0/boost

  # Path to GetFEM main folder
   export mkGetFEMHome = /opt/getfem

  # Path to Eigen includes
   export mkEigenInc = /u/geo2/sw/Packages/libs/eigen/3/include/eigen3


BEWARE: 
Recall to add the library path to LD_LIBRARY_PATH. Example:
```
$ export LD_LIBRARY_PATH=/home/...path/to.../getfem/lib
```

======================

### Installation
Build the whole project with:
``` 
$ make
``` 
It first builds the (static) library "libproblem3d1d" by calling
the Makefile in `MixedFormulation/`:
...
$ make -C MixedFormulation/
...
 which will call
... 
$ make -C include/
``` 
Then, it calls the inner makefiles provided for all examples.

As second it builds the (static) library "libc_problem3d1d" by calling:
...
$ make -C CurvedFormulation/
...
 which will call
... 
$ make -C include/
``` 
Then, it calls the inner makefiles provided for all examples.

As last it builds the (static) library "libgraphgenerator" by calling:
...
$ make -C GraphGenerator/
...
 which will call
... 
$ make -C include/
``` 
Then, it calls the inner makefiles provided for all examples.


If you don't need the GraphGenerator part just type:
...
$ make MANwork
...

which will compile only the "MixedFormulation" and "CurvedFormulation" parts.




BEWARE: 
If you want non-optimized program type:
``` 
$ make DEBUG=yes 
``` 
By defaul DEBUG=no.

The following macro are defined and exported
Recall that any macro may be overled by specifying it when calling 
make. Example: 
``` 
$ make CXXFLAGS+=-DSOMETHING OPTFLAGS=-g
``` 

======================

### Documentation
The documentation is produced by doxygen. The file Doxyfile contains 
the common doxygen configuration for all examples.
Build the code documentation with:
``` 
$ make pdf
``` 
It first fills doc/ with code documentation ($ make doc) and then compile
the .tex files to produce a portable file ($ pdflatex doc/latex/refman.tex).
You can visualize the documentation with
``` 
$ evince doc/latex/refman.pdf
``` 

## MAKE OPTIONS
All examples are provided with a Makefile which accepts the following
options:
-  all       : makes the example
-  clean     : as it says
-  distclean : clean and also deletes temporary file and local doc directory
Being "all" the first target of the makefile, to compile the examples is
sufficient to type make. 
In addition the external Makefile (./Makefile) has the following options:
-  doc       : produces the documentation (html, tex)
-  pdf       : produces a guide in portable format
- library    : build the library from files in include/

## RUN EXAMPLES
To run a specif example, go to the related subdirectory
``` 
$ cd ./MixedFormulation/src/3_Yshaped
``` 
Build the program
``` 
$ make
``` 
Execute the program with specific input
``` 
$ ./M3D1D input.param
``` 
Each program contains the file input.param that specifies 

- Some flags to identify the particular example
  -  TEST_PARAM = 1  # import parameters in dimensionless form
  -  VTK_EXPORT = 1  # export results in vtk format
  -  ...

- The mesh
  - For the 3D mesh you can either provide instruction to build a simple
  regular mesh (TEST_GEOMETRY = 1) or the absolute path to import a mesh
  pre-built with Gmsh (.gmsh)
  - For the 1D mesh specify the path to the file of points (.pts). All
  examples come with a possible pts file

- GetFEM++ descriptors (FEM, ...)

- Problem parameters (dimensional or dimensionless)

- Boundary conditions. You can choose an arbitrary combination of
  Dirichlet-type conditions on pt and/or Robin-type conditions
  on the flux, namely:

  % Faces:   x=0  x=L  y=0  y=L  z=0  z=L

  % BC labels (DIR / MIX)

  BClabel = 'DIR  DIR  DIR  DIR  DIR  DIR'

  % BC values

  BCvalue = '0.0  0.0  0.0  0.0  0.0  0.0'
  
  
BEWARE: All paths in file param must be ABSOLUTE

##  DEV ENVIRONMENT
OS         : Ubuntu 14.04 LTS 64-bit

Processor  : Intel® Core™ i5-2410M CPU @ 2.30GHz × 4 

Compiler   : g++-4.8

GetFEM lib : 5.0

