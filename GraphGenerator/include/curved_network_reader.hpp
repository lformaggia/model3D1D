/* -*- c++ -*- (enableMbars emacs c++ mode) */ 
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Curved Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2017 Giorgio Raimondi
======================================================================*/
/*! 
  @file   curved_network_reader.hpp
  @author Giorgio Raimondi <giorgio3.raimondi@mail.polimi.it>
  @date   May 2017.
  @brief  Function definition for Reading generic ASCII format for bspline_geometry.
 */    
 /*! @defgroup graph reader routines */
#ifndef HH_CURV_NET_READ_HH
#define HH_CURV_NET_READ_HH

#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <boost/graph/adjacency_list.hpp>

#include "netdiff_graph_properties.hpp"
#include "reader_netdiff.hpp"
#include "spline_reader.hpp"
#include "writer_pts.hpp"
#include <getfem/getfem_export.h>   //For ftool::md_param


using namespace boost;
using namespace BGLgeom;
using namespace std;


namespace NetDiff{

/*!
	This function reads branch by branch the input value for generate the graph, 
	construct with the bspline_geometry each branch and return a pts file with 
	all the parameters needed by the c_problem class in order to compute the flow
	in the network.
	This fuction need as input the name of the input.param file from which is 
	possible to import some information for the construction of the file.

	In the input.param file is mandatory to specify the following information:

		1) The name and path of the input ascii file.
			example: INPUT_ASCII="./input"

		2) The name and the location for the output pts file.
			example: OUTPUT_PTS="../src/graph.pts"

	More over is possible to choose which method use for generate the bspline.
	If you set INTERP=1, the function will generate bspline using interpolation
	method.
	Otherwise, even without specify INTERP, is used the approximation method.



	The  ASCII Input file has to be of the following format:
	
		Each branch is written in the same format and consecutively to the previous one.
		They need:
		1) Boundaries conditions for the source and the target point.
		2) The number of point insert for the branch.
		3) The size of the mesh we want in the output for that branch.
		4) The source of the branch.
		5) The target of the branch.
		6) All the other nodes in the correct order.

	Examples:
	
		DIR 1.0 INT 0.0 		(first branch)
		3 100					
		0.0 0.0 0.0				(sorce node)
		1.0 1.0 1.0				(target node)
		0.5 0.5 0.5				(internal node)
		INT 0.0 DIR 2.0 		(second branch)
		5 50
		1.0 1.0 1.0
		3.0 3.0 3.0
		1.5 1.6 1.5
		2.0 2.0 2.0
		2.4 2.3 2.4
*/
void curved_reader(int argc, char *argv[]);


}
#endif