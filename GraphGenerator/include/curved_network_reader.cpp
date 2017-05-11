/* -*- c++ -*- (enableMbars emacs c++ mode) */ 
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Curved Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2017 Giorgio Raimondi
======================================================================*/
/*! 
  @file   curved_network_reader.cpp
  @author Giorgio Raimondi <giorgio3.raimondi@mail.polimi.it>
  @date   May 2017.
  @brief  Function declaration for Reading generic ASCII format for bspline_geometry.
 */    
 /*! @defgroup graph reader routines */
#include "curved_network_reader.hpp"

namespace NetDiff{



void curved_reader(int argc, char *argv[]){

	using Graph = adjacency_list<vecS, vecS, directedS, Vertex_prop, B_Edge_prop>;

	ftool::md_param PARAM;
 	PARAM.read_command_line(argc, argv);	
 	
	std::string in_filename=PARAM.string_value("INPUT_ASCII","Pass the input file's path for the graph initialization");
	bool Interpolation=PARAM.int_value("INTERP");
	BSP_type _type;
	if(Interpolation==1)
		_type=BSP_type::Interp;
	else
		_type=BSP_type::Approx;

	Graph G;
	B_reader_netdiff R(in_filename);
	
	// Utilities to read the data and build the graph
	Vertex_prop src_prop, tgt_prop;
	//Edge_prop e_prop;
	Edge_desc<Graph> e;
	Vertex_desc<Graph> src, tgt;
	vect_pts CP;
	vector<unsigned int> Mesh_dim;
	unsigned int count = 0;
	
	// Reading the file
	while(!R.is_eof()){
		// Reading data
		R.get_data();
		// Creating edge and setting it up
		src_prop = R.get_source_data();
		tgt_prop = R.get_target_data();
		src = new_vertex(src_prop, G);
		tgt = new_vertex(tgt_prop, G);
		CP = R.get_control_points();
		Mesh_dim.push_back(R.get_mesh_dimension());
		e = new_bspline_edge<Graph,3>(src, tgt,CP,_type, G);
		G[e].index = count;
		count++;
	}
	// Creating a mesh on every edge
	Edge_iter<Graph> e_it, e_end;
	int edge_index;
	for(std::tie(e_it, e_end) = edges(G); e_it != e_end; ++e_it){
		edge_index=G[*e_it].index;
		G[*e_it].make_uniform_mesh(Mesh_dim[edge_index]);
	}
	
	// Writing on a pts output
	std::string out_pts_filename=PARAM.string_value("OUTPUT_PTS","Pass the output file name for the graph");
	writer_pts<Graph,3> W(out_pts_filename);
	W.export_pts(G);	
}



}