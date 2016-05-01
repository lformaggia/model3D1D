/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
       Copyright (C) 2015 Domenico Notaro, 2014 Laura Cattaneo
======================================================================*/
/*! 
  @file   mesh3d.hpp
  @author Laura Cattaneo <laura.cattaneo1@polimi.it>
  @date   April 2015.
  @brief  Miscelleanous handlers for 1D network mesh.
 */

#ifndef M3D1D_MESH_3D_HPP_
#define M3D1D_MESH_3D_HPP_

#include <iostream>
#include <iomanip>
#include <fstream>


namespace getfem {


/*! 
	Import the 3D mesh files produced by GMSH (.msh format).
	\ingroup geom
 */
//! \note It works for TETRAHEDRAL meshes only!
static 
void 
import_msh_file(std::string filename, mesh& m) 
{
	m.clear();
	std::ifstream f(filename.c_str());
	GMM_ASSERT1(f.good(), "Can't open file " << filename);

	//char tmp[1024];
	
	bgeot::read_until(f, "$Nodes");

	size_type nb_node; 

	f >> nb_node;

	dal::dynamic_tree_sorted<size_type> msh_node_2_getfem_node;

	for (size_type node_cnt=0; node_cnt < nb_node; ++node_cnt) {
		size_type node_id;
		base_node n(3); n[0]=0.0; n[1]=0.1; n[2]=1e30;
		f >> node_id >> n[0] >> n[1] >> n[2];
		msh_node_2_getfem_node.add_to_index(m.add_point(n), node_id);
	}
	bgeot::read_until(f, "$EndNodes");
	bgeot::read_until(f, "$Elements");

	size_type nb_cv;
	f >> nb_cv;

	std::vector<size_type> cv_nodes;

	for (size_type cv=0; cv < nb_cv; ++cv) {
		size_type cv_id, cv_type, cv_region_id, cv_elm_id, cv_partition, cv_nb_nodes ;
		f >> cv_id >> cv_type >> cv_region_id >> cv_elm_id >> cv_partition;

		if (cv_type == size_type(4)) cv_nb_nodes = size_type(4);
		else GMM_ASSERT1(0, "No thetraedra mesh");
	
		cv_id--; 
		cv_nodes.resize(cv_nb_nodes);

		for (size_type i=0; i < cv_nb_nodes; ++i) {
			size_type j;
			f >> j;
			cv_nodes[i] = msh_node_2_getfem_node.search(j);
			GMM_ASSERT1(cv_nodes[i] != size_type(-1),
				"Invalid node ID " << j << " in gmsh convex " << cv_id);
		}

		m.add_simplex(3,cv_nodes.begin());

	}
	bgeot::read_until(f, "$EndElements");
	f.close();

}


} /* end of namespace */

#endif
