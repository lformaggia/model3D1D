/* -*- c++ -*- (enableMbars emacs c++ mode) */ 
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Curved Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2017 Giorgio Raimondi
======================================================================*/
/*! 
  @file   spline_reader.hpp
  @author Giorgio Raimondi <giorgio3.raimondi@mail.polimi.it>
  @date   May 2017.
  @brief  Class definition for Reading generic ASCII format for bspline_geometry.
 */    
 /*! @defgroup graph reader routines */
#ifndef HH_NETDIFF_GRAPH_BSPLINE_HH
#define HH_NETDIFF_GRAPH_BSPLINE_HH



#include "base_properties.hpp"
#include "bspline_geometry.hpp"
#include "netdiff_graph_properties.hpp"
#include <functional>

namespace NetDiff{


//! Base properties for the edges (no more properties required)
using B_Edge_prop = BGLgeom::Edge_base_property< BGLgeom::bspline_geometry<3>, 3>;
using points = BGLgeom::point<3>;
using vect_pts = std::vector<points>;

//!Basic Class used to read ASCII file format
class B_reader_netdiff : public BGLgeom::reader_ASCII	<NetDiff::Vertex_prop,
													 NetDiff::B_Edge_prop> {
	public:
		//! Constructor
		B_reader_netdiff(std::string _filename) : BGLgeom::reader_ASCII	<NetDiff::Vertex_prop,
																		 NetDiff::B_Edge_prop>(_filename),
												SRC() , TGT() , BC_src(),BC_tgt(),dim_cp(0),dim_mesh(0),CP() {}
												
		//! Reading form input
		void get_data(){
			this->in_file >> BC_src >> BC_tgt>>dim_cp>>dim_mesh>> SRC >> TGT;
			points PTS;
			CP.push_back(SRC);
			for(unsigned int i=2; i<dim_cp;i++){
				this->in_file>>PTS;
				CP.push_back(PTS);
			}
			CP.push_back(TGT);
		}
		
		//!Returning data on the source
		NetDiff::Vertex_prop get_source_data(){
			return NetDiff::Vertex_prop(SRC, {BC_src});
		}
		
		//! Returning data on the target
		NetDiff::Vertex_prop get_target_data(){
			return NetDiff::Vertex_prop(TGT, {BC_tgt});
		}

		//! Returning the data dimension
		unsigned int get_data_dimension(){
			return dim_cp;
		}

		//! Returning the dimension of the mesh required
		unsigned int get_mesh_dimension(){
			return dim_mesh;
		}

		//! Returning control points for the bspline method
		vect_pts get_control_points(){
			return CP;
		}
		/*!
			@brief Returning data on the edge
			
			We don't have particular properties on the edges in this application,
			so we return a default edge property. Anyway, there will be no need 
			to use this member, since creating a new edge already means to 
			default construct all its properties
		*/
		NetDiff::B_Edge_prop get_edge_data(){
			return NetDiff::B_Edge_prop();	//volendo qua si può restituire vuoto, tanto non c'è nessuna property particolare
		}
		
		/*!
			@brief Returning topological data
			
			No need of topological data, leaving it blank using the empty
			struct provided in the BGLgeom
		*/
		BGLgeom::no_topological_data get_topological_data(){
			return BGLgeom::no_topological_data();
		}
			
	private:
		//! Coordinates of source and target
		points SRC,TGT;
		BGLgeom::boundary_condition BC_src, BC_tgt;
		unsigned int dim_cp;
		unsigned int dim_mesh;
		vect_pts CP;

};	//reader_netdiff

}	//NetDiff


#endif	//HH_NETDIFF_GRAPH_BSPLINE_HH