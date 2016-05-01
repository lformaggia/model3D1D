/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   node.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   November 2015.
  @brief  Definition of the class node.
 */
 
#ifndef M3D1D_NODE_HPP_
#define M3D1D_NODE_HPP_

namespace getfem{

//! Class to handle the boundary and junction nodes
struct node {

	//! Label ('INT','DIR','MIX','JUN')
	std::string label; 
	//! Boundary value (useless for 'JUN')
	scalar_type value;
	//! Global index
	size_type   idx;
	//! Associated mesh region
	size_type   rg;
	//! Possible list of intersecting vessel branches
	std::vector<long signed int> branches;
	
	//! Constructor
	node(const std::string & label_="", 
		 const scalar_type & value_=0, 
		 const size_type & idx_=0, 
		 const size_type & rg_=0) 
		 : label(label_), value(value_), idx(idx_), rg(rg_)
	{}
	//! Overloading of the output operator
	friend std::ostream & operator << (
		std::ostream & out, const node & N
		)
	{ 
		out << "('" << N.label << "'," 
			<< N.value << "," 
			<< N.idx   << ","
			<< N.rg    << ")";

		return out;            
	}

};

}

#endif
