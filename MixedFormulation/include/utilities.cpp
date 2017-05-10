/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   utilities.cpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   January 2016.
  @brief  Definition of some miscelleanous auxiliary routines.
 */
 
#include <vector>
#include <string>
#include <sstream>
#include <utilities.hpp>

namespace getfem {

// Aux function to compute the diameter of an element
scalar_type 
estimate_h(const mesh & mesh, const size_type i) 
{
	std::vector<size_type> cpt = mesh.ind_points_of_convex(i);
	std::vector<size_type>::const_iterator icpt, jcpt;
	scalar_type d=0.0;
	for (icpt = cpt.begin(); icpt != cpt.end(); icpt++ ) {
		for (jcpt = cpt.begin(); jcpt != cpt.end(); jcpt++ ) {
			d = std::max(d, 
				gmm::vect_norm2(mesh.points()[*icpt]-mesh.points()[*jcpt]));
		}
	}
	return d;
}

// Read an array of string split by a delim. 
// Store the results in a pre-constructed vector
std::vector<std::string> &
split(const std::string & s, 
	  char delim, 
	  std::vector<std::string> & elems
	  ) 
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) 
			elems.push_back(item);
    }
    return elems;
}

// Read an array of string split by a delim. 
// Return a new vector
std::vector<std::string>
split(const std::string & s, 
	  char delim
	  ) 
{
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}



//! Exact pressure 
double sol_pt(const bgeot::base_node & x){ 
	return sin(2.0*pi/Lx*x[0])*sin(2.0*pi/Ly*x[1])*sin(2*pi/Lz*x[2]);
}
//! Exact x velocity
double sol_utx(const bgeot::base_node & x){
	return -2.0*pi*kappat/Lx*cos(2.0*pi/Lx*x[0])*sin(2.0*pi/Ly*x[1])*sin(2.0*pi/Lz*x[2]);
}
//! Exact y velocity
double sol_uty(const bgeot::base_node & x){
	return -2.0*pi*kappat/Ly*sin(2.0*pi/Lx*x[0])*cos(2.0*pi/Ly*x[1])*sin(2.0*pi/Lz*x[2]);
} 
//! Exact z velocity
double sol_utz(const bgeot::base_node & x){
	return -2.0*pi*kappat/Lz*sin(2.0*pi/Lx*x[0])*sin(2.0*pi/Ly*x[1])*cos(2.0*pi/Lz*x[2]);
}
//! Exact velocity magnitude
double sol_utm(const bgeot::base_node & x){
	return sqrt(sol_utx(x)*sol_utx(x)+sol_uty(x)*sol_uty(x)+sol_utz(x)*sol_utz(x));
}
//! Exact vectorial velocity
std::vector<double> sol_ut(const bgeot::base_node & x){
	std::vector<double> out(3);
	out[0] = sol_utx(x); out[1] = sol_uty(x); out[2] = sol_utz(x);
	return out;
}
//! Exact rhs
double sol_gt(const bgeot::base_node & x){
	return 4.0*pi*pi*kappat*(1.0/(Lx*Lx)+1.0/(Ly*Ly)+1.0/(Lz*Lz))*sin(2.0*pi/Lx*x[0])*sin(2.0*pi/Ly*x[1])*sin(2.0*pi/Lz*x[2]);
}


}
