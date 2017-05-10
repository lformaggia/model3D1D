/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   utilities.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   January 2016.
  @brief  Declaration of some miscelleanous auxiliary routines.
 */
#ifndef M3D1D_UTILITIES_HPP_
#define M3D1D_UTILITIES_HPP_

#include <getfem/getfem_assembling.h> 
#include <getfem/getfem_mesh.h>
#include <gmm/gmm.h>
#include <defines.hpp>

namespace getfem {

//! Aux function to compute the diameter of an element
scalar_type 
estimate_h(const mesh & mesh, const size_type i); 

//! Aux function to read an array of string split by a delim. 
//! Store the results in a pre-constructed vector
std::vector<std::string> &
split(const std::string & s, 
	  char delim, 
	  std::vector<std::string> & elems
	  ) ;

//! Aux function to read an array of string split by a delim. 
//! Return a new vector
std::vector<std::string>
split(const std::string & s, 
	  char delim
	  ) ;

//! Build the integral of FE base functions
//! @f$ \Phi = \int_{\Omega} \phi(x)~dx @f$
/*!
 @param V    The array of computed values
 @param mim  The integration method to be used
 @param mf   The finite element method
 @param rg   The region where to integrate
 */ 
template
<typename VEC>
void 
asm_basis_function
	(VEC & V,
	 const mesh_im & mim,
	 const mesh_fem & mf,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf.get_qdim() == 1, 
		"invalid data mesh fem for pressure (Qdim=1 required)");
	generic_assembly 
	assem("V$1(#1)+=comp(Base(#1));");
	assem.push_mi(mim);
	assem.push_mf(mf);
	assem.push_vec(V);
	assem.assembly(rg);
}

//! Aux function to extract the radius of the ith branch, R[i] 
template
<typename VEC>
scalar_type
compute_radius
	(const mesh_im & mim,
	 const mesh_fem & mf_coef,
	 const VEC & R,
	 const size_type & rg
	 ) 
{
	VEC rbasis(mf_coef.nb_dof());
	asm_basis_function(rbasis, mim, mf_coef, rg);
	size_type firstcv = mf_coef.linked_mesh().region(rg).index().first_true();
	size_type rg_size = mf_coef.linked_mesh().region(rg).size();
	gmm::scale(rbasis, 1.0/estimate_h(mf_coef.linked_mesh(), firstcv));
	return gmm::vect_sp(rbasis, R)/rg_size;
}

//! Compute the integral of the solution
template
<typename VEC>
scalar_type 
asm_mean(const mesh_fem &mf, const mesh_im &mim, const VEC &U) 
{
	getfem::generic_assembly assem;
	assem.set("u=data(#1); V()+=u(i).comp(Base(#1))(i)");
	assem.push_mi(mim);
	assem.push_mf(mf);
	assem.push_data(U);
	std::vector<scalar_type> v(1);
	assem.push_vec(v);
	assem.assembly();
	
	return v[0];
}

//! Aux function to save a gmm::col_vector as a plain list of values
template 
<typename L> 
void 
write_col_vector(std::ostream &o, const L &l) 
{
    for (size_type i = 0; i < gmm::mat_nrows(l); ++i) {

	if (gmm::mat_ncols(l) != 0) o << ' ' << l(i, 0);
	for (size_type j = 1; j < gmm::mat_ncols(l); ++j) o << ", " << l(i, j); 
      o << " \n";
    }
}


//! Parameteres for exact solution (1_uncoupled)
/*! 
	\todo Read from param Lx, Ly, Lz and Kt
 */
const double Lx = 1.0, Ly = 1.0, Lz = 1.0, kappat = 1.0; 
	//! Exact pressure 
double sol_pt(const bgeot::base_node & x);
	//! Exact x velocity
double sol_utx(const bgeot::base_node & x);
	//! Exact y velocity
double sol_uty(const bgeot::base_node & x);
	//! Exact z velocity
double sol_utz(const bgeot::base_node & x);
	//! Exact velocity magnitude
double sol_utm(const bgeot::base_node & x);
	//! Exact vectorial velocity
std::vector<double> sol_ut(const bgeot::base_node & x);
	//! Exact rhs
double sol_gt(const bgeot::base_node & x);




} /* end of namespace */

#endif
