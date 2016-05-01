/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   assembling3d1d.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   January 2016.
  @brief  Miscelleanous assembly routines for the 3D/1D coupling.
 */
#ifndef M3D1D_ASSEMBLING_3D1D_HPP_
#define M3D1D_ASSEMBLING_3D1D_HPP_
#include <defines.hpp>
#include <utilities.hpp>

namespace getfem {

/*!
	Build the averaging matrix @f$\bar{\Pi}_{tv}@f$ and the interpolation
	matrix @f${\Pi}_{tv}@f$. 
	Used to build the exchange matrices @f$B_{ii}@f$.
	@ingroup asm
 */
template<typename MAT, typename VEC>
void 
asm_exchange_aux_mat
	(MAT & Mbar, MAT & Mlin,
	 const getfem::mesh_im & mim,
	 const getfem::mesh_fem & mf_t,
	 const getfem::mesh_fem & mf_v,
	 const VEC & RADIUS,
	 const size_type NInt
	 ) 
{
	gmm::clear(Mbar); gmm::clear(Mlin);
	// Aux params
	const scalar_type Pi = 2*acos(0.0);
	size_type nb_dof_t     = mf_t.nb_dof();
	size_type nb_dof_v     = mf_v.nb_dof();

	// Linear interpolation map mf_t --> mf_v
	getfem::interpolation(mf_t, mf_v, Mlin);  

	// Aux vectors for local interpolation
	std::vector<scalar_type> Pbari(NInt); 
	std::vector<scalar_type> Pt(nb_dof_t); 
	size_type counter = 0;
	for (size_type i = 0; i < nb_dof_v; i++){
		counter++;
		if (counter*100 >= nb_dof_v) {
			counter = 0; 
			#ifdef M3D1D_VERBOSE_
			cout << "*"; cout.flush();
			#endif
		}
		// We need the following interpolation tool <getfem_interpolation.h>      
		getfem::mesh_trans_inv mti(mf_t.linked_mesh());
		// Build the list of point on the i-th circle:
		// ... first consider an orthonormal system v0, v1, v2:
		base_node v0;
		if (i==0) {
			v0 = mf_v.point_of_basic_dof(i+1) - mf_v.point_of_basic_dof(i);			
		} else {
			v0 = mf_v.point_of_basic_dof(i) - mf_v.point_of_basic_dof(i-1);
		}
		base_node v1(0.0, -v0[2], v0[1]);
		base_node v2(v0[1]*v0[1] +v0[2]*v0[2], -v0[0]*v0[1], -v0[0]*v0[2]);
		if (gmm::vect_norm2(v2) < 1.0e-8 * gmm::vect_norm2(v0)) {
			v1[0] = -v0[1]; v1[1] = v0[0]; v1[2] = 0.0;
			v2[0] = -v0[0]*v0[2]; v2[1] = -v0[1]*v0[2]; v2[2] = v0[0]*v0[0] +v0[1]*v0[1];
		}
		v1 = v1 / gmm::vect_norm2(v1);
		v2 = v2 / gmm::vect_norm2(v2);
		// ... then parametrize the circle:
		for (size_type j = 0; j < NInt; j++){ 	
			mti.add_point( mf_v.point_of_basic_dof(i) + 
						   RADIUS[i]*(cos(2*Pi*j/NInt)*v1 + sin(2*Pi*j/NInt)*v2) );
		}
		// Get the local interpolation matrix Mbari
		MAT Mbari(NInt, nb_dof_t); gmm::clear(Mbari);
		interpolation(mf_t, mti, Pt, Pbari, Mbari, 1);
		scalar_type sum_row = 0.0;
		for (size_type j=0; j < NInt; ++j) {
			typename gmm::linalg_traits<MAT>::const_sub_row_type 
				row = mat_const_row(Mbari,j);
			gmm::linalg_traits< gmm::rsvector<scalar_type> >::const_iterator 
				it_nz = vect_const_begin(row);
			gmm::linalg_traits< gmm::rsvector<scalar_type> >::const_iterator 
				ite_nz = vect_const_end(row);
			for (; it_nz != ite_nz ; ++it_nz) {
				Mbar(i, it_nz.index()) += (*it_nz);
				sum_row += (*it_nz);
			}
		}
		typename gmm::linalg_traits<MAT>::sub_row_type 
			row = mat_row(Mbar,i);
		gmm::linalg_traits< gmm::rsvector<scalar_type> >::iterator
			it_nz = vect_begin(row);
		gmm::linalg_traits< gmm::rsvector<scalar_type> >::iterator
			ite_nz = vect_end(row);
		for (; it_nz != ite_nz ; ++it_nz) {
			(*it_nz)/=sum_row;
		}

	} /* end of outer for loop */
	cout << endl;

} /* end of build_aux_matrices */


/*!
	Build the exchange matrices
	@f$B_{tt} = \Pi^T_{tv} M_{vv} \bar{\Pi}_{tv}@f$,
	@f$B_{tv} = \Pi^T_{tv} M_{vv}@f$,
	@f$B_{vt} = M_{vv} \bar{\Pi}_{tv}@f$,
	@f$B_{vv} = M_{vv}@f$.
	If ALT_FORM==true we substitute @f${\Pi}_{tv}@f$ with
	@f$\bar{\Pi}_{tv}@f$.
	@ingroup asm
 */
template<typename MAT, typename VEC>
void 
asm_exchange_mat
	(MAT & Btt, MAT & Btv, MAT & Bvt, MAT & Bvv, 
	 const getfem::mesh_im & mim,
	 const getfem::mesh_fem & mf_v, 
	 const getfem::mesh_fem & mf_coefv,
	 const MAT & Mbar, const MAT & Mlin,
	 const VEC & Q,
	 const bool ALT_FORM
	 ) 
{
	#ifdef M3D1D_VERBOSE_
	cout << "    Assembling Bvv ..." << endl;  
	#endif
	getfem::asm_mass_matrix_param(Bvv, mim, mf_v, mf_coefv, Q); 
	#ifdef M3D1D_VERBOSE_
	cout << "    Assembling Bvt ..." << endl;
	#endif
	gmm::mult(Bvv, Mbar, Bvt);
	if (ALT_FORM){
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btv (alternative form) ..." << endl;
		#endif
		gmm::copy(gmm::transposed(Bvt), Btv); 
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btt (alternative form) ..." << endl;
		#endif
		gmm::mult(gmm::transposed(Mbar), Bvt, Btt); 
	}
	else{
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btv ..." << endl;
		#endif
		gmm::mult(gmm::transposed(Mlin), Bvv, Btv); 
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btt ..." << endl;
		#endif
		gmm::mult(gmm::transposed(Mlin), Bvt, Btt); 
	}
	
} /* end of build_exchange_matrices */


} /* end of namespace */

#endif
