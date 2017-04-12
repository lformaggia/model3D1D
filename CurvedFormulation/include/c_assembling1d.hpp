/* -*- c++ -*- (enableMbars emacs c++ mode) */ 
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Curved Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2017 Giorgio Raimondi
======================================================================*/
/*! 
  @file   c_assembling.hpp
  @author Giorgio Raimondi <giorgio3.raimondi@mail.polimi.it>
  @date   May 2017.
  @brief Miscelleaneous assembly routines for the 1D network preblem for non Linear terms.
 */    
 /*! @defgroup asm Assembly routines */
#ifndef C_M3D1D_ASSEMBLING_1D_HPP_
#define C_M3D1D_ASSEMBLING_1D_HPP_
#include <defines.hpp>
#include <node.hpp>
#include <utilities.hpp>

namespace getfem {


//! Build the advective non linear matrix for the 1D Navier-Stokes curved preblem,
//! @f$ NLM = \int_{\Lambda} c~u~uold~\nabla v \cdot \mathbf{\lambda}\,~ds @f$
/*!
	@param NLM       Computed non linear matrix term
	@param mim       The integration method to be used
	@param mf_u      The finite element method for the velocity @f$ \mathbf{u} @f$
	@param mf_data   The finite element method for the tangent versor on @f$ \Lambda @f$
	@param coef      The coefficient for NLM
	@param lambdax   First cartesian component of the tangent versor  @f$ \mathbf{\lambda} @f$
	@param lambday   Second cartesian component of the tangent versor @f$ \mathbf{\lambda} @f$
	@param lambdaz   Third cartesian component of the tangent versor @f$ \mathbf{\lambda} @f$
	@param u_old     The velocity at the previous iteration
	@param rg        The region where to integrate

	@ingroup asm
 */
template<typename MAT, typename VEC>
void 
asm_network_nonlinear
	(MAT & NLM,
	 const mesh_im & mim,
	 const mesh_fem & mf_u,
	 const mesh_fem & mf_data,
	 const VEC & coef,
	 const VEC & lambdax, const VEC & lambday, const VEC & lambdaz,
	 const VEC & u_old,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf_u.get_qdim() == 1,
		"invalid data mesh fem (Qdim=1 required)");
	// Build the local matrix for the non linear term NLMvvi
	generic_assembly 
	assem("l1=data$1(#2); l2=data$2(#2); l3=data$3(#2); u=data$4(#1); c=data$5(#2);"
		  "t=comp(Base(#1).Base(#1).Grad(#1).Base(#2).Base(#2));"
		  "M$1(#1,#1)+=	 t(j,:,:,1,m,k).l1(m).u(j).c(k)+t(j,:,:,2,m,k).l2(m).u(j).c(k)+t(j,:,:,3,m,k).l3(m).u(j).c(k);");
	assem.push_mi(mim);
	assem.push_mf(mf_u);
	assem.push_mf(mf_data);
	assem.push_data(lambdax);
	assem.push_data(lambday);
	assem.push_data(lambdaz);
	assem.push_data(u_old);
	assem.push_data(coef);
	assem.push_mat(NLM);
	assem.assembly(rg);
}

/*!
	Compute the network junction matrix non linear term @f$J=\langle[u],[v]\rangle_{\Lambda}@f$.
	@ingroup asm
 */
template<typename MAT, typename VEC, typename VEC_I>
void
asm_network_nonlinear_junctions
	(MAT & NLJ,
	 const mesh_im & mim,
	 const std::vector<mesh_fem> & mf_ui,
	 const mesh_fem & mf_data,
	 const scalar_type Gam,
	 const VEC & radius,
	 const VEC & G,
	 const VEC_I & k,
	 const VEC & u_old,
	 const bool Is_Newton
	) 
{

	GMM_ASSERT1 (mf_ui[0].get_qdim() == 1, 
		"invalid data mesh fem for velocity (Qdim=1 required)");

	gmm::clear(NLJ);
	scalar_type alfa  = (Gam+2.0)/(Gam+1.0);
	scalar_type beta  = (Gam+2.0)/(Gam+4.0);
	scalar_type gamma = (Gam+2.0)*(Gam+2.0)/3/(Gam+3.0)/(Gam+6.0);
	size_type shift=0;
	size_type dofi=mf_ui[0].nb_dof();

	for (size_type i=0; i<mf_ui.size(); ++i){ /* branch loop */
		if(i>0) shift+=dofi;

		dofi=mf_ui[i].nb_dof();
		scalar_type Ri = compute_radius(mim, mf_data, radius, i);
		scalar_type Gi = G[0]; // Dynamic pressure gain term
		scalar_type ki_s= k[i][0]; //Curvature at the beginning of the branch
		scalar_type ki_e= k[i].back(); //Curvature at the end of the branch
		scalar_type ui_s= u_old[shift]; //Velocity at the previous Iteration at the beginning of the branch
		scalar_type ui_e= u_old[shift+dofi-1]; //Velocity at the previous Iteration at the end of the branch
		scalar_type Method=1.0; // Kind of method used

		if(Is_Newton) //Using newton method
			Method*=2.0; 

		//Element at the beginning of the branch
		NLJ(shift, shift)-=Method*pi*Ri*Ri*Gi*ui_s*
						( 	alfa  + 
						  	beta  * ki_s * ki_s * Ri * Ri +   
					  		gamma * ki_s * ki_s * ki_s *ki_s * Ri * Ri * Ri * Ri
						);				

		//Element at the end of the branch
		NLJ(shift+dofi-1, shift+dofi-1)+= Method*pi*Ri*Ri*Gi*ui_e*
									 ( 	alfa  +
								  	 	beta  * ki_e * ki_e * Ri * Ri +   
								  	 	gamma * ki_e * ki_e * ki_e * ki_e * Ri * Ri * Ri * Ri
									  );	

	}
}


} /* end of namespace */

#endif
