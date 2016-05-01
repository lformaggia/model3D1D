/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   assembling3d.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   January 2016.
  @brief  Miscelleanous assembly routines for the 3D tissue problem.
 */
#ifndef M3D1D_ASSEMBLING_3D_HPP_
#define M3D1D_ASSEMBLING_3D_HPP_
#include <defines.hpp>
#include <utilities.hpp>

namespace getfem {

//! Build the mass and divergence matrices for the 3D Darcy's problem,
//! @f$ M = \int_{\Omega} \frac{1}{\kappa}\,\mathbf{u}\cdot\mathbf{v}~dx @f$ and
//! @f$ D = \int_{\Omega} div(\mathbf{u})\,p~dx @f$
/*!
	@param M     Darcy's mass matrix
	@param D     Darcy's divergence matrix
	@param mim   The integration method to be used
	@param mf_u  The finite element method for the velocity @f$ \mathbf{u} @f$
	@param mf_p  The finite element method for the pressure @f$ p @f$
	@param rg    The region where to integrate

	@ingroup asm
 */ 
template<typename MAT>
void 
asm_tissue_darcy
	(MAT & M, MAT & D,
	 const mesh_im & mim,
	 const mesh_fem & mf_u,
	 const mesh_fem & mf_p,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf_p.get_qdim() == 1, 
		"invalid data mesh fem for pressure (Qdim=1 required)");
	GMM_ASSERT1(mf_u.get_qdim() > 1, 
		"invalid data mesh fem for velocity (Qdim>1 required)");
	// Build the mass matrix Mtt
	getfem::asm_mass_matrix(M, mim, mf_u, rg);
	// Build the divergence matrix Dtt
	generic_assembly 
	assem("M$1(#2,#1)+=comp(Base(#2).vGrad(#1))(:,:,i,i);");
	assem.push_mi(mim);
	assem.push_mf(mf_u);
	assem.push_mf(mf_p);
	assem.push_mat(D);
	assem.assembly(rg);
}

/*! Build the mixed boundary conditions for Darcy's problem
    @f$ M = \int_{\Gamma_u} \frac{1}{\beta}\,(u.n)(v.n)~d\sigma
    @f$ F = - \int_{\Gamma_u} p0\,(v.n)~d\sigma - \int_{\Gamma_p} g\,(v.n)~d\sigma
 */
/*!
	@param M        BC contribution to Darcy's mass matrix
	@param F        BC contribution to Darcy's rhs
	@param mim      The integration method to be used
	@param mf_u     The finite element method for the velocity @f$ \mathbf{u} @f$
	@param mf_data  The finite element method for the coefficients
	@param BC       Array of values of tissue boundary conditions
	@param G        Array of values of the boundary datum @f$g@f$
	@param P0       Array of values of the external pressure @f$p_0@f$
	@param rg       The region where to integrate

	@ingroup asm
 */ 
template<typename MAT, typename VEC>
void
asm_tissue_bc
	(MAT & M, VEC & F,
	 const mesh_im & mim,
	 const mesh_fem & mf_u,
	 const mesh_fem & mf_data,
	 const std::vector<getfem::node> & BC,
	 const VEC & P0,
	 const VEC & coef
	 )
{
	GMM_ASSERT1(mf_u.get_qdim()>1,  "invalid data mesh fem (Qdim>1 required)");
	GMM_ASSERT1(mf_data.get_qdim()==1, "invalid data mesh fem (Qdim=1 required)");

	std::vector<scalar_type> G(mf_data.nb_dof());
	std::vector<scalar_type> ones(mf_data.nb_dof(), 1.0);

	// Define assembly for velocity bc (\Gamma_u)
	generic_assembly 
	assemU("g=data$1(#2);" "c=data$2(#2);" 
		   "V$1(#1)+=-g(i).comp(Base(#2).vBase(#1).Normal())(i,:,k,k);"
		   "M$1(#1,#1)+=c(i).comp(Base(#2).vBase(#1).Normal().vBase(#1).Normal())(i,:,j,j,:,k,k);");
	assemU.push_mi(mim);
	assemU.push_mf(mf_u);
	assemU.push_mf(mf_data);
	assemU.push_data(P0);
	assemU.push_data(coef);
	assemU.push_vec(F);
	assemU.push_mat(M);
	// Define assembly for pressure bc (\Gamma_p)
	generic_assembly 
	assemP("p=data$1(#2);" 
		   "V$1(#1)+=-p(i).comp(Base(#2).vBase(#1).Normal())(i,:,k,k);");
	assemP.push_mi(mim);
	assemP.push_mf(mf_u);
	assemP.push_mf(mf_data);
	assemP.push_data(G);
	assemP.push_vec(F);

	for (size_type f=0; f < BC.size(); ++f) {

		GMM_ASSERT1(mf_u.linked_mesh().has_region(f), 
				"missed mesh region" << f);
		if (BC[f].label=="DIR") { // Dirichlet BC
			gmm::copy(gmm::scaled(ones, BC[f].value), G);	
			assemP.assembly(mf_u.linked_mesh().region(BC[f].rg));
		} 
		else if (BC[f].label=="MIX") { // Robin BC
			assemU.assembly(mf_u.linked_mesh().region(BC[f].rg));
		}
		else if (BC[f].label=="INT") { // Internal Node
			DAL_WARNING1("internal node passed as boundary.");
		}
		else if (BC[f].label=="JUN") { // Junction Node
			DAL_WARNING1("junction node passed as boundary.");
		}
		else {
			GMM_ASSERT1(0, "Unknown Boundary Condition " << BC[f].label << endl);
		}
	}

} /* end of asm_tissue_bc */

} /* end of namespace */

#endif
