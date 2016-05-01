/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   assembling1d.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   January 2016.
  @brief  Miscelleanous assembly routines for the 1D network problem.
 */
/** @defgroup asm Assembly routines */

#ifndef M3D1D_ASSEMBLING_1D_HPP_
#define M3D1D_ASSEMBLING_1D_HPP_
#include <defines.hpp>
#include <node.hpp>
#include <utilities.hpp>

namespace getfem {

//! Build the mass and divergence matrices for the 1D Poiseuille's problem,
//! @f$ M = \int_{\Lambda} c~u~v~ds @f$ and
//! @f$ D = \int_{\Lambda} \nabla u \cdot \mathbf{\lambda}\,p~ds @f$
/*!
	@param M         Computed mass matrix
	@param D         Computed divergence matrix
	@param mim       The integration method to be used
	@param mf_u      The finite element method for the velocity @f$ \mathbf{u} @f$
	@param mf_p      The finite element method for the pressure @f$ p @f$
	@param mf_data   The finite element method for the tangent versor on @f$ \Lambda @f$
	@param coef      The coefficient for M
	@param lambdax   First cartesian component of the tangent versor  @f$ \mathbf{\lambda} @f$
	@param lambday   Second cartesian component of the tangent versor @f$ \mathbf{\lambda} @f$
	@param lambdaz   Third cartesian component of the tangent versor @f$ \mathbf{\lambda} @f$
	@param rg        The region where to integrate

	@ingroup asm
 */ 
template<typename MAT, typename VEC>
void 
asm_network_poiseuille
	(MAT & M, MAT & D,
	 const mesh_im & mim,
	 const mesh_fem & mf_u,
	 const mesh_fem & mf_p,
	 const mesh_fem & mf_data,
	 const VEC & coef,
	 const VEC & lambdax, const VEC & lambday, const VEC & lambdaz,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf_p.get_qdim() == 1 && mf_u.get_qdim() == 1,
		"invalid data mesh fem (Qdim=1 required)");
	// Build the local mass matrix Mvvi
	getfem::asm_mass_matrix_param(M, mim, mf_u, mf_data, coef, rg);
	// Build the local divergence matrix Dvvi
	generic_assembly 
	assem("l1=data$1(#3); l2=data$2(#3); l3=data$3(#3);"
		  "t=comp(Base(#2).Grad(#1).Base(#3));"
		  "M$1(#2,#1)+=t(:,:,1,i).l1(i)+t(:,:,2,i).l2(i)+t(:,:,3,i).l3(i);");
	assem.push_mi(mim);
	assem.push_mf(mf_u);
	assem.push_mf(mf_p);
	assem.push_mf(mf_data);
	assem.push_data(lambdax);
	assem.push_data(lambday);
	assem.push_data(lambdaz);
	assem.push_mat(D);
	assem.assembly(rg);
}


/*! Build the mixed boundary conditions for Poiseuille's problem
    @f$ M=\int_{\mathcal{E}_u} \frac{1}{\beta}\,u\,v~d\sigma@f$ and
    @f$ F=-\int_{\mathcal{E}_u} p0\,v~d\sigma-\int_{\mathcal{E}_p} g\,v~d\sigma@f$
 */
/*!
	@param M        BC contribution to Poiseuille's mass matrix
	@param F        BC contribution to Poiseuille's rhs
	@param mim      The integration method to be used
	@param mf_u     The finite element method for the velocity @f$u@f$
	@param mf_data  The finite element method for the coefficients
	@param BC       Array of values of network boundary points
	@param P0       Array of values of the external pressure @f$p_0@f$
	@param R        Network radii
	@param rg       The region where to integrate
	@ingroup asm
 */ 
template<typename MAT, typename VEC>
void
asm_network_bc
	(MAT & M, VEC & F,
	 const mesh_im & mim,
	 const std::vector<mesh_fem> & mf_u,
	 const mesh_fem & mf_data,
	 const std::vector<getfem::node> & BC,
	 const VEC & P0,
	 const VEC & R,
	 const scalar_type beta
	 ) 
{
	// Aux data
	std::vector<scalar_type> ones(mf_data.nb_dof(), 1.0);

	for (size_type bc=0; bc < BC.size(); bc++) {

		size_type i = abs(BC[bc].branches[0]);
		size_type start = i*mf_u[i].nb_dof();
		scalar_type Ri = compute_radius(mim, mf_data, R, i);

		if (BC[bc].label=="DIR") { // Dirichlet BC
			// Add gv contribution to Fv
			scalar_type BCVal = BC[bc].value*pi*Ri*Ri;
			getfem::asm_source_term(gmm::sub_vector(F, gmm::sub_interval(start,mf_u[i].nb_dof())), 
				mim, mf_u[i], mf_data, gmm::scaled(ones, BCVal), BC[bc].rg);
		} 
		else if (BC[bc].label=="MIX") { // Robin BC
			// Add correction to Mvv
			MAT Mi(mf_u[i].nb_dof(), mf_u[i].nb_dof());
			getfem::asm_mass_matrix(Mi,
				mim, mf_u[i], BC[bc].rg);
			gmm::scale(Mi, pi*pi*Ri*Ri*Ri*Ri/beta);				
			gmm::add(gmm::scaled(Mi, -1.0), 
				gmm::sub_matrix(M,
					gmm::sub_interval(start, mf_u[i].nb_dof()),
					gmm::sub_interval(start, mf_u[i].nb_dof())));
			gmm::clear(Mi);	
			// Add p0 contribution to Fv
			getfem::asm_source_term(gmm::sub_vector(F, gmm::sub_interval(start,mf_u[i].nb_dof())), 
				mim, mf_u[i], mf_data, gmm::scaled(P0, pi*Ri*Ri), BC[bc].rg);			
		}
		else if (BC[bc].label=="INT") { // Internal Node
			DAL_WARNING1("internal node passed as boundary.");
		}
		else if (BC[bc].label=="JUN") { // Junction Node
			DAL_WARNING1("junction node passed as boundary.");
		}
		else {
			GMM_ASSERT1(0, "Unknown Boundary Condition"<< BC[bc].label << endl);
		}
	}

}


/*!
	Compute the network junction matrix @f$J=\langle[u],p\rangle_{\Lambda}@f$.
	@ingroup asm
 */
template<typename MAT, typename VEC>
void
asm_network_junctions
	(MAT & J,
	 const mesh_im & mim,
	 const std::vector<mesh_fem> & mf_u,
	 const mesh_fem & mf_p,
	 const mesh_fem & mf_data,
	 const std::vector<getfem::node> & J_data,
	 const VEC & radius
	 ) 
{
	GMM_ASSERT1 (mf_p.get_qdim() == 1, 
		"invalid data mesh fem for pressure (Qdim=1 required)");
	GMM_ASSERT1 (mf_u[0].get_qdim() == 1, 
		"invalid data mesh fem for velocity (Qdim=1 required)");
	GMM_ASSERT1 (getfem::name_of_fem(mf_p.fem_of_element(0)) != "FEM_PK(1,0)" &&
		getfem::name_of_fem(mf_p.fem_of_element(0)) != "FEM_PK_DISCONTINUOUS(1,0)",
		"invalid data mesh fem for pressure (k>0 required)");
	
	for (size_type i=0; i<mf_u.size(); ++i){ /* branch loop */

		scalar_type Ri = compute_radius(mim, mf_data, radius, i);
		//cout << "Region " << i << " : radius=" << Ri << endl;

		for (size_type j=0; j<J_data.size(); ++j){

			// Identify pressure dof corresponding to junction node
			VEC psi(mf_p.nb_dof());
			asm_basis_function(psi, mim, mf_p, J_data[j].rg);
			size_type row = 0;
			bool found = false;
			while (!found && row<mf_p.nb_dof()){
				found = (1.0 - psi[row] < 1.0E-06);
				if (!found) row++;
			}
			GMM_ASSERT1 (row!=0 && found,  // No junction in first point
				"Error in assembling pressure basis function");
			std::vector<long signed int>::const_iterator bb = J_data[j].branches.begin();
			std::vector<long signed int>::const_iterator be = J_data[j].branches.end();
			// Outflow branch contribution
			if (std::find(bb, be, i) != be){
				J(row, (i+1)*mf_u[i].nb_dof()-1) -= pi*Ri*Ri; //col to be generalized!
			}
			// Inflow branch contribution
			if (i!=0 && std::find(bb, be, -i) != be){
				J(row, i*mf_u[i].nb_dof()) += pi*Ri*Ri;	//col to be generalized!
			}
		}
	}

} /* end of asm_junctions */


} /* end of namespace */

#endif
