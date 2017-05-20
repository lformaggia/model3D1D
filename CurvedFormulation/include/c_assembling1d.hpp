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
		scalar_type Ri = simple_compute_radius(mim, mf_data, radius, i);
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

template<typename MAT, typename VEC, typename PARAM>
void
asm_network_junctions
	(MAT & J,
	 const mesh_im & mim,
	 const std::vector<mesh_fem> & mf_u,
	 const mesh_fem & mf_p,
	 const mesh_fem & mf_data,
	 const std::vector<getfem::node> & J_data,
	 const VEC & radius,
	 PARAM & P
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

		scalar_type Ri = simple_compute_radius(mim, mf_data, radius, i);
		
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
			size_type last_, first_;
			vector_type dof_enum;
			size_type fine=0;
			for(getfem::mr_visitor mrv(mf_u[i].linked_mesh().region(i)); !mrv.finished(); ++mrv){
				for(auto b: mf_u[i].ind_basic_dof_of_element(mrv.cv())){
					dof_enum.emplace_back(b);
					fine++;
				}
			}
			first_=dof_enum[0];
			last_ =dof_enum[fine-1];
			dof_enum.clear();
			
			if (std::find(bb, be, i) != be){
				J(row, i*mf_u[i].nb_dof()+last_) -= pi*Ri*Ri; //col to be generalized!
			}
			// Inflow branch contribution
			if (i!=0 && std::find(bb, be, -i) != be){
				J(row, i*mf_u[i].nb_dof()+first_) += pi*Ri*Ri;	//col to be generalized!
			}
		}
	}
} /* end of asm_junctions */

//! Aux function to extract the radius of the ith branch, R[i] 
template
<typename VEC>
scalar_type
simple_compute_radius
	(const mesh_im & mim,
	 const mesh_fem & mf_coef,
	 const VEC & R,
	 const size_type & rg
	 ) 
{
	vector_type dof_enum;
	size_type fine=0;
	for(getfem::mr_visitor mrv(mf_coef.linked_mesh().region(rg)); !mrv.finished(); ++mrv){
				for(auto b: mf_coef.ind_basic_dof_of_element(mrv.cv())){
					dof_enum.emplace_back(b);
					fine++;
				}
			}
	size_type first_=dof_enum[0];
	return R[first_];
}


} /* end of namespace */

#endif
