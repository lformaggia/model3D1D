
#ifndef C_M3D1D_ASSEMBLING_1D_HPP_
#define C_M3D1D_ASSEMBLING_1D_HPP_
#include <defines.hpp>
#include <node.hpp>
#include <utilities.hpp>

namespace getfem {

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
		scalar_type Gi = G[0];
		scalar_type ki_s= k[i][0];
		scalar_type ki_e= k[i].back();
		scalar_type ui_s= u_old[shift];
		scalar_type ui_e= u_old[shift+dofi-1];

		if(Is_Newton)
			Gi=Gi*2;

		NLJ(shift, shift)-=pi*Ri*Ri*Gi*ui_s*
						( 	alfa  + 
						  	beta  * ki_s * ki_s * Ri * Ri +   
					  		gamma * ki_s * ki_s * ki_s *ki_s * Ri * Ri * Ri * Ri
						);				
		NLJ(shift+dofi-1, shift+dofi-1)+= pi*Ri*Ri*Gi*ui_e*
									 ( 	alfa  +
								  	 	beta  * ki_e * ki_e * Ri * Ri +   
								  	 	gamma * ki_e * ki_e * ki_e * ki_e * Ri * Ri * Ri * Ri
									  );	

	}
}



template<typename VEC_getfem, typename VEC, typename MESH_FEM>
scalar_type
euclidean_norm
(
	VEC_getfem & Unew,
	VEC & Uold,
	MESH_FEM & mf
) 
{
	scalar_type dist=0;
	size_type dofi=0;
	size_type shift;
	for(size_type branch=0; branch<mf.size();branch++){
		dofi=mf[branch].nb_dof();
		if(branch>0) shift=shift+dofi;
		for(size_type el=0; el<dofi; el++){
			dist=dist+abs(Unew[shift+el]-Uold[shift+el]);
		}
	}
	return dist;
}


} /* end of namespace */

#endif
