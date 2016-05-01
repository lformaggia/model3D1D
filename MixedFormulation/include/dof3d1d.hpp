/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   dof3d1d.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   January 2016.
  @brief  Definition of the aux class for the number of degrees of freedom.
 */
#ifndef M3D1D_DOF3D1D_HPP_
#define M3D1D_DOF3D1D_HPP_

namespace getfem {

//! Class to store the number of degrees of freedom of used FEMs
struct dof3d1d {

	//! Number of dof of the interstitial velocity FEM mf_Ut
	size_type Ut_;
	//! Number of dof of the interstitial pressure FEM mf_Pt
	size_type Pt_;
	//! Number of dof of the interstitial coefficients FEM mf_coeft
	size_type ct_;
	//! Number of dof of the vessel network velocity
	//! It is the sum of local vessel branch dof
	size_type Uv_;
	//! Number of dof of the vessel pressure FEM mf_Pv
	size_type Pv_;
	//! Number of dof of the vessel network coefficients FEM mf_coefv
	//! It is NOT the sum of local vessel branch dof
	size_type cv_;
	//! Total number of dof of the tissue problem
	size_type tissue_;
	//! Total number of dof of the vessel problem
	size_type vessel_;
	//! Total number of dof of the coupled problem
	size_type tot_;
	
	//! Compute the number of dof of given FEM
	void set (
			const getfem::mesh_fem & mf_Ut, const getfem::mesh_fem & mf_Pt,
			const std::vector<getfem::mesh_fem> & mf_Uv, const getfem::mesh_fem & mf_Pv,
			const getfem::mesh_fem & mf_coeft, const getfem::mesh_fem & mf_coefv
			)
	{
		Ut_ = mf_Ut.nb_dof(); 
		Pt_ = mf_Pt.nb_dof();
		ct_ = mf_coeft.nb_dof();
		Uv_ = 0;
		for (size_type i = 0; i<mf_Uv.size(); ++i) Uv_ += mf_Uv[i].nb_dof();
		Pv_ = mf_Pv.nb_dof();
		cv_ = mf_coefv.nb_dof();
		
		tissue_ = Pt_ + Ut_; 
		vessel_ = Pv_ + Uv_;
		tot_    = tissue_ + vessel_;
	}
	//! Accessor to the number of dof of mf_Ut
	inline size_type Ut (void) { return Ut_; } const
	//! Accessor to the number of dof of mf_Pt
	inline size_type Pt (void) { return Pt_; } const
	//! Accessor to the number of dof of mf_Uv
	inline size_type Uv (void) { return Uv_; } const
	//! Accessor to the number of dof of mf_Pv
	inline size_type Pv (void) { return Pv_; } const
	//! Accessor to the number of dof of mf_coeft
	inline size_type coeft (void) { return ct_; } const
	//! Accessor to the number of dof of mf_coefv
	inline size_type coefv (void) { return cv_; } const
	//! Accessor to the number of dof of tissue problem
	inline size_type tissue (void) { return tissue_; } const
	//! Accessor to the number of dof of vessel problem
	inline size_type vessel (void) { return vessel_; } const
	//! Accessor to the number of dof of coupled problem
	inline size_type tot (void) { return tot_; } const

	//! Overloading of the output operator
	friend std::ostream & operator << (
			std::ostream & out, const dof3d1d & dof
			)
	{ 
		out << "--- DEGREES OF FREEDOM --- " << endl;
		out << "  nb_dof_Pt     : " 		 << dof.Pt_  << endl;
		out << "  nb_dof_Ut     : " 		 << dof.Ut_  << endl;
		out << "  nb_dof_coeft  : " 		 << dof.ct_  << endl;
		out << "  nb_dof_Pv     : " 		 << dof.Pv_  << endl;
		out << "  nb_dof_Uv     : " 		 << dof.Uv_  << endl;
		out << "  nb_dof_coefv  : " 		 << dof.cv_  << endl;
		out << "  nb_dof_tot    : " 		 << dof.tot_ << endl;
		out << "-------------------------- " << endl;
		return out;            
	}

};

} /* end of namespace */

#endif
