/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   problem3d1d.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   January 2016.
  @brief  Declaration of the main class for the 3D/1D coupled problem.
 */

#ifndef M3D1D_PROBLEM3D1D_HPP_
#define M3D1D_PROBLEM3D1D_HPP_

// GetFem++ libraries
#include <getfem/getfem_assembling.h> 
#include <getfem/getfem_import.h>
#include <getfem/getfem_export.h>   
#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_mesh.h>
#include <getfem/getfem_derivatives.h>
#include <getfem/getfem_superlu.h>
#include <getfem/getfem_partial_mesh_fem.h>
#include <getfem/getfem_interpolated_fem.h>
#include <gmm/gmm.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm_iter_solvers.h>
// Standard libraries
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string>
// Project headers
#include <defines.hpp>
#include <mesh3d.hpp>       
#include <mesh1d.hpp>
#include <utilities.hpp>
#include <assembling1d.hpp>          
#include <assembling3d.hpp>        
#include <assembling3d1d.hpp>
#include <node.hpp>
#include <dof3d1d.hpp>
#include <descr3d1d.hpp>
#include <param3d1d.hpp>
//#include <defines.hpp>

namespace getfem {

//!	Main class defining the coupled 3D/1D fluid problem.
class problem3d1d {

public:

	//! Constructor
	/*! 
		It links integration methods and finite element methods to the meshes
	 */
	problem3d1d(void) : 
		mimt(mesht),  mimv(meshv),
		mf_Pt(mesht), mf_coeft(mesht), mf_Ut(mesht),
		mf_Pv(meshv), mf_coefv(meshv)
	{} 
	//! Initialize the problem
	/*!
		1. Read the .param filename from standard input
		2. Import problem descriptors (file paths, GetFEM types, ...)
		3. Import mesh for tissue (3D) and vessel network (1D)
		4. Set finite elements and integration methods
		5. Build problem parameters
		6. Build the list of tissue boundary data
		7. Build the list of vessel boundary (and junction) data
	 */
	void init (int argc, char *argv[]);
	//! Assemble the problem
	/*!
		1. Initialize problem matrices and vectors
		2. Build the monolithic matrix AM
		3. Build the monolithic rhs FM
	 */
	void assembly (void);
	//! Solve the problem
	/*!
		Solve the monolithic system AM*UM=FM (direct or iterative)
	 */
	bool solve (void);
	//! Solve the problem with arterial-venous network
	/*!
		Merge arterial and venous networks
		Solve the monolithic system AM*UM=FM (direct or iterative)
	 */
	friend bool merge_and_solve (problem3d1d & Pa, problem3d1d & Pv);
	//! Export results into vtk files
	/*!
		Export solutions Ut, Pt, Uv, Pv from the monolithic array UM
	 */
	void export_vtk (const string & suff = "");
	//! Compute mean tissue pressure
	inline scalar_type mean_pt (void){ 
		return asm_mean(mf_Pt, mimt, 
			gmm::sub_vector(UM, gmm::sub_interval(dof.Ut(), dof.Pt()))); 
	}
	//! Compute mean vessel pressure
	inline scalar_type mean_pv (void){ 
		return asm_mean(mf_Pv, mimv, 
			gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()))); 
	}
	//! Compute total flow rate (network to tissue) - pressure
	inline scalar_type flow_rate (void) { return TFR; };

protected:

	//! Mesh for the interstitial tissue @f$\Omega@f$ (3D)
	mesh mesht;
	//! Mesh for the vessel network @f$\Lambda@f$ (1D)
	mesh meshv;
	//! Intergration Method for the interstitial tissue @f$\Omega@f$
	mesh_im mimt;
	//! Intergration Method for the vessel network @f$\Lambda@f$
	mesh_im mimv;
	//! Finite Element Method for the interstitial velocity @f$\mathbf{u}_t@f$
	mesh_fem mf_Ut;
	//! Finite Element Method for the interstitial pressure @f$p_t@f$
	mesh_fem mf_Pt;
	//! Finite Element Method for PDE coefficients defined on the interstitial volume
	mesh_fem mf_coeft;  
	//! Finite Element Method for the vessel velocity @f$u_v@f$
	//! \note Array of local FEMs on vessel branches  @f$\Lambda_i@f$ (@f$i=1,\dots,N@f$)
	vector<mesh_fem> mf_Uvi;
	//! Finite Element Method for the vessel pressure @f$p_v@f$
	mesh_fem mf_Pv; 
	//! Finite Element Method for PDE coefficients defined on the vessel branches
	//! \note Array of local FEMs on vessel branches  @f$\Lambda_i@f$ (@f$i=1,\dots,N@f$)
	vector<mesh_fem> mf_coefvi;
	//! Finite Element Method for PDE coefficients defined on the network
	mesh_fem mf_coefv;

	////////////////////////////////////////////////////////////////////
	
	//! Input file
	ftool::md_param PARAM;
	//! Algorithm description strings (mesh files, FEM types, solver info, ...) 
	descr3d1d descr;
	//! Physical parameters (dimensionless)
	param3d1d param;
	//! Dimension of the tissue domain (3)
	size_type DIMT;
	//! Number of vertices per branch in the vessel network
	vector_size_type nb_vertices;
	//! Number of branches in the vessel network
	size_type nb_branches;
	//! Number of extrema of the vessel network
	size_type nb_extrema;
	//! Number of junctions of the vessel network
	size_type nb_junctions;
	//! Number of degrees of freedom
	dof3d1d dof;
	//! Total flow rate from network to tissue
	scalar_type TFR;

	////////////////////////////////////////////////////////////////////

	//! List of BC nodes on the tissue
	vector< node > BCt;
	//! List of BC nodes of the network
	vector< node > BCv;
	//! List of junction nodes of the network
	vector< node > Jv;
	
	////////////////////////////////////////////////////////////////////

	//! Monolithic matrix for the coupled problem
	sparse_matrix_type AM;
	//! Monolithic array of unknowns for the coupled problem
	vector_type        UM;
	//! Monolithic right hand side for the coupled problem
	vector_type        FM;

	////////////////////////////////////////////////////////////////////
	
	// Aux methods for init
	//! Import algorithm specifications
	void import_data(void);
	//! Import mesh for tissue (3D) and vessel (1D)  
	void build_mesh(void); 
	//! Set finite elements methods and integration methods 
	void set_im_and_fem(void);
	//! Build problem parameters
	void build_param(void);
	//! Build the list of tissue boundary data 
	/*!	Face numbering:
		  0 : {x = 0 }  "back"
		  1 : {x = Lx}  "front"
		  2 : {y = 0 }  "left"
		  3 : {y = Ly}  "right"
		  4 : {z = 0 }  "bottom"
		  5 : {z = Lz}  "top"
	 */
	void build_tissue_boundary(void);
	//! Build the list of vessel boundary (and junctions) data 
	void build_vessel_boundary(void);
	// Aux methods for assembly
	//! Build the monolithic matrix AM by blocks
	void assembly_mat(void);
	//! Build the monolithic rhs FM by blocks
	void assembly_rhs(void);
	//! Assemble RHS source term for stand-alone tissue problem
	void assembly_tissue_test_rhs(void);

}; /* end of class problem3d1d */


} /* end of namespace */

#endif
