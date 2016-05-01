/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   descr3d1d.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   January 2016.
  @brief  Definition of the aux class for algorithm description strings.
 */
/** @defgroup input User-defined parameters  */

#ifndef M3D1D_DESCR3D1D_HPP_
#define M3D1D_DESCR3D1D_HPP_

#include <string>

namespace getfem {

//! Class to import the descriptors of the coupled 3D/1D solver
/*!
	\ingroup input
 */
struct descr3d1d {

	// GetFEM specifications
	//! Absolute path to the tissue mesh file
	std::string MESH_FILET;
	//! Absolute path to the vessel mesh file
	std::string MESH_FILEV;
	//! Identifier of tissue mesh tipe
	std::string MESH_TYPET;
	//! Identifier of vessel mesh tipe
	std::string MESH_TYPEV;
	//! Identifier of tissue velocity's FEM type
	std::string FEM_TYPET;
	//! Identifier of tissue pressure's FEM type
	std::string FEM_TYPET_P;
	//! Identifier of tissue coefficients' FEM type
	std::string FEM_TYPET_DATA;
	//! Identifier of tissue velocity's FEM type
	std::string FEM_TYPEV;
	//! Identifier of vessel pressure's FEM type
	std::string FEM_TYPEV_P;
	//! Identifier of vessel coefficients' FEM type
	std::string FEM_TYPEV_DATA;
	//! Identifier of tissue integration method type
	std::string IM_TYPET;
	//! Identifier of vessel integration method type
	std::string IM_TYPEV;
	//! Output directory
	std::string OUTPUT;	
	// Solver information
	//! Identifief of the monolithic solver
	std::string SOLVE_METHOD;
	//! Maximum number of iterations (iterative solvers)
	size_type   MAXITER;
	//! Mamimum residual (iterative solvers)
	scalar_type RES; 
	//! Number of target points for the tissue-to-vessel average
	size_type   NInt;
	// Utils
	//! File .param
	ftool::md_param FILE_;
	//! Import algorithm specifications from file .param
	void import(ftool::md_param & fname) 
	{
		FILE_ = fname;
		MESH_FILET  = FILE_.string_value("MESH_FILET","3D mesh file");
		MESH_FILEV  = FILE_.string_value("MESH_FILEV","1D points file");
		MESH_TYPET  = FILE_.string_value("MESH_TYPET","3D mesh type");
		MESH_TYPEV  = FILE_.string_value("MESH_TYPEV","1D mesh type");
		FEM_TYPET   = FILE_.string_value("FEM_TYPET","FEM 3D tissue - velocity");
		FEM_TYPET_P = FILE_.string_value("FEM_TYPET_P","FEM 3D tissue - pressure");
		FEM_TYPET_DATA = FILE_.string_value("FEM_TYPET_DATA");
		FEM_TYPEV   = FILE_.string_value("FEM_TYPEV","FEM 1D vessel - velocity");
		FEM_TYPEV_P = FILE_.string_value("FEM_TYPEV_P","FEM 1D vessel - pressure");
		FEM_TYPEV_DATA = FILE_.string_value("FEM_TYPEV_DATA");
		IM_TYPET 	= FILE_.string_value("IM_TYPET","Name of integration method");
		IM_TYPEV 	= FILE_.string_value("IM_TYPEV","Name of integration method");
		if(FEM_TYPET_DATA=="") FEM_TYPET_DATA = "FEM_PK(3,0)";
		if(FEM_TYPEV_DATA=="") FEM_TYPEV_DATA = "FEM_PK(1,0)";

		SOLVE_METHOD = FILE_.string_value("SOLVE_METHOD", "Monolithic Solver"); 
		if (SOLVE_METHOD != "SuperLU") { // iterative solver
			MAXITER  = FILE_.int_value("MAXITER", "Max number of sub-iterations");
			RES = FILE_.real_value("RES"); if (RES == 0.) RES = 2.0e-10;
		}
		NInt = size_type(FILE_.int_value("NInt", "Node numbers on the circle for the nonlocal term"));  
		OUTPUT = FILE_.string_value("OUTPUT","Output Directory");
	}

	//! Overloading of the output operator
	friend std::ostream & operator << (
		std::ostream & out, const descr3d1d & descr
		)
	{ 
		cout << "---- PROBLEM DESCRIPTORS--------------------------" << endl;
		cout << " MESH FILE 3D problem      : " << descr.MESH_FILET  << endl;
		cout << " MESH FILE 1D problem      : " << descr.MESH_FILEV  << endl;
		cout << " IM  TYPE  3D problem      : " << descr.IM_TYPET    << endl;
		cout << " IM  TYPE  1D problem      : " << descr.IM_TYPEV	<< endl;
		cout << " FEM TYPE  3D velocity     : " << descr.FEM_TYPET   << endl;
		cout << " FEM TYPE  3D pressure     : " << descr.FEM_TYPET_P << endl;
		cout << " FEM TYPE  3D coefficients : " << descr.FEM_TYPET_DATA << endl;
		cout << " FEM TYPE  1D velocity     : " << descr.FEM_TYPEV   << endl;
		cout << " FEM TYPE  1D pressure     : " << descr.FEM_TYPEV_P << endl;
		cout << " FEM TYPE  1D coefficients : " << descr.FEM_TYPEV_DATA << endl;
		cout << "--------------------------------------------------" << endl;

		return out;            
	}

}; /* end of class */

} /* end of namespace */

#endif
