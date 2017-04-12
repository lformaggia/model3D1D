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
  @brief  Definition of the aux class for algorithm description strings for curved problem.
 */
/*! @defgroup input User-defined parameters  */

#ifndef C_M3D1D_DESCR3D1D_HPP_
#define C_M3D1D_DESCR3D1D_HPP_

#include <string>

namespace getfem {

//! Class to import the descriptors of the coupled 3D/1D solver
/*!
	\ingroup input
 */
struct c_descr3d1d {

	// GetFEM specifications
	//! Absolute path to the tissue mesh file
	std::string MESH_CURVE;

	std::string FEM_TYPECURVE;

	std::string FIXED_POINT_METHOD;

	std::string INCREMENT_VELOCITY_TYPE;

	scalar_type MAXITER;

	scalar_type MINERROR;
	// Utils
	//! File .param
	ftool::md_param FILE_;
	//! Import algorithm specifications from file .param
	void import(ftool::md_param & fname) 
	{
		FILE_ = fname;

		MESH_CURVE = FILE_.string_value("CURVE_FILE","Path of the mesh of the curve parameters");
		FEM_TYPECURVE = FILE_.string_value("FEM_TYPEV_DATA","Name of Finite element on parameters");

		if(FILE_.int_value("Newton","Bool variable which say if you want to use Newton Method")==1)
			FIXED_POINT_METHOD="Newton Method";
		else
			FIXED_POINT_METHOD="Fixed Point Method";

		MAXITER = FILE_.real_value("Max_iter","Maximum number of iteration of Fixed Point Method");

		MINERROR = FILE_.real_value("minERR","Minimum value of increment for Convergence");
	}

	//! Overloading of the output operator
	friend std::ostream & operator << (
		std::ostream & out, const c_descr3d1d & descr
		)
	{ 
		out << "---- CURVE PROBLEM DESCRIPTORS--------------------" << endl;
		out << " MESH CURVE 1D problem          : " << descr.MESH_CURVE  << endl;
		out << " FEM TYPE   1D curve parameters : " << descr.FEM_TYPECURVE<< endl;
		out << " FIXED POINT Method             : " << descr.FIXED_POINT_METHOD  << endl;
		out << " MAX number of Iterations       : " << descr.MAXITER <<endl;
		out << " MIN increment accepted         : " << descr.MINERROR <<endl;
		out << "--------------------------------------------------" << endl;

		return out;            
	}

}; /* end of class */

} /* end of namespace */

#endif
