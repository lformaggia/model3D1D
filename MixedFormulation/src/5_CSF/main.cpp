/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   main.cpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   April 2015.
  @brief  Main program for CSF flow simulations.
  @detail 
    We solve the coupled 3D/1D problem of cerebrospinal fluid flow for 
    a schematic arterial-venous network \Lambda immersed in a 3D
    interstitial volume \Omega
    
    *****************************************
      Application : CSF flow 
      Mixed finite elements approximation
      Monolithic resolution by SuperLU 3.0
    *****************************************
    
	See Section "Cerebrospinal Fluid Flow in Human Brain"
 */

#include <iostream>
#include <problem3d1d.hpp>

//! main program
int main(int argc, char *argv[]) 
{

	GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
	FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

	try {   
		// Declare problems
		getfem::problem3d1d pba, pbv;
		// Initialize the problem
		pba.init(argc-1, argv);
		pbv.init(argc, argv);
		// Build the monolithic systems
		pba.assembly();
		pbv.assembly();
		// Solve the problem 
		merge_and_solve(pba, pbv);
		// Save results in .vtk format
		pbv.export_vtk("v");  
		pba.export_vtk("a"); // overwrite Ut, Pt
		// Display some global results: mean pressures, total flow rate
		std::cout << "--- FINAL RESULTS --------------------------" << std::endl; 
		std::cout << "  Pt average            = " << pba.mean_pt()   << std::endl;
		std::cout << "  Pav average           = " << pba.mean_pv()   << std::endl;
		std::cout << "  Network-to-Tissue TFR = " << pba.flow_rate() << std::endl;
		std::cout << "  Pvv average           = " << pbv.mean_pv()   << std::endl;
		std::cout << "  Network-to-Tissue TFR = " << pbv.flow_rate() << std::endl;
		std::cout << "--------------------------------------------" << std::endl; 	
	}

	GMM_STANDARD_CATCH_ERROR;

	return 0; 
	
} /* end of main program */
