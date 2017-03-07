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
  @brief  Main program for test simulations.
  @details
    We solve the coupled 3D/1D problem of fluid exchange between a 1D 
    network \Lambda and the 3D interstitial tissue \Omega
    
    *****************************************
      Benchmark : single-vessel network 
      Mixed finite elements approximation
      Monolithic resolution by SuperLU 3.0
    *****************************************
    
	See Section "Code verification: test-cases"
 */
#include <iostream>
#include <getfem/bgeot_config.h> // for FE_ENABLE_EXCEPT
#include <problem3d1d.hpp>
#include <c_problem3d1d.hpp>

using namespace getfem;

int main(int argc, char * argv[]){


	GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
	FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

	try {  
		c_problem3d1d prob;

		prob.init(argc, argv);
		prob.assembly();
		prob.solve();
		prob.export_vtk();

		// Display some global results: mean pressures, total flow rate
		std::cout << "--- FINAL RESULTS -------------------------" << std::endl; 
		std::cout << "  Pt average            = " << prob.mean_pt()   << std::endl;
		std::cout << "  Pv average            = " << prob.mean_pv()   << std::endl;
		std::cout << "  Network-to-Tissue TFR = " << prob.flow_rate() << std::endl;
		std::cout << "-------------------------------------------" << std::endl; 
	}

	GMM_STANDARD_CATCH_ERROR;

	return 0;
}

