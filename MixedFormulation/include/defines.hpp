/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   defines.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   January 2016.
  @brief  Miscellaneous definitions for the 3D/1D coupling.
 */
#ifndef M3D1D_DEFINES_HPP_
#define M3D1D_DEFINES_HPP_

namespace getfem {

// Some useful abbreviations
using std::cout; 
using std::cin; 
using std::cerr; 
using std::ends; 
using std::endl; 
using std::vector;
using std::string;

// Useful type definitions (built using the predefined types in Gmm++)
//! Special class for small (dim < 16) vectors 
using bgeot::base_small_vector;  
//! Geometrical nodes (derived from base_small_vector)
using bgeot::base_node;   	 	 
//! Double-precision FP numbers 
using bgeot::scalar_type; 	 	
//! Unsigned long integers 
using bgeot::size_type;   	 	 
//! Short integers 
using bgeot::short_type;         
//! Type for vector of integers 
typedef std::vector<size_type> vector_size_type;
//! Type for dense vector
typedef std::vector<scalar_type> vector_type;
//! Type for sparse vector (std::vector) 
typedef gmm::rsvector<scalar_type> sparse_vector_type;
//! Type for dense matrix
typedef gmm::row_matrix<vector_type> matrix_type;
//! Type for sparse matrix
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;    


//! Definition of \pi
const scalar_type pi = std::atan(1)*4; 

} /* end of namespace */

#endif
