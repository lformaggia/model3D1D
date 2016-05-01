/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
       Copyright (C) 2015 Domenico Notaro, 2014 Laura Cattaneo
======================================================================*/
/*! 
  @file   mesh1d.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   April 2015.
  @brief  Miscelleanous handlers for 1D network mesh.
  @details
  Create the edge sequence and build the related 1D mesh. \n
  The input stream ist is used to read a file with the following format (gambit-like): \n

		BEGIN_LIST
		BEGIN_ARC 
		BC KEYWA [VALUES] 
		BC KEYWB [VALUES]
		 106       0.4421      -1.6311       2.5101		start
		 107       0.4421      -1.6311       7.5101		end  
		 108       0.3546      -1.6524       2.5539		point
		 109       0.2621      -1.6695       2.5880		point
		... 
		END_ARC 
		... 
		BEGIN_ARC  
		... 
		END_ARC 
		... 
		END_LIST 

  1. The list of points of the IS ordered as follows:
     - first we have the coordinates of TWO ENDS (A,B) (i.e. A=start and B=end)
     - then we have all the remaining nodes of the arc, from A to B
  2. If a node is shared by more arcs, the arcs will be CONNECTED in the resulting 1D mesh.
  3. BC KEYWA [VALUES] and BC KEYWB [VALUES] are keywords/values related to boundary conditions
     to be imposed at nodes A, B. Each KEYW [VALUES] entry can be one of the following: \n
     - BC DIR 1.1
     - BC MIX 
     - BC INT \n
     Correspondingly, each node will be marked with the associated boundary condition, that are:
     - Dirichlet node (pressure = 1.1)
     - Mixed     node (flux = coef*(pressure - p0))
     - Internal  node
     At this stage, this is only meant to assign such BC labels to each node.
     If one end is INTERNAL, the corresponding BC will be ignored 
     (for clarity, please write the INT keyword).     
*/
/*!
	\defgroup geom Problem geometry
 */

#ifndef M3D1D_MESH_1D_HPP_
#define M3D1D_MESH_1D_HPP_

#include <node.hpp>

namespace getfem {

/*!
	Import the network points from the file of points (pts) and build the mesh.

	\ingroup geom
 */
//! \note It also build vessel mesh regions (#=0 for branch 0, #=1 for branch 1, ...).
template<typename VEC>
void 
import_pts_file(
		std::istream & ist, 
		getfem::mesh & mh1D, 
		std::vector<getfem::node> &  BCList,
		VEC & Nn,
		const std::string & MESH_TYPE
		) 
{

	size_type Nb = 0; // nb of branches
	Nn.resize(0); Nn.clear();
	mh1D.clear();
	
	ist.precision(16);
	ist.seekg(0); ist.clear();
	GMM_ASSERT1(bgeot::read_until(ist, "BEGIN_LIST"), 
		"This seems not to be a data file");

	size_type globalBoundaries = 0;

	while (bgeot::read_until(ist, "BEGIN_ARC")) {
	
		Nb++;
		Nn.emplace_back(0);

		std::vector<base_node> lpoints; 

		dal::dynamic_array<scalar_type> tmpv;
		std::string tmp, BCtype, value;
		bool thend = false; 
		size_type bcflag = 0;
		size_type bcintI = 0, bcintF = 0;
		node BCA, BCB;

		// Read an arc from data file and write to lpoints
		while (!thend) {
			bgeot::get_token(ist, tmp, 1023);
			if (tmp.compare("END_ARC") == 0) { 
				thend = true;
			}
			else if (ist.eof()) {
				GMM_ASSERT1(0, "Unexpected end of stream");
			}
			else if (tmp.compare("BC") == 0) { 
				bcflag++;
				bgeot::get_token(ist, BCtype, 4);
				if (BCtype.compare("DIR") == 0) {
					bgeot::get_token(ist, value, 1023);
					if (bcflag == 1) {
						BCA.label = BCtype; 
						BCA.value = stof(value); 
						//BCA.ind = globalBoundaries;
						globalBoundaries++;
					}
					else if (bcflag == 2) {
						BCB.label = BCtype; 
						BCB.value = stof(value); 
						//BCB.ind = globalBoundaries;
						globalBoundaries++;
					}
					else
						GMM_ASSERT1(0, "More than 2 BC found on this arc!");
				}
				else if (BCtype.compare("MIX") == 0) {
					bgeot::get_token(ist, value, 1023);
					if (bcflag == 1) {
						BCA.label = BCtype; 
						BCA.value = stof(value); 
						//BCA.ind = globalBoundaries;
						globalBoundaries++;
					}
					else if (bcflag == 2) {
						BCB.label = BCtype; 
						BCB.value = stof(value); 
						//BCB.ind = globalBoundaries;
						globalBoundaries++;
					}
				}
				else if (BCtype.compare("INT") == 0) {
					if (bcflag == 1) {
						bcintI++;
						BCA.label = BCtype; 
						//BCA.value = stof(value); //Error: no number to read
					}
					else if (bcflag == 2) {
						bcintF++;
						BCB.label = BCtype; 
						//BCB.value = stof(value); //Error: no number to read
					}
					else
						GMM_ASSERT1(0, "More than 2 BC found on this arc!");
				}
				else
					GMM_ASSERT1(0, "Unknown Boundary condition");	  
			
			} /* end of "BC" case */
			else if (tmp.size() == 0) {
				GMM_ASSERT1(0, "Syntax error in file, at token '" 
							 << tmp << "', pos=" << std::streamoff(ist.tellg()));
			} 
			else { /* "point" case */
				Nn[Nb-1]++;
				int d = 0;
				while ( (isdigit(tmp[0]) != 0) || tmp[0] == '-' || tmp[0] == '+' || tmp[0] == '.'){ 
					tmpv[d++] = stof(tmp); 
					bgeot::get_token(ist, tmp, 1023); 
				}
				if (d != 4) GMM_ASSERT1(0, "Points must have 3 coordinates");
				base_node tmpn(tmpv[1], tmpv[2], tmpv[3]);
				lpoints.push_back(tmpn);
				if (tmp.compare("END_ARC") == 0) { thend = true; Nn[Nb-1]--; }
			} 
						
		} /* end of inner while */
		
		// Insert the arc into the 1D mesh and build a new region for the corresponding branch
		// Check validity of branch region
		GMM_ASSERT1(mh1D.has_region(Nb-1)==0, "Overload in meshv region assembling!");
	
		for (size_type i=0; i<lpoints.size()-1; ++i ) {
			std::vector<size_type> ind(2);
			size_type ii = (i==0) ? 0 : i+1;
			size_type jj;
			
			if (ii == lpoints.size()-1) jj = 1;
			else if (ii == 0) jj = 2;
			else jj = ii+1;
			
			ind[0] = mh1D.add_point(lpoints[ii]);
			ind[1] = mh1D.add_point(lpoints[jj]);
			size_type cv;
			cv = mh1D.add_convex(bgeot::geometric_trans_descriptor(MESH_TYPE), ind.begin());
		
			// Build branch regions
			mh1D.region(Nb-1).add(cv);
			
			if ((bcflag>0) && (ii==0)&& (bcintI==0)) {
				BCA.idx = ind[0];
				BCList.push_back(BCA);
			}
			if ((bcflag>1) && (jj==1) && (bcintF==0)) {
				BCB.idx = ind[1];
				BCList.push_back(BCB);
			}

		} /* end of inner for */
		
	} /* end of outer while */	
		
} /* end of import_pts_file */


/*!
	Compute the cartesian components of the network tangent versor 
	@f$\mathbf{\lambda}@f$.

	\ingroup geom
 */
template<typename VEC>
void 
asm_tangent_versor(
		std::istream & ist, 
		VEC & lambdax, 
		VEC & lambday, 
		VEC & lambdaz
		) 
{
	lambdax.resize(0); lambdax.clear();
	lambday.resize(0); lambday.clear();
	lambdaz.resize(0); lambdaz.clear();
	
	ist.precision(16);
	ist.seekg(0); ist.clear();
	GMM_ASSERT1(bgeot::read_until(ist, "BEGIN_LIST"), "This seems not to be a data file");

	int Nb = 0;
	
	// Read a branch
	while (bgeot::read_until(ist, "BEGIN_ARC")) {
		
		Nb++;
		std::vector<base_node> lpoints; 
		std::vector<scalar_type> tmpv(4);
		bool thend=false;
		std::string tmp;
		int d=0;
	
		while (!thend) {
			bgeot::get_token(ist, tmp, 1023);
			if (tmp.compare("END_ARC") == 0) { 
				thend = true;
			}
			else if(tmp.find("BC")!=std::string::npos){
				bgeot::get_token(ist, tmp, 1023);
				if(tmp.find("DIR")!=std::string::npos){
					bgeot::get_token(ist, tmp, 1023);
				}
				else if(tmp.find("MIX")!=std::string::npos){
					bgeot::get_token(ist, tmp, 1023);
				}
				else if(tmp.find("INT")!=std::string::npos){
					//
				}
				else GMM_ASSERT1(0, "Syntax error in file, at token '" << tmp << "', pos=" << std::streamoff(ist.tellg()));
			}
			else if (tmp.size() == 0) {
				GMM_ASSERT1(0, "Syntax error in file, at token '" << tmp << "', pos=" << std::streamoff(ist.tellg()));
			} 
			else { /* "point" case */
				int i=d++%5;
				if(i==1 || i==2 || i==3) tmpv[i] = stof(tmp);
				if(i==3){
					base_node tmpn(tmpv[1], tmpv[2], tmpv[3]);
					lpoints.push_back(tmpn);
				}
				if (tmp.compare("END_ARC") == 0) thend = true;
			} 
		}
		// Compute the tangent versor
		// It is assumed to be constant over each branch
		scalar_type lx = lpoints[1][0]-lpoints[0][0];
		scalar_type ly = lpoints[1][1]-lpoints[0][1];
		scalar_type lz = lpoints[1][2]-lpoints[0][2];
		scalar_type lnorm = sqrt(lx*lx+ly*ly+lz*lz);
		
		lambdax.emplace_back(lx/lnorm);
		lambday.emplace_back(ly/lnorm);
		lambdaz.emplace_back(lz/lnorm);

	}
	
}

/*!
	Import the network radius @f$R(s)=\sum_{i=1}^N R_i~\delta_{\Lambda_i}(s)@f$.
	\ingroup geom
 */
//! \note It assumes that a proper file of points (pts) does exist
template<typename VEC>
void
import_network_radius
	(VEC & Radius,
 	 std::istream & ist, 
 	 const mesh_fem & mf_data
 	 ) 
{
	// Try to read data from pts file
	vector_type Rdata;
	ist.precision(16);
	ist.seekg(0); ist.clear();
	GMM_ASSERT1(bgeot::read_until(ist, "BEGIN_LIST"), "This seems not to be a data file");
	std::string line;
	size_type nb_branches = 0;
	bool thend = false; 
	while (!thend){
		bgeot::get_token(ist, line, 1023);
		thend = (line=="END_LIST");
		if (!thend){
			Rdata.emplace_back(stof(line));
			nb_branches++;
		} 
	}	
	// Assemble the POdisc radius
	gmm::resize(Radius, mf_data.nb_dof()); 
	gmm::clear(Radius);
	for (size_type b=0; b<nb_branches; ++b)
		for (getfem::mr_visitor mrv(mf_data.linked_mesh().region(b)); !mrv.finished(); ++mrv)
			for (auto i : mf_data.ind_basic_dof_of_element(mrv.cv()))
				Radius[i] = Rdata[b];
}


} /* end of namespace */
#endif
