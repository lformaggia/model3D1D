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
  @brief Miscelleaneous handlers for 1D curved network mesh and curved parameters.
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


     The input stream istc is used to read a file with the following format (gambit-like): \n

		BEGIN_LIST
		BEGIN_ARC 
		BC KEYWA [VALUES] 
		BC KEYWB [VALUES]
		 106       0.4421      -1.6311       2.5101		0.1221		0.3454		0.2943		0.7561 start
		 107       0.4421      -1.6311       7.5101		0.1221		0.3454		0.2943		0.7561 end  
		 108       0.3546      -1.6524       2.5539		0.2235		0.2541		0.3012		0.7214 point
		 109       0.2621      -1.6695       2.5880		0.2173		0.2497		0.2886		0.7493 point
		... 
		END_ARC 
		... 
		BEGIN_ARC  
		... 
		END_ARC 
		... 
		END_LIST 

  1. The list of points of the IS ordered as follows:
     - first we have the parameters for the TWO ENDS coordinate (A,B) (i.e. A=start and B=end)
     - then we have all the remaining parameters for the nodes of the arc, from A to B
  2. If a node is shared by more arcs, the arcs will be CONNECTED in the resulting 1D mesh.
  3. BC KEYWA [VALUES] and BC KEYWB [VALUES] are keywords/values related to boundary conditions
     to be imposed at nodes A, B. For this kind of problem is not required any boundary condition
     fot the curved parameters, because all the BC are imposed with the mesh file input. If the 
     reader would find any BC, it will simply avoid to read it.

*/
/*!
	\defgroup geom Problem geometry
 */    

#ifndef CURVED_M3D1D_MESH_1D_HPP_
#define CURVED_M3D1D_MESH_1D_HPP_

#include <node.hpp>
#include <defines.hpp>
#include <cmath>

namespace getfem {

/*!
	Import the network points and the curved paremeters from the file pts and build the 
	mesh.

	\ingroup geom
 */
//! \note It also build vessel mesh regions (#=0 for branch 0, #=1 for branch 1, ...).
template<typename VEC_SIZE, typename PARAM>
void 
import_pts_file(
		std::istream & ist, 
		std::istream & istc,
		getfem::mesh & mh1D, 
		std::vector<getfem::node> &  BCList,
		VEC_SIZE & Nn,
		const std::string & MESH_TYPE,
		PARAM & param
		) 
{
	size_type Nold=0;
	size_type Nnode = 0; // Total number of node
	size_type Nb = 0; // nb of branches
	Nn.resize(0); Nn.clear();
	mh1D.clear();
	
	//Vectors Storing the curved parameters
	vector<vector_type> lx;
	vector<vector_type> ly;
	vector<vector_type> lz;
	vector<vector_type> Curv;

	//Opening ifstream 
	ist.precision(16);
	ist.seekg(0); ist.clear();
	istc.precision(16);
	istc.seekg(0); istc.clear();

	GMM_ASSERT1(bgeot::read_until(ist, "BEGIN_LIST"), 
		"This seems not to be a data file");
	GMM_ASSERT1(bgeot::read_until(istc, "BEGIN_LIST"), 
		"This seems not to be a data file");

	size_type globalBoundaries = 0;

	while (bgeot::read_until(ist, "BEGIN_ARC") && bgeot::read_until(istc, "BEGIN_ARC")) { //Reading an arc for time
		Nb++;
		Nn.emplace_back(0);
		//temporary vectors for read point and curved parameters for the branch
		std::vector<base_node> lpoints;
		vector_type Curv_b;
		vector_type lx_b;
		vector_type ly_b;
		vector_type lz_b;
		//temporary value used to compute euclidean norm
		scalar_type lx_t;
		scalar_type ly_t;
		scalar_type lz_t;
		scalar_type lmod;
		
		dal::dynamic_array<scalar_type> tmpv;

		std::string tmp,BCtype, value;
		bool thend = false; //End of reading an arc for points
		bool thendc= false; //End of reading an arc for curved paramenters
		size_type bcflag = 0;
		size_type bcintI = 0, bcintF = 0; //Boundary condition
		node BCA, BCB; //Boundary node

		// Read an arc from data file and write to lpoints		
		while(!thend){
		
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
		}
		
		///////////////////////////////////////////////////////////////FINE LETTURA PUNTI MESH
		
		//////////////////////////////////////////////////////////////LETTURA VERSORE TANGENTE E CURVATURA
		while (!thendc){
				bgeot::get_token(istc, tmp, 1023);
				if (tmp.compare("END_ARC") == 0) { 
					thendc = true;
				}
				else if (ist.eof()) {
					GMM_ASSERT1(0, "Unexpected end of stream");
				}
				else if (tmp.compare("BC") == 0) { 
					bcflag++;
					bgeot::get_token(istc, BCtype, 4);
					if (BCtype.compare("DIR") == 0) {
						bgeot::get_token(istc, value, 1023);
					}
					else if (BCtype.compare("MIX") == 0) {
						bgeot::get_token(istc, value, 1023);
					}
					else if (BCtype.compare("INT") == 0) {
					}
					else{
						GMM_ASSERT1(0, "Unknown Boundary condition");	  
					}
				} /* end of "BC" case */
				else if (tmp.size() == 0) {
					GMM_ASSERT1(0, "Syntax error in file, at token '" 
								 << tmp << "', pos=" << std::streamoff(ist.tellg()));
				} 
				else { /* "curvature" case */
					int d = 0;
					while ( (isdigit(tmp[0]) != 0) || tmp[0] == '-' || tmp[0] == '+' || tmp[0] == '.'){ 
						tmpv[d++] = stof(tmp); 
						bgeot::get_token(istc, tmp, 1023);
					}
					if (d != 8) GMM_ASSERT1(0, "Points must have 7 coordinates: 3 for tangent versor, 3 for normal vector, 1 for curvature");
					Curv_b.push_back(tmpv[7]);
					if (tmp.compare("END_ARC") == 0) { thendc = true; }
				} 			
		}

		GMM_ASSERT1(lpoints.size() == Curv_b.size(), 
			"The point file contain different number of element then the normal versor file on the branch "<<Nb);	
		// Insert the arc into the 1D mesh and build a new region for the corresponding branch
		// Check validity of branch region
		GMM_ASSERT1(mh1D.has_region(Nb-1)==0, "Overload in meshv region assembling!");

		//Generating mesh for the branch
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


			lx_t=lpoints[jj][0]-lpoints[ii][0];
			ly_t=lpoints[jj][1]-lpoints[ii][1];
			lz_t=lpoints[jj][2]-lpoints[ii][2];
			lmod=sqrt(lx_t*lx_t+ly_t*ly_t+lz_t*lz_t);

			lx_b.push_back(lx_t/lmod);
			ly_b.push_back(ly_t/lmod);
			lz_b.push_back(lz_t/lmod);

		} /* end of inner for */
		Curv.push_back(Curv_b);
		lx.push_back(lx_b);
		ly.push_back(ly_b);
		lz.push_back(lz_b);

	} /* end of outer while */	
	param.get_curve(Curv,lx,ly,lz);//Saving the parameters
} 

/*!
	Reassemble the curved parameters on the data finite elements.
	The parameters are given in the coordinates of the mesh, so it is
	assumed that the parameters are imported using polynom of degree 1.

	/ingroup geom
*/
template<typename VEC, typename VEC_FEM>
void rasm_curve_parameter(
		VEC_FEM & mf_Coefi,
		VEC & Curv,
		VEC & lx,
		VEC & ly,
		VEC & lz
	)
{
	VEC Curv_temp=Curv;
	VEC lx_temp=lx;
	VEC ly_temp=ly;
	VEC lz_temp=lz;
	size_type Nb=Curv_temp.size(); //Number of branch
	size_type Ni=0; //Number of dof of the coefficient mesh at branch i
	pfem fd_curv = fem_descriptor("FEM_PK(1,1)"); //Fem descriptor for curvature
	pfem fd_tgv  = fem_descriptor("FEM_PK(1,0)"); //Fem descriptor for tangent versor
	scalar_type ct;

	for(size_type b=0; b<Nb;++b){
		Ni=Curv_temp[b].size();
	
		//Reodering the value of the parameters
		ct=Curv[b][1];
		   //Shifting all elements
		for(size_type i=2; i<Curv_temp[b].size(); i++){
			Curv_temp[b][i-1]=Curv[b][i];
		}
		Curv_temp[b].back()=ct;

		//Adapting parameters to tbe finite element interpolation
		Ni=mf_Coefi[b].nb_dof();
		Curv[b].resize(Ni); gmm::clear(Curv[b]);
		 lx[b].resize(Ni);  gmm::clear(lx[b]);
		 ly[b].resize(Ni);  gmm::clear(lx[b]);
		 lz[b].resize(Ni);  gmm::clear(lx[b]);
		
		//Generating the P1 finite element for curvature
		getfem::mesh_fem mf_curv(mf_Coefi[b].linked_mesh());
		mf_curv.set_finite_element(mf_Coefi[b].linked_mesh().region(b).index(), fd_curv);

		//Generating the P0 finite element for tangent versor
		getfem::mesh_fem mf_tgv(mf_Coefi[b].linked_mesh());
		mf_tgv.set_finite_element(mf_Coefi[b].linked_mesh().region(b).index(), fd_tgv);

		//Interpolating v curvature on a different mesh fem 
		gmm::row_matrix<vector_type> M_curv(mf_Coefi[b].nb_dof(),mf_curv.nb_dof());
		getfem::interpolation(mf_curv,mf_Coefi[b],M_curv);

		gmm::mult(M_curv,Curv_temp[b],Curv[b]);



		//Interpolating v curvature on a different mesh fem 
		gmm::row_matrix<vector_type> M_tgv(mf_Coefi[b].nb_dof(),mf_tgv.nb_dof());
		getfem::interpolation(mf_tgv,mf_Coefi[b],M_tgv);

		gmm::mult(M_tgv,lx_temp[b],lx[b]);
		gmm::mult(M_tgv,ly_temp[b],ly[b]);
		gmm::mult(M_tgv,lz_temp[b],lz[b]);

		gmm::clear(M_tgv);
		gmm::clear(M_curv);
	}
}

/*! 
	Fuction to Compute the curvature given 3 nodes
	K=4*Apar/(d(x,y)*d(y,z)*d(z,x))
*/
template<typename Point3d>
scalar_type
curvature3d(
	Point3d & X,
	Point3d & Y,
	Point3d & Z
	)
{

	scalar_type dXY,dYZ, dZX , K;
	scalar_type A=0.0;
	vector_type XY(3);
	vector_type XZ(3);
	
	//compute distance of Points
	dXY=sqrt((X[1]-Y[1])*(X[1]-Y[1])+(X[2]-Y[2])*(X[2]-Y[2])+(X[3]-Y[3])*(X[3]-Y[3]));
	dYZ=sqrt((Z[1]-Y[1])*(Z[1]-Y[1])+(Z[2]-Y[2])*(Z[2]-Y[2])+(Z[3]-Y[3])*(Z[3]-Y[3]));
	dZX=sqrt((X[1]-Z[1])*(X[1]-Z[1])+(X[2]-Z[2])*(X[2]-Z[2])+(X[3]-Z[3])*(X[3]-Z[3]));

	//Compute Area of the parallelogram
	for(int i=0;i<3;i++){
		XY[i]=Y[i]-X[i];
		XZ[i]=Z[i]-X[i];
	}

	for(int i=0; i<3;i++){
		A=A+(XY[i]*XZ[(i+1)%3]-XY[(i+1)%3]*XZ[i]);
	}
	A=A/2.0;

	//Compute Curvature
	K=4.0*A/dXY/dYZ/dZX;
	return K;
}


/*!
	Import the network points from the file of points(pts), compute the curve parameters and build the mesh.

	\ingroup geom
 */
//! \note It also build vessel mesh regions (#=0 for branch 0, #=1 for branch 1, ...).
template<typename VEC, typename PARAM>
void 
import_pts_file(
		std::istream & ist, 
		getfem::mesh & mh1D, 
		std::vector<getfem::node> &  BCList,
		VEC & Nn,
		const std::string & MESH_TYPE,
		PARAM & param
		) 
{

	size_type Nb = 0; // nb of branches
	Nn.resize(0); Nn.clear();
	mh1D.clear();
	
	ist.precision(16);
	ist.seekg(0); ist.clear();
	GMM_ASSERT1(bgeot::read_until(ist, "BEGIN_LIST"), 
		"This seems not to be a data file");


	vector<vector_type> lx;
	vector<vector_type> ly;
	vector<vector_type> lz;
	vector<vector_type> Curv;

	size_type globalBoundaries = 0;

	while (bgeot::read_until(ist, "BEGIN_ARC")) {
	
		Nb++;
		Nn.emplace_back(0);

		std::vector<base_node> lpoints; 
		vector_type Curv_b;
		vector_type lx_b;
		vector_type ly_b;
		vector_type lz_b;
		//temporary value used to compute euclidean norm
		scalar_type lx_t;
		scalar_type ly_t;
		scalar_type lz_t;
		scalar_type lmod;

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
		scalar_type KK=0;
		for (size_type i=0; i<lpoints.size()-1; ++i ) {
			std::vector<size_type> ind(2);
			size_type ii = (i==0) ? 0 : i+1;
			size_type jj;
			size_type ll= (i==1) ? 0 : i;

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

			//Compute Curvature
			if(ii!=0 ){
				KK=curvature3d(lpoints[ll],lpoints[ii],lpoints[jj]);
				if(ii==2){
					Curv_b.push_back(KK);
				}
				Curv_b.push_back(KK);
				if(ii==lpoints.size()-1){
					Curv_b.push_back(KK);
				}
			}

			//Compute Tanget Versor
			lx_t=lpoints[jj][0]-lpoints[ii][0];
			ly_t=lpoints[jj][1]-lpoints[ii][1];
			lz_t=lpoints[jj][2]-lpoints[ii][2];
			lmod=sqrt(lx_t*lx_t+ly_t*ly_t+lz_t*lz_t);

			lx_b.push_back(lx_t/lmod);
			ly_b.push_back(ly_t/lmod);
			lz_b.push_back(lz_t/lmod);
		} /* end of inner for */
		Curv.push_back(Curv_b);
		lx.push_back(lx_b);
		ly.push_back(ly_b);
		lz.push_back(lz_b);
	} /* end of outer while */	
	param.get_curve(Curv,lx,ly,lz);//Saving the parameters
} /* end of import_pts_file */

} /* end of namespace */
#endif