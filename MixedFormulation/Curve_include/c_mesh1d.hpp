#ifndef CURVED_M3D1D_MESH_1D_HPP_
#define CURVED_M3D1D_MESH_1D_HPP_

#include <node.hpp>
#include <defines.hpp>

namespace getfem {

template<typename VEC, typename PARAM>
void 
import_pts_file(
		std::istream & ist, 
		std::istream & istc,
		std::istream & istt,
		getfem::mesh & mh1D, 
		std::vector<getfem::node> &  BCList,
		VEC & Nn,
		const std::string & MESH_TYPE,
		PARAM & param
		) 
{
	size_type Nold=0;
	size_type Nnode = 0; // Total number of node
	size_type Nb = 0; // nb of branches
	Nn.resize(0); Nn.clear();
	mh1D.clear();
	
	vector<vector_type> lx;
	vector<vector_type> ly;
	vector<vector_type> lz;
	vector<vector_type> Curv;

  //////////////////////////////////////////////////////////////////// IFSTREAM
	ist.precision(16);
	ist.seekg(0); ist.clear();
	istc.precision(16);
	istc.seekg(0); istc.clear();
	istt.precision(16);
	istt.seekg(0); istt.clear();
	////////////////////////////////////////////////////////////////// END IFSTREAM


	GMM_ASSERT1(bgeot::read_until(ist, "BEGIN_LIST"), 
		"This seems not to be a data file");
	GMM_ASSERT1(bgeot::read_until(istc, "BEGIN_LIST"), 
		"This seems not to be a data file");
	GMM_ASSERT1(bgeot::read_until(istt, "BEGIN_LIST"), 
		"This seems not to be a data file");


	size_type globalBoundaries = 0;

	while (bgeot::read_until(ist, "BEGIN_ARC") && bgeot::read_until(istc, "BEGIN_ARC") && bgeot::read_until(istt, "BEGIN_ARC")) {
	    ///////////////////////////////////////////////////////////////ADATTAZIONE VARIABILI
		Nb++;
		Nn.emplace_back(0);
		std::vector<base_node> lpoints;
		vector_type Curv_b;
		vector_type lx_b;
		vector_type ly_b;
		vector_type lz_b;
		
		dal::dynamic_array<scalar_type> tmpv;

		std::string tmp,BCtype, value;
		bool thend = false; 
		bool thendc= false;
		bool thendt= false;
		size_type bcflag = 0;
		size_type bcintI = 0, bcintF = 0;
		node BCA, BCB;
		///////////////////////////////////////////////////////////////FINE ADATTAZIONE VARIABILI

		///////////////////////////////////////////////////////////////LETTURA DATI PER IL RAMO
		// Read an arc from data file and write to lpoints
		
		/////////////////////////////////////////////////////////////LETTURA PUNTI MESH
		
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
		
		//////////////////////////////////////////////////////////////LETTURA PUNTI VETTORE NORMALE
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
						cout<<"\n\n"<<tmp<<"\n\n";
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
					if (d != 4) GMM_ASSERT1(0, "Points must have 3 coordinates");
					base_node tmpn(tmpv[1], tmpv[2], tmpv[3]);
					Curv_b.push_back(gmm::vect_norm2(tmpn));
					if (tmp.compare("END_ARC") == 0) { thendc = true; }
				} 			
		}
		/////////////////////////////////////////////////////////////FINE LETTURA PUNTI VETTORE NORMALE
		/////////////////////////////////////////////////////////////LETTURA PUNTI VERSORE TANGENTE
		while(!thendt){
				bgeot::get_token(istt, tmp, 1023);
				if (tmp.compare("END_ARC") == 0) { 
					thendt = true;
				}
				else if (ist.eof()) {
					GMM_ASSERT1(0, "Unexpected end of stream");
				}
				else if (tmp.compare("BC") == 0) { 
					bcflag++;
					bgeot::get_token(istt, BCtype, 4);
					if (BCtype.compare("DIR") == 0) {
						bgeot::get_token(istt, value, 1023);
					}
					else if (BCtype.compare("MIX") == 0) {
						bgeot::get_token(istt, value, 1023);
					}
					else if (BCtype.compare("INT") == 0) {
					}
					else{
						cout<<"\n\n"<<tmp<<"\n\n";
						GMM_ASSERT1(0, "Unknown Boundary condition");	  
					}
				} /* end of "BC" case */
				else if (tmp.size() == 0) {
					GMM_ASSERT1(0, "Syntax error in file, at token '" 
								 << tmp << "', pos=" << std::streamoff(ist.tellg()));
				} 
				else { /* "tangent versor" case */
					int d = 0;
					while ( (isdigit(tmp[0]) != 0) || tmp[0] == '-' || tmp[0] == '+' || tmp[0] == '.'){ 
						tmpv[d++] = stof(tmp); 
						bgeot::get_token(istt, tmp, 1023); 
					}
					if (d != 4) GMM_ASSERT1(0, "Points must have 3 coordinates");
					lx_b.push_back(tmpv[1]);
					ly_b.push_back(tmpv[2]);
					lz_b.push_back(tmpv[3]);
					if (tmp.compare("END_ARC") == 0) { thendc = true; }
				} 			
		}
		/////////////////////////////////////////////////////////////FINE LETTURA PUNTI VERSORE TANGENTE

		GMM_ASSERT1(lpoints.size() == Curv_b.size(), 
			"The point file contain different number of element then the normal versor file on the branch "<<Nb);	
		GMM_ASSERT1(lpoints.size() == lx_b.size(), 
			"The point file contain different number of element then the tangential velocity versor file on the branch "<<Nb);		
		///////////////////////////////////////////////////////////////FINE LETTURA DATI PER IL RAMO

		// Insert the arc into the 1D mesh and build a new region for the corresponding branch
		// Check validity of branch region
		GMM_ASSERT1(mh1D.has_region(Nb-1)==0, "Overload in meshv region assembling!");
		Curv.push_back(Curv_b);
		lx.push_back(lx_b);
		ly.push_back(ly_b);
		lz.push_back(lz_b);

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
       
			/////////////////////////////////////////////////////////////////////////////////BOUNDARY
			if ((bcflag>0) && (ii==0)&& (bcintI==0)) {
				BCA.idx = ind[0];
				BCList.push_back(BCA);
			}
			if ((bcflag>1) && (jj==1) && (bcintF==0)) {
				BCB.idx = ind[1];
				BCList.push_back(BCB);
			}
			/////////////////////////////////////////////////////////////////////////////////END BOUNDARY

		} /* end of inner for */
		if(0){
			cout<<"Al branch "<<Nn.size()<<" Curv=[";
			for(size_type i=0; i<Curv[Nb-1].size();++i) cout<<" "<<Curv[Nb-1][i];
			cout<<" ];"<<endl;
		}
	} /* end of outer while */	
	//cout<<endl;
	param.get_curve(Curv,lx,ly,lz);
} 

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
	size_type Nb=Curv_temp.size();
	size_type Ni=0;
	pfem pf_coefvi = fem_descriptor("FEM_PK(1,1)");

	//Adattando per ogni ramo i parametri curvi alla mesh in ogni ramo
	for(size_type b=0; b<Nb;++b){

		getfem::mesh_fem mf_tmp(mf_Coefi[b].linked_mesh());
		mf_tmp.set_finite_element(mf_Coefi[b].linked_mesh().region(b).index(), pf_coefvi);
		Ni=mf_Coefi[b].nb_dof();
		Curv[b].resize(Ni); lx[b].resize(Ni); ly[b].resize(Ni); lz[b].resize(Ni); 
		gmm::row_matrix<vector_type> M(mf_Coefi[b].nb_dof(),mf_tmp.nb_dof());
		getfem::interpolation(mf_tmp,mf_Coefi[b],M);

		gmm::mult(M,Curv_temp[b],Curv[b]);
		gmm::mult(M,lx_temp[b],lx[b]);
		gmm::mult(M,ly_temp[b],ly[b]);
		gmm::mult(M,lz_temp[b],lz[b]);
	}
	if(0){
		for(size_type b=0; b<Nb;++b){
			cout<<"Curv_new["<<b+1<<"] = [";
			for(size_type i=0; i<Curv[b].size();++i){	cout<<" "<<Curv[b][i];}
			cout<<" ]\n\n";
		}
		for(size_type b=0; b<Nb;++b){
			cout<<"Curv_old["<<b+1<<"] = [";
			for(size_type i=0; i<Curv_temp[b].size();++i){	cout<<" "<<Curv_temp[b][i];}
			cout<<" ]\n\n";
		}
	}
}



} /* end of namespace */
#endif