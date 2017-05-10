/* -*- c++ -*- (enableMbars emacs c++ mode) */ 
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   problem3d1d.cpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   January 2016.
  @brief  Definition of the main class for the 3D/1D coupled problem.
 */

#include <problem3d1d.hpp>

namespace getfem {




/////////// Initialize the problem ///////////////////////////////////// 
void 
problem3d1d::init(int argc, char *argv[])
{
	//1. Read the .param filename from standard input
	PARAM.read_command_line(argc, argv);
	//2. Import data (algorithm specifications, boundary conditions, ...)
	import_data();
	//3. Import mesh for tissue (3D) and vessel network (1D)
	build_mesh();
	//4. Set finite elements and integration methods
	set_im_and_fem();
	//5. Build problem parameters
	build_param();
	//6. Build the list of tissue boundary data
	build_tissue_boundary();
	//7. Build the list of tissue boundary (and junction) data
	build_vessel_boundary();
}

void
problem3d1d::import_data(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Importing descriptors for tissue and vessel problems ..." << endl;
	#endif
	descr.import(PARAM);
	#ifdef M3D1D_VERBOSE_
	cout << descr;
	#endif
}


void
problem3d1d::build_mesh(void)
{
	bool test = 0;
	test = PARAM.int_value("TEST_GEOMETRY");
	if(test==0){
		#ifdef M3D1D_VERBOSE_
		cout << "Importing the 3D mesh for the tissue ...  "   << endl;
		#endif
		 import_msh_file(descr.MESH_FILET, mesht);
	}else{
		#ifdef M3D1D_VERBOSE_
		cout << "Building the regular 3D mesh for the tissue ...  "   << endl;
		#endif
		string st("GT='" + PARAM.string_value("GT_T") + "'; " +
					   "NSUBDIV=" + PARAM.string_value("NSUBDIV_T") + "; " +  
					   "ORG=" + PARAM.string_value("ORG_T") + "; " +  
					   "SIZES=" + PARAM.string_value("SIZES_T") + "; " +  
					   "NOISED=" + PARAM.string_value("NOISED_T")); 
		cout << "mesht description: " << st << endl;
		regular_mesh(mesht, st);
	}

	#ifdef M3D1D_VERBOSE_
	cout << "Importing the 1D mesh for the vessel ... "   << endl;
	#endif
	std::ifstream ifs(descr.MESH_FILEV);
	GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr.MESH_FILEV);
	import_pts_file(ifs, meshv, BCv, nb_vertices, descr.MESH_TYPEV);
	nb_branches = nb_vertices.size();
	ifs.close();
}

void
problem3d1d::set_im_and_fem(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs for tissue and vessel problems ..." << endl;
	#endif
	pintegration_method pim_t = int_method_descriptor(descr.IM_TYPET);
	pintegration_method pim_v = int_method_descriptor(descr.IM_TYPEV);
	mimt.set_integration_method(mesht.convex_index(), pim_t);
	mimv.set_integration_method(meshv.convex_index(), pim_v);
	
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEMs for tissue and vessel problems ..." << endl;
	#endif
	bgeot::pgeometric_trans pgt_t = bgeot::geometric_trans_descriptor(descr.MESH_TYPET);
	bgeot::pgeometric_trans pgt_v = bgeot::geometric_trans_descriptor(descr.MESH_TYPEV);
	pfem pf_Ut = fem_descriptor(descr.FEM_TYPET);
	pfem pf_Pt = fem_descriptor(descr.FEM_TYPET_P);
	pfem pf_Uv = fem_descriptor(descr.FEM_TYPEV);
	pfem pf_Pv = fem_descriptor(descr.FEM_TYPEV_P);
	pfem pf_coeft = fem_descriptor(descr.FEM_TYPET_DATA);
	pfem pf_coefv = fem_descriptor(descr.FEM_TYPEV_DATA);
	DIMT = pgt_t->dim();	//DIMV = 1;
	mf_Ut.set_qdim(bgeot::dim_type(DIMT)); 
	
	mf_Ut.set_finite_element(mesht.convex_index(), pf_Ut);
	GMM_ASSERT1(mf_Ut.get_qdim() == mf_Ut.fem_of_element(0)->target_dim(), 
		"Intrinsic vectorial FEM used"); // RT0 IS INTRINSIC VECTORIAL!!!
	mf_Pt.set_finite_element(mesht.convex_index(), pf_Pt);
	mf_coeft.set_finite_element(mesht.convex_index(), pf_coeft); 
	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for vessel branches ..." << endl;
	#endif
	mf_Uvi.reserve(nb_branches);
	mf_coefvi.reserve(nb_branches);
	for(size_type i=0; i<nb_branches; ++i){
		
		mesh_fem mf_tmp(meshv);
		mf_tmp.set_finite_element(meshv.region(i).index(), pf_coefv);
		mf_coefvi.emplace_back(mf_tmp);
		mf_tmp.clear();
		
		mf_tmp.set_finite_element(meshv.region(i).index(), pf_Uv);
		mf_Uvi.emplace_back(mf_tmp);
		mf_tmp.clear();
	}
	mf_Pv.set_finite_element(meshv.convex_index(), pf_Pv);
	mf_coefv.set_finite_element(meshv.convex_index(), pf_coefv);
	
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEM dimensions for tissue and vessel problems ..." << endl;
	#endif
	dof.set(mf_Ut, mf_Pt, mf_Uvi, mf_Pv, mf_coeft, mf_coefv);
	#ifdef M3D1D_VERBOSE_
	cout << std::scientific << dof;
	#endif
}

void
problem3d1d::build_param(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Building parameters for tissue and vessel problems ..." << endl;
	#endif
	param.build(PARAM, mf_coeft, mf_coefv);
	#ifdef M3D1D_VERBOSE_
	cout << param ;
	#endif
}

void
problem3d1d::build_tissue_boundary (void) 
{
	#ifdef M3D1D_VERBOSE_
	cout << "Building tissue boundary ..." << endl;
	#endif
	BCt.clear();
	BCt.reserve(2*DIMT);
	// Parse BC data
	string label_in = PARAM.string_value("BClabel", "Array of tissue boundary labels");
	string value_in = PARAM.string_value("BCvalue", "Array of tissue boundary values");
	vector<string> labels = split(label_in, ' ');
	vector<string> values = split(value_in, ' ');
	GMM_ASSERT1(labels.size()==2*DIMT, "wrong number of BC labels");
	GMM_ASSERT1(values.size()==2*DIMT, "wrong number of BC values");
	for (unsigned f=0; f<2*DIMT; ++f) {
		BCt.emplace_back(labels[f], std::stof(values[f]), 0, f);
		#ifdef M3D1D_VERBOSE_
		cout << "  face " << f << " : " << BCt.back() << endl;
		#endif
	}
	// Build mesht regions
	mesh_region border_faces;
	outer_faces_of_mesh(mesht, border_faces);

	for (mr_visitor i(border_faces); !i.finished(); ++i) {

		assert(i.is_face());

		// Unit outward normal : used to identify faces
		//! \todo Use getfem 5.0's function select_faces_of_normal?
		base_node un = mesht.normal_of_face_of_convex(i.cv(), i.f());
		un /= gmm::vect_norm2(un);

		if (gmm::abs(un[0] + 1.0) < 1.0E-7)      // back
			mesht.region(0).add(i.cv(), i.f());
		else if (gmm::abs(un[0] - 1.0) < 1.0E-7) // front
			mesht.region(1).add(i.cv(), i.f());
		else if (gmm::abs(un[1] + 1.0) < 1.0E-7) // left
			mesht.region(2).add(i.cv(), i.f());
		else if (gmm::abs(un[1] - 1.0) < 1.0E-7) // right
			mesht.region(3).add(i.cv(), i.f());
		else if (gmm::abs(un[2] + 1.0) < 1.0E-7) // bottom
			mesht.region(4).add(i.cv(), i.f());
		else if (gmm::abs(un[2] - 1.0) < 1.0E-7) // top
			mesht.region(5).add(i.cv(), i.f());
		
	} /* end of border_faces loop */
	// Export an indicator function for BC regions
	if (PARAM.int_value("VTK_EXPORT")){

		vector_type ones(dof.coeft(), 1.0);
		vector_type indicator(dof.coeft());
		for (unsigned f=0; f<2*DIMT; ++f) {
			asm_source_term(indicator, mimt, mf_coeft, mf_coeft, 
				gmm::scaled(ones, BCt[f].value), mesht.region(BCt[f].rg));
		}
		vtk_export rgvtk(descr.OUTPUT+"mesht_boundary.vtk");
		rgvtk.exporting(mf_coeft);
		rgvtk.write_mesh();
		rgvtk.write_point_data(mf_coeft, indicator, "1t");
	}
}

void 
problem3d1d::build_vessel_boundary(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Building vessel boundary ..." << endl;
	#endif
try {

	dal::bit_vector junctions; // global idx of junctions vertices in meshv
	dal::bit_vector extrema;   // global idx of extreme vertices in meshv

	Jv.clear();
	nb_extrema=0; 
	nb_junctions=0;
	
	size_type fer = nb_branches; // first empty region
	GMM_ASSERT1(meshv.has_region(fer)==0, 
		"Overload in meshv region assembling!");
	
	// List all the convexes
	dal::bit_vector nn = meshv.convex_index();
	bgeot::size_type cv;
	for (cv << nn; cv != bgeot::size_type(-1); cv << nn) {
		
		bgeot::pconvex_structure cvs = meshv.structure_of_convex(cv);
		if (cvs->nb_points()>2) 
			cerr << "Error: convex #" << cv << "has more than 2 vertices!" << endl;
		if (cvs->nb_faces()>2)  
			cerr << "Error: convex #" << cv << "has more than 2 faces!" << endl;

		// Build regions for BCs and junctions
		// Global idx of mesh vertices
		size_type i0 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(1)[0]];
		size_type i1 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(0)[0]];
		// Identify vertex type
		if (meshv.convex_to_point(i0).size()==1){ /* inflow extremum */
			// Update information
			extrema.add(i0);
			nb_extrema++;
			// Build a new region made by a single face
			GMM_ASSERT1(meshv.has_region(fer)==0, 
				"Overload in meshv region assembling!");
			meshv.region(fer).add(cv, 1);
			// Store the current index and then update it
			size_type bc = 0; 
			bool found = false;
			while (!found && (bc<BCv.size())) {
				found = (i0 == BCv[bc].idx);
				if (!found) bc++;
			}
			GMM_ASSERT1(found=true, "Miss a boundary node in BCv list!");
			BCv[bc].rg = fer; 
			fer++;
			// Store the containing branch index
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			BCv[bc].branches.emplace_back(branch); 
		}
		else if (meshv.convex_to_point(i0).size()==2){ /* trivial inflow junction */
			// DO NOTHING
		}
		else if (meshv.convex_to_point(i0).size()>=2){ /* non-trivial inflow junction */
			// Check if jucntion has been already stored, 
			// if not add to the junction list (J) and build a new region
			dal::bit_vector tmp; tmp.add(i0);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i0);
				nb_junctions++;
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshv.region(fer).add(cv, 1); // single-face region
				// Create a new junction node
				Jv.emplace_back("JUN", 0, i0, fer);
				fer++;
			}
			// Search for index of containing branch (\mathcal{P}^{in}_j)
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			// Add the inflow branch (to the right junction node)
			size_type jj = 0;
			bool found = false;
			while (!found && jj < nb_junctions){
				found = (i0 == Jv[jj].idx);
				if (!found) jj++;
			}
			//cout << "Branch -" << branch << " added to junction " << jj << endl;
			Jv[jj].value += param.R(mimv, branch);
			Jv[jj].branches.emplace_back(-branch);
			GMM_ASSERT1(branch>0, 
				"Error in network labeling: -0 makes no sense");
		}
		
		if (meshv.convex_to_point(i1).size()==1){ 
			size_type bc = 0; 
			bool found = false;
			while (!found && (bc<BCv.size())) {
				found = (i1 == BCv[bc].idx);
				if (!found) bc++;
			}
			if (found){ /* outlow extremum */
				extrema.add(i1); 
				nb_extrema++; 
				// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 0);
				// Store the current index and then update it
				BCv[bc].value *= -1.0;
				BCv[bc].rg = fer; 
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv[bc].branches.emplace_back(branch); 
			}
			else { /* interior -> Mixed point */
				// "MIX" label via post-processing
				// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 0);
				BCv.emplace_back("MIX", 0.0, i1, fer);
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv.back().branches.emplace_back(branch); 
			}
		}
		else if (meshv.convex_to_point(i1).size()==2){ /* trivial outflow junction */

			// Search for index of first containing branch (\mathcal{P}^{out}_j)
			size_type firstbranch = 0; 
			bool contained = false;
			while (!contained && firstbranch<nb_branches ) {
				contained = meshv.region(firstbranch).is_in(cv);
				if (!contained) firstbranch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i1!");

			// Check if i1 is a trivial junction (or a INT point)
			size_type cv1 = meshv.convex_to_point(i1)[0];
			size_type cv2 = meshv.convex_to_point(i1)[1];
			bool is_junc = (meshv.region(firstbranch).is_in(cv1) < 1 ||
							meshv.region(firstbranch).is_in(cv2) < 1 );
							
			if (is_junc){
				cout << "Found a trivial junction at i1 = " << i1 << endl;
				// Check if jucntion has been already stored, 
				// if not add to the junction list (J) and build a new region
				dal::bit_vector tmp; tmp.add(i1);
				if(!junctions.contains(tmp)){
					// Store the junction vertex
					junctions.add(i1);
					nb_junctions++;
					GMM_ASSERT1(meshv.has_region(fer)==0, 
						"Overload in meshv region assembling!");
					// Build a new region with idx "first empty region"
					meshv.region(fer).add(cv, 0);
					// Create a new junction node
					Jv.emplace_back("JUN", 0, i1, fer);
					fer++;
					//} QUI C'ERA LA CHIUSURA DI JUNCTION.CONTAINS(TMP)
					// Search for index of second containing branch (\mathcal{P}^{out}_j)
					size_type secondbranch = 0; 
					size_type secondcv = (( cv1 == cv) ? cv2 : cv1);
					size_type firstcv  = (( cv1 != cv) ? cv2 : cv1);
					contained = false;
					while (!contained && secondbranch<nb_branches ) {
						if(secondbranch != firstbranch)
							contained = meshv.region(secondbranch).is_in(secondcv);
						if (!contained) secondbranch++;
					}
					GMM_ASSERT1(contained=true, "No branch region contains node i1!");
					// Add the two branches
					scalar_type in=0;
					if(meshv.ind_points_of_convex(firstcv)[0]==i1) in=-1;
					else if (meshv.ind_points_of_convex(firstcv)[1]==i1) in=+1;
					GMM_ASSERT1(in!=0, "There is something wrong in firstbranch convex index");
					Jv.back().branches.emplace_back(in*firstbranch);

					in=0;
					if(meshv.ind_points_of_convex(secondcv)[0]==i1) in=-1;
					else if (meshv.ind_points_of_convex(secondcv)[1]==i1) in=+1;
					GMM_ASSERT1(in!=0, "There is something wrong in secondbranch convex index");
					Jv.back().branches.emplace_back(in*secondbranch);
					Jv.back().value += param.R(mimv, firstbranch);
					Jv.back().value += param.R(mimv, secondbranch);
				}
			}
		}
		else if (meshv.convex_to_point(i1).size()>=2){ /* non-trivial outflow junction */

			// Search for index of containing branch (\mathcal{P}^{out}_j)
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");

			// Check if jucntion has been already stored, 
			// if not add to the junction list (J) and build a new region
			dal::bit_vector tmp; tmp.add(i1);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i1);
				nb_junctions++;
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshv.region(fer).add(cv, 0);
				// Create a new junction node
				Jv.emplace_back("JUN", 0, i1, fer);
				// Add the outflow branch
				Jv.back().branches.emplace_back(+branch);
				Jv.back().value += param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << i1 << endl;
				fer++;
			}
			else {
				// Add the outflow branch (to the right junction node)
				size_type jj = 0;
				bool found = false;
				while (!found && jj < nb_junctions){
					found = (i1 == Jv[jj].idx);
					if (!found) jj++;
				}
				Jv[jj].branches.emplace_back(+branch);
				Jv[jj].value += param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << jj << endl;
			}
		}

	} /* end of convexes loop */
	
	// Ckeck network assembly
	#ifdef M3D1D_VERBOSE_
	cout << "--- NETWORK ASSEMBLY ------------------ "   << endl;
	cout << "  Branches:   " << nb_branches << endl
		 << "  Vertices:   " << nn.size()+1 << endl;
	cout << "  Extrema:    " << extrema << endl;	  
	for (size_type i=0; i<BCv.size(); ++i)
		cout << "    -  label=" << BCv[i].label 
			 << ", value=" << BCv[i].value << ", ind=" << BCv[i].idx 
			 << ", rg=" << BCv[i].rg << ", branches=" << BCv[i].branches << endl; 
	cout << "  Junctions: " << junctions << endl;
	for (size_type i=0; i<Jv.size(); ++i)
		cout << "    -  label=" << Jv[i].label 
			 << ", value=" << Jv[i].value << ", ind=" << Jv[i].idx 
			 << ", rg=" << Jv[i].rg << ", branches=" << Jv[i].branches << endl; 
	cout << "---------------------------------------- "   << endl;
	#endif

} 
GMM_STANDARD_CATCH_ERROR; // catches standard errors

} /* end of build_vessel_boundary */

//////// Assemble the problem ////////////////////////////////////////// 
void
problem3d1d::assembly(void)
{
	//1 Build the monolithic matrix AM
	assembly_mat();
	//2 Build the monolithic rhs FM
	assembly_rhs();
}


void 
problem3d1d::assembly_mat(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Allocating AM, UM, FM ..." << endl;
	#endif
	gmm::resize(AM, dof.tot(), dof.tot()); gmm::clear(AM);
	gmm::resize(UM, dof.tot()); gmm::clear(UM);
	gmm::resize(FM, dof.tot()); gmm::clear(FM);
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic matrix AM ..." << endl;
	#endif
	// Mass matrix for the interstitial problem
	sparse_matrix_type Mtt(dof.Ut(), dof.Ut());
	// Divergence matrix for the interstitial problem
	sparse_matrix_type Dtt(dof.Pt(), dof.Ut());
	// Junction compatibility matrix for the network problem
	sparse_matrix_type Jvv(dof.Pv(), dof.Uv());
	// Tissue-to-tissue exchange matrix
	sparse_matrix_type Btt(dof.Pt(), dof.Pt());
	// Vessel-to-tissue exchange matrix
	sparse_matrix_type Btv(dof.Pt(), dof.Pv());
	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Bvt(dof.Pv(), dof.Pt());
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Bvv(dof.Pv(), dof.Pv());
	// Aux tissue-to-vessel averaging matrix
	sparse_matrix_type Mbar(dof.Pv(), dof.Pt());
	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type Mlin(dof.Pv(), dof.Pt());
	
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Mtt and Dtt ..." << endl;
	#endif
	asm_tissue_darcy(Mtt, Dtt, mimt, mf_Ut, mf_Pt);
	gmm::scale(Mtt, 1.0/param.kt(0)); // kt scalar
	// Copy Mtt
	gmm::add(Mtt, 
			  gmm::sub_matrix(AM, 
					gmm::sub_interval(0, dof.Ut()), 
					gmm::sub_interval(0, dof.Ut()))); 
	// Copy -Dtt^T
	gmm::add(gmm::scaled(gmm::transposed(Dtt), -1.0),  
			  gmm::sub_matrix(AM, 
					gmm::sub_interval(0, dof.Ut()),
					gmm::sub_interval(dof.Ut(), dof.Pt()))); 
	// Copy Dtt
	gmm::add(Dtt,
			  gmm::sub_matrix(AM, 
					gmm::sub_interval(dof.Ut(), dof.Pt()),
					gmm::sub_interval(0, dof.Ut()))); 

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling the tangent versor ..." << endl;
	#endif
	vector_type lambdax; // tangent versor: x component
	vector_type lambday; // tangent versor: y component
	vector_type lambdaz; // tangent versor: z component
	std::ifstream ifs(descr.MESH_FILEV);
	GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr.MESH_FILEV);
	asm_tangent_versor(ifs, lambdax, lambday, lambdaz);
	ifs.close();

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Mvv and Dvv ..." << endl;
	#endif
	// Local matrices
	size_type shift = 0;
	for(size_type i=0; i<nb_branches; ++i){

		if(i>0) shift += mf_Uvi[i-1].nb_dof();
		scalar_type Ri = param.R(mimv, i);
		// Coefficient \pi^2*Ri'^4/\kappa_v
		vector_type ci(mf_coefvi[i].nb_dof(), pi*pi*Ri*Ri*Ri*Ri/param.kv(mimv,i));
		// Allocate temp local matrices
		sparse_matrix_type Mvvi(mf_Uvi[i].nb_dof(), mf_Uvi[i].nb_dof());
		sparse_matrix_type Dvvi(dof.Pv(), mf_Uvi[i].nb_dof());
		// Allocate temp local tangent versor
		vector_type lambdax_K, lambday_K, lambdaz_K;
		for(size_type j=0; j<mf_coefvi[i].nb_dof(); j++){
			lambdax_K.emplace_back(lambdax[i]);
			lambday_K.emplace_back(lambday[i]);
			lambdaz_K.emplace_back(lambdaz[i]);
		}
		// Build Mvvi and Dvvi
		asm_network_poiseuille(Mvvi, Dvvi, 
			mimv, mf_Uvi[i], mf_Pv, mf_coefvi[i],
			ci, lambdax_K, lambday_K, lambdaz_K, meshv.region(i));
		gmm::scale(Dvvi, pi*Ri*Ri);
		// Copy Mvvi and Dvvi
		gmm::add(Mvvi, 
			gmm::sub_matrix(AM, 
				gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof()), 
				gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof()))); 
		gmm::add(gmm::scaled(gmm::transposed(Dvvi), -1.0),
			gmm::sub_matrix(AM, 
				gmm::sub_interval(dof.Ut()+dof.Pt()+shift,    mf_Uvi[i].nb_dof()),
				gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()))); 
		gmm::add(Dvvi, 
			gmm::sub_matrix(AM, 
				gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()),
				gmm::sub_interval(dof.Ut()+dof.Pt()+shift,     mf_Uvi[i].nb_dof()))); 
		gmm::clear(Mvvi); 
		gmm::clear(Dvvi);
		
	} /* end of branches loop */
	
   if (nb_junctions > 0){
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Jvv" << " ..." << endl;
	#endif
	asm_network_junctions(Jvv, mimv, mf_Uvi, mf_Pv, mf_coefv, 
		Jv, param.R());
	#ifdef M3D1D_VERBOSE_
	cout << "  Copying -Jvv^T" << " ..." << endl;
	#endif		
	gmm::add(gmm::scaled(gmm::transposed(Jvv), -1.0),
		gmm::sub_matrix(AM,
			 gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv()), 
			 gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())));
	#ifdef M3D1D_VERBOSE_
	cout << "  Copying Jvv" << " ..." << endl;
	#endif		
	gmm::add(Jvv,
		gmm::sub_matrix(AM,
			 gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()),
			 gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));
  }
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling aux exchange matrices Mbar and Mlin ..." << endl;
	#endif
	asm_exchange_aux_mat(Mbar, Mlin, 
			mimv, mf_Pt, mf_Pv, param.R(), descr.NInt);
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling exchange matrices ..." << endl;
	#endif
	bool NEWFORM = PARAM.int_value("NEW_FORMULATION");
	asm_exchange_mat(Btt, Btv, Bvt, Bvv,
			mimv, mf_Pv, mf_coefv, Mbar, Mlin, param.Q(), NEWFORM);

	// Copying Btt
	gmm::add(Btt, 
			  gmm::sub_matrix(AM, 
					gmm::sub_interval(dof.Ut(), dof.Pt()), 
					gmm::sub_interval(dof.Ut(), dof.Pt()))); 
	// Copying -Btv
	gmm::add(gmm::scaled(Btv, -1),
	 		  gmm::sub_matrix(AM, 
					gmm::sub_interval(dof.Ut(), dof.Pt()),
					gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()))); 
	// Copying -Bvt
	gmm::add(gmm::scaled(Bvt,-1), 
			  gmm::sub_matrix(AM, 
			  		gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv()	, dof.Pv()),
					gmm::sub_interval(dof.Ut(), dof.Pt()))); 
	// Copying Bvv
	gmm::add(Bvv, 
			  gmm::sub_matrix(AM, 
					gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()), 
					gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()))); 

	// De-allocate memory
	gmm::clear(Mtt);  gmm::clear(Dtt); 
	gmm::clear(Mbar); gmm::clear(Mlin);
	gmm::clear(Btt);  gmm::clear(Btv);
	gmm::clear(Bvt);  gmm::clear(Bvv);

}

void 
problem3d1d::assembly_rhs(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic rhs FM ... " << endl;
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Initializing RHS for FM ..." << endl;
	#endif
	// Right Hand Side for the interstitial problem 
	vector_type Ft(dof.Ut());
	// Right Hand Side for the network problem 
	vector_type Fv(dof.Uv());

	// Coefficients for tissue BCs
	scalar_type bcoef  = PARAM.real_value("BETA", "Coefficient for mixed BC");
	scalar_type p0coef = PARAM.real_value("P0"); // default: 0

	#ifdef M3D1D_VERBOSE_
	cout << "  Building tissue boundary term ..." << endl;
	#endif
	vector_type beta(dof.coeft(), 1.0/bcoef);
	vector_type P0(dof.coeft(), p0coef);
	vector_type P0_vel(dof.coefv(),p0coef);
	
	if (PARAM.int_value("TEST_RHS")) {
		#ifdef M3D1D_VERBOSE_
		cout << "  ... as the divergence of exact velocity ... " << endl;
		#endif
		assembly_tissue_test_rhs();
	}
	else {
		sparse_matrix_type Mtt(dof.Ut(), dof.Ut());
		asm_tissue_bc(Mtt, Ft, mimt, mf_Ut, mf_coeft, BCt, P0, beta);
		gmm::add(Mtt, 
			gmm::sub_matrix(AM,
				gmm::sub_interval(0, dof.Ut()),
				gmm::sub_interval(0, dof.Ut())));
		gmm::add(Ft, gmm::sub_vector(FM, gmm::sub_interval(0, dof.Ut())));
		// De-allocate memory
		gmm::clear(Mtt); 
	}
	#ifdef M3D1D_VERBOSE_
	cout << "  Building vessel boundary term ..." << endl;
	#endif
	sparse_matrix_type Mvv(dof.Uv(), dof.Uv());
	asm_network_bc(Mvv, Fv, 
			mimv, mf_Uvi, mf_coefv, BCv, P0_vel, param.R(), bcoef);
	gmm::add(Mvv, 
		gmm::sub_matrix(AM,
			gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv()),
			gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));
	gmm::add(Fv, gmm::sub_vector(FM, gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));
	// De-allocate memory
	gmm::clear(Ft); gmm::clear(Fv); gmm::clear(Mvv);
	
}

void
problem3d1d::assembly_tissue_test_rhs(void)
{
	// Exact rhs (see defines.hpp)
	vector_type sol_Gt(dof.coeft());
	interpolation_function(mf_coeft, sol_Gt, sol_gt);
	#ifdef M3D1D_VERBOSE_
	cout << "    Assemble divergence source term ... " << endl;
	#endif
	vector_type Gt(dof.Pt());
	asm_source_term(Gt, mimt, mf_Pt, mf_coeft, sol_Gt);	
	gmm::add(Gt, gmm::sub_vector(FM, gmm::sub_interval(dof.Ut(), dof.Pt())));

  if (PARAM.int_value("VTK_EXPORT")) {

	#ifdef M3D1D_VERBOSE_
	cout << "    Compute theoretical expectations ... " << endl;
	#endif
	// FE spaces for exact velocity
	mesh_fem mf_data(mesht), mf_data_vec(mesht);
	mf_data_vec.set_qdim(DIMT);
	bgeot::pgeometric_trans pgt_t = bgeot::geometric_trans_descriptor(descr.MESH_TYPET);
	mf_data.set_classical_discontinuous_finite_element(1);
	mf_data_vec.set_classical_discontinuous_finite_element(1);

	GMM_ASSERT1(mf_data.nb_dof()*DIMT==mf_data_vec.nb_dof(), 
		"Wrong size of mf_data_vec"); 
	// Exact pressure (see defines.hpp)
	vector_type sol_Pt(dof.coeft());
	interpolation_function(mf_coeft, sol_Pt, sol_pt);
	// Exact x velocity (see defines.hpp)
	vector_type sol_Utx(mf_data.nb_dof());
	interpolation_function(mf_data, sol_Utx, sol_utx);
	// Exact y velocity (see defines.hpp)
	vector_type sol_Uty(mf_data.nb_dof());
	interpolation_function(mf_data, sol_Uty, sol_uty);
	// Exact z velocity (see defines.hpp)
	vector_type sol_Utz(mf_data.nb_dof());
	interpolation_function(mf_data, sol_Utz, sol_utz);
	// Exact velocity magnitude (see defines.hpp)
	vector_type sol_Utm(mf_data.nb_dof());
	interpolation_function(mf_data, sol_Utm, sol_utm);
	// Exact vectorial velocity (see defines.hpp)
	vector_type sol_Ut; sol_Ut.reserve(mf_data_vec.nb_dof());
	for( size_type i=0; i<mf_data_vec.nb_dof()/DIMT; ++i ){
		sol_Ut.emplace_back(sol_Utx[i]);
		sol_Ut.emplace_back(sol_Uty[i]);
		sol_Ut.emplace_back(sol_Utz[i]);
	}

	#ifdef M3D1D_VERBOSE_
	cout << "    Export theoretical expectations ... " << endl;
	#endif
	vtk_export vtk_sol_Gt(descr.OUTPUT+"sol_Gt.vtk");
	vtk_sol_Gt.exporting(mf_coeft);
	vtk_sol_Gt.write_mesh();
	vtk_sol_Gt.write_point_data(mf_coeft, sol_Gt, "sol_Gt");

	vtk_export vtk_sol_Pt(descr.OUTPUT+"sol_Pt.vtk");
	vtk_sol_Pt.exporting(mf_coeft);
	vtk_sol_Pt.write_mesh();
	vtk_sol_Pt.write_point_data(mf_coeft, sol_Pt, "sol_Pt");

	vtk_export vtk_sol_Utm(descr.OUTPUT+"sol_Utm.vtk");
	vtk_sol_Utm.exporting(mf_data);
	vtk_sol_Utm.write_mesh();
	vtk_sol_Utm.write_point_data(mf_data, sol_Utm, "sol_Utm");

	vtk_export vtk_sol_Ut(descr.OUTPUT+"sol_Ut.vtk");
	vtk_sol_Ut.exporting(mf_data_vec);
	vtk_sol_Ut.write_mesh();
	vtk_sol_Ut.write_point_data(mf_data_vec, sol_Ut, "sol_Ut");

  }

}

////////// Solve the problem ///////////////////////////////////////////    
bool
problem3d1d::solve(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Solving the monolithic system ... " << endl;
	#endif
	gmm::csc_matrix<scalar_type> A;
	gmm::clean(AM, 1E-12);
	gmm::copy(AM, A);
	gmm::clear(AM); // to be postponed for preconditioner
	double time = gmm::uclock_sec();	
	
	if ( descr.SOLVE_METHOD == "SuperLU" ) { // direct solver //
		#ifdef M3D1D_VERBOSE_
		cout << "  Applying the SuperLU method ... " << endl;
		#endif
		scalar_type cond;
		gmm::SuperLU_solve(A, UM, FM, cond);
		cout << "  Condition number : " << cond << endl;
	}
	else { // Iterative solver //

		// Iterations
		gmm::iteration iter(descr.RES);  // iteration object with the max residu
		iter.set_noisy(1);               // output of iterations (2: sub-iteration)
		iter.set_maxiter(descr.MAXITER); // maximum number of iterations

		// Preconditioners
		//! \todo Add preconditioner choice to param file
		// See \link http://download.gna.org/getfem/html/homepage/gmm/iter.html
		gmm::identity_matrix PM; // no precond
		//gmm::diagonal_precond<sparse_matrix_type> PM(AM); // diagonal preocond
		//gmm::ilu_precond<sparse_matrix_type> PM(AM);
		// ...
		//gmm::clear(AM);
		// See <http://download.gna.org/getfem/doc/gmmuser.pdf>, pag 15
	
		if ( descr.SOLVE_METHOD == "CG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Conjugate Gradient method ... " << endl;
			#endif
			gmm::identity_matrix PS;  // optional scalar product
			gmm::cg(AM, UM, FM, PS, PM, iter);
		}
		else if ( descr.SOLVE_METHOD == "BiCGstab" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the BiConjugate Gradient Stabilized method ... " << endl;
			#endif
			gmm::bicgstab(AM, UM, FM, PM, iter);
		}
		else if ( descr.SOLVE_METHOD == "GMRES" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Generalized Minimum Residual method ... " << endl;
			#endif
			size_type restart = 50;
			gmm::gmres(A, UM, FM, PM, restart, iter);
		}
		else if ( descr.SOLVE_METHOD == "QMR" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Quasi-Minimal Residual method ... " << endl;
			#endif
			gmm::qmr(AM, UM, FM, PM, iter);
		}
		else if ( descr.SOLVE_METHOD == "LSCG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the unpreconditionned Least Square CG method ... " << endl;
			#endif
			gmm::least_squares_cg(AM, UM, FM, iter);
		}
		// Check
		if (iter.converged())
			cout << "  ... converged in " << iter.get_iteration() << " iterations." << endl;
		else if (iter.get_iteration() == descr.MAXITER)
			cerr << "  ... reached the maximum number of iterations!" << endl;

	}
	cout << "... time to solve : " << gmm::uclock_sec() - time << " seconds\n";

	#ifdef M3D1D_VERBOSE_
	cout << "Compute the total flow rate ... " << endl;
	#endif
	// Aux vector
	vector_type Uphi(dof.Pv()); 
	// Extracting matrices Bvt, Bvv
	sparse_matrix_type Bvt(dof.Pv(), dof.Pt());
	sparse_matrix_type Bvv(dof.Pv(), dof.Pv());
	gmm::copy(gmm::sub_matrix(A, 
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv()	, dof.Pv()),
			gmm::sub_interval(dof.Ut(), dof.Pt())),
				Bvt); 
	gmm::copy(gmm::sub_matrix(A, 
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()), 
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())),
				Bvv); 
	// Extracting solutions Pt, Pv 
	vector_type Pt(dof.Pt()); 
	vector_type Pv(dof.Pv()); 
	gmm::copy(gmm::sub_vector(UM, 
		gmm::sub_interval(dof.Ut(), dof.Pt())), Pt);
	gmm::copy(gmm::sub_vector(UM, 
		gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())), Pv);
	// Computing Bvv*Pv - Bvt*Pt
	gmm::mult(Bvt, Pt, Uphi);
	gmm::mult_add(Bvv, Pv, Uphi);
	TFR = std::accumulate(Uphi.begin(), Uphi.end(), 0.0);

	// De-allocate memory
	gmm::clear(Bvt); gmm::clear(Bvv);
	gmm::clear(Pt);  gmm::clear(Pv);  
	gmm::clear(Uphi);

		size_type shift=dof.Ut()+dof.Pt();
		size_type Cb=0;
		cout<<"U problem3d1d:\n\n";
		for(size_type b=0;b<nb_branches;++b){
			cout<<" b"<<b<<" [";
			if(b>0) shift+= Cb;
			Cb=mf_Uvi[b].nb_dof();
			for(size_type i=0; i<Cb; ++i)
				cout<<" "<<UM[shift+i];
			cout<<" ]\n"; 
		}
		cout<<"\n";

	return true;
}


////////// Export results into vtk files ///////////////////////////////
void 
problem3d1d::export_vtk(const string & suff)
{
  if (PARAM.int_value("VTK_EXPORT"))
  {
	#ifdef M3D1D_VERBOSE_
	cout << "Exporting the solution (vtk format) to " << descr.OUTPUT << " ..." << endl;
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Saving the results from the monolithic unknown vector ... " << endl;
	#endif
	// Array of unknown dof of the interstitial velocity
	vector_type Ut(dof.Ut()); 
	// Array of unknown dof of the interstitial pressure
	vector_type Pt(dof.Pt()); 
	// Array of unknown dof of the network velocity
	vector_type Uv(dof.Uv()); 
	// Array of unknown dof of the network pressure
	vector_type Pv(dof.Pv()); 
	gmm::copy(gmm::sub_vector(UM, 
		gmm::sub_interval(0, dof.Ut())), Ut);
	gmm::copy(gmm::sub_vector(UM, 
		gmm::sub_interval(dof.Ut(), dof.Pt())), Pt);
	gmm::copy(gmm::sub_vector(UM, 
		gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())), Uv);
	gmm::copy(gmm::sub_vector(UM, 
		gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())), Pv);

	#ifdef M3D1D_VERBOSE_
	// Save vessel solution for test-cases
	if (nb_branches==1){
		std::ofstream outUv("Uv.txt");
		outUv << gmm::col_vector(Uv);
		outUv.close();
		std::ofstream outPv("Pv.txt");
		outPv << gmm::col_vector(Pv);
		outPv.close();
	}
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Ut ..." << endl;
	#endif
	pfem pf_Ut = fem_descriptor(descr.FEM_TYPET);
	if(pf_Ut->is_lagrange()==0){ 
		/*
			There is no built-in export for non-lagrangian FEM.
			If this is the case, we need to project before exporting.
		 */
		#ifdef M3D1D_VERBOSE_
		cout << "    Projecting Ut on P1 ..." << endl;
		#endif
		mesh_fem mf_P1(mesht);
		mf_P1.set_qdim(bgeot::dim_type(DIMT)); 
		mf_P1.set_classical_finite_element(1);
		sparse_matrix_type M_RT0_P1(mf_P1.nb_dof(), dof.Ut());
		sparse_matrix_type M_P1_P1(mf_P1.nb_dof(), mf_P1.nb_dof());
		vector_type Ut_P1(mf_P1.nb_dof());
		asm_mass_matrix(M_RT0_P1, mimt, mf_P1, mf_Ut);
		asm_mass_matrix(M_P1_P1,  mimt, mf_P1, mf_P1);
		
		vector_type Utt(mf_P1.nb_dof());
		gmm::mult(M_RT0_P1, Ut, Utt);
		double cond1;
		gmm::SuperLU_solve(M_P1_P1, Ut_P1, Utt, cond1);

		vtk_export exp1(descr.OUTPUT+"Ut.vtk");
		exp1.exporting(mf_P1);
		exp1.write_mesh();
		exp1.write_point_data(mf_P1, Ut_P1, "Ut");
	}	
	else {
		vtk_export exp_Ut(descr.OUTPUT+"Ut.vtk");
		exp_Ut.exporting(mf_Ut);
		exp_Ut.write_mesh();
		exp_Ut.write_point_data(mf_Ut, Ut, "Ut");	 
	}
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Pt ..." << endl;
	#endif
	vtk_export exp_Pt(descr.OUTPUT+"Pt.vtk");
	exp_Pt.exporting(mf_Pt);
	exp_Pt.write_mesh();
	exp_Pt.write_point_data(mf_Pt, Pt, "Pt");

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Uv ..." << endl;
	#endif
	size_type start = 0;
	size_type length = 0;
	for(size_type i=0; i<nb_branches; ++i){
		if(i>0) start += mf_Uvi[i-1].nb_dof();
		length = mf_Uvi[i].nb_dof();
		vtk_export exp_Uv(descr.OUTPUT+"Uv"+suff+std::to_string(i)+".vtk");
		exp_Uv.exporting(mf_Uvi[i]);
		exp_Uv.write_mesh();
		exp_Uv.write_point_data(mf_Uvi[i], 
			gmm::sub_vector(Uv, gmm::sub_interval(start, length)), "Uv"); 
	}

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Pv ..." << endl;
	#endif
	vtk_export exp_Pv(descr.OUTPUT+"Pv"+suff+".vtk");
	exp_Pv.exporting(mf_Pv);
	exp_Pv.write_mesh();
	exp_Pv.write_point_data(mf_Pv, Pv, "Pv");

	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
	#endif
  }

} /* end of export_vtk */

bool 
merge_and_solve(problem3d1d & Pba, problem3d1d & Pbv)
{

	//! \todo Improve code re-utilisation. Derived class?
	
	#ifdef M3D1D_VERBOSE_
	cout << "Merge arterial and venous monolithic matrix  ..." << endl;
	#endif
	// Dim assumption
	GMM_ASSERT1(Pba.dof.tot() == Pbv.dof.tot(),
		"arterial and venous problem must have same dimension");
	// Dimensions	
	size_type dof_t   = Pba.dof.Ut() + Pba.dof.Pt();
	size_type dof_v   = Pba.dof.Uv() + Pba.dof.Pv();
	size_type dof_ut  = Pba.dof.Ut();
	size_type dof_pt  = Pba.dof.Pt();
	size_type dof_uv  = Pba.dof.Uv();
	size_type dof_pv  = Pba.dof.Pv();
	size_type dof_tot = dof_t + 2*dof_v;
	// Temp arterial-venous monolithic matrix 
	sparse_matrix_type AMav(dof_tot, dof_tot);
	#ifdef M3D1D_VERBOSE_
	cout << "  Copying tissue and arterial matrices ..." << endl;
	#endif
	gmm::copy(Pba.AM, 
		gmm::sub_matrix(AMav,
			gmm::sub_interval(0, dof_t + dof_v),
			gmm::sub_interval(0, dof_t + dof_v)));
	#ifdef M3D1D_VERBOSE_
	cout << "  Copying venous matrices ..." << endl;
	#endif
	// Mvv, -Dvv^T, Dvv, Bvv
	gmm::copy(
		gmm::sub_matrix(Pbv.AM, 
			gmm::sub_interval(dof_t, dof_v), 
			gmm::sub_interval(dof_t, dof_v)),
		gmm::sub_matrix(AMav,
			gmm::sub_interval(dof_t+dof_v, dof_v),
			gmm::sub_interval(dof_t+dof_v, dof_v)));
	// Btt
	gmm::add(
		gmm::sub_matrix(Pbv.AM, 
			gmm::sub_interval(dof_ut, dof_pt), 
			gmm::sub_interval(dof_ut, dof_pt)),
		gmm::sub_matrix(AMav,
			gmm::sub_interval(dof_ut, dof_pt), 
			gmm::sub_interval(dof_ut, dof_pt)));
	// Btv
	gmm::copy(
		gmm::sub_matrix(Pbv.AM, 
			gmm::sub_interval(dof_ut, dof_pt), 
			gmm::sub_interval(dof_t+dof_uv, dof_pv)),
		gmm::sub_matrix(AMav,
			gmm::sub_interval(dof_ut, dof_pt), 
			gmm::sub_interval(dof_t+dof_v+dof_uv, dof_pv)));
	// Bvt
	gmm::copy(
		gmm::sub_matrix(Pbv.AM, 
			gmm::sub_interval(dof_t+dof_uv, dof_pv),
			gmm::sub_interval(dof_ut, dof_pt)), 
		gmm::sub_matrix(AMav,
			gmm::sub_interval(dof_t+dof_v+dof_uv, dof_pv),
			gmm::sub_interval(dof_ut, dof_pt))); 
	// De-allocate memory
	gmm::clear(Pba.AM); gmm::clear(Pbv.AM);

	#ifdef M3D1D_VERBOSE_
	cout << "Merge arterial and venous monolithic rhs  ..." << endl;
	#endif
	vector_type FMav(dof_tot);
	#ifdef M3D1D_VERBOSE_
	cout << "  Copying tissue and arterial RHS ..." << endl;
	#endif
	gmm::copy(Pba.FM,
		gmm::sub_vector(FMav, gmm::sub_interval(0, dof_t+dof_v)));	
	#ifdef M3D1D_VERBOSE_
	cout << "  Copying venous RHS ..." << endl;
	#endif
	gmm::copy(
		gmm::sub_vector(Pbv.FM,
			gmm::sub_interval(dof_t, dof_v)),
		gmm::sub_vector(FMav, 
			gmm::sub_interval(dof_t+dof_v, dof_v)));
	// De-allocate memory
	gmm::clear(Pba.FM); gmm::clear(Pbv.FM);

	#ifdef M3D1D_VERBOSE_
	cout << "Solving the extended arterial-venous monolithic rhs  ..." << endl;
	#endif
	vector_type UMav(dof_tot);
	double time = gmm::uclock_sec();	
	gmm::csc_matrix<scalar_type> Aav;
	
	if ( Pba.descr.SOLVE_METHOD == "SuperLU" ) { // direct solver //
		gmm::clean(AMav, 1E-12);
		gmm::copy(AMav, Aav);
		gmm::clear(AMav);
		#ifdef M3D1D_VERBOSE_
		cout << "  Applying the SuperLU method ... " << endl;
		#endif
		scalar_type cond;
		gmm::SuperLU_solve(Aav, UMav, FMav, cond);
		cout << "  Condition number : " << cond << endl;
	}
	else { // Iterative solver //

		// Iterations
		gmm::iteration iter(Pba.descr.RES);  // iteration object with the max residu
		iter.set_noisy(1);               // output of iterations (2: sub-iteration)
		iter.set_maxiter(Pba.descr.MAXITER); // maximum number of iterations

		// Preconditioners
		gmm::identity_matrix PMav; // no precond
	
		if ( Pba.descr.SOLVE_METHOD == "CG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Conjugate Gradient method ... " << endl;
			#endif
			gmm::identity_matrix PSav;  // optional scalar product
			gmm::cg(AMav, UMav, FMav, PSav, PMav, iter);
		}
		else if ( Pba.descr.SOLVE_METHOD == "BiCGstab" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the BiConjugate Gradient Stabilized method ... " << endl;
			#endif
			gmm::bicgstab(AMav, UMav, FMav, PMav, iter);
		}
		else if ( Pba.descr.SOLVE_METHOD == "GMRES" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Generalized Minimum Residual method ... " << endl;
			#endif
			size_type restart = 50;
			gmm::gmres(AMav, UMav, FMav, PMav, restart, iter);
		}
		else if ( Pba.descr.SOLVE_METHOD == "QMR" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Quasi-Minimal Residual method ... " << endl;
			#endif
			gmm::qmr(AMav, UMav, FMav, PMav, iter);
		}
		else if ( Pba.descr.SOLVE_METHOD == "LSCG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the unpreconditionned Least Square CG method ... " << endl;
			#endif
			gmm::least_squares_cg(AMav, UMav, FMav, iter);
		}
		// Check
		if (iter.converged())
			cout << "  ... converged in " << iter.get_iteration() << " iterations." << endl;
		else if (iter.get_iteration() == Pba.descr.MAXITER)
			cerr << "  ... reached the maximum number of iterations!" << endl;

	}
	cout << "... time to solve : " << gmm::uclock_sec() - time << " seconds\n";

	#ifdef M3D1D_VERBOSE_
	cout << "Saving results of venous problem ... " << endl;
	#endif
	gmm::copy(
		gmm::sub_vector(UMav, 
			gmm::sub_interval(0, dof_t+dof_v)), 
		gmm::sub_vector(Pba.UM, 
			gmm::sub_interval(0, dof_t+dof_v)));
	gmm::copy(
		gmm::sub_vector(UMav, 
			gmm::sub_interval(dof_t+dof_v, dof_v)), 
		gmm::sub_vector(Pbv.UM, 
			gmm::sub_interval(dof_t, dof_v)));

	#ifdef M3D1D_VERBOSE_
	cout << "Compute the total flow rate ... " << endl;
	#endif
	// Aux vector
	vector_type Uphi(dof_pv); 
	// Extracting matrices Bvt, Bvv
	sparse_matrix_type Bvt(dof_pv, dof_pt);
	sparse_matrix_type Bvv(dof_pv, dof_pv);
	gmm::copy(gmm::sub_matrix(Aav, 
			gmm::sub_interval(dof_t+dof_uv, dof_pv),
			gmm::sub_interval(dof_ut, dof_pt)),
				Bvt); 
	gmm::copy(gmm::sub_matrix(Aav, 
			gmm::sub_interval(dof_t+dof_uv, dof_pv), 
			gmm::sub_interval(dof_t+dof_uv, dof_pv)),
				Bvv); 
	// Extracting solutions Pt, Pv 
	vector_type Pt(dof_pt); 
	vector_type Pv(dof_pv); 
	gmm::copy(gmm::sub_vector(UMav, 
		gmm::sub_interval(dof_ut, dof_pt)), Pt);
	gmm::copy(gmm::sub_vector(UMav, 
		gmm::sub_interval(dof_t+dof_uv, dof_pv)), Pv);
	// Computing Bvv*Pv - Bvt*Pt
	gmm::mult(Bvt, Pt, Uphi);
	gmm::mult_add(Bvv, Pv, Uphi);
	Pba.TFR = std::accumulate(Uphi.begin(), Uphi.end(), 0.0);

	// Extracting matrices Bvt, Bvv
	gmm::copy(gmm::sub_matrix(Aav, 
			gmm::sub_interval(dof_t+dof_v+dof_uv, dof_pv),
			gmm::sub_interval(dof_ut, dof_pt)),
				Bvt); 
	gmm::copy(gmm::sub_matrix(Aav, 
			gmm::sub_interval(dof_t+dof_v+dof_uv, dof_pv), 
			gmm::sub_interval(dof_t+dof_v+dof_uv, dof_pv)),
				Bvv); 
	// Extracting solutions Pt, Pv 
	gmm::copy(gmm::sub_vector(UMav, 
		gmm::sub_interval(dof_t+dof_v+dof_uv, dof_pv)), Pv);
	// Computing Bvv*Pv - Bvt*Pt
	gmm::clear(Uphi);
	gmm::mult(Bvt, Pt, Uphi);
	gmm::mult_add(Bvv, Pv, Uphi);
	Pbv.TFR = std::accumulate(Uphi.begin(), Uphi.end(), 0.0);

	// De-allocate memory
	gmm::clear(AMav); gmm::clear(UMav); gmm::clear(FMav);
	gmm::clear(Bvt); gmm::clear(Bvv);
	gmm::clear(Pt);  gmm::clear(Pv);  
	gmm::clear(Uphi);

	return true;

} /* end of merge_and_solve */


} /* end of namespace */
