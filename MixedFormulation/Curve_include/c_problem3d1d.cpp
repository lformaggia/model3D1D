
#include <c_problem3d1d.hpp>

namespace getfem {

	//! Parameteres for exact solution (1_uncoupled)
	/*! 
		\todo Read from param Lx, Ly, Lz and Kt
	 */
	double Lx1 = 1.0, Ly1 = 1.0, Lz1 = 1.0, kappat1 = 1.0; 
	//! Exact pressure 
	double sol_pt1(const bgeot::base_node & x){
		return sin(2.0*pi/Lx1*x[0])*sin(2.0*pi/Ly1*x[1])*sin(2*pi/Lz1*x[2]);
	}
	//! Exact x velocity
	double sol_utx1(const bgeot::base_node & x){
		return -2.0*pi*kappat1/Lx1*cos(2.0*pi/Lx1*x[0])*sin(2.0*pi/Ly1*x[1])*sin(2.0*pi/Lz1*x[2]);
	}
	//! Exact y velocity
	double sol_uty1(const bgeot::base_node & x){
		return -2.0*pi*kappat1/Ly1*sin(2.0*pi/Lx1*x[0])*cos(2.0*pi/Ly1*x[1])*sin(2.0*pi/Lz1*x[2]);
	}
	//! Exact z velocity
	double sol_utz1(const bgeot::base_node & x){
		return -2.0*pi*kappat1/Lz1*sin(2.0*pi/Lx1*x[0])*sin(2.0*pi/Ly1*x[1])*cos(2.0*pi/Lz1*x[2]);
	}
	//! Exact velocity magnitude
	double sol_utm1(const bgeot::base_node & x){
		return sqrt(sol_utx1(x)*sol_utx1(x)+sol_uty1(x)*sol_uty1(x)+sol_utz1(x)*sol_utz1(x));
	}
	//! Exact vectorial velocity
	std::vector<double> sol_ut1(const bgeot::base_node & x){
		std::vector<double> out(3);
		out[0] = sol_utx1(x); out[1] = sol_uty1(x); out[2] = sol_utz1(x);
		return out;
	}
	//! Exact rhs
	double sol_gt1(const bgeot::base_node & x){
		return 4.0*pi*pi*kappat1*(1.0/(Lx1*Lx1)+1.0/(Ly1*Ly1)+1.0/(Lz1*Lz1))*sin(2.0*pi/Lx1*x[0])*sin(2.0*pi/Ly1*x[1])*sin(2.0*pi/Lz1*x[2]);
	}



/////////// Initialize the problem ///////////////////////////////////// 
void 
c_problem3d1d::init(int argc, char *argv[])
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
c_problem3d1d::import_data(void)
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
c_problem3d1d::build_mesh(void)
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
	bool Import=PARAM.int_value("IMPORT_CURVE");

	if(Import){
		std::ifstream ifc(PARAM.string_value("CURVATURE_FILE","curvature file"));
		std::ifstream ift(PARAM.string_value("VELOCITY_VERSOR_FILE","velocity versor file"));
		GMM_ASSERT1(ifc.good(), "impossible to read from file " << PARAM.string_value("CURVATURE_FILE","curvature file"));
		GMM_ASSERT1(ift.good(), "impossible to read from file " << PARAM.string_value("VELOCITY_VERSOR_FILE","velocity versor file"));

		import_pts_file(ifs,ifc,ift, meshv, BCv, nb_vertices, descr.MESH_TYPEV, c_param);

		ifc.close();
		ift.close();
	} else{
		import_pts_file(ifs, meshv, BCv, nb_vertices, descr.MESH_TYPEV);
	}

	nb_branches = nb_vertices.size();
	ifs.close();
}

void
c_problem3d1d::set_im_and_fem(void)
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
	mf_Pt.set_finite_element(mesht.convex_index(),  pf_Pt);
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
c_problem3d1d::build_param(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Building parameters for tissue and vessel problems ..." << endl;
	#endif
	c_param.build(PARAM, mf_coeft, mf_coefv,mf_coefvi);
	#ifdef M3D1D_VERBOSE_
	cout << c_param ;
	#endif
}

void
c_problem3d1d::build_tissue_boundary (void) 
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
c_problem3d1d::build_vessel_boundary(void)
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
			Jv[jj].value += c_param.R(mimv, branch);
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
				}
				// Search for index of second containing branch (\mathcal{P}^{out}_j)
				size_type secondbranch = firstbranch+1; 
				size_type secondcv = (( cv1 == cv) ? cv2 : cv1);
				contained = false;
				while (!contained && secondbranch<nb_branches ) {
					contained = meshv.region(secondbranch).is_in(secondcv);
					if (!contained) secondbranch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				// Add the two branches
				Jv.back().branches.emplace_back(+firstbranch);
				Jv.back().branches.emplace_back(-secondbranch);
				Jv.back().value += c_param.R(mimv, firstbranch);
				Jv.back().value += c_param.R(mimv, secondbranch);
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
				Jv.back().value += c_param.R(mimv, branch);
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
				Jv[jj].value += c_param.R(mimv, branch);
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
c_problem3d1d::assembly(void)
{
	//1 Build the monolithic matrix AM
	assembly_mat();
	//2 Build the monolithic rhs FM
	assembly_rhs();
}

void 
c_problem3d1d::assembly_mat(void)
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
	gmm::scale(Mtt, 1.0/c_param.kt(0)); // kt scalar
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
 
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Mvv and Dvv ..." << endl;
	#endif
	// Local matrices
	size_type shift = 0;
	for(size_type i=0; i<nb_branches; ++i){

		if(i>0) shift += mf_Uvi[i-1].nb_dof();
		scalar_type Ri = c_param.R(mimv, i);
		// Coefficient \pi^2*Ri'^4/\kappa_v
		vector_type ci(mf_coefvi[i].nb_dof());
		for(size_type j=0; j<mf_coefvi[i].nb_dof(); ++j){
			ci[j]=pi*pi*Ri*Ri*Ri*Ri/c_param.kv(i)*(1+(c_param.Curv())[i][j]*(c_param.Curv())[i][j]*Ri*Ri);
		}
		// Allocate temp local matrices
		sparse_matrix_type Mvvi(mf_Uvi[i].nb_dof(), mf_Uvi[i].nb_dof());
		sparse_matrix_type Dvvi(dof.Pv(), mf_Uvi[i].nb_dof());
 

		// Build Mvvi and Dvvi
		asm_network_poiseuille(Mvvi, Dvvi, 
			mimv, mf_Uvi[i], mf_Pv, mf_coefvi[i],
			ci, (c_param.lambdax())[i], (c_param.lambday())[i], (c_param.lambdaz())[i], meshv.region(i));

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
			Jv, c_param.R());
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
			mimv, mf_Pt, mf_Pv, c_param.R(), descr.NInt);
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling exchange matrices ..." << endl;
	#endif
	bool NEWFORM = PARAM.int_value("NEW_FORMULATION");
	asm_exchange_mat(Btt, Btv, Bvt, Bvv,
			mimv, mf_Pv, mf_coefv, Mbar, Mlin, c_param.Q(), NEWFORM);

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
c_problem3d1d::assembly_rhs(void)
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
			mimv, mf_Uvi, mf_coefv, BCv, P0_vel, c_param.R(), bcoef);
	gmm::add(Mvv, 
		gmm::sub_matrix(AM,
			gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv()),
			gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));
	gmm::add(Fv, gmm::sub_vector(FM, gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));
	// De-allocate memory
	gmm::clear(Ft); gmm::clear(Fv); gmm::clear(Mvv);	
}

void
c_problem3d1d::assembly_tissue_test_rhs(void)
{
	// Exact rhs (see defines.hpp)
	vector_type sol_Gt(dof.coeft());
	interpolation_function(mf_coeft, sol_Gt, sol_gt1);
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
	vector_type sol_Pt1(dof.coeft());
	interpolation_function(mf_coeft, sol_Pt1, sol_pt1);
	// Exact x velocity (see defines.hpp)
	vector_type sol_Utx1(mf_data.nb_dof());
	interpolation_function(mf_data, sol_Utx1, sol_utx1);
	// Exact y velocity (see defines.hpp)
	vector_type sol_Uty1(mf_data.nb_dof());
	interpolation_function(mf_data, sol_Uty1, sol_uty1);
	// Exact z velocity (see defines.hpp)
	vector_type sol_Utz1(mf_data.nb_dof());
	interpolation_function(mf_data, sol_Utz1, sol_utz1);
	// Exact velocity magnitude (see defines.hpp)
	vector_type sol_Utm1(mf_data.nb_dof());
	interpolation_function(mf_data, sol_Utm1, sol_utm1);
	// Exact vectorial velocity (see defines.hpp)
	vector_type sol_Ut1; sol_Ut1.reserve(mf_data_vec.nb_dof());
	for( size_type i=0; i<mf_data_vec.nb_dof()/DIMT; ++i ){
		sol_Ut1.emplace_back(sol_Utx1[i]);
		sol_Ut1.emplace_back(sol_Uty1[i]);
		sol_Ut1.emplace_back(sol_Utz1[i]);
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
	vtk_sol_Pt.write_point_data(mf_coeft, sol_Pt1, "sol_Pt");

	vtk_export vtk_sol_Utm(descr.OUTPUT+"sol_Utm.vtk");
	vtk_sol_Utm.exporting(mf_data);
	vtk_sol_Utm.write_mesh();
	vtk_sol_Utm.write_point_data(mf_data, sol_Utm1, "sol_Utm");

	vtk_export vtk_sol_Ut(descr.OUTPUT+"sol_Ut.vtk");
	vtk_sol_Ut.exporting(mf_data_vec);
	vtk_sol_Ut.write_mesh();
	vtk_sol_Ut.write_point_data(mf_data_vec, sol_Ut1, "sol_Ut");

  }
}

void
c_problem3d1d::assembly_nonlinear_mat(void)
{
	gmm::clear(NL);
	scalar_type alfa  = (c_param.Gamma()+2.0)/(c_param.Gamma()+1.0);
	scalar_type beta  = (c_param.Gamma()+2.0)/(c_param.Gamma()+4.0);
	scalar_type gamma = (c_param.Gamma()+2.0)*(c_param.Gamma()+2.0)/3.0/(c_param.Gamma()+3.0)/(c_param.Gamma()+6.0);
	size_type shift=0;
	size_type dofi=mf_Uvi[0].nb_dof();
	bool newton=PARAM.int_value("Newton");

	for(size_type i=0; i<nb_branches; ++i){ 
		if(i>0) shift += dofi;
		scalar_type Ri = c_param.R(mimv, i);
		scalar_type Gi= c_param.G(0);
		dofi=mf_Uvi[i].nb_dof();
		scalar_type Method=1.0;
		if(newton)
			Method=2.0;

		vector_type ci(mf_coefvi[i].nb_dof());
		for(size_type j=0; j<mf_coefvi[i].nb_dof(); ++j){
			ci[j]=Method*pi*Ri*Ri*Gi*
				(	alfa  +
					beta  * (c_param.Curv())[i][j] * (c_param.Curv())[i][j] * Ri * Ri +
					gamma * (c_param.Curv())[i][j] * (c_param.Curv())[i][j] * (c_param.Curv())[i][j] * (c_param.Curv())[i][j] * Ri * Ri * Ri * Ri );
		}

		// Allocate temp local matrices
		sparse_matrix_type NMvvi(dofi,dofi);
		vector_type Uv_oldi(dofi);
		gmm::copy(gmm::sub_vector(Uv_old,gmm::sub_interval(shift,dofi)),Uv_oldi);
 
		// Build NLMvvi 
		asm_network_nonlinear(NMvvi, 
			mimv, mf_Uvi[i], mf_coefvi[i],
			ci, (c_param.lambdax())[i], (c_param.lambday())[i], (c_param.lambdaz())[i], Uv_oldi, meshv.region(i));

		//l'ho inserito trasposto
		gmm::add(gmm::scaled(gmm::transposed(NMvvi), -1.0), 
			gmm::sub_matrix(NLMvv,  
				gmm::sub_interval(shift, dofi), 
				gmm::sub_interval(shift, dofi))); 
	
		gmm::clear(NMvvi); 
	} /* end of branches loop */

	asm_network_nonlinear_junctions(NLJvv,  
		mimv, mf_Uvi, mf_coefv,
		c_param.Gamma(),c_param.R(), c_param.G(), (c_param.Curv()),Uv_old,newton); 

	gmm::add(NLJvv, NL); 
	gmm::add(NLMvv, NL); 

	gmm::clear(NLJvv);
	gmm::clear(NLMvv);
}

////////// Solve the problem ///////////////////////////////////////////
bool
c_problem3d1d::solve_pass(size_type iter)
{
	#ifdef M3D1D_VERBOSE_
	cout << "  Solving the monolithic system ... " << endl;
	#endif
	
		gmm::copy(AM,CM);
		gmm::copy(FM,CFM);
		if(iter>0){
			gmm::add( NL, 
				gmm::sub_matrix(CM, 
					gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv()), 
					gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv()))); 

			if(PARAM.int_value("Newton")){
				gmm::mult(NL,Uv_old,NLF); //NLF va diviso per 2
				gmm::add(gmm::scaled(NLF,1.0/2.0),
					gmm::sub_vector(CFM,
						gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));

			} 
		}

		gmm::csc_matrix<scalar_type> A;
		gmm::clean(CM, 1E-12);
		gmm::copy(CM, A);
		gmm::clear(CM);


	//double time = gmm::uclock_sec();	
	
	if ( descr.SOLVE_METHOD == "SuperLU" ) { // direct solver //
		#ifdef M3D1D_VERBOSE_
		//cout << "  Applying the SuperLU method ... " << endl;
		#endif
		scalar_type cond;
		gmm::SuperLU_solve(A, UM, CFM, cond);
		//cout << "  Condition number : " << cond << endl;
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
	//cout << "... time to solve : " << gmm::uclock_sec() - time << " seconds\n";

	#ifdef M3D1D_VERBOSE_
	//cout << "Compute the total flow rate ... " << endl;
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
	
	return true;
}

bool  
c_problem3d1d::solve(void)
{
	cout<<" Solving the system with a Fixed Point Method..."<<endl;
 
	gmm::resize(NLMvv, dof.Uv(),dof.Uv()); gmm::clear(NLMvv);
	gmm::resize(NLJvv, dof.Uv(),dof.Uv()); gmm::clear(NLJvv);
	gmm::resize(Uv_old,dof.Uv());          gmm::clear(Uv_old);
	gmm::resize(NL, dof.Uv(),dof.Uv());    gmm::clear(NL);
	gmm::resize(CM, dof.tot(),dof.tot());  gmm::clear(CM);
	gmm::resize(CFM, dof.tot());           gmm::clear(CFM);
	gmm::resize(NLF, dof.Uv());            gmm::clear(NLF);

	size_type N_iter = 0;
	ii=N_iter;

	size_type Max_iter = PARAM.int_value("Max_iter");
	scalar_type minERR = PARAM.real_value("minERR");
	bool Plot= PARAM.int_value("Plot");
	bool ERRH1= PARAM.int_value("ERRH1");
	
	scalar_type Error  = minERR+1.0;
	
	double time = gmm::uclock_sec();	
	size_type shift=0;
	size_type dofi=mf_Uvi[0].nb_dof();

	#ifdef M3D1D_VERBOSE_
	cout<<"\n MaxIter "<<Max_iter<<"\n";
	cout<<"Solving iter=0 (Stokes Problem)\n";
	#endif
	
	if(!(solve_pass(N_iter))) 
		GMM_ASSERT1(0," At the step "<<N_iter<<" the solver stopped");
	cout<<endl;
	N_iter++;
	if(Plot) Uplot();
	gmm::copy(gmm::sub_vector(UM,
					gmm::sub_interval( dof.Ut()+dof.Pt() , dof.Uv() ) ), 
				Uv_old);

	while( Error>minERR  && N_iter<Max_iter ) {   
		#ifdef M3D1D_VERBOSE_
		cout<<"Solving iter "<<N_iter<<"  with minIncremental Error "<<minERR <<endl;
		cout<<"  Assembling nonlinear term at iter "<<N_iter<<endl;
		#endif

		assembly_nonlinear_mat();

		//Risoluzione al passo N_iter
		if(!(solve_pass(N_iter))) 
			GMM_ASSERT1(0," At the step "<<N_iter<<" the solver stopped");
		
		//Calcolo Incremento
		Error=0.0;
		shift=0;
		dofi=mf_Uvi[0].nb_dof();
		for(size_type b=0; b<nb_branches; b++){
			if(b>0) shift+=dofi;
			dofi=mf_Uvi[b].nb_dof();
			scalar_type ErrorBranch;
			if(ERRH1==0)
				ErrorBranch=asm_L2_dist(mimv, 
							mf_Uvi[b] , gmm::sub_vector(UM    ,gmm::sub_interval(dof.Ut()+dof.Pt()+shift, dofi)), 
							mf_Uvi[b] , gmm::sub_vector(Uv_old,gmm::sub_interval(shift, dofi)), 
							meshv.region(b) );
			else if(ERRH1==1)
				ErrorBranch=asm_H1_dist(mimv, 
							mf_Uvi[b] , gmm::sub_vector(UM    ,gmm::sub_interval(dof.Ut()+dof.Pt()+shift, dofi)), 
							mf_Uvi[b] , gmm::sub_vector(Uv_old,gmm::sub_interval(shift, dofi)), 
							meshv.region(b) );
			Error=Error+ErrorBranch;
		}
		#ifdef M3D1D_VERBOSE_
		cout<<"Increment at iter "<<N_iter<<": "<<Error<<endl<<endl;
		#endif
		Errore.push_back(Error);
		if(Plot)	Uplot();

		//Aggiornamento
		N_iter++;
		gmm::clear(Uv_old);
		gmm::copy(gmm::sub_vector(UM,gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())), Uv_old);
	}

	double time_end=gmm::uclock_sec() - time;

	if((Error>minERR )&&(Max_iter>1 ))
	{
		#ifdef M3D1D_VERBOSE_
		cout<<"*---------------------------------------------------------------------------*\n";
		cout<<" The method doesn't converged for minumum Increment Error="<<minERR<<endl;
		cout<<"    Increment Error="<<Error<<endl;
		cout<<"*---------------------------------------------------------------------------*\n\n";
		#endif
		return false;
	} else{
		#ifdef M3D1D_VERBOSE_
		cout<<"*---------------------------------------------------------------------------*\n";
		cout<<"* The ";
		if(PARAM.int_value("Newton"))
			cout<<"Newton ";
		else
			cout<<"Fpoint ";
		cout<<"Method converged in "<<N_iter<<" iterations with ";
		if(!ERRH1)
			cout<<"L2:";
		else
			cout<<"H1:";

		cout<<" error "<<Error<<"   *"<<endl;
		cout << "* Time to solve the Problem : " << time_end << " seconds                          *"<<endl;
		cout<<"*---------------------------------------------------------------------------*\n\n";
		#endif
		return true;
	}


}

void 
c_problem3d1d::export_vtk(const string & suff)
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
} // end of vtk

} /* end of namespace */
