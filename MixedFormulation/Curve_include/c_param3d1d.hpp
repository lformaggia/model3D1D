#ifndef CURVED_M3D1D_PARAM3D1D_HPP_
#define CURVED_M3D1D_PARAM3D1D_HPP_

#include <mesh1d.hpp>    // import_network_radius
#include <utilities.hpp> // compute_radius
#include <param3d1d.hpp>
#include <c_mesh1d.hpp>

namespace getfem {

//! Class to handle the physical parameter of the coupled 3D/1D model
struct c_param3d1d : public param3d1d {
	//My Variable
	vector_type G_;
	scalar_type Gamma_;
	vector<vector_type> lambdax_;
	vector<vector_type> lambday_;
	vector<vector_type> lambdaz_;	
	vector<vector_type> Curv_;


	// Methods
	//! Build the arrays of dimensionless parameters
	void build(ftool::md_param & fname, 
			const getfem::mesh_fem & mf_datat,
			const getfem::mesh_fem & mf_datav,
			const vector<getfem::mesh_fem> & mf_datavi
			) 
	{
		FILE_ = fname;
		mf_datat_ = mf_datat;
		mf_datav_ = mf_datav;
		size_type dof_datat = mf_datat_.nb_dof();
		size_type dof_datav = mf_datav_.nb_dof();
		size_type n_branch= mf_datavi.size();
		 
		bool IMPORT_RADIUS = FILE_.int_value("IMPORT_RADIUS");
		bool NONDIM_PARAM  = FILE_.int_value("TEST_PARAM");
		bool EXPORT_PARAM  = FILE_.int_value("EXPORT_PARAM");
		bool IMPORT_CURVE = FILE_.int_value("IMPORT_CURVE");

		
		// Check
		if (IMPORT_RADIUS)
			GMM_ASSERT1(NONDIM_PARAM == 0,
				"try to import non constant (dimensionless) radius:" 
				"please insert dimensional parameters");
		#ifdef M3D1D_VERBOSE_
		cout << "  Assembling dimensionless radius R'... "   << endl;
		#endif
		if (!IMPORT_RADIUS) { 	/* case R' = const */
			if (NONDIM_PARAM) // we assume it is already non-dimensional
				Rav_ = FILE_.real_value("RADIUS", "Vessel average radius");
			else // to be non-dimensionalized
				Rav_ = FILE_.real_value("RADIUS", "Vessel average radius")/FILE_.real_value("d");
			R_.assign(dof_datav, Rav_);
		} else { 				/* case R' = R'(s) */
			std::string RFILE = FILE_.string_value("RFILE"); 
			cout << "  Importing radius values from file " << RFILE << " ..." << endl;
			std::ifstream ist(RFILE);
			if (!ist) cerr << "impossible to read from file " << RFILE << endl;
			import_network_radius(R_, ist, mf_datav_);
		}

		if(!IMPORT_CURVE){
			for(size_type b=0;b<n_branch;++b){
				Curv_[b].resize(mf_datavi[b].nb_dof()); Curv_.clear();
				Curv_[b].assign(mf_datavi[b].nb_dof(), 0.0);
				
				scalar_type lx_temp=lambdax_[b][0];
				scalar_type ly_temp=lambday_[b][0];
				scalar_type lz_temp=lambdaz_[b][0];
					
				lambdax_[b].resize(mf_datavi[b].nb_dof()); lambdax_[b].clear();
				lambday_[b].resize(mf_datavi[b].nb_dof()); lambday_[b].clear();
				lambdaz_[b].resize(mf_datavi[b].nb_dof()); lambdaz_[b].clear();

				lambdax_[b].assign(mf_datavi[b].nb_dof(),lx_temp);
				lambday_[b].assign(mf_datavi[b].nb_dof(),ly_temp);
				lambdaz_[b].assign(mf_datavi[b].nb_dof(),lz_temp);
			}
		} else {
			rasm_curve_parameter(mf_datavi,Curv_,lambdax_,lambday_,lambdaz_);
			for(size_type b=0;b<n_branch;++b)
				gmm::scaled(Curv_[b],1.0/FILE_.real_value("d"));
		}

		#ifdef M3D1D_VERBOSE_
		cout << "  Assembling dimensionless permeabilities kt, Q, kv ... "   << endl;
		#endif
		if (NONDIM_PARAM) {
			// Import dimensionless params from FILE_
			scalar_type ktval = FILE_.real_value("Kt"); 
			scalar_type Qval  = FILE_.real_value("Q"); 
			scalar_type kvval = FILE_.real_value("Kv");
			scalar_type Gval  = FILE_.real_value("G");
			Gamma_ = FILE_.real_value("Gamma");
			// Fill the data arrays
			kt_.assign(dof_datat, ktval);
			kv_.assign(dof_datav, kvval);
			Q_.assign(dof_datav,  Qval);
			G_.assign(dof_datav,Gval);
		} 
		else {
			// Import dimensional params from FILE_
			P_  = FILE_.real_value("P", "average interstitial pressure [Pa]"); 
			U_  = FILE_.real_value("U", "characteristic flow speed in the capillary bed [m/s]"); 
			d_  = FILE_.real_value("d", "characteristic length of the problem [m]"); 
			k_  = FILE_.real_value("k", "permeability of the interstitium [m^2]"); 
			mu_ = FILE_.real_value("mu", "fluid viscosity [kg/ms]"); 
			Lp_ = FILE_.real_value("Lp", "permeability of the vessel walls [m^2 s/kg]"); 
			Gamma_= FILE_.real_value("Gamma");
			// Compute the dimenless params
			kt_.assign(dof_datat, k_/mu_*P_/U_/d_);
			G_.assign(dof_datav, U_*U_/P_);
			for (auto r : R_){ // C++11-only!
				kv_.emplace_back(pi/2.0/(Gamma_+2.0)/mu_*P_*d_/U_*r*r*r*r);
				Q_.emplace_back(2.0*pi*Lp_*P_/U_*r);
			}
		}
		// Check values
		GMM_ASSERT1(kt_[0] != 0, "wrong tissue conductivity (kt>0 required)"); 
		GMM_ASSERT1(kv_[0] != 0, "wrong vessel bed conductivity (kv>0 required)");
		if (Q_[0] == 0) cout << "Warning: uncoupled problem (Q=0)" << endl;
		if (G_[0] == 0) cout << "Warning: creating linear problem (G=0)" << endl;
		
		if (EXPORT_PARAM){
			std::string ODIR = FILE_.string_value("OutputDir","OutputDirectory");
			getfem::vtk_export exp(ODIR+"radius.vtk");
			exp.exporting(mf_datav_);
			exp.write_mesh();
			exp.write_point_data(mf_datav_, R_, "R");
			getfem::vtk_export expQ(ODIR+"conductivity.vtk");
			expQ.exporting(mf_datav_);
			expQ.write_mesh();
			expQ.write_point_data(mf_datav_, Q_, "Q");
		}

	}
    /////////////////////////////////////////////////////////  MIO USATO IN MESH1D
	void get_curve(
		vector<vector_type> & Curv, 
		vector<vector_type> & lambdax, 
		vector<vector_type> & lambday, 
		vector<vector_type> & lambdaz
	)
	{
		Curv_=Curv;
		lambdax_=lambdax;
		lambday_=lambday;
		lambdaz_=lambdaz;
	}


	inline scalar_type G  (size_type i) { return G_[i];  } const
	inline scalar_type Gamma (void) { return Gamma_;} const
	vector_type G(void){ return G_;}

	vector<vector_type> & lambdax (void) { return lambdax_; }
	vector<vector_type> & lambday (void) { return lambday_; }
	vector<vector_type> & lambdaz (void) { return lambdaz_; }
	vector<vector_type> & Curv (void) { return Curv_; }
 
}; /* end of class */

} /* end of namespace */

#endif