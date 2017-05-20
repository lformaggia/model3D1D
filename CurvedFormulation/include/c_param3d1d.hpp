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
	@brief  Definition of the aux class for physical parameters.
	@details 
	Assemble the dimensionless parameters of the coupled 3D/1D model:
	- Radius @f$R'(s)@f$,
	- Tissue permeability @f$\kappa_t@f$,
	- Vessel wall permeability @f$Q(s)@f$,
	- Vessel bed permeability @f$\kappa_v(s)@f$.
	- Vessel dynamic pressure gain @f$G(s)@f$.
	- Curvature of the vessel @f$\kappa'(s)@f$.
	- Tangent versor of the vessel @f$\Lambda'(s)@f$.


	being @f$s\in\Lambda@f$ the arc-lenght over the vessel network.
	\note @f$\kappa_t@f$ is assumed to be constant.

	\ingroup input
 */

#ifndef CURVED_M3D1D_PARAM3D1D_HPP_
#define CURVED_M3D1D_PARAM3D1D_HPP_

#include <mesh1d.hpp>    // import_network_radius
#include <utilities.hpp> // compute_radius
#include <param3d1d.hpp>
#include <c_mesh1d.hpp> //rasm_curve_parameter 
#include <c_assembling1d.hpp> //simple_compute_radius

namespace getfem {

//! Class to handle the physical parameter of the coupled 3D/1D model
struct c_param3d1d : public param3d1d {
	// Dynamic Pressure Gain
	vector_type G_;
	// Order of Velocity Profile
	scalar_type Gamma_;
	// Density of the fluid
	scalar_type rho_;
	// Trasport Coefficents for curved problem
	scalar_type alfa;
	scalar_type beta;
	scalar_type gamma;
	// Mesh tangent versor
	vector<vector_type> lambdax_;
	vector<vector_type> lambday_;
	vector<vector_type> lambdaz_;	
	// Mesh curvature
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
			gmm::scaled(R_,1.0/FILE_.real_value("d"));
		}

		if(!IMPORT_CURVE){
			#ifdef M3D1D_VERBOSE_
			cout<<"CURVE NOT IMPORTED, THE PROBLEM IS CONSIDERED LINEAR FOR ALL BRANCHES\n\n";
			#endif

			vector_type lx_temp,ly_temp,lz_temp;
			std::ifstream ifst(FILE_.string_value("MESH_FILEV","1D points file"));
			GMM_ASSERT1(ifst.good(), "impossible to read from file " << FILE_.string_value("MESH_FILEV","1D points file"));
			asm_tangent_versor(ifst, lx_temp,ly_temp,lz_temp);
			Curv_.resize(n_branch);
			lambdax_.resize(n_branch);
			lambday_.resize(n_branch);
			lambdaz_.resize(n_branch); 

			for(size_type b=0;b<n_branch;++b){
				size_type dofi=mf_datavi[b].nb_dof();
				Curv_[b].resize(dofi); Curv_[b].clear();
				Curv_[b].assign(dofi, 0.0);
				
				gmm::resize(lambdax_[b],dofi);
				gmm::resize(lambday_[b],dofi);
				gmm::resize(lambdaz_[b],dofi);
				
				lambdax_[b].assign(dofi,lx_temp[b]);
				lambday_[b].assign(dofi,ly_temp[b]);
				lambdaz_[b].assign(dofi,lz_temp[b]);
				
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
			rho_= FILE_.real_value("rho","density of the fluid");
			// Compute the dimenless params
			kt_.assign(dof_datat, k_/mu_*P_/U_/d_);
			G_.assign(dof_datav, U_*U_*rho_/P_);
			for (auto r : R_){ // C++11-only!
				kv_.emplace_back(pi/2.0/(Gamma_+2.0)/mu_*P_*d_/U_*r*r*r*r);
				Q_.emplace_back(2.0*pi*Lp_*P_/U_*r);
			}
		}
		alfa  = (Gamma_+2.0)/(Gamma_+1.0);
		beta  = (Gamma_+2.0)/(Gamma_+4.0);
		gamma = (Gamma_+2.0)*(Gamma_+2.0)/3.0/(Gamma_+3.0)/(Gamma_+6.0);


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


	//! Saving the curved parameters during the initialisation
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
	//! Get the vessel dynamic pressure gain coefficent element i
	inline scalar_type G  (size_type i) { return G_[i];  } const
	//! Get the esponent of velocity profile
	inline scalar_type Gamma (void) { return Gamma_;} const
	//! Get the vessel dynamic pressure gain coefficent vector
	vector_type G(void){ return G_;}
	//! Get vessel tangent versor x component
	vector<vector_type> & lambdax (void) { return lambdax_; }
	//! Get vessel tangent versor y component
	vector<vector_type> & lambday (void) { return lambday_; }
	//! Get vessel tangent versor z component
	vector<vector_type> & lambdaz (void) { return lambdaz_; }
	//! Get vessel curvature
	vector<vector_type> & Curv (void) { return Curv_; }
	//! Get the vessel Radius at a given dof
	scalar_type R  (const getfem::mesh_im & mim, const size_type rg) { 
		return simple_compute_radius(mim, mf_datav_, R_, rg);  
	}
	//! Get the vessel wall permeability at a given dof
	scalar_type kv  (const getfem::mesh_im & mim, const size_type rg) { 
		return simple_compute_radius(mim, mf_datav_, kv_, rg);  
	}
	//! Get the radius
	vector_type & R (void) { return R_; }
 
}; /* end of class */

} /* end of namespace */

#endif 