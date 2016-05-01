  
namespace getfem {

  // struct to handle list of Boundary Conditions

  class BCData {
  public :
    char label[4];     // DIR, NEU or MIX
    size_type ind;     // An index
    scalar_type Val1;  // A Value
    scalar_type Val2;  // Another Value (optional, used for mixed BC)
    size_type cv;      // Convex index (optional, used for networks)
  };



  //Compute the integral of the solution

  template<typename VEC>
  scalar_type asm_mean(const mesh_fem &mf, const mesh_im &mim, const VEC &U) {
    getfem::generic_assembly assem;
    assem.set("u=data(#1); V()+=u(i).comp(Base(#1))(i)");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(U);
    std::vector<scalar_type> v(1);
    assem.push_vec(v);
    assem.assembly();
    return v[0];
  }


  scalar_type element_estimate_h(const mesh &mesh, const size_type i) {
    std::vector<size_type> cpt = mesh.ind_points_of_convex(i);
    std::vector<size_type>::const_iterator icpt, jcpt;
    scalar_type d=0.0;
    for (icpt = cpt.begin(); icpt != cpt.end(); icpt++ ) {
      for (jcpt = cpt.begin(); jcpt != cpt.end(); jcpt++ ) {
	d = std::max(d, gmm::vect_norm2(mesh.points()[*icpt] - mesh.points()[*jcpt]));
      }
    }
    return d;
  }


  template<typename MAT>
  void build_exchange_matrices(MAT &Mbar, MAT &Mlin, MAT &M0, 
			       const getfem::mesh_fem &mf_ut, 
			       const getfem::mesh_fem &mf_coeft, 
			       const getfem::mesh_fem &mf_uv, 
			       const getfem::mesh_fem &mf_coefv,
                               const getfem::mesh_im &mim,
			       const scalar_type RADIUS,
			       const size_type NInt) {

    // Useful constants:
    // Pi
    const scalar_type Pi = 2* acos(0.0);

    size_type nb_dof_t = mf_ut.nb_dof();
    size_type nb_dof_coeft = mf_coeft.nb_dof();
    size_type nb_dof_v = mf_uv.nb_dof();
    size_type nb_dof_coefv = mf_coefv.nb_dof();


    gmm::clear(Mbar);
    gmm::clear(Mlin);
    gmm::clear(M0); 

    std::vector<scalar_type> Ubari(NInt); 
    std::vector<scalar_type> Ut(nb_dof_t); 
    std::vector<scalar_type> Uv(nb_dof_v); 


    // Compute Mlin --> We need the following interpolation tool <getfem_interpolation.h>   
   
    interpolation(mf_ut, mf_uv,Mlin);

    // Compute Mbar.
    cout << "1...5...10...15...20...25...30...35...40...45...50"
	 << "...55...60...65...70...75...80...85...90...95..100" << endl;
    size_type counter = 0;
    for (size_type i = 0; i < nb_dof_coefv; i++){
      counter++;
      if ( counter*100 >= nb_dof_coefv ) {
	counter = 0; cout << "*"; cout.flush();
      }
      // We need the following interpolation tool <getfem_interpolation.h>      
      getfem::mesh_trans_inv mti(mf_ut.linked_mesh());
      // Build the list of point on the i-th circle:
      // ... first consider an orthonormal system v0, v1, v2:
      size_type cv_id = mf_coefv.first_convex_of_dof(i);
      base_node v0 = mf_coefv.linked_mesh().points_of_convex(cv_id)[0] - 
	mf_coefv.linked_mesh().points_of_convex(cv_id)[1];
      base_node v1(0.0, -v0[2], v0[1]);
      base_node v2(v0[1]*v0[1] +v0[2]*v0[2], -v0[0]*v0[1], -v0[0]*v0[2]);
      if ( gmm::vect_norm2(v2) < 1.0e-8 * gmm::vect_norm2(v0) ) {
	v1[0] = -v0[1]; v1[1] = v0[0]; v1[2] = 0.0;
	v2[0] = -v0[0]*v0[2]; v2[1] = -v0[1]*v0[2]; v2[2] = v0[0]*v0[0] +v0[1]*v0[1];
      }
      v1 = v1 / gmm::vect_norm2(v1);
      v2 = v2 / gmm::vect_norm2(v2);
      // ... then parametrize the circle:
      for(size_type j = 0; j < NInt; j++) 
	mti.add_point( mf_coefv.point_of_dof(i) + 
		       RADIUS*(cos(2*Pi*j/NInt)*v1 + sin(2*Pi*j/NInt)*v2) );

      // Get the local interpolation matrix Mbari, Mbar0:
      MAT Mbari(NInt, nb_dof_t); gmm::clear(Mbari);
      // Remark: here the last argument value determines wether
      // we write on Ubari (0) or Mbari (1). In the latter case,
      // vector Ut is obviously not used by the routine.
     interpolation(mf_ut, mti, Ut, Ubari, Mbari, 1);
     scalar_type sum_row = 0.0;
     for (size_type j=0; j < NInt; ++j) {
	typename gmm::linalg_traits<MAT>::const_sub_row_type
	  row = mat_const_row(Mbari,j);
	gmm::linalg_traits< gmm::rsvector<scalar_type> >::const_iterator
	  it_nz = vect_const_begin(row);
	gmm::linalg_traits< gmm::rsvector<scalar_type> >::const_iterator
	  ite_nz = vect_const_end(row);
	for (; it_nz != ite_nz ; ++it_nz) {
	  Mbar(i, it_nz.index()) += (*it_nz);
	  sum_row += (*it_nz);
	}
      }
      typename gmm::linalg_traits<MAT>::sub_row_type
	row = mat_row(Mbar,i);
      gmm::linalg_traits< gmm::rsvector<scalar_type> >::iterator
	it_nz = vect_begin(row);
      gmm::linalg_traits< gmm::rsvector<scalar_type> >::iterator
	ite_nz = vect_end(row);
      for (; it_nz != ite_nz ; ++it_nz) {
	(*it_nz)/=sum_row;
	}
      }


    // We will need this Pk-P0 mixed mass matrix:
    getfem::asm_mass_matrix(M0, mim, mf_uv, mf_coefv);     
    
  }



  template<typename MAT>
  void coefficient_exchange_matrices(MAT &Mtt, MAT &Mtv, MAT &Mvt, MAT &Mvv, 
			       MAT &Mbar, MAT &Mlin, MAT &M0, 
			       const getfem::mesh_fem &mf_ut, 
			       const getfem::mesh_fem &mf_coeft, 
			       const getfem::mesh_fem &mf_uv, 
			       const getfem::mesh_fem &mf_coefv,
                               const getfem::mesh_im &mimt,
                               const getfem::mesh_im &mimv,
			       std::vector<scalar_type> &Cvold,
			       std::vector<scalar_type> &Ctold, 
			       const scalar_type OMEGA1,
			       const scalar_type OMEGA2,
			       const scalar_type GAMMA,
			       const scalar_type EPSILON,
			       const scalar_type ALPHA,	
			       const scalar_type Cm,		       
			       const scalar_type RADIUS,
			       const size_type NInt) {



    cout << "Sono entrato" << endl;
    size_type nb_dof_t = mf_ut.nb_dof();
    size_type nb_dof_coeft = mf_coeft.nb_dof();
    size_type nb_dof_v = mf_uv.nb_dof();
    size_type nb_dof_coefv = mf_coefv.nb_dof();

    // Structure to store the new coefficient

    // We will need this PO-Pk mixed mass matrix, to compute w1*ccold


    MAT MCoef(nb_dof_coefv, nb_dof_v);
    std::vector<scalar_type> omega1_v(nb_dof_coefv, OMEGA1);
    getfem::asm_mass_matrix(MCoef, mimv, mf_coefv, mf_uv); 
    gmm::scale(MCoef,OMEGA1);

    MAT MbarCoef(nb_dof_coefv, nb_dof_t); gmm::clear(MbarCoef);
    gmm::copy(Mbar, MbarCoef);
    gmm::scale(MbarCoef,OMEGA2);

    std::vector<scalar_type> coldparam;
    gmm::resize(coldparam, nb_dof_coefv); gmm::clear(coldparam); 

    std::vector<scalar_type> epsilon_v(nb_dof_coefv, EPSILON);
    std::vector<scalar_type> cm_v(nb_dof_coefv, Cm);

    gmm::mult_add(MCoef, Cvold, coldparam);

    gmm::mult_add(MbarCoef, Ctold, coldparam);

    gmm::scale(coldparam, ALPHA);

    gmm::scale(cm_v, 1.0-ALPHA);

    gmm::add(cm_v, coldparam);

    gmm::scale(coldparam, GAMMA);

    gmm::add(epsilon_v, coldparam );

    MAT Mbartemp(nb_dof_coefv, nb_dof_t); clear(Mbartemp);
    gmm::copy(Mbar, Mbartemp);

    for (size_type i = 0; i < nb_dof_coefv; i++){
      typename gmm::linalg_traits<MAT>::sub_row_type
	row = mat_row(Mbartemp,i);
      gmm::linalg_traits< gmm::rsvector<scalar_type> >::iterator
	it_nz = vect_begin(row);
      gmm::linalg_traits< gmm::rsvector<scalar_type> >::iterator
	ite_nz = vect_end(row);
      for (; it_nz != ite_nz ; ++it_nz) {
	(*it_nz) *=coldparam[i];
	}
        }


   // Build Mvv:
    scalar_type time = gmm::uclock_sec();  

    getfem::asm_mass_matrix_param(Mvv, mimv, mf_uv, mf_coefv, coldparam); 

    cout << "-o Mvv assembled in " << gmm::uclock_sec() - time 
	 << " seconds." << endl;  

    // Now we build Mvt

    time = gmm::uclock_sec();
    cout << "-o Assembling Exchange 1D-3D terms..." << endl;
    cout << endl << "-o Creating Mvt" << endl;

    gmm::mult(M0, Mbartemp, Mvt);

    cout << "-o Creating Mtt" << endl;

    MAT Mtemp(nb_dof_t, nb_dof_coefv);

    gmm::mult(gmm::transposed(Mlin), M0, Mtemp);

    gmm::mult(Mtemp, Mbartemp, Mtt);

    gmm::clear(Mtemp); // this DOES free the allocated memory

    cout << "-o Creating Mtv" << endl;

    gmm::mult(gmm::transposed(Mlin), Mvv, Mtv);

    gmm::clear(Mbartemp);

}


  template<typename MAT, typename VEC>
  void build_boundary_terms(MAT &M, 
			    const getfem::mesh_fem &mf, 
			    const getfem::mesh_fem &mf_data, 
                            const getfem::mesh_im &mim,
			    const VEC &K,
			    VEC &B,
			    const scalar_type GAMMA, 
			    const std::vector< BCData > &BCDataVect ) {
    
    
    // Aux data
    std::vector<scalar_type> inv_h_mesh(mf_data.nb_dof());
    std::vector<scalar_type> ones(mf_data.nb_dof(), 1.0);
    std::vector<scalar_type> onesu(mf.nb_dof(), 1.0);
    
    for(size_type i=0; i<mf_data.nb_dof(); i++) {
      inv_h_mesh[i] = GAMMA/element_estimate_h(mf.linked_mesh(),mf_data.first_convex_of_dof(i)) ; ///mf_data.linked_mesh().convex_radius_estimate(mf_data.first_convex_of_dof(i));
     }
    
    // Surface integration, consistency + Dirichlet BC (penalization)
    getfem::generic_assembly
      AssemD("wa=data$1(#2);"
	     "wb=data$2(#2);"
	     "a=comp(Base(#1).Grad(#1).Normal().Base(#2));" 
	     "at=comp(Grad(#1).Base(#1).Normal().Base(#2));"
	     "b=comp(Base(#1).Base(#1).Base(#2));"
	     "M(#1,#1)+=(-a(:, :,i, i, k).wa(k) - at(:,i, :, i, k).wa(k) + b(:,:,h).wb(h));"
	     );
    AssemD.push_mi(mim);
    AssemD.push_mf(mf);
    AssemD.push_mf(mf_data);
    AssemD.push_data(K);
    AssemD.push_data(inv_h_mesh);
    AssemD.push_mat(M);
    
    for (size_type i=0; i < BCDataVect.size(); i++) {
      if (strcmp(BCDataVect[i].label, "DIR")==0) { // Dirichlet BC
	size_type BoundaryLabel =  BCDataVect[i].ind;
	scalar_type Pressure =  BCDataVect[i].Val1;
	AssemD.assembly(BoundaryLabel);
	getfem::asm_source_term(B, mim, mf, mf_data, 
				gmm::scaled(inv_h_mesh, Pressure), BoundaryLabel);
	getfem::generic_assembly
	  AssemF("wa=data$1(#2);"
		 "u=data$2(#1);"
		 "at=comp(Grad(#1).Base(#1).Normal().Base(#2));"
		 "V(#1)+=( -at(:,i, j, i, k).u(j).wa(k) );"
		 );
        AssemF.push_mi(mim);
	AssemF.push_mf(mf);
	AssemF.push_mf(mf_data);
	AssemF.push_data(K);
	AssemF.push_data(gmm::scaled(onesu, Pressure));
	AssemF.push_vec(B);
	AssemF.assembly(BoundaryLabel);
      } 
      else if (strcmp(BCDataVect[i].label, "NEU")==0) { // Neumann BC
	size_type BoundaryLabel =  BCDataVect[i].ind;
	scalar_type Flux =  BCDataVect[i].Val1;
	getfem::asm_source_term(B, mim, mf, mf_data, 
				gmm::scaled(ones, Flux), BoundaryLabel);
      }
      else if (strcmp(BCDataVect[i].label, "MIX")==0) { // Mixed BC
	size_type BoundaryLabel =  BCDataVect[i].ind;
	scalar_type Coef     =  BCDataVect[i].Val1;
	scalar_type Pressure =  BCDataVect[i].Val2;
	getfem::asm_mass_matrix_param(M, mim, mf, mf_data, gmm::scaled(ones, Coef), 
				      BoundaryLabel);
	getfem::asm_source_term(B, mim, mf, mf_data, gmm::scaled(ones, Coef*Pressure), 
				BoundaryLabel);
      }
      else if (strcmp(BCDataVect[i].label, "INT")==0) { // Internal Node
	// Do nothing
	cout << "Warning: an internal node has been passed as boundary." << endl;
      }
      else {
	cout << "Unknown BC: " << BCDataVect[i].label << endl;
	DAL_THROW(dal::failure_error, "Unknown Boundary Condition");
      }
    }

  }
  
  // Compute average quantities on different boundaries.
  // Flow rates \int (-K dp/dn) and mean pressures \int p are computed.
  template <typename VEC>
  void integrate_boundary_quantities(const getfem::mesh_fem &mf, 
				     const getfem::mesh_fem &mf_data, 
                                     const getfem::mesh_im &mim,
				     const VEC &U,
				     const VEC &K,
				     const std::vector< BCData > &BCDataVect ) {
    
    
    // Aux data
    std::vector<scalar_type> v(1);    
    
    // Surface integration for Flow rate:
    getfem::generic_assembly
      AssemF("u=data$1(#1);"
	     "w=data$2(#2);"
	     "a=comp(Grad(#1).Normal().Base(#2));"
	     "V()+=-a(i,k, k, h).u(i).w(h);"
	     );
    AssemF.push_mi(mim);
    AssemF.push_mf(mf);
    AssemF.push_mf(mf_data);
    AssemF.push_data(U);
    AssemF.push_data(K);
    AssemF.push_vec(v);

   
    // Surface integration for integrating pressure:
    getfem::generic_assembly
          AssemM("u=data(#1); V()+=u(i).comp(Base(#1))(i)");
    AssemM.push_mi(mim);
    AssemM.push_mf(mf);
    AssemM.push_data(U);
    AssemM.push_vec(v);



    // Surface integration for computing areas:
    std::vector<scalar_type> onesu(mf.nb_dof(), 1.0);
    getfem::generic_assembly
      AssemA("u=data(#1); V()+=u(i).comp(Base(#1))(i)");
    AssemA.push_mi(mim);
    AssemA.push_mf(mf);
    AssemA.push_data(onesu);
    AssemA.push_vec(v);


    
    for (size_type i=0; i < BCDataVect.size(); i++) {
      if (strcmp(BCDataVect[i].label, "DIR")==0) { // Dirichlet BC
	size_type BoundaryLabel =  BCDataVect[i].ind;
	scalar_type Pressure =  BCDataVect[i].Val1;
        gmm::clear(v);
        AssemA.assembly(BoundaryLabel);
	scalar_type Area = v[0];
	gmm::clear(v);
	AssemF.assembly(BoundaryLabel);
	cout << "Boundary label " << BoundaryLabel << "\t"
	     << "[DIRICHLET]" << "\t"
	     << "Area : " << Area << "\t"
	     << "Flow rate: " << v[0] << "\t";
	gmm::clear(v);
	AssemM.assembly(BoundaryLabel); 
	cout << "Mean pressure: " << v[0] / Area << "\t"
	     << "Imposed Pressure was: " << Pressure << endl;
      } 
      else if (strcmp(BCDataVect[i].label, "NEU")==0) { // Neumann BC
	size_type BoundaryLabel =  BCDataVect[i].ind;
	scalar_type Flux =  BCDataVect[i].Val1;
	gmm::clear(v);
	AssemA.assembly(BoundaryLabel);
	scalar_type Area = v[0];
	gmm::clear(v);
	AssemF.assembly(BoundaryLabel);
	cout << "Boundary label " << BoundaryLabel << "\t"
	     << "[NEUMANN  ]" << "\t" 
	     << "Area : " << Area << "\t"
	     << "Flow rate: " << v[0] << "\t";
	gmm::clear(v);
	AssemM.assembly(BoundaryLabel); 
	cout << "Mean pressure: " << v[0] / Area << "\t"
	     << "Imposed Flow Rate was: " << Flux*Area << endl;
      }
      else if (strcmp(BCDataVect[i].label, "MIX")==0) { // Mixed BC
	size_type BoundaryLabel =  BCDataVect[i].ind;
	scalar_type Coef     =  BCDataVect[i].Val1;
	scalar_type Pressure =  BCDataVect[i].Val2;
	gmm::clear(v);
	AssemA.assembly(BoundaryLabel);
	scalar_type Area = v[0];
	gmm::clear(v);
	AssemF.assembly(BoundaryLabel);
	cout << "Boundary label " << BoundaryLabel << "\t"
	     << "[ROBIN    ]" << "\t" 
	     << "Area : " << Area << "\t"
	     << "Flow rate: " << v[0] << "\t";
	gmm::clear(v);
	AssemM.assembly(BoundaryLabel); 
	cout << "Mean pressure: " << v[0] / Area << "\t";
	//
	VEC FU(U.size(), 0.0);
	for(size_type i=0; i<FU.size(); ++i)
	  FU[i] = Coef*(Pressure - U[i]);
	// Surface integration for Flow rate:
	getfem::generic_assembly
	  AssemF1("u=data$1(#1);"
		  "a=comp(Base(#1));"
		  "V()+=-a(i).u(i);"
		  );
        AssemF1.push_mi(mim);
	AssemF1.push_mf(mf);
	AssemF1.push_mf(mf_data);
	AssemF1.push_data(FU);
	AssemF1.push_vec(v);
	//
	gmm::clear(v);
	//
	AssemF1.assembly(BoundaryLabel); 
	cout << "Flow rate based on Lp: " << v[0] << endl;
	//


      }
      else if (strcmp(BCDataVect[i].label, "INT")==0) { // Internal Node
	// Do nothing
      }
      else
	DAL_THROW(dal::failure_error, "Unknown Boundary Condition");
    }
  }





}
