/*===========================================================================
 
 Copyright (C) 2002-2012 Yves Renard, Julien Pommier.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/



/**
 * We implement a IBM for a coupled problem.
 * We consider the evolution equation
 *
 * FLUID FLOW PROBLEM

 * - div(K p_t \nabla p_t) - f_f \delta = 0                             in \Omega, the tissue
 *
 *                    - (K p_c')' + f_f = 0                             on \Lambda, the vessel
 *
 *
 * DRUG TRANSPORT PROBLEM
 *
 *   Ct dc_t/dt - div(Kt \nabla c_t) + lambda c_t - f_c \delta = 0      in \Omega, the tissue
 *
 *               Cv dc_v /dt - (Kv c_v')' + f_c = 0                     on \Lambda, the vessel
 *
 * 
 * where \delta is the Dirac measure on the vessel \Lambda, 
 * and f_f, f_c are:: 
 * 
 * f_f(s) =  2\pi R' Q(p_c - \overline{p}_t)
 *
 * f_c(s) = 2\pi R' f'(c_c, \overline{c}_t) [Q(p_c - \overline{p}_t) -\Pi (c_c - \overline{c}_t)] + 2\pi R' P(c_c - \overline{c}_t)
 * 
 * being \overline{u}(s) = (2\pi)^{-1}\int_0^{2\pi} u(s,\rho(s),\theta) d\theta.
 * 
 * In general, \Lambda is a network; \Lambda and \Omega are not
 * conforming.
*/



#include <getfem/getfem_assembling.h> /* assembly methods (and comp. of norms) */
#include <getfem/getfem_import.h>
#include <getfem/getfem_export.h>   /* export functions (save solutions in a file) */

#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_mesh.h>
#include <getfem/getfem_derivatives.h>
#include <getfem/getfem_superlu.h>
#include <gmm/gmm.h>
#include <gmm/gmm_inoutput.h>
#include <getfem/getfem_derivatives.h>


#include "utilities.cc"          // finite element operators
#include "import_gmshv2.cc"         // the .neu mesh reader
#include "import_network.cc"     // the importer for the network

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

/* 888888888888888888888888888888888888888888888888888888888888888888888888888888 */

/* funzioni da aggiungere

#include "interpolation.C"      // interpolation in a FE space

//#include "utilities.C"          // finite element operators
#include "utilities_nonlinear3.C"          // finite element operators
//#include "utilities-modified.C"          // finite element operators
 
 
#include "ensight.C"            // the Ensight7 exporter (useful for time-dependent problems) 

*/


/* some Getfem++ types that we will be using */
using bgeot::base_small_vector;  /* special class for small (dim < 16) vectors */
using bgeot::base_node;   /* geometrical nodes (derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::short_type;


/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef gmm::rsvector<scalar_type> sparse_vector_type;
typedef gmm::wsvector<scalar_type> wsparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef gmm::row_matrix<wsparse_vector_type> wsparse_matrix_type;
typedef gmm::col_matrix<sparse_vector_type> col_sparse_matrix_type;


struct myProblem {

  //*********** Meshes: ****************************************************

  // For the 3D domain: we use subscript 't' to denote the tissue.
  // For the 1D domain (line): we use subscript 'v' to denote the vessel.

  getfem::mesh mesht;  // the mesh for the tissue (3D)
  getfem::mesh meshv;  // the mesh for the capillaries network (1D)


  // ************** FLUID PROBLEM *******************************************

  // ************** Setting the Finite Elements: ****************************

  getfem::mesh_im mimtF;      			// the integration methods.
  getfem::mesh_fem mf_utF;     			// the main mesh_fem, pressure
  getfem::mesh_im mim_coeftF;
  getfem::mesh_fem mf_coeftF;  			// the mesh_fem to represent pde coefficients 
  getfem::mesh_im  mim_gradienttF;    
  getfem::mesh_fem mf_gradienttF; 		// mesh_fem to represent P0 gradients
  getfem::mesh_im  mim_veltF; 
  getfem::mesh_fem mf_veltF; 			//mesh_fem to represent velocity P0, Q=3;

  getfem::mesh_im mimvF;      			// the integration methods. 
  getfem::mesh_fem mf_uvF;    			// the main mesh_fem, for the solution on the line  
  getfem::mesh_im mim_coefvF;
  getfem::mesh_fem mf_coefvF;  			// the mesh_fem to represent pde coefficients
  getfem::mesh_im  mim_gradientvF;   
  getfem::mesh_fem mf_gradientvF; 		// mesh_fem to represent P0 gradients
  getfem::mesh_im  mim_velvF; 
  getfem::mesh_fem mf_velvF; 			//mesh_fem to represent velocity P0, Q=3;

 // **************Setting vectors and matrices: **************************************************

  // Preconditioner, 3D (t) and 1D (v) problems
  gmm::ilu_precond<sparse_matrix_type> tPreconditionerF;
  gmm::ilu_precond<sparse_matrix_type> vPreconditionerF;

  sparse_matrix_type AttF;             		// stiffness matrix for fluid problem
  std::vector<scalar_type> Pt, Ptold, BtF;     	// main unknown, and right hand side 
  std::vector<scalar_type> PPt;         	// 3D pressure of fluid problem
  std::vector<scalar_type> gradPt;		//pressure gradient
  std::vector<scalar_type> velUt; 		//velocity Ut
  std::vector<scalar_type> velUtD;      	//velocity Ut to pass to drug problem


  sparse_matrix_type AvvF;             		// stiffness matrix 
  std::vector<scalar_type> Pv, Pvold, BvF;    	// main unknown, and right hand side 
  std::vector<scalar_type> PPv;	       		// 1D pressure of fluid problem
  std::vector<scalar_type> gradPv;	       	// pressure gradient
  std::vector<scalar_type> velUv;       	// velocity Ut
  std::vector<scalar_type> velUvD;     		//velocity Uv to pass to drug problem

  sparse_matrix_type MvvF;             		// vessel matrix
  sparse_matrix_type MttF;             		// tissue matrix
  sparse_matrix_type MvtF;             		// mixed vessel-tissue matrix
  sparse_matrix_type MtvF;             		// mixed tissue-vessel matrix

  sparse_matrix_type KKvvF;             	// vessel matrix (rhs term)
  sparse_matrix_type KKttF;            		// tissue matrix (rhs term)
  sparse_matrix_type KKvtF;             	// mixed vessel-tissue matrix (rhs term)
  sparse_matrix_type KKtvF;             	// mixed tissue-vessel matrix (rhs term)

  std::vector<scalar_type> Cvold, Ctold;


 // **************Setting vectors and matrices for the monolitic solver: ******************************

  sparse_matrix_type AM;
  std::vector<scalar_type> PM, BM;
   
  //*************************************************************************************************
  //*************************************************************************************************

  sparse_matrix_type M_tissueF;        		// unit mass matrix tissue
  sparse_matrix_type M_vesselF;        		// unit mass matrix vessel

  //*************************************************************************************************

  // matrices for the exchange terms

  sparse_matrix_type Mbar;        
  sparse_matrix_type M0;        
  sparse_matrix_type Mlin;  

  //*************************************************************************************************     
  // ************** Setting SCALAR VARIABLES: **************************************************

  scalar_type residu;        // max residual for the iterative solvers 

  scalar_type meanPt; 
  scalar_type meanPv;

  scalar_type FlowRateP;


  scalar_type       RADIUS;        // actual vessel radius
  scalar_type       GAMMA;         // Penalization parameter 
  size_type         NInt;          // Node number (integration on circles)

  const char       *SOLVE_METHOD_T;// 'GS', 'CG', 'SUPERLU', 'GMRES', etc.;
  const char       *SOLVE_METHOD_V;// 'GS', 'CG', 'SUPERLU', 'GMRES', etc.;
  const char       *SOLVE_METHOD_M;// 'GS', 'CG', 'SUPERLU', 'GMRES', etc.;


  size_type         MAX_SUBITERF;
  scalar_type       OMEGA1;         	// weight for the non linear term
  scalar_type       OMEGA2;         	// weigth for the non linear term
  scalar_type wt;
  scalar_type wc;

  scalar_type       Peclet;


  scalar_type       ApresF;          	// exchange constant, fluid problem
  scalar_type       BpresF;          	// exchange constant (Kedem - Katchalsky), fluid problem
  scalar_type       AconcF;          	// exchange constant, fluid problem
  scalar_type       BconcF;          	// exchange constant (Kedem - Katchalsky), fluid problem

  scalar_type       alphaF;         	// non linear term selection
  scalar_type       Cm;             	// reference concentration 

  scalar_type       TFINAL;        	// Final time
  scalar_type       DT;            	// dt

  scalar_type VOLT;
  scalar_type VOLV;
  scalar_type velT;
  scalar_type velV;

  //*************************************************************************************************     
  // ************** FILE reading and OUTPUT tool **************************************************

  ftool::md_param PARAM;     

  std::string OutputDir;    		// Output Directory



  //*************************************************************************************************     
  // ************** METHODS **************************************************



  myProblem(void) : mimtF(mesht), mf_utF(mesht), mim_coeftF(mesht), mf_coeftF(mesht), mim_gradienttF(mesht), mf_gradienttF(mesht), mim_veltF(mesht), mf_veltF(mesht), 
		    mimvF(meshv), mf_uvF(meshv), mim_coefvF(meshv), mf_coefvF(meshv), mim_gradientvF(meshv), mf_gradientvF(meshv), mim_velvF(meshv), mf_velvF(meshv) {}

  void Pressureaverage(void);

  void Computew(void);

  void init(void);

  void assemblyFluid(void);

  void compute_velocity_Ut(void);

  void compute_velocity_Uv(void);

  void assemblyRHSFluid(void);

  bool time_advanceFluid(const scalar_type& time, const scalar_type& dt);

  bool time_advanceFluidMONO(void);
  


  // ************** FLUID TRANSPORT PROBLEM: *******************************************************

  // TISSUE
  // Compliance
  scalar_type CompltF;  

  // Permeability
  scalar_type KtValF;
  scalar_type KtFuncF(const base_node &x) {
    return KtValF;
  }
  std::vector<scalar_type> KtF;  

  // Mass term
  scalar_type MtValF;
  scalar_type MtFuncF(const base_node &x) {
    return MtValF;
  }
  std::vector<scalar_type> MtF;        

  // reaction term
  scalar_type RtValF;
  scalar_type RtFuncF(const base_node &x) {
    return RtValF;
  }
  std::vector<scalar_type> RtF;

  // VESSEL
  // Compliance
  scalar_type ComplvF;  
  
  // Permeability
  scalar_type KvValF;
  scalar_type KvFuncF(const base_node &x) { 
    return KvValF; 
  }
  std::vector<scalar_type> KvF;        

  // Mass term
  scalar_type MvValF;
  scalar_type MvFuncF(const base_node &x) {
    return MvValF;
  }
  std::vector<scalar_type> MvF;

  // Reaction term
  scalar_type RvValF;
  scalar_type RvFuncF(const base_node &x) {
    return RvValF;
  }
  std::vector<scalar_type> RvF;


  // **************  STORE THE BCs: *******************************************************

  size_type         BCList_SizeF;
  std::vector< getfem::BCData > BCListF;
  std::vector< getfem::BCData > networkBCListF;

};

/**************************************************************************/
/* Methods definition.                                                         */
/**************************************************************************/


  void myProblem::Pressureaverage(void){
  meanPt = getfem::asm_mean(mf_utF, mimtF, Pt);
  meanPv = getfem::asm_mean(mf_uvF, mimvF, Pv);
  }


  void myProblem::Computew(void){
  Pressureaverage(); 
  scalar_type meanPvV =  meanPv/VOLV;
  scalar_type meanPtV =  meanPt/VOLT;
  scalar_type Pe = Peclet*abs(meanPvV - meanPtV);
  wt = 1.0/Pe- 1.0/(exp(Pe)-1.0);
  cout << "Valore di wt " << wt << endl;
  wc = exp(Pe)/(exp(Pe)-1)-1.0/Pe;
  cout << "Valore di wc " << wc << endl;
  }




/* ***** Read parameters from the .param file, build the mesh, set finite element
         and integration methods and select the boundaries. */

 void myProblem::init(void)
 {

 std::ofstream timetxtSFE( "./Vtk/Secomb93b/timeSFElc025.txt" ); 
 double timeSFE = gmm::uclock_sec();


 std::string MESH_FILE = PARAM.string_value("MESH_FILE",".MSH Mesh file");
 std::string MESH_TYPET = PARAM.string_value("MESH_TYPET","3D mesh type");
 std::string MESH_TYPEV = PARAM.string_value("MESH_TYPEV","1D mesh type");
 std::string FEM_TYPET  = PARAM.string_value("FEM_TYPET","FEM 3D tissue");
 std::string FEM_TYPEV  = PARAM.string_value("FEM_TYPEV","FEM 1D vessels");
 std::string INTEGRATION_T = PARAM.string_value("INTEGRATION_T","Name of integration method");
 std::string INTEGRATION_V = PARAM.string_value("INTEGRATION_V","Name of integration method");
 const char *LINE_FILEF = PARAM.string_value("LINE_FILEF","File with line points for fluid").c_str();


  RADIUS = PARAM.real_value("RADIUS", "Line radius");
  GAMMA  = PARAM.real_value("GAMMA", "Penalty Coefficient");
  OMEGA1 = PARAM.real_value("OMEGA1", "weight coefficient");      // weight for the non linear term
  OMEGA2 = PARAM.real_value("OMEGA2", "weight coefficient");      // weigth for the non linear term
  VOLT = PARAM.real_value("VOLT", "Tissue volume");
  VOLV  = PARAM.real_value("VOLV", "Total capillaries length");
  velT = PARAM.real_value("velT", "velocity coefficient T");
  velV  = PARAM.real_value("velV", "velocity coefficient V");

  Peclet = PARAM.real_value("PECLET", "To compute the Peclet numeber"); 

  ApresF = PARAM.real_value("ApresF", "exchange constant, fluid problem, pressure term");          // exchange constant, fluid problem
  BpresF = PARAM.real_value("BpresF", "exchange constant, fluid problem, pressure term");          // exchange constant fluid problem
  AconcF = PARAM.real_value("AconcF", "exchange constant, fluid problem, concentration term");        // exchange constant, fluid problem
  BconcF = PARAM.real_value("BconcF", "exchange constant, fluid problem, concentration term");          // exchange constant, fluid problem

  alphaF = PARAM.real_value("alphaF", "non linear term selection");         // non linear term selection
  Cm =  PARAM.real_value("Cm", "reference concentration");            // reference concentration


  // ********************************* Algorithm information *********************************

  TFINAL = PARAM.real_value("TFINAL", "Final time");
  DT     = PARAM.real_value("DT", "dt");

  NInt = size_type(PARAM.int_value("NInt",
				   "Node numbers on the circle for the nonlocal term"));  
  SOLVE_METHOD_T = PARAM.string_value("SOLVE_METHOD_T", "Solver (3D)").c_str();  
  SOLVE_METHOD_V = PARAM.string_value("SOLVE_METHOD_V", "Solver (1D)").c_str();
  SOLVE_METHOD_M = PARAM.string_value("SOLVE_METHOD_M", "Monolitic Solver").c_str();
  MAX_SUBITERF = PARAM.int_value("MAX_SUBITERF", "Max number of sub-iterations");
  residu = PARAM.real_value("RESIDU"); if (residu == 0.) residu = 1e-10;

  // ********************************* FLUID TRANSPORT PARAMETERS *********************************

  CompltF = PARAM.real_value("CtF", "Compliance, tissue, FLUID PROBLEM");
  KtValF  = PARAM.real_value("KtF", "Diffusivity, tissue, FLUID PROBLEM");
  MtValF  = PARAM.real_value("MtF", "Mass term, tissue, FLUID PROBLEM");
  RtValF  = PARAM.real_value("RtF", "Reaction term, tissue, FLUID PROBLEM");
  ComplvF = PARAM.real_value("CvF", "Compliance, tissue, FLUID PROBLEM");
  KvValF  = PARAM.real_value("KvF", "Diffusivity, vessel, FLUID PROBLEM");
  MvValF  = PARAM.real_value("MvF", "Mass term, vessel, FLUID PROBLEM");
  RvValF  = PARAM.real_value("RvF", "Reaction term, vessel, FLUID PROBLEM");


  // ********************************* BOUNDARIES *********************************

  BCList_SizeF 
    = (size_type)(PARAM.int_value("NumberOfBCF",
				  "Number of Boundary Conditions (tissue)"));


  // ********************************* Output directory *********************************

  std::string ODIR = PARAM.string_value("OutputDir","OutputDirectory");
  OutputDir = ODIR;

  // ********************************* Check output *********************************

  cout << "MESH_FILE         = " << MESH_FILE << endl;
  cout << "FEM TYPE tissue   = " << FEM_TYPET << endl;
  cout << "FEM TYPE vessels  = " << FEM_TYPEV << endl;
  cout << "LINE_FILEF        = " << LINE_FILEF << endl;
  cout << "INTEGRATION_T     = " << INTEGRATION_T << endl;
  cout << "INTEGRATION_V     = " << INTEGRATION_V << endl << endl;
  cout << "RADIUS            = " << RADIUS  << endl;
  cout << "ApresF            = " << ApresF << endl;
  cout << "BpresF            = " << BpresF << endl;
  cout << "AconcF            = " << AconcF << endl;
  cout << "BconcF            = " << BconcF << endl;
  cout << "NInt              = " << NInt << endl;
  cout << "SOLVE_METHOD_T    = " << SOLVE_METHOD_T << endl;
  cout << "SOLVE_METHOD_V    = " << SOLVE_METHOD_V << endl;
  cout << "MAX_SUBITERF      = " << MAX_SUBITERF << endl;
  cout << "TFINAL            = " << TFINAL       << endl;
  cout << "DT                = " << DT           << endl;


  // ***************** FLUID PROBLEM: read and save the BCs for the 3D mesh:********************

  for (size_type i=0; i<BCList_SizeF; i++) {
    getfem::BCData BC;
    char PNUM[1024]; 
    sprintf(PNUM, "BCF%dType", i); 
    std::string stringBCTYPE = PARAM.string_value(PNUM, "BC Type");
    const char *BCTYPE = stringBCTYPE.c_str();
    sprintf(PNUM, "BCF%dLabel", i);
    size_type pLabel = (size_type)PARAM.int_value(PNUM, "Boundary Label");
    if (strcmp(BCTYPE, "DIR")==0) { // Dirichlet BC
      sprintf(PNUM, "BCF%dPressure", i);
      strcpy(BC.label, BCTYPE);
      BC.ind = pLabel;
      BC.Val1  = PARAM.real_value(PNUM, "Pressure");
    } 
    else if (strcmp(BCTYPE, "NEU")==0) { // Neumann BC
      sprintf(PNUM, "BCF%dFlux", i);
      strcpy(BC.label, BCTYPE);
      BC.ind = pLabel;
      BC.Val1  = PARAM.real_value(PNUM, "Flux");      
    }
    else if (strcmp(BCTYPE, "MIX")==0) { // Mixed BC
      sprintf(PNUM, "BCF%dCoef", i);
      strcpy(BC.label, BCTYPE);
      BC.ind = pLabel;
      BC.Val1  = PARAM.real_value(PNUM, "Coef");
      sprintf(PNUM, "BCF%dPressure", i);
      BC.Val2  = PARAM.real_value(PNUM, "Pressure");      
    }
    else
      DAL_THROW(dal::failure_error, "Unknown Boundary Condition");
    BCListF.push_back(BC);
  }


  // ***************** Importing the 3D mesh *****************
  cout << "Importing the 3D mesh  "   << endl;
  getfem::import_mshv2_file(MESH_FILE, mesht);
  cout << "Imported mesh from file " << MESH_FILE << " [OK]" << endl;  

  // ***************** Importing the vessel network *****************
  cout << "Importing the vessel network "   << endl;
  std::ifstream ifsF(LINE_FILEF);

  if (!ifsF)
    DAL_THROW(dal::failure_error, "impossible to read from file " << LINE_FILEF);
  cout << " Chiamo import_pts_file su networkBCListF " << endl;
  getfem::import_pts_file(ifsF, meshv, MESH_TYPEV, networkBCListF);


  // ***************** set the finite element on the mf_utF and mf_utD *****************

  cout << "set the finite element on the mf_utF"   << endl;
  getfem::pfem pf_ut = getfem::fem_descriptor(FEM_TYPET);
  getfem::pintegration_method ppi_t = getfem::int_method_descriptor(INTEGRATION_T);
  mimtF.set_integration_method(mesht.convex_index(), ppi_t);
  mf_utF.set_finite_element(mesht.convex_index(), pf_ut);

  // set the finite element on the mf_uvF and mf_uvD
  cout << "set the finite element on the mf_uvF"   << endl;
  getfem::pfem pf_uv = getfem::fem_descriptor(FEM_TYPEV);
  getfem::pintegration_method ppi_v = getfem::int_method_descriptor(INTEGRATION_V);
  mimvF.set_integration_method(meshv.convex_index(), ppi_v);
  mf_uvF.set_finite_element(meshv.convex_index(), pf_uv);



  bgeot::pgeometric_trans pgt_t = bgeot::geometric_trans_descriptor(MESH_TYPET);
  size_type N = pgt_t->dim();

  // set the finite element on mf_coeft. Here we use a P0 (piecewise constant) elements
  cout << "set the finite element on the mf_coeftF"   << endl;
  mim_coeftF.set_integration_method(mesht.convex_index(), ppi_t);
  mf_coeftF.set_finite_element(mesht.convex_index(),
			     getfem::classical_fem(pgt_t,0));


 
  bgeot::pgeometric_trans pgt_v = bgeot::geometric_trans_descriptor(MESH_TYPEV);

  cout << "set the finite element on the mf_coefvF"   << endl;
  mim_coefvF.set_integration_method(meshv.convex_index(), ppi_v);
  mf_coefvF.set_finite_element(meshv.convex_index(),
			      getfem::classical_fem(pgt_v,0));

  // set the finite element on mf_gradient and mf_vel. We use P0 elements

  cout << "set the finite element on the mf_gradientvF, mf_velvF and mf_veltF"   << endl;
  mim_gradienttF.set_integration_method(mesht.convex_index(), ppi_t);
  mf_gradienttF.set_finite_element(mesht.convex_index(), getfem::classical_fem(pgt_t,0));

  mim_veltF.set_integration_method(mesht.convex_index(), ppi_t);
  mf_veltF.set_finite_element(mesht.convex_index(), getfem::classical_fem(pgt_t,0)); 
  mf_veltF.set_qdim(3);
 
  mim_gradientvF.set_integration_method(meshv.convex_index(), ppi_v);
  mf_gradientvF.set_finite_element(meshv.convex_index(), getfem::classical_fem(pgt_v,0));

  mim_velvF.set_integration_method(meshv.convex_index(), ppi_v);
  mf_velvF.set_finite_element(meshv.convex_index(), getfem::classical_fem(pgt_v,0)); 
  mf_velvF.set_qdim(3);


  //******************** FLUID PROBLEM: SET THE BCs: ************************


  cout << "Selecting BCs for the fluid problem, tissue:" << endl;
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesht, border_faces);
  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
    assert(i.is_face());
    base_node un = mesht.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);
    if (gmm::abs(un[N-1] - 1.0) < 1.0E-7) { 
      mesht.region(0).add(i.cv(), i.f());
    }else if (gmm::abs(un[N-1] + 1.0) < 1.0E-7)
	mesht.region(0).add(i.cv(), i.f());
	else {
      mesht.region(1).add(i.cv(), i.f());
    	}
  }



  cout << "Selecting BCs for the fluid problem, vessel " << endl;
  getfem::convex_face_ct border_facesv;
  extrema_of_oneDmesh(meshv, meshv.convex_index(), border_facesv);
  for (getfem::convex_face_ct::const_iterator it = border_facesv.begin(); it != border_facesv.end(); ++it) {   //
    size_type i=0; bool found = false;
    while (!found && (i<networkBCListF.size())) {
      found = (it->cv == networkBCListF[i].cv);
      if (!found) i++;
    }
    if (found && (strcmp(networkBCListF[i].label,"INT") != 0)) {
    meshv.region(networkBCListF[i].ind).add(it->cv, it->f);
    }
  }


  timetxtSFE << " Reading parameters, set finite elements and BCs, done in  = " << gmm::uclock_sec() - timeSFE << " seconds" << endl;
  cout << " done in  = " << gmm::uclock_sec() - timeSFE << " seconds" << endl;
  timetxtSFE.close();

  //******************** Initialization exchange matrices ************************


  std::ofstream timetxtEM( "./Vtk/Secomb93b/timeEMlc025.txt" ); 
  double timeF = gmm::uclock_sec();

  size_type nb_dof_tF = mf_utF.nb_dof();
  size_type nb_dof_coeftF = mf_coeftF.nb_dof();
  size_type nb_dof_vF = mf_uvF.nb_dof();
  size_type nb_dof_coefvF = mf_coefvF.nb_dof();


  cout<< "Dimensioni del fluid problem: "<< endl;
  cout<< " nb_dof_tF = " << nb_dof_tF << endl;
  cout<< " nb_dof_coeftF = " << nb_dof_coeftF << endl;
  cout<< " nb_dof_vF = " << nb_dof_vF << endl;
  cout<< " nb_dof_coefvF = " << nb_dof_coefvF << endl;
  
  gmm::resize(Mbar, nb_dof_coefvF, nb_dof_tF); gmm::clear(Mbar);
  gmm::resize(Mlin, nb_dof_vF, nb_dof_tF); gmm::clear(Mlin);
  gmm::resize(M0, nb_dof_vF, nb_dof_coefvF); gmm::clear(M0);


  cout << " Construction of exchange matrices, fluid problem" << endl;

  getfem::build_exchange_matrices(Mbar, Mlin, M0, mf_utF, mf_coeftF, mf_uvF, mf_coefvF, mimvF, RADIUS, NInt);

  timetxtEM << " Construction of exchange matrices, fluid problem, done in  = " << gmm::uclock_sec() - timeF << " seconds" << endl;
  cout << " done in  = " << gmm::uclock_sec() - timeF << " seconds" << endl;
  timetxtEM.close();

 }


/*
 * Discrete problem: solve the following system
 *
 *      [ Att + Mtt    -Mtv    ] [ u ] = [ bu ]
 *      [  -Mvt      Avv + Mvv ] [ v ] = [ bv ]
 *
 * where A is the stiffness matrix. When the nonlocal Robin term is neglected
 * M is the mass matrix on \Lambda, and subscripts
 * indicate the mixed integration between tissue (t, 3D) and vessel (v, 1D) 
 * finite element test functions. Note that Mtv = transpose(Mvt).
 * When the nonlocal term is taken into account, only Mvv is a mass matrix,
 * wheras Mtt and Mtv correspond to the bilinear form 
 *                r(u,v) = \int_\Lambda( \bar{u}(s) v(s) ) ds
 *
 */


void myProblem::assemblyFluid(void)
 { 

  // STAGE 1. INITIALIZATION ===============================================
  cout << "Initialization " << endl;
  
  size_type nb_dof_t = mf_utF.nb_dof();
  size_type nb_dof_coeft = mf_coeftF.nb_dof();
  size_type nb_dof_gradt = mf_gradienttF.nb_dof();  
  size_type nb_dof_velt  = mf_veltF.nb_dof();

  size_type nb_dof_v = mf_uvF.nb_dof();
  size_type nb_dof_coefv = mf_coefvF.nb_dof();
  size_type nb_dof_gradv = mf_gradientvF.nb_dof();
  size_type nb_dof_velv  = mf_velvF.nb_dof();


  // FLUID components *******************************************************
 
  gmm::resize(BtF, nb_dof_t); gmm::clear(BtF);
  gmm::resize(Pt, nb_dof_t); gmm::clear(Pt); 
  gmm::resize(Ptold, nb_dof_t); gmm::clear(Ptold); 
  gmm::resize(PPt, nb_dof_t); gmm::clear(PPt); 
  gmm::resize(gradPt, 3*nb_dof_gradt); gmm::clear(gradPt); 
  gmm::resize(velUt, 3*nb_dof_gradt); gmm::clear(velUt); 
  gmm::resize(velUtD, nb_dof_velt); gmm::clear(velUtD); 

  gmm::resize(AttF, nb_dof_t, nb_dof_t); gmm::clear(AttF);

  gmm::resize(BvF, nb_dof_v); gmm::clear(BvF);
  gmm::resize(Pv, nb_dof_v); //gmm::clear(Pv); 
  gmm::resize(Pvold, nb_dof_v); gmm::clear(Pvold); 
  gmm::resize(PPv, nb_dof_v); gmm::clear(PPv);
  gmm::resize(gradPv, 3*nb_dof_gradv); gmm::clear(gradPv);
  gmm::resize(velUv, 3*nb_dof_gradv); gmm::clear(velUv);
  gmm::resize(velUvD, nb_dof_velv); gmm::clear(velUvD);

  gmm::resize(AvvF, nb_dof_v, nb_dof_v); gmm::clear(AvvF);
  
  gmm::resize(MvtF, nb_dof_v, nb_dof_t); gmm::clear(MvtF);
  gmm::resize(MtvF, nb_dof_t, nb_dof_v); gmm::clear(MtvF);
  gmm::resize(MvvF, nb_dof_v, nb_dof_v); gmm::clear(MvvF);
  gmm::resize(MttF, nb_dof_t, nb_dof_t); gmm::clear(MttF);

  gmm::resize(KKvtF, nb_dof_v, nb_dof_t); gmm::clear(KKvtF);
  gmm::resize(KKtvF, nb_dof_t, nb_dof_v); gmm::clear(KKtvF);
  gmm::resize(KKvvF, nb_dof_v, nb_dof_v); gmm::clear(KKvvF);
  gmm::resize(KKttF, nb_dof_t, nb_dof_t); gmm::clear(KKttF);

  gmm::resize(M_tissueF, nb_dof_t, nb_dof_t); gmm::clear(M_tissueF);
  gmm::resize(M_vesselF, nb_dof_v, nb_dof_v); gmm::clear(M_vesselF);

  gmm::resize(AM, nb_dof_t+nb_dof_v, nb_dof_t+nb_dof_v); gmm::clear(AM);
  gmm::resize(PM, nb_dof_t+nb_dof_v); gmm::clear(PM);
  gmm::resize(BM, nb_dof_t+nb_dof_v); gmm::clear(BM);
  

  // FLUID PARAMETERS *******************************************************

  gmm::resize(KtF, nb_dof_coeft);
  gmm::resize(MtF, nb_dof_coeft);
  gmm::resize(RtF, nb_dof_coeft);
  for (size_type i = 0; i < nb_dof_coeft; ++i) {
    KtF[i] = KtFuncF(mf_coeftF.point_of_dof(i));
    MtF[i] = MtFuncF(mf_coeftF.point_of_dof(i));
    RtF[i] = RtFuncF(mf_coeftF.point_of_dof(i));
  }
  gmm::resize(KvF, nb_dof_coefv);
  gmm::resize(MvF, nb_dof_coefv);
  gmm::resize(RvF, nb_dof_coefv);
  for (size_type i = 0; i < nb_dof_coefv; ++i) {
    KvF[i] = KvFuncF(mf_coefvF.point_of_dof(i));
    MvF[i] = MvFuncF(mf_coefvF.point_of_dof(i));
    RvF[i] = RvFuncF(mf_coefvF.point_of_dof(i));
  }

  // STAGE 2. ASSEMBLING AttF, AvvF ***************************************************

  cout << "Assembling stage" << endl;
  double time; 
  std::ofstream timetxtAss( "./Vtk/Secomb93b/timeAsslc025.txt" );

  // Build unit lumped mass matrices
  cout << "Assembling mass matrices for fluid problem..." << endl;
  getfem::asm_mass_matrix(M_tissueF, mimtF, mf_utF, mf_utF);
  //cout << "mass lumping..." << endl;
  //getfem::masslumping(M_tissueF);
  cout << "Assembling mass matrices for fluid problem vessels..." << endl;
  getfem::asm_mass_matrix(M_vesselF, mimvF, mf_uvF, mf_uvF);
  //cout << "mass lumping..." << endl;
  //getfem::masslumping(M_vesselF);

  cout << "...[OK]" << endl;

  
  // Build AttF as a stiffness matrix:
  time = gmm::uclock_sec();
  // Step 1tt: add a Kt stiffness matrix
  getfem::asm_stiffness_matrix_for_laplacian(AttF, mimtF, mf_utF, mf_coeftF, KtF);
  // Step 2tt: add a Mt mass matrix 
  getfem::asm_mass_matrix_param(AttF, mimtF, mf_utF, mf_coeftF, MtF);
  // Step 3tt: add time dependent mass term
  gmm::add( gmm::scaled(M_tissueF, CompltF/DT), AttF );
  //
  timetxtAss << "AttF assembled in  = " << gmm::uclock_sec() - time << " seconds" << endl;

  time = gmm::uclock_sec(); 

  // Build AvvF as a stiffness matrix:
  time = gmm::uclock_sec();
  // Step 1vv: add a Kv stiffness matrix
  getfem::asm_stiffness_matrix_for_laplacian(AvvF, mimvF, mf_uvF, mf_coefvF, KvF);
  // Step 2vv: add a Mv mass matrix 
  getfem::asm_mass_matrix_param(AvvF, mimvF, mf_uvF, mf_coefvF, MvF);
  // Step 3vv: add time dependent mass term
  gmm::add( gmm::scaled(M_vesselF, ComplvF/DT), AvvF );
  //

  timetxtAss << "AvvF assembled in  = " << gmm::uclock_sec() - time << " seconds" << endl;

  time = gmm::uclock_sec(); 

  // STAGE 4. Boundary terms ================================================

  // Fluid Tissue 
  getfem::build_boundary_terms(AttF, mf_utF, mf_coeftF,
			       mimtF,
			       KtF, BtF,
			       GAMMA, 
			       BCListF);

  // Fluid Vessel
  getfem::build_boundary_terms(AvvF, mf_uvF, mf_coefvF,
			       mimvF,
			       KvF, BvF,
			       GAMMA,
			       networkBCListF);


  cout << "   Computing the Preconditioners for the fluid problem" << endl;
  tPreconditionerF.build_with(AttF);
  vPreconditionerF.build_with(AvvF);

  //cout << "...done ("
  //    << gmm::uclock_sec() - time << " seconds)" << endl;
  timetxtAss << "Computing the Boundary terms and Preconditioners for the fluid problem, done in  = " << gmm::uclock_sec() - time << " seconds" << endl;
  timetxtAss.close();
 }


  //******************** Compute the velocities ************************

void myProblem::compute_velocity_Ut(void){
  getfem::compute_gradient(mf_utF, mf_gradienttF, Pt, gradPt); 
  gmm::copy(gradPt,velUt);
  gmm::scale(velUt, velT);
  gmm::copy(velUt, velUtD);
  
}

void myProblem::compute_velocity_Uv(void){
  getfem::compute_gradient(mf_uvF, mf_gradientvF, Pv, gradPv); 
  gmm::copy(gradPv,velUv);
  gmm::scale(velUv, velV);
  gmm::copy(velUv, velUvD);
}


  //******************** Assembly the RHS ************************

void myProblem::assemblyRHSFluid(void)
{

  size_type nb_dof_t = mf_utF.nb_dof();
  size_type nb_dof_v = mf_uvF.nb_dof(); 

  gmm::resize(Cvold, nb_dof_v); gmm::clear(Cvold); 
  gmm::resize(Ctold, nb_dof_t); gmm::clear(Ctold); 
 
  
  getfem::coefficient_exchange_matrices(KKttF, KKtvF, KKvtF, KKvvF, Mbar, Mlin, M0, 
				  mf_utF, mf_coeftF, mf_uvF, mf_coefvF, mimtF, mimvF,
				  Cvold, Ctold, 
				  OMEGA1, OMEGA2,
				  AconcF, BconcF,
				  alphaF,  Cm,
				  RADIUS, NInt);

  cout << "Costruisco BtF "<<endl;
  gmm::mult_add(KKtvF, PPv, BtF);
  gmm::mult_add(KKttF, PPt, BtF);

  cout << "Costruisco BvF "<<endl;
  gmm::mult_add(KKvtF, PPt, BvF);
  gmm::mult_add(KKvvF, PPv, BvF);

  cout << "Esco da assemblyRHSFluid "<<endl;

}


 //********************************** Time advance for the fluid problem **************************


bool myProblem::time_advanceFluid(const scalar_type& time, const scalar_type& dt) {


  double timeS = gmm::uclock_sec();
  std::ofstream timetxtSI( "./Vtk/Secomb93b/timeSIlc025.txt" );

  size_type nb_dof_t = mf_utF.nb_dof();
  size_type nb_dof_coeft = mf_coeftF.nb_dof();
  size_type nb_dof_v = mf_uvF.nb_dof();


  cout << "TIME t = " << time
       << "      ======================================================"
       << endl;

  //std::vector<scalar_type> Pt_old(nb_dof_t, 0.0);
  std::vector<scalar_type> BBt(nb_dof_t, 0.0);
  std::vector<scalar_type> dPt(nb_dof_t, 0.0);

  //std::vector<scalar_type> Pv_old(nb_dof_v, 0.0); 
  std::vector<scalar_type> BBv(nb_dof_v, 0.0);
  std::vector<scalar_type> dPv(nb_dof_v, 0.0);
      
  // Iterative method tissue-vessel

  scalar_type error = 1.0, tolerance = 1.0e-16;
  size_type it_num = 0;

  cout << "Copio la p vecchia" << endl;

  gmm::copy(Pt, Ptold);
  gmm::copy(Pv, Pvold);


  // Exchange terms


  cout << " entro in build exchange 3 " << endl;

  getfem::coefficient_exchange_matrices(MttF, MtvF, MvtF, MvvF, Mbar, Mlin, M0, 
				  mf_utF, mf_coeftF, mf_uvF, mf_coefvF,
				  mimtF, mimvF,
		 	  	  Pvold, Ptold, 
				  OMEGA1, OMEGA2,
				  ApresF, BpresF,
				  alphaF, Cm,
				  RADIUS, NInt);


  timetxtSI << "AB_{} matrices done in " << gmm::uclock_sec() - timeS << " seconds" << endl;
  cout << "B_{} matrices done in  " << gmm::uclock_sec() - timeS << "seconds" << endl;
 

  // Add MvvF to AvvF:
  gmm::add(MvvF,AvvF);
  
  // Add MttF to AttF
  gmm::add(MttF, AttF); 
 
  timeS = gmm::uclock_sec();

  assemblyRHSFluid();

  timetxtSI << "fluid RHS done in " << gmm::uclock_sec() - timeS << " seconds" << endl;
  cout << " fluid RHS done in " << gmm::uclock_sec() - timeS << "seconds" << endl;

  timeS = gmm::uclock_sec();

  

  while ( (error > tolerance) && (it_num < MAX_SUBITERF) ) {
    
    it_num++;
    cout << "SUB-ITERATION # " << it_num
         << " ______________________________________________________"
         << endl;

    gmm::copy(Pt, dPt);
    gmm::copy(Pv, dPv);

    gmm::copy(BtF, BBt);
    gmm::copy(BvF, BBv);
    gmm::mult_add( M_tissueF, gmm::scaled(Ptold, CompltF/DT), BBt);
    gmm::mult_add( M_vesselF, gmm::scaled(Pvold, ComplvF/DT), BBv);
    gmm::mult_add( MtvF, Pv, BBt);
    gmm::mult_add( MvtF, Pt, BBv);


    gmm::iteration iter(residu, 0, 10000);


    if (strcmp(SOLVE_METHOD_T,"CG")==0) {
      gmm::cg(AttF, Pt, BBt, tPreconditionerF, iter);
    }
    else if (strcmp(SOLVE_METHOD_T,"GMRES")==0) {
      gmm::gmres(AttF, Pt, BBt, tPreconditionerF, 50, iter);
    }
    else if (strcmp(SOLVE_METHOD_T,"BICGSTAB")==0) {
      gmm::bicgstab(AttF, Pt, BBt, tPreconditionerF, iter);
    }
    else if (strcmp(SOLVE_METHOD_T,"QMR")==0) {
      gmm::qmr(AttF, Pt, BBt, tPreconditionerF, iter);
    }
    else {
      scalar_type condest;
      gmm::SuperLU_solve(AttF, Pt, BBt, condest);
      cout << " -o 3D problem, condition number: " << 1.0/condest << endl;
    }
 
    cout << " -o 3D solve OK, #iterations = " << iter.get_iteration() << endl;

    iter.init();

    if (strcmp(SOLVE_METHOD_V,"CG")==0) {
      gmm::cg(AvvF, Pv, BBv, vPreconditionerF, iter);
    }
    else if (strcmp(SOLVE_METHOD_V,"GMRES")==0) {
      gmm::gmres(AvvF, Pv, BBv, vPreconditionerF, 50, iter);
    }
    else if (strcmp(SOLVE_METHOD_V,"BICGSTAB")==0) {
      gmm::bicgstab(AvvF, Pv, BBv, vPreconditionerF, iter);
    }
    else if (strcmp(SOLVE_METHOD_V,"QMR")==0) {
      gmm::qmr(AvvF, Pv, BBv, vPreconditionerF, iter);
    }
    else {
      scalar_type condest;
      gmm::SuperLU_solve(AvvF, Pv, BBv, condest);
      cout << " -o 1D problem, condition number: " << 1.0/condest << endl;
    }
    cout << " -o 1D solve OK, #iterations = " << iter.get_iteration() << endl;

    gmm::add(gmm::scaled(Pt, -1.0), dPt);
    gmm::add(gmm::scaled(Pv, -1.0), dPv);
    
    scalar_type tError = getfem::asm_L2_norm(mimtF, mf_utF, dPt);  
    scalar_type vError = getfem::asm_L2_norm(mimvF, mf_uvF, dPv);
    
    cout << " -o Error 3D: " << tError << endl;
    cout << " -o Error 1D: " << vError << endl;
    error = tError + vError;
  }
   

  cout << "--- TISSUE Boundary report  --------" << endl;
  integrate_boundary_quantities(mf_utF, mf_coeftF, mimtF, Pt, KtF, BCListF);

  cout << "--- NETWORK Boundary report --------" << endl;
  integrate_boundary_quantities(mf_uvF, mf_coefvF, mimvF, Pv, KvF, networkBCListF);

  cout << "--- NETWORK TO TISSUE flow rate ----" << endl;
  std::vector<scalar_type> Uphi(nb_dof_v); gmm::clear(Uphi);
  gmm::mult(MvtF, gmm::scaled(Pt, -1.0), Uphi);
  gmm::mult_add(MvvF, Pv, Uphi);
  FlowRateP = 0.0;
  for (size_type i = 0; i < Uphi.size(); ++i)
    FlowRateP += Uphi[i];
  cout << "Total Flow Rate Network -> Tissue pressure: "
       << FlowRateP << endl;
  timetxtSI << "Iterative scheme done in " << gmm::uclock_sec() - timeS << " seconds" << endl;
  timetxtSI.close();
  cout << "Iterative scheme done in " << gmm::uclock_sec() - timeS << " seconds" << endl;
  
}


 //********************************** Time advance for the fluid problem **************************


bool myProblem::time_advanceFluidMONO(void) {

  double timeS = gmm::uclock_sec();
  std::ofstream timetxtSI( "./Vtk/Secomb93b/timeSIlc025.txt" );

  size_type nb_dof_t = mf_utF.nb_dof();
  size_type nb_dof_coeft = mf_coeftF.nb_dof();
  size_type nb_dof_v = mf_uvF.nb_dof();


  std::vector<scalar_type> BBt(nb_dof_t, 0.0);
  std::vector<scalar_type> BBv(nb_dof_v, 0.0);

  // Exchange terms


  cout << " entro in build exchange 3 " << endl;

  getfem::coefficient_exchange_matrices(MttF, MtvF, MvtF, MvvF, Mbar, Mlin, M0, 
				  mf_utF, mf_coeftF, mf_uvF, mf_coefvF,
				  mimtF, mimvF,
		 	  	  Pvold, Ptold, 
				  OMEGA1, OMEGA2,
				  ApresF, BpresF,
				  alphaF, Cm,
				  RADIUS, NInt);


  // Add MvvF to AvvF:
  gmm::add(MvvF, AvvF);
  
  // Add MttF to AttF
  gmm::add(MttF, AttF); 

  timetxtSI << "BB_{} matrices done in " << gmm::uclock_sec() - timeS << " seconds" << endl;
  timeS = gmm::uclock_sec();
 
  //assemblyRHSFluid();

  gmm::copy(BtF, BBt);

  gmm::copy(BvF, BBv);

  // construction of AM 
  
  gmm::copy(AttF, gmm::sub_matrix(AM, gmm::sub_interval(0, nb_dof_t), gmm::sub_interval(0, nb_dof_t))); 

  gmm::copy(scaled(MtvF, -1), gmm::sub_matrix(AM, gmm::sub_interval(0, nb_dof_t), gmm::sub_interval(nb_dof_t, nb_dof_v))); 

  gmm::copy(scaled(MvtF,-1), gmm::sub_matrix(AM, gmm::sub_interval(nb_dof_t, nb_dof_v), gmm::sub_interval(0, nb_dof_t))); 

  gmm::copy(AvvF, gmm::sub_matrix(AM, gmm::sub_interval(nb_dof_t, nb_dof_v), gmm::sub_interval(nb_dof_t, nb_dof_v))); 


  // construction of BM

  gmm::copy(BBt, gmm::sub_vector(BM, gmm::sub_interval(0, nb_dof_t)));

  gmm::copy(BBv, gmm::sub_vector(BM, gmm::sub_interval(nb_dof_t, nb_dof_v)));

  gmm::ilu_precond<sparse_matrix_type> P(AM);

  //gmm::ilut_precond<sparse_matrix_type> P(AM, 50, 1E-9);

  gmm::iteration iter(residu, 0, 10000);

  if (strcmp(SOLVE_METHOD_V,"CG")==0) {

      gmm::cg(AM, PM, BM, P, iter);
 
      }

  else if (strcmp(SOLVE_METHOD_T,"GMRES")==0) {

  		gmm::gmres(AM, PM, BM, P, 50, iter);

                }

  gmm::copy(gmm::sub_vector(PM, gmm::sub_interval(0, nb_dof_t)), Pt);
  
  gmm::copy(gmm::sub_vector(PM, gmm::sub_interval(nb_dof_t, nb_dof_v)), Pv);

  timetxtSI << "Monolitic system solved in " << gmm::uclock_sec() - timeS << " seconds" << endl;
  timeS = gmm::uclock_sec();

  cout << "--- TISSUE Boundary report  --------" << endl;
  integrate_boundary_quantities(mf_utF, mf_coeftF, mimtF, Pt, KtF, BCListF);

  cout << "--- NETWORK Boundary report --------" << endl;
  integrate_boundary_quantities(mf_uvF, mf_coefvF, mimvF, Pv, KvF, networkBCListF);

  cout << "--- NETWORK TO TISSUE flow rate ----" << endl;
  std::vector<scalar_type> Uphi(nb_dof_v); gmm::clear(Uphi);
  gmm::mult(MvtF, gmm::scaled(Pt, -1.0), Uphi);
  gmm::mult_add(MvvF, Pv, Uphi);
  FlowRateP = 0.0;
  for (size_type i = 0; i < Uphi.size(); ++i)
    FlowRateP += Uphi[i];
  cout << "Total Flow Rate Network -> Tissue pressure: "
       << FlowRateP << endl;


  timetxtSI << "Computed flow in" << gmm::uclock_sec() - timeS << " seconds" << endl;
  timetxtSI.close();


}


/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/


int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.


 try {   
 
    myProblem p;
    p.PARAM.read_command_line(argc, argv);
    p.init(); 
    p.assemblyFluid();

    std::ofstream statistics( "./Vtk/rat98/statisticslc025.txt" ); 


    scalar_type time = 0.0;
    scalar_type dt = p.DT;
    size_type timeStep = 0;
 
    p.time_advanceFluidMONO();

    
    cout << " ... Computing the velocity ... " << endl;

    p.compute_velocity_Ut();
    p.compute_velocity_Uv();


    char str_tF[50], str_vF[50], str_gradPt[50], str_Ut[50], str_gradPv[50], str_Uv[50];

    sprintf(str_tF, "./Vtk/rat98/Pt.vtk");
    sprintf(str_vF, "./Vtk/rat98/Pv.vtk");
    sprintf(str_gradPt, "./Vtk/rat98/gradPt.vtk");
    sprintf(str_gradPv, "./Vtk/rat98/gradPv.vtk");
    sprintf(str_Ut, "./Vtk/rat98/Ut.vtk");
    sprintf(str_Uv, "./Vtk/rat98/Uv.vtk");


	{
        cout << "esporto pt" << endl;
	getfem::vtk_export vtk_tF(str_tF);
	vtk_tF.exporting(p.mf_utF);
	vtk_tF.write_mesh();
	vtk_tF.write_point_data(p.mf_utF, p.Pt, "Pt");
        cout << "esporto pv" << endl;
	getfem::vtk_export vtk_vF(str_vF);
	vtk_vF.exporting(p.mf_uvF);
	vtk_vF.write_mesh();
	vtk_vF.write_point_data(p.mf_uvF, p.Pv, "Pv");

        cout << "esporto gradPt" << endl;
	getfem::vtk_export vtk_gradPt(str_gradPt);
	vtk_gradPt.exporting(p.mf_veltF);
	vtk_gradPt.write_mesh();
	vtk_gradPt.write_point_data(p.mf_veltF, p.gradPt, "gradPt");

        cout << "esporto gradPv" << endl;
	getfem::vtk_export vtk_gradPv(str_gradPv);
	vtk_gradPv.exporting(p.mf_velvF);
	vtk_gradPv.write_mesh();
	vtk_gradPv.write_point_data(p.mf_velvF, p.gradPv, "gradPv");


        cout << "esporto velUt" << endl;
	getfem::vtk_export vtk_Ut(str_Ut);
	vtk_Ut.exporting(p.mf_veltF);
	vtk_Ut.write_mesh();
	vtk_Ut.write_point_data(p.mf_veltF, p.velUt, "Ut");

        cout << "esporto velUv" << endl;
	getfem::vtk_export vtk_Uv(str_Uv);
	vtk_Uv.exporting(p.mf_velvF);
	vtk_Uv.write_mesh();
	vtk_Uv.write_point_data(p.mf_velvF, p.velUv, "Uv");

	}

    p.Pressureaverage();

    cout << " integrale di Pt = " << p.meanPt << endl;
    cout << " integrale di Pv = " << p.meanPv << endl;


    statistics << " integrale di Pt = " << p.meanPt << endl;
    statistics << " integrale di Pv = " << p.meanPv << endl;

    statistics << "Total Flow Rate Network -> Tissue, pressure: " << p.FlowRateP << endl;

    statistics.close();

     }

  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}


