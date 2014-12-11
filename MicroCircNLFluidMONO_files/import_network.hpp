
namespace getfem {

/*
  Create the edge sequence in mh3D and build the related 1D mesh mh1D.
  The input stream ist is used to read a file with the following format (gambit-like):
  
  BEGIN_LIST
  BEGIN_ARC
  BC KEYWA [VALUES] 
  BC KEYWB [VALUES]
     106       0.4421      -1.6311       2.5101 vertex.78
     107       0.4421      -1.6311       7.5101 vertex.90
     108       0.3546      -1.6524       2.5539 edge.30
     109       0.2621      -1.6695       2.5880 edge.30
  ...
  END_ARC

  BEGIN_ARC 
  ...
  END_ARC
  END_LIST

  IMPORTANT: 
  1. The list of points of the IS ordered as follows:
     - first we have the coordinates of TWO ENDS (A,B) (i.e. A=vertex.78 and B=vertex.90)
     - then we have all the remaining nodes of the arc, from A to B
  2. If a node is shared by more arcs, the arcs will be CONNECTED in the resulting 1D mesh.
  3. BC KEYWA [VALUES] and BC KEYWB [VALUES] are keywords/values related to boundary conditions
     to be imposed at nodes A, B. Each KEYW [VALUES] entry can be one of the following:

     BC DIR 1.1
     BC NEU 12.3
     BC MIX 1.2 3.4
     BC INT

     Corrrespondingly, each node will be marked with the associated boundary condition,
     that are:

       Dirichlet node (pressure = 1.1)
       Neumann   node (flux = 12.3)
       Mixed     node (flux + 1.2 * (pressure - 3.4) = 0)
       Internal  node

     At this stage, this is only meant to assign such BC labels to each node.
     If one end is INTERNAL, the corresponding BC will be ignored 
     (for clarity, please write the INTERNAL keyword).
*/

void import_pts_file(std::istream &ist, 
		     getfem:: mesh &mh1D, 
		     const std::string &MESH_TYPE,
		     std::vector< BCData > &BCList
		     ) {
  
  // IMPORTING DATA FILE --------------------------

  ist.precision(16);
  ist.seekg(0); ist.clear();
  if (!bgeot::read_until(ist, "BEGIN_LIST"))
    DAL_THROW(dal::failure_error, "This seems not to be a data file");

  size_type globalBoundaries = 0;
  
  while (bgeot::read_until(ist, "BEGIN_ARC")) {
    
    std::vector<base_node> lpoints; 
    
    dal::dynamic_array<scalar_type> tmpv;
    std::string tmp, BCtype, Val1, Val2;;
    //char  tmp[1024], BCtype[4], Val1[1024], Val2[1024];   
    bool thend = false; 
    size_type bcflag = 0;
    size_type bcintI=0, bcintF = 0;
    BCData BCA, BCB;

    while (!thend) {
      bgeot::get_token(ist, tmp, 1023);
      if (strcmp(tmp.c_str(), "END_ARC") == 0) { 
	thend = true;
      }
      else if (ist.eof()) {
	DAL_THROW(dal::failure_error, "Unexpected end of stream");
      }
      else if (strcmp(tmp.c_str(), "BC") == 0) { 
	bcflag++;
	bgeot::get_token(ist, BCtype, 4);
	if (strcmp(BCtype.c_str(), "DIR") == 0) {
	  bgeot::get_token(ist, Val1, 1023);
	  if (bcflag == 1) {
	    strcpy(BCA.label, BCtype.c_str()); BCA.Val1 = atof(Val1.c_str()); BCA.Val2 = 0.0;
	    BCA.ind = globalBoundaries;
	    globalBoundaries++;
	  }
	  else if (bcflag == 2) {
	    strcpy(BCB.label, BCtype.c_str()); BCB.Val1 = atof(Val1.c_str()); BCB.Val2 = 0.0;	    
	    BCB.ind = globalBoundaries;
	    globalBoundaries++;
	  }
	  else
	    DAL_THROW(dal::failure_error, "More than 2 BC found on this arc!");
	}
	else if (strcmp(BCtype.c_str(), "NEU") == 0) {
	  bgeot::get_token(ist, Val1, 1023);
	  if (bcflag == 1) {
	    strcpy(BCA.label, BCtype.c_str()); BCA.Val1 = atof(Val1.c_str()); BCA.Val2 = 0.0;	    
	    BCA.ind = globalBoundaries;
	    globalBoundaries++;
	  }
	  else if (bcflag == 2) {
	    strcpy(BCB.label, BCtype.c_str()); BCB.Val1 = atof(Val1.c_str()); BCB.Val2 = 0.0;	    	    
	    BCB.ind = globalBoundaries;
	    globalBoundaries++;
	  }
	}
	else if (strcmp(BCtype.c_str(), "MIX") == 0) {
	  bgeot::get_token(ist, Val1, 1023);
	  bgeot::get_token(ist, Val2, 1023);
	  if (bcflag == 1) {
	    strcpy(BCA.label, BCtype.c_str()); BCA.Val1 = atof(Val1.c_str()); BCA.Val2 = atof(Val2.c_str());
	    BCA.ind = globalBoundaries;
	    globalBoundaries++;
	  }
	  else if (bcflag == 2) {
	    strcpy(BCB.label, BCtype.c_str()); BCB.Val1 = atof(Val1.c_str()); BCB.Val2 = atof(Val2.c_str());
	    BCB.ind = globalBoundaries;
	    globalBoundaries++;
	  }
	  else
	    DAL_THROW(dal::failure_error, "More than 2 BC found on this arc!");	  
	}
	else if (strcmp(BCtype.c_str(), "INT") == 0) {
	  if (bcflag == 1) {
            bcintI++;
	    strcpy(BCA.label, BCtype.c_str()); BCA.Val1 = atof(Val1.c_str()); BCA.Val2 = 0.0;
	  }
	  else if (bcflag == 2) {
            bcintF++;
	    strcpy(BCB.label, BCtype.c_str()); BCB.Val1 = atof(Val1.c_str()); BCB.Val2 = 0.0;
	  }
	  else
	    DAL_THROW(dal::failure_error, "More than 2 BC found on this arc!");
	}
	else
	  DAL_THROW(dal::failure_error, "Unknown Boundary condition");	  
      }
      else if (strlen(tmp.c_str()) == 0) {
	DAL_THROW(dal::failure_error, "Syntax error in file, at token '" << tmp
		  << "', pos=" << std::streamoff(ist.tellg()));
      } 
      else {
	int d = 0;
	while ( (isdigit(tmp[0]) != 0) || tmp[0] == '-' || tmp[0] == '+'
		|| tmp[0] == '.')
	  { tmpv[d++] = atof(tmp.c_str()); bgeot::get_token(ist, tmp, 1023); }
	if (d != 4)
	  DAL_THROW(dal::failure_error, "Points must have 3 coordinates");
	base_node tmpn(tmpv[1], tmpv[2], tmpv[3]);
	lpoints.push_back(tmpn);
	if (strcmp(tmp.c_str(), "END_ARC") == 0) thend = true;
      } 
    };
    
    cout << " Found an arc with " << lpoints.size()
	 << " points." << endl;

    // Insert the arc into the 1D mesh
    for (size_type i=0; i<lpoints.size()-1; ++i ) {
      std::vector<size_type> ind(2);
      size_type ii = (i==0) ? 0 : i+1;
      size_type jj;
      if (ii == lpoints.size()-1) jj = 1;
      else if (ii == 0) jj = 2;
      else jj = ii+1;
      ind[0] = mh1D.add_point(lpoints[ii]);
      ind[1] = mh1D.add_point(lpoints[jj]);
      //std::string MESH_TYPE1 = "GT_PK(1," + poly_deg + ")";
      size_type cv;
      cv = mh1D.add_convex(bgeot::geometric_trans_descriptor(MESH_TYPE), 
			   ind.begin());
      if ((bcflag>0) && (ii==0)&& (bcintI==0)) {
	BCA.cv = cv;
	BCList.push_back(BCA);
	cout << " - Boundary node: " << BCA.label << ", assigned flag " << BCA.ind << endl;
      }
      if ((bcflag>1) && (jj==1)&& (bcintF==0)) {
	BCB.cv = cv;
	BCList.push_back(BCB);
	cout << " - Boundary node: " << BCB.label << ", assigned flag " << BCB.ind << endl;
      }

    };
    
  };

  // SAVE the 1D MESH
  //mh1D.write_to_file("mesh_1D.mesh");
  //cout << "saved mesh in mesh_1D.mesh" << endl; 
  
}
  

// Get the boudary nodes (modified from getfem_mesh::outer_faces_of_mesh) 
void extrema_of_oneDmesh(const getfem::mesh &m,
                         const dal::bit_vector& cvlst,
                         getfem::convex_face_ct& flist) {
  for (dal::bv_visitor ic(cvlst); !ic.finished(); ++ic) {
    for (size_type f = 0; f < m.structure_of_convex(ic)->nb_faces(); f++) {
      if (m.neighbour_of_convex(ic,f)==-1) {
        flist.push_back(getfem::convex_face(ic,f));
      }
    }
  }
}

}
