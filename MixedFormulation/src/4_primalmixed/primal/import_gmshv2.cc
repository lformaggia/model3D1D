
#include <iostream>
#include <iomanip>
#include <fstream>


namespace getfem {

  //   Import TETRAHEDRAL mesh file from GMSH via the .msh format.

  static void import_mshv2_file(std::string filename, mesh& m) {

    m.clear();
    std::ifstream f(filename.c_str());
    if (!f.good()) 
      DAL_THROW(failure_error, "Can't open file " << filename);

    char tmp[1024];

    bgeot::read_until(f, "$Nodes");

    size_type nb_node; 

    f >> nb_node;

    dal::dynamic_tree_sorted<size_type> msh_node_2_getfem_node;

    for (size_type node_cnt=0; node_cnt < nb_node; ++node_cnt) {
      size_type node_id;
      base_node n(3); n[0]=0.0; n[1]=0.1; n[2]=1e30;
      f >> node_id >> n[0] >> n[1] >> n[2];
      msh_node_2_getfem_node.add_to_index(m.add_point(n), node_id);
    }

    bgeot::read_until(f, "$EndNodes");
    bgeot::read_until(f, "$Elements");

    size_type nb_cv;
    f >> nb_cv;

    std::vector<size_type> cv_nodes;

    for (size_type cv=0; cv < nb_cv; ++cv) {
      size_type cv_id, cv_type, cv_region_id, cv_elm_id, cv_partition, cv_nb_nodes ;
      f >> cv_id >> cv_type >> cv_region_id >> cv_elm_id >> cv_partition;

      if (cv_type == size_type(4))
      cv_nb_nodes = size_type(4);
      else
      DAL_THROW(failure_error, "No thetraedra mesh");
      cv_id--; 
      cv_nodes.resize(cv_nb_nodes);
      for (size_type i=0; i < cv_nb_nodes; ++i) {
	size_type j;
	f >> j;
	cv_nodes[i] = msh_node_2_getfem_node.search(j);
	if (cv_nodes[i] == size_type(-1)) 
	  DAL_THROW(failure_error, "Invalid node ID " << j 
		    << " in gmsh convex " << cv_id);
      }

	m.add_simplex(3,cv_nodes.begin());

   }

   bgeot::read_until(f, "$EndElements");
   f.close();


}
}
