#ifndef HH_CURV_NET_READ_HH
#define HH_CURV_NET_READ_HH

#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <boost/graph/adjacency_list.hpp>

#include "netdiff_graph_properties.hpp"
#include "reader_netdiff.hpp"
#include "spline_reader.hpp"
#include "writer_pts.hpp"
#include "graph_builder.hpp"
#include <getfem/getfem_export.h>   //For ftool::md_param


using namespace boost;
using namespace BGLgeom;
using namespace std;


namespace NetDiff{

void curved_reader(int argc, char *argv[]);


}
#endif