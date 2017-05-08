#include <curved_network_reader.hpp>
#include <c_problem3d1d.hpp>



int main(int argc, char* argv[]){
 //Creating the input mesh file
 NetDiff::curved_reader(argc, argv);
 //Creating, solving and exporting the problem results
 
 getfem::c_problem3d1d prob;
 prob.init(argc, argv);
 prob.assembly();
 prob.solve();
 prob.export_vtk();
 
 return 0;

}