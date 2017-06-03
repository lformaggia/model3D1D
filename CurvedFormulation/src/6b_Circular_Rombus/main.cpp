#include <c_problem3d1d.hpp>

using namespace getfem;


int main(int argc, char * argv[]){
 	c_problem3d1d prob;
 	prob.init(argc, argv);
 	prob.assembly();
	prob.solve();
	prob.export_vtk();


	std::cout << "--- FINAL RESULTS -------------------------" << std::endl; 
	std::cout << "  Pt average            = " << prob.mean_pt()   << std::endl;
	std::cout << "  Pv average            = " << prob.mean_pv()   << std::endl;
	std::cout << "  Uv average            = " << prob.mean_uv()   << std::endl;
	std::cout << "  Network-to-Tissue TFR = " << prob.flow_rate() << std::endl;
	std::cout << "  Dynamic Pressure Gain = " << prob.G() << std::endl;
	std::cout << "-------------------------------------------" << std::endl;

}