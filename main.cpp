#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include "code/mathelpers.h"
#include "class/Param_State.h"
#include "code/utils.h"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "libs/writeOBJ.h"
#include "code/OptimProblemIsoDist.h"

using namespace std;
using namespace cv;
using namespace Eigen;

int main()
{
	//Example of gecko deformation
	int num_iter = 2500;
	double TolX = 1e-10;
	double TolFun = 1e-6;

	Param_State mesh;

//Seccion of read a mesh
    read_mesh_2D("data_gecko/V.csv", "data_gecko/F.csv","data_gecko/eq_lhs.csv", "data_gecko/eq_rhs.csv", mesh);

    /*cout << "F: " << mesh.F.rows() << " - " << mesh.F.cols() << endl;
    cout << "V: " << mesh.V.rows() << " - " << mesh.V.cols() << endl;
    cout << "eq_lhs: " << mesh.eq_lhs.rows() << " - " << mesh.eq_lhs.cols() << endl;
    cout << "eq_rhs: " << mesh.eq_rhs.rows() << " - " << mesh.eq_rhs.cols() << endl;*/

    update_F(mesh.F);

    MatrixXd V0; 
    OptimProblemIsoDist optimProblem(mesh, V0, 25);

    cout << optimProblem.T.rows() << " - "<< optimProblem.T.cols()<<endl;
    cout << optimProblem.areas.rows()<< " - "<< optimProblem.areas.cols() << endl;

    
   



    /*MatrixXd nV(mesh.V.rows(), mesh.V.cols()+1);
    create_column_zeros(mesh.V, nV);
    mesh.V = nV;

    igl::writeOBJ("prueba.obj", mesh.V, mesh.F);*/
   
	return 0;
}