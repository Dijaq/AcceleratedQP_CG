#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "../class/Param_State.h"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"

using namespace std;
using namespace cv;
using namespace Eigen;

void loadMatrix(std::istream& inStream, MatrixXd &matrix)
{
    string         line;
    stringstream   lineStream;
    string         value;

    int wrows = 0;
    int wcols = 0;
    while(getline(inStream, line) )
    {
        lineStream.clear();
        lineStream.str(line);
        //read cells
        //cout << "row=" << wrows << " lineStream.str() = " << lineStream.str() << std::endl;
        while(getline(lineStream, value, ','))
        {
            matrix(wrows, wcols) = stof(value);
            wcols++;
            //cout << "cell=" << value << std::endl;
        }
        wrows++;
        wcols = 0;
    }
}

void load_SparseMatrix(std::istream& inStream, SparseMatrix<double> &matrix)
{
    string         line;
    stringstream   lineStream;
    string         value;

    int wrows = 0;
    int wcols = 0;
    while(getline(inStream, line) )
    {
        lineStream.clear();
        lineStream.str(line);
        //read cells
        //cout << "row=" << wrows << " lineStream.str() = " << lineStream.str() << std::endl;
        while(getline(lineStream, value, ','))
        {
            if(stof(value) != 0)
                matrix.insert(wrows, wcols) = stof(value);
            wcols++;
            //cout << "cell=" << value << std::endl;
        }
        wrows++;
        wcols = 0;
    }
}

void export_mat_to_excel(MatrixXd matrix, string name)
{
    ofstream myfile;
    myfile.open(name+".csv");
    for(int i=0; i<matrix.rows(); i++)
    {
        for(int j=0; j<matrix.cols(); j++)
        {
            myfile << matrix(i,j);
            myfile << ",";
        }
        myfile << "\n";
    }
    myfile.close();
}


void export_sparsemat_to_excel(SparseMatrix<double> matrix, string name)
{
    ofstream myfile;
    myfile.open(name+".csv");
    
    for(int i=0; i<matrix.outerSize(); i++)
    {
        for(int j=0; j< matrix.rows(); j++)
        {
            for(SparseMatrix<double>::InnerIterator it(matrix, i); it; ++it)
            {
                if(it.row() == j)
                    myfile << it.col() << "," << it.row() << "," << it.value() << "\n";
                /*cout << "value: " << it.value() << endl;
                cout << "row: "<< it.row() << endl;
                cout << "col: "<< it.col() << endl;
                cout << it.index() << endl;*/
            }
        }
    }


    myfile.close();
}


void countRowsColsCsv(istream& inStream, int &rows, int &cols)
{
    string         line;
    stringstream   lineStream;
    string         value;
    //read lines
    while(getline(inStream, line) )
    {
        lineStream.clear();
        lineStream.str(line);
        //cout << "row=" << rows << " lineStream.str() = " << lineStream.str() << std::endl;
        //read cells
        cols=0;
        rows++;
        while(getline(lineStream, value, ','))
        {
            cols++;
            //cout << "cell=" << value << std::endl;
        }
    }

}

void read_mesh_2D(string PATH_V, string PATH_F, string PATH_eq_lhs, string PATH_eq_rhs, Param_State &mesh)
{
    //read faces
    ifstream infile_F(PATH_F);
    int rows = 0;
    int cols = 0;
    countRowsColsCsv(infile_F, rows, cols);
    infile_F.close();

    ifstream infile_F_read(PATH_F);
    MatrixXd F(rows, cols);
    loadMatrix(infile_F_read, F);
    infile_F_read.close();

//read vertices
    ifstream infile_V(PATH_V);
    rows = 0;
    cols = 0;
    countRowsColsCsv(infile_V, rows, cols);
    infile_V.close();

    ifstream infile_V_read(PATH_V);
    MatrixXd V(rows, cols);
    loadMatrix(infile_V_read, V);
    infile_V_read.close();

//read eq_lhs

    ifstream infile_eq_lhs(PATH_eq_lhs);
    rows = 0;
    cols = 0;
    countRowsColsCsv(infile_eq_lhs, rows, cols);
    infile_eq_lhs.close();

    ifstream infile_eq_lhs_read(PATH_eq_lhs);
    SparseMatrix<double> eq_lhs(rows, cols);
    load_SparseMatrix(infile_eq_lhs_read, eq_lhs);
    infile_eq_lhs_read.close();

//read eq_rhs

    ifstream infile_eq_rhs(PATH_eq_rhs);
    rows = 0;
    cols = 0;
    countRowsColsCsv(infile_eq_rhs, rows, cols);
    infile_eq_rhs.close();

    ifstream infile_eq_rhs_read(PATH_eq_rhs);
    MatrixXd eq_rhs(rows, cols);
    loadMatrix(infile_eq_rhs_read, eq_rhs);
    infile_eq_rhs_read.close();

    mesh.F = F;
    mesh.V = V;
    mesh.eq_lhs = eq_lhs;
    mesh.eq_rhs = eq_rhs;

}

void MatrixXd_to_Mat(MatrixXd matrix, Mat dst)
{
    for(int i=0; i<matrix.rows(); i++)
    {
        for(int j=0; j<matrix.cols(); j++)
        {
            dst.at<float>(i,j) = matrix(i,j);
        }
    }
}

void create_column_zeros(MatrixXd V, MatrixXd &new_V)
{
    for(int i=0; i<V.rows(); i++)
    {
        for(int j=0; j<V.cols()+1; j++)
        {
            if(j != 2)
                new_V(i,j) = V(i,j);
            else
                new_V(i,j) = 0;
        }
    }

}

//This is because the relation between F with V is since the number 1 to the final of the list
//But in MatrixXd "V" is saved since the coordinate 0
void update_F(MatrixXd &F)
{
    //cout << "Filas: " << F.rows() << " columnas: " << F.cols() << endl;
    for(int i=0; i<F.rows(); i++)
    {
        for(int j=0; j<F.cols(); j++)
        {
           // cout << "i: " << F(i,j) << endl;
            F(i,j) = F(i,j)-1;
        }
    }
}

//Matrix should be a stach for for easy reshape
void matrix_reshape(MatrixXd &matrix, int rows, int cols)
{
    MatrixXd reshapeM(rows, cols);
    for(int j=0; j<cols; j++)
    {
        for(int i=0; i<rows; i++)
        {
            reshapeM(i,j) = matrix(i+j*rows,0);
        }
    }

    matrix = reshapeM;

}

void print_dimensions(string etiqueta, MatrixXd x)
{
    cout << etiqueta <<" rows: " << x.rows() << ", cols: " << x.cols() << endl;
}

#endif