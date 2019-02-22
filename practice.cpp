#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "code/mathelpers.h"
#include "code/utils.h"

using namespace std;
using namespace Eigen;

class CSVMatrix
{
    public:
        void loadMatrix(std::istream& inStream)
        {
            string         line;
            stringstream   lineStream;
            string         cell;

            int row=0;

            //read lines
            while(getline(inStream, line) )
            {
                lineStream.clear();
                lineStream.str(line);
                cout << "row=" << row++
                    << " lineStream.str() = " << lineStream.str() << std::endl;

                //read cells
                while(getline(lineStream, cell, ','))
                {
                    cout << "cell=" << cell << std::endl;
                }
            }
        }
};

int main()
{
	string PATH_F = "data_gecko/F.csv";
    ifstream infile_count(PATH_F);
    int rows = 0;
    int cols = 0;
    countRowsColsCsv(infile_count, rows, cols);

    ifstream infile_read(PATH_F);
    MatrixXd F(rows, cols);
    loadMatrix(infile_read, F);

    
    //CSVMatrix matrix;

    //matrix.loadMatrix(infile);
    //MatrixXd F = loadMatrix(infile);

    //cout << F.rows() << "-" << F.cols() << endl;

    for(int i=0; i<F.rows(); i++)
    {
    	for(int j=0; j<F.cols(); j++)
    	{
    		cout << F(i,j) << "-";
    	}
    	cout << endl;
    }

    return 0;
}

/*int main()
{
	fstream fin;
	fin.open("data_gecko.V", ios::in);
	int rollnum, rool2, count=0;
	cout << "Enter";
	cin >> rollnum;

	vector<string> row;
	string word;
	string line, temp;

	while(fin >> temp)
	{
		row.clear();
		getline(fin, line);
		stringstream s(line);

		while(s.getline(s, word, ', '))
		{
			row.push_back(word);
		}
	}

	cout << "size: " << row.size() << endl;

	return 0;
}*/

/*int main()
{
	MatrixXd m(2,3);
	m(0,0) = 1;
	m(0,1) = 3;
	m(0,2) = 3.5;
	m(1,0) = 2;
	m(1,1) = 4;
	m(1,2) = 4.5;

	cout << m << endl;
	MatrixXd ms = colStack(m);
	cout << "colStack" << endl;
	cout << ms << endl;

	return 0;
}*/