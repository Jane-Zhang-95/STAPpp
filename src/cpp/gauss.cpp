
#include "gauss.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//CT3

//!	Constructor
CT3Gauss::CT3Gauss(){}
    
//!	Desconstructor
CT3Gauss::~CT3Gauss(){
        if (w)
            delete [] w;
        if (csi)
            delete [] csi;
}

//	Read ngp data from stream Input
bool CT3Gauss::Read(ifstream& Input)
{
	Input >> ngp;	// Number of property set

	return true;
}

//	Write ngp data to Stream 输出ngp
void CT3Gauss::Write(COutputter& output)
{
	output << setw(16) << ngp << endl;
}

void CT3Gauss::T3_GPoint()
{
    csi,w = new double[3];
    switch(ngp)
    {
        case 1:
            csi[0] = 0.333333333;
			csi[1] = 0.333333333;
			csi[2] = 0.333333333;
            w[0] = 0.5;
			w[1] = 0.5;
			w[2] = 0.5;
            break;
        case 2:
            csi[0] = 0.666666667;
			csi[1] = 0.166666667;
			csi[2] = 0.166666667;
            w[0] = 0.166666667;
			w[1] = 0.166666667;
			w[2] = 0.166666667;
            break;
        case 3:
            csi[0] = 0.5;
			csi[1] = 0.5;
			csi[2] = 0.0;
            w[0] = 0.166666667;
			w[1] = 0.166666667;
			w[2] = 0.166666667;
            break;
        default:
            std::cerr << "ngp " << ngp << " not available.." << std::endl;
            exit(5);
    }
}