#include "H8.h"
#include "Gauss.h"


#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//  Constructor
CH8::CH8()
{
	NEN_ = 8; 
	nodes_ = new CNode * [NEN_];

	ND_ = 24; 
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//  Desconstructor
CH8::~CH8()
{
}

//  Read element data from stream Input
bool CH8::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int Mset;
	unsigned int N1, N2, N3, N4, N5, N6, N7, N8;

	Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> Mset;
	ElementMaterial_ = dynamic_cast<C3DMaterial*>(MaterialSets) + Mset - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
    	nodes_[4] = &NodeList[N5 - 1];
	nodes_[5] = &NodeList[N6 - 1];
	nodes_[6] = &NodeList[N7 - 1];
	nodes_[7] = &NodeList[N8 - 1];

	return true;
}

//  Write element data to stream
void CH8::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(9) << nodes_[2]->NodeNumber
		   << setw(9) << nodes_[3]->NodeNumber
           << setw(9) << nodes_[4]->NodeNumber
		   << setw(9) << nodes_[5]->NodeNumber
		   << setw(9) << nodes_[6]->NodeNumber
           << setw(9) << nodes_[7]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
    
}

//  Calculate element stiffness matrix
//  Upper triangular matrix stored as an array column by colum starting from the diagonal element
void CH8::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	C3DMaterial* material_ = dynamic_cast<C3DMaterial*>(ElementMaterial_);

	double E = material_->E;
	double mu = material_->mu;
    	double lambda = E * mu/((1 + mu) * (1 - 2 * mu));
    	double miu = E/(2*(1 + mu));
	
    	double D[6][6] = {0};
    	D[0][0] = lambda + 2*miu;
    	D[1][1] = lambda + 2*miu;
    	D[2][2] = lambda + 2*miu;
    	D[3][3] = miu;
    	D[4][4] = miu;
    	D[5][5] = miu;
    	D[0][1] = lambda;
    	D[0][2] = lambda;
    	D[1][0] = lambda;
    	D[1][2] = lambda;
    	D[2][0] = lambda;
    	D[2][1] = lambda;

	double K[24][24] = {0};

    	double w[2] = {0};
    	double gp[2] = {0};
   	GaussianQuadrature H8gauss;
   	unsigned int ngp = H8gauss.ngp;
    	H8gauss.Gauss_Calculate();
    	for(unsigned int i = 0; i<ngp; i++)
    	{
        	w[i] = H8gauss.g_w[i];
        	gp[i] = H8gauss.g_p[i];
    	}

	for (unsigned int i = 0; i < 2; i++) {
		for (unsigned int j = 0; j < 2; j++) {
            		for (unsigned int k = 0; k < 2; k++) {
		                // Gaussian quadrature weights and points
		                double w_i = w[i];
		                double w_j = w[j];
		                double w_k = w[k];
		                double p_i = gp[i];
		                double p_j = gp[j];
		                double p_k = gp[k];
		
		                double n[8][3] = {0};
		
		                n[0][0] = -1.0/8 * (1-p_j)*(1-p_k);
		                n[0][1] = -1.0/8 * (1-p_i)*(1-p_k);
		                n[0][2] = -1.0/8 * (1-p_i)*(1-p_j);
		
		                n[1][0] = 1.0/8 * (1-p_j)*(1-p_k);
		                n[1][1] = -1.0/8 * (1+p_i)*(1-p_k);
		                n[1][2] = -1.0/8 * (1+p_i)*(1-p_j);
		
		                n[2][0] = 1.0/8 * (1+p_j)*(1-p_k);
		                n[2][1] = 1.0/8 * (1+p_i)*(1-p_k);
		                n[2][2] = -1.0/8 * (1+p_i)*(1+p_j);
		
		                n[3][0] = -1.0/8 * (1+p_j)*(1-p_k);
		                n[3][1] = 1.0/8 * (1-p_i)*(1-p_k);
		                n[3][2] = -1.0/8 * (1-p_i)*(1+p_j);
		
		                n[4][0] = -1.0/8 * (1-p_j)*(1+p_k);
		                n[4][1] = -1.0/8 * (1-p_i)*(1+p_k);
		                n[4][2] = 1.0/8 * (1-p_i)*(1-p_j);
		
		                n[5][0] = 1.0/8 * (1-p_j)*(1+p_k);
		                n[5][1] = -1.0/8 * (1+p_i)*(1+p_k);
		                n[5][2] = 1.0/8 * (1+p_i)*(1-p_j);
		
		                n[6][0] = 1.0/8 * (1+p_j)*(1+p_k);
		                n[6][1] = 1.0/8 * (1+p_i)*(1+p_k);
		                n[6][2] = 1.0/8 * (1+p_i)*(1+p_j);
		
		                n[7][0] = -1.0/8 * (1-p_j)*(1-p_k);
		                n[7][1] = 1.0/8 * (1-p_i)*(1+p_k);
		                n[7][2] = 1.0/8 * (1-p_i)*(1+p_j);
		
		                // Jacobian matrix and det
		                double J_11 = 1.0/8 * ((p_j-1)*(1-p_k)*(nodes_[0]->XYZ[0]) + 
		                                    (1-p_j)*(1-p_k)*(nodes_[1]->XYZ[0]) + 
		                                    (1+p_j)*(1-p_k)*(nodes_[2]->XYZ[0]) + 
		                                    (-p_j-1)*(1-p_k)*(nodes_[3]->XYZ[0]) + 
		                                    (p_j-1)*(1+p_k)*(nodes_[4]->XYZ[0]) + 
		                                    (1-p_j)*(1+p_k)*(nodes_[5]->XYZ[0]) + 
		                                    (1+p_j)*(1+p_k)*(nodes_[6]->XYZ[0]) + 
		                                    (-1-p_j)*(1+p_k)*(nodes_[7]->XYZ[0]));
		                // printf("J_11:%f\n",J_11);
		                double J_12 = 1.0/8 * ((p_i-1)*(1-p_k)*(nodes_[0]->XYZ[0]) + 
		                                    (1+p_i)*(p_k-1)*(nodes_[1]->XYZ[0]) + 
		                                    (1+p_i)*(1-p_k)*(nodes_[2]->XYZ[0]) + 
		                                    (1-p_i)*(1-p_k)*(nodes_[3]->XYZ[0]) + 
		                                    (p_i-1)*(1+p_k)*(nodes_[4]->XYZ[0]) + 
		                                    (-1-p_i)*(1+p_k)*(nodes_[5]->XYZ[0]) + 
		                                    (1+p_i)*(1+p_k)*(nodes_[6]->XYZ[0]) + 
		                                    (1-p_i)*(1+p_k)*(nodes_[7]->XYZ[0]));
		                // printf("J_12:%f\n",J_12);
		                double J_13 = 1.0/8 * ((p_i-1)*(1-p_j)*(nodes_[0]->XYZ[0]) + 
		                                    (1+p_i)*(p_j-1)*(nodes_[1]->XYZ[0]) + 
		                                    (-1-p_i)*(1+p_j)*(nodes_[2]->XYZ[0]) + 
		                                    (p_i-1)*(1+p_j)*(nodes_[3]->XYZ[0]) + 
		                                    (1-p_i)*(1-p_j)*(nodes_[4]->XYZ[0]) + 
		                                    (1+p_i)*(1-p_j)*(nodes_[5]->XYZ[0]) + 
		                                    (1+p_i)*(1+p_j)*(nodes_[6]->XYZ[0]) + 
		                                    (1-p_i)*(1+p_j)*(nodes_[7]->XYZ[0]));
		                // printf("J_13:%f\n",J_13);
		                double J_21 = 1.0/8 * ((p_j-1)*(1-p_k)*(nodes_[0]->XYZ[1]) + 
		                                    (1-p_j)*(1-p_k)*(nodes_[1]->XYZ[1]) + 
		                                    (1+p_j)*(1-p_k)*(nodes_[2]->XYZ[1]) + 
		                                    (-p_j-1)*(1-p_k)*(nodes_[3]->XYZ[1]) + 
		                                    (p_j-1)*(1+p_k)*(nodes_[4]->XYZ[1]) + 
		                                    (1-p_j)*(1+p_k)*(nodes_[5]->XYZ[1]) + 
		                                    (1+p_j)*(1+p_k)*(nodes_[6]->XYZ[1]) + 
		                                    (-1-p_j)*(1+p_k)*(nodes_[7]->XYZ[1]));
		                // printf("J_21:%f\n",J_21);
		                double J_22 = 1.0/8 * ((p_i-1)*(1-p_k)*(nodes_[0]->XYZ[1]) + 
		                                    (1+p_i)*(p_k-1)*(nodes_[1]->XYZ[1]) + 
		                                    (1+p_i)*(1-p_k)*(nodes_[2]->XYZ[1]) + 
		                                    (1-p_i)*(1-p_k)*(nodes_[3]->XYZ[1]) + 
		                                    (p_i-1)*(1+p_k)*(nodes_[4]->XYZ[1]) + 
		                                    (-1-p_i)*(1+p_k)*(nodes_[5]->XYZ[1]) + 
		                                    (1+p_i)*(1+p_k)*(nodes_[6]->XYZ[1]) + 
		                                    (1-p_i)*(1+p_k)*(nodes_[7]->XYZ[1]));
		                // printf("J_22:%f\n",J_22);
		                double J_23 = 1.0/8 * ((p_i-1)*(1-p_j)*(nodes_[0]->XYZ[1]) + 
		                                    (1+p_i)*(p_j-1)*(nodes_[1]->XYZ[1]) + 
		                                    (-1-p_i)*(1+p_j)*(nodes_[2]->XYZ[1]) + 
		                                    (p_i-1)*(1+p_j)*(nodes_[3]->XYZ[1]) + 
		                                    (1-p_i)*(1-p_j)*(nodes_[4]->XYZ[1]) + 
		                                    (1+p_i)*(1-p_j)*(nodes_[5]->XYZ[1]) + 
		                                    (1+p_i)*(1+p_j)*(nodes_[6]->XYZ[1]) + 
		                                    (1-p_i)*(1+p_j)*(nodes_[7]->XYZ[1]));
		                // printf("J_23:%f\n",J_23);
		                double J_31 = 1.0/8 * ((p_j-1)*(1-p_k)*(nodes_[0]->XYZ[2]) + 
		                                    (1-p_j)*(1-p_k)*(nodes_[1]->XYZ[2]) + 
		                                    (1+p_j)*(1-p_k)*(nodes_[2]->XYZ[2]) + 
		                                    (-p_j-1)*(1-p_k)*(nodes_[3]->XYZ[2]) + 
		                                    (p_j-1)*(1+p_k)*(nodes_[4]->XYZ[2]) + 
		                                    (1-p_j)*(1+p_k)*(nodes_[5]->XYZ[2]) + 
		                                    (1+p_j)*(1+p_k)*(nodes_[6]->XYZ[2]) + 
		                                    (-1-p_j)*(1+p_k)*(nodes_[7]->XYZ[2]));
		                // printf("J_31:%f\n",J_31);
		                double J_32 = 1.0/8 * ((p_i-1)*(1-p_k)*(nodes_[0]->XYZ[2]) + 
		                                    (1+p_i)*(p_k-1)*(nodes_[1]->XYZ[2]) + 
		                                    (1+p_i)*(1-p_k)*(nodes_[2]->XYZ[2]) + 
		                                    (1-p_i)*(1-p_k)*(nodes_[3]->XYZ[2]) + 
		                                    (p_i-1)*(1+p_k)*(nodes_[4]->XYZ[2]) + 
		                                    (-1-p_i)*(1+p_k)*(nodes_[5]->XYZ[2]) + 
		                                    (1+p_i)*(1+p_k)*(nodes_[6]->XYZ[2]) + 
		                                    (1-p_i)*(1+p_k)*(nodes_[7]->XYZ[2]));
		                // printf("J_32:%f\n",J_32);
		                double J_33 = 1.0/8 * ((p_i-1)*(1-p_j)*(nodes_[0]->XYZ[2]) + 
		                                    (1+p_i)*(p_j-1)*(nodes_[1]->XYZ[2]) + 
		                                    (-1-p_i)*(1+p_j)*(nodes_[2]->XYZ[2]) + 
		                                    (p_i-1)*(1+p_j)*(nodes_[3]->XYZ[2]) + 
		                                    (1-p_i)*(1-p_j)*(nodes_[4]->XYZ[2]) + 
		                                    (1+p_i)*(1-p_j)*(nodes_[5]->XYZ[2]) + 
		                                    (1+p_i)*(1+p_j)*(nodes_[6]->XYZ[2]) + 
		                                    (1-p_i)*(1+p_j)*(nodes_[7]->XYZ[2]));
		                // printf("J_33:%f\n",J_33);
		                // double J[3][3] = {0};
		                // for (unsigned int ii = 0; ii < 3; ++ii) {
		                //     for (unsigned int jj = 0; jj < 3; ++jj) {
		                //         for (unsigned int kk = 0; kk < 8; ++kk) {
		                //             J[ii][jj] += n[kk][ii] * nodes_[kk]->XYZ[jj];
		                //         }
		                //     }
		                // }
		    
		                double det_J = J_11*(J_22*J_33-J_23*J_32)-J_12*(J_21*J_33-J_23*J_31)+J_13*(J_21*J_32-J_31*J_22);
		                // double det_J = J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[2][0]*J[1][1]);
		
		                // Inv of J
		                double J_inv_11 = (J_22*J_33-J_23*J_32) / det_J;
		                double J_inv_21 = - (J_21*J_33-J_23*J_31) / det_J;
		                double J_inv_31 = (J_21*J_32-J_31*J_22) / det_J;
		                double J_inv_12 = - (J_12*J_33-J_32*J_13) / det_J;
		                double J_inv_22 = (J_11*J_33-J_13*J_31) / det_J;
		                double J_inv_32 = - (J_11*J_32-J_31*J_12) / det_J;
		                double J_inv_13 = (J_12*J_23-J_22*J_13) / det_J;
		                double J_inv_23 = - (J_11*J_23-J_21*J_13) / det_J;
		                double J_inv_33 = (J_11*J_22-J_21*J_12) / det_J;
				double J_inv[3][3] = {0};
				J_inv[0][0] = J_inv_11;
		                J_inv[0][1] = J_inv_12;
		                J_inv[0][2] = J_inv_13;
		                J_inv[1][0] = J_inv_21;
		                J_inv[1][1] = J_inv_22;
		                J_inv[1][2] = J_inv_23;
		                J_inv[2][0] = J_inv_31;
		                J_inv[2][1] = J_inv_32;
		                J_inv[2][2] = J_inv_33;
		
		                // Grad of shape function
		                double N1_x = 1.0/8 * (J_inv_11 * (p_j-1)*(1-p_k) + J_inv_12 * (p_i-1)*(1-p_k) + J_inv_13 * (p_i-1)*(1-p_j));
		                double N2_x = 1.0/8 * (J_inv_11 * (1-p_j)*(1-p_k) + J_inv_12 * (1+p_i)*(p_k-1) + J_inv_13 * (1+p_i)*(p_j-1));
		                double N3_x = 1.0/8 * (J_inv_11 * (1+p_j)*(1-p_k) + J_inv_12 * (1+p_i)*(1-p_k) + J_inv_13 * (-1-p_i)*(1+p_j));
		                double N4_x = 1.0/8 * (J_inv_11 * (1+p_j)*(p_k-1) + J_inv_12 * (1-p_i)*(1-p_k) + J_inv_13 * (p_i-1)*(1+p_j));
		                double N5_x = 1.0/8 * (J_inv_11 * (p_j-1)*(1+p_k) + J_inv_12 * (p_i-1)*(1+p_k) + J_inv_13 * (1-p_i)*(1-p_j));
		                double N6_x = 1.0/8 * (J_inv_11 * (1-p_j)*(1+p_k) + J_inv_12 * (-1-p_i)*(1+p_k) + J_inv_13 * (1+p_i)*(1-p_j));
		                double N7_x = 1.0/8 * (J_inv_11 * (1+p_j)*(1+p_k) + J_inv_12 * (1+p_i)*(1+p_k) + J_inv_13 * (1+p_i)*(1+p_j));
		                double N8_x = 1.0/8 * (J_inv_11 * (-1-p_j)*(1+p_k) + J_inv_12 * (1-p_i)*(1+p_k) + J_inv_13 * (1-p_i)*(1+p_j));
		
		                double N1_y = 1.0/8 * (J_inv_21 * (p_j-1)*(1-p_k) + J_inv_22 * (p_i-1)*(1-p_k) + J_inv_23 * (p_i-1)*(1-p_j));
		                double N2_y = 1.0/8 * (J_inv_21 * (1-p_j)*(1-p_k) + J_inv_22 * (1+p_i)*(p_k-1) + J_inv_23 * (1+p_i)*(p_j-1));
		                double N3_y = 1.0/8 * (J_inv_21 * (1+p_j)*(1-p_k) + J_inv_22 * (1+p_i)*(1-p_k) + J_inv_23 * (-1-p_i)*(1+p_j));
		                double N4_y = 1.0/8 * (J_inv_21 * (1+p_j)*(p_k-1) + J_inv_22 * (1-p_i)*(1-p_k) + J_inv_23 * (p_i-1)*(1+p_j));
		                double N5_y = 1.0/8 * (J_inv_21 * (p_j-1)*(1+p_k) + J_inv_22 * (p_i-1)*(1+p_k) + J_inv_23 * (1-p_i)*(1-p_j));
		                double N6_y = 1.0/8 * (J_inv_21 * (1-p_j)*(1+p_k) + J_inv_22 * (-1-p_i)*(1+p_k) + J_inv_23 * (1+p_i)*(1-p_j));
		                double N7_y = 1.0/8 * (J_inv_21 * (1+p_j)*(1+p_k) + J_inv_22 * (1+p_i)*(1+p_k) + J_inv_23 * (1+p_i)*(1+p_j));
		                double N8_y = 1.0/8 * (J_inv_21 * (-1-p_j)*(1+p_k) + J_inv_22 * (1-p_i)*(1+p_k) + J_inv_23 * (1-p_i)*(1+p_j));
		
		                double N1_z = 1.0/8 * (J_inv_31 * (p_j-1)*(1-p_k) + J_inv_32 * (p_i-1)*(1-p_k) + J_inv_33 * (p_i-1)*(1-p_j));
		                double N2_z = 1.0/8 * (J_inv_31 * (1-p_j)*(1-p_k) + J_inv_32 * (1+p_i)*(p_k-1) + J_inv_33 * (1+p_i)*(p_j-1));
		                double N3_z = 1.0/8 * (J_inv_31 * (1+p_j)*(1-p_k) + J_inv_32 * (1+p_i)*(1-p_k) + J_inv_33 * (-1-p_i)*(1+p_j));
		                double N4_z = 1.0/8 * (J_inv_31 * (1+p_j)*(p_k-1) + J_inv_32 * (1-p_i)*(1-p_k) + J_inv_33 * (p_i-1)*(1+p_j));
		                double N5_z = 1.0/8 * (J_inv_31 * (p_j-1)*(1+p_k) + J_inv_32 * (p_i-1)*(1+p_k) + J_inv_33 * (1-p_i)*(1-p_j));
		                double N6_z = 1.0/8 * (J_inv_31 * (1-p_j)*(1+p_k) + J_inv_32 * (-1-p_i)*(1+p_k) + J_inv_33 * (1+p_i)*(1-p_j));
		                double N7_z = 1.0/8 * (J_inv_31 * (1+p_j)*(1+p_k) + J_inv_32 * (1+p_i)*(1+p_k) + J_inv_33 * (1+p_i)*(1+p_j));
		                double N8_z = 1.0/8 * (J_inv_31 * (-1-p_j)*(1+p_k) + J_inv_32 * (1-p_i)*(1+p_k) + J_inv_33 * (1-p_i)*(1+p_j));
		
				double N_xyz[8][3] = {0};
		                for (unsigned int ii = 0; ii < 8; ++ii) {
		                    for (unsigned int jj = 0; jj < 3; ++jj) {
		                        for (unsigned int kk = 0; kk < 3; ++kk) {
		                            N_xyz[ii][jj] += n[ii][kk]*J_inv[jj][kk];
		                        }
		                    }
		                }
		
				double B[6][24] = {0};
		                for (unsigned int ii = 0; ii < 8; ++ii) {
		                    B[0][0+3*ii] = N_xyz[ii][0];
		                    B[1][1+3*ii] = N_xyz[ii][1];
		                    B[2][2+3*ii] = N_xyz[ii][2];
		                    B[3][0+3*ii] = N_xyz[ii][1];
		                    B[3][1+3*ii] = N_xyz[ii][0];
		                    B[4][1+3*ii] = N_xyz[ii][2];
		                    B[4][2+3*ii] = N_xyz[ii][1];
		                    B[5][0+3*ii] = N_xyz[ii][2];
		                    B[5][2+3*ii] = N_xyz[ii][0];
		                }
		                printf("B[2][20]:%f\n",B[2][20]);
		
		                double BTD[24][6] = {0};
		                for (unsigned int ii = 0; ii < 24; ++ii) {
		                    for (unsigned int jj = 0; jj < 6; ++jj) {
		                        for (unsigned int kk = 0; kk < 6; ++kk) {
		                            BTD[ii][jj] += B[kk][ii]*D[kk][jj];
		                        }
		                    }
		                }
		
		                for (unsigned int ii = 0; ii < 24; ++ii) {
		                    for (unsigned int jj = 0; jj < 24; ++jj) {
		                        for (unsigned int kk = 0; kk < 6; ++kk) {
		                            K[ii][jj] += w_i*w_j*w_k* BTD[ii][kk]*B[kk][jj] * det_J;
		                        }
		                    }
		                }
		
		                // // L1
		                // Matrix[0] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_x*N1_x + miu*N1_y*N1_y + miu*N1_z*N1_z) * det_J;
		
		                // // L2
		                // Matrix[1] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_y*N1_y + miu*N1_x*N1_x + miu*N1_z*N1_z) * det_J;
		                // Matrix[2] += w_i*w_j*w_k* (lambda*N1_x*N1_y + miu*N1_x*N1_y) * det_J;
		
		                // // L3
		                // Matrix[3] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_z*N1_z + miu*N1_x*N1_x + miu*N1_y*N1_y) * det_J;
		                // Matrix[4] += w_i*w_j*w_k* (lambda*N1_y*N1_z + miu*N1_y*N1_z) * det_J;
		                // Matrix[5] += w_i*w_j*w_k* (lambda*N1_x*N1_z + miu*N1_x*N1_z) * det_J;
		
		                // // L4
		                // Matrix[6] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_x*N2_x + miu*N2_y*N2_y + miu*N2_z*N2_z) * det_J;
		                // Matrix[7] += w_i*w_j*w_k* (lambda*N1_z*N2_x + miu*N1_x*N2_z) * det_J;
		                // Matrix[8] += w_i*w_j*w_k* (lambda*N1_y*N2_x + miu*N1_x*N2_y) * det_J;
		                // Matrix[9] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_x*N2_x + miu*N1_y*N2_y + miu*N1_z*N2_z) * det_J;
		
		                // // L5
		                // Matrix[10] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_y*N2_y + miu*N2_x*N2_x + miu*N2_z*N2_z) * det_J;
		                // Matrix[11] += w_i*w_j*w_k* (lambda*N2_y*N2_x + miu*N2_x*N2_y) * det_J;
		                // Matrix[12] += w_i*w_j*w_k* (lambda*N2_y*N1_z + miu*N1_y*N2_z) * det_J;
		                // Matrix[13] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_y*N2_y + miu*N1_x*N2_x + miu*N1_z*N2_z) * det_J;
		                // Matrix[14] += w_i*w_j*w_k* (lambda*N1_x*N2_y + miu*N1_y*N2_x) * det_J;
		                
		                // // L6
		                // Matrix[15] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_z*N2_z + miu*N2_x*N2_x + miu*N2_y*N2_y) * det_J;
		                // Matrix[16] += w_i*w_j*w_k* (lambda*N2_y*N2_z + miu*N2_z*N2_y) * det_J;
		                // Matrix[17] += w_i*w_j*w_k* (lambda*N2_x*N2_z + miu*N2_z*N2_x) * det_J;
		                // Matrix[18] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_z*N2_z + miu*N1_x*N2_x + miu*N1_y*N2_y) * det_J;
		                // Matrix[19] += w_i*w_j*w_k* (lambda*N1_y*N2_z + miu*N1_z*N2_y) * det_J;
		                // Matrix[20] += w_i*w_j*w_k* (lambda*N1_x*N2_z + miu*N1_z*N2_x) * det_J;
		                
		                // // L7
		                // Matrix[21] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_x*N3_x + miu*N3_y*N3_y + miu*N3_z*N3_z) * det_J;
		                // Matrix[22] += w_i*w_j*w_k* (lambda*N2_z*N3_x + miu*N2_x*N3_z) * det_J;
		                // Matrix[23] += w_i*w_j*w_k* (lambda*N2_y*N3_x + miu*N2_x*N3_y) * det_J;
		                // Matrix[24] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_x*N3_x + miu*N2_y*N3_y + miu*N2_z*N3_z) * det_J;
		                // Matrix[25] += w_i*w_j*w_k* (lambda*N1_z*N3_x + miu*N1_x*N3_z) * det_J;
		                // Matrix[26] += w_i*w_j*w_k* (lambda*N1_y*N3_x + miu*N1_x*N3_y) * det_J;
		                // Matrix[27] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_x*N3_x + miu*N1_y*N3_y + miu*N1_z*N3_z) * det_J;
		
		                // // L8
		                // Matrix[28] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_y*N3_y + miu*N3_x*N3_x + miu*N3_z*N3_z) * det_J;
		                // Matrix[29] += w_i*w_j*w_k* (lambda*N3_x*N3_y + miu*N3_y*N3_x) * det_J;
		                // Matrix[30] += w_i*w_j*w_k* (lambda*N2_z*N3_y + miu*N2_y*N3_z) * det_J;
		                // Matrix[31] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_y*N3_y + miu*N2_x*N3_x + miu*N2_z*N3_z) * det_J;
		                // Matrix[32] += w_i*w_j*w_k* (lambda*N2_x*N3_y + miu*N2_y*N3_x) * det_J;
		                // Matrix[33] += w_i*w_j*w_k* (lambda*N1_z*N3_y + miu*N1_y*N3_z) * det_J;
		                // Matrix[34] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_y*N3_y + miu*N1_x*N3_x + miu*N1_z*N3_z) * det_J;
		                // Matrix[35] += w_i*w_j*w_k* (lambda*N1_x*N3_y + miu*N1_y*N3_x) * det_J;
		
		                // // L9
		                // Matrix[36] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_z*N3_z + miu*N3_x*N3_x + miu*N3_y*N3_y) * det_J;
		                // Matrix[37] += w_i*w_j*w_k* (lambda*N3_y*N3_z + miu*N3_z*N3_y) * det_J;
		                // Matrix[38] += w_i*w_j*w_k* (lambda*N3_x*N3_z + miu*N3_z*N3_x) * det_J;
		                // Matrix[39] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_z*N3_z + miu*N2_x*N3_x + miu*N2_y*N3_y) * det_J;
		                // Matrix[40] += w_i*w_j*w_k* (lambda*N2_y*N3_z + miu*N2_z*N3_y) * det_J;
		                // Matrix[41] += w_i*w_j*w_k* (lambda*N2_x*N3_z + miu*N2_z*N3_x) * det_J;
		                // Matrix[42] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_z*N3_z + miu*N1_x*N3_x + miu*N1_y*N3_y) * det_J;
		                // Matrix[43] += w_i*w_j*w_k* (lambda*N1_y*N3_z + miu*N1_z*N3_y) * det_J;
		                // Matrix[44] += w_i*w_j*w_k* (lambda*N1_x*N3_z + miu*N1_z*N3_x) * det_J;
		
		                // // L10
		                // Matrix[45] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_x*N4_x + miu*N4_y*N4_y + miu*N4_z*N4_z) * det_J;
		                // Matrix[46] += w_i*w_j*w_k* (lambda*N3_z*N4_x + miu*N3_x*N4_z) * det_J;
		                // Matrix[47] += w_i*w_j*w_k* (lambda*N3_y*N4_x + miu*N3_x*N4_y) * det_J;
		                // Matrix[48] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_x*N4_x + miu*N3_y*N4_y + miu*N3_z*N4_z) * det_J;
		                // Matrix[49] += w_i*w_j*w_k* (lambda*N2_z*N4_x + miu*N2_x*N4_z) * det_J;
		                // Matrix[50] += w_i*w_j*w_k* (lambda*N2_y*N4_x + miu*N2_x*N4_y) * det_J;
		                // Matrix[51] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_x*N4_x + miu*N2_y*N4_y + miu*N2_z*N4_z) * det_J;
		                // Matrix[52] += w_i*w_j*w_k* (lambda*N1_z*N4_x + miu*N1_x*N4_z) * det_J;
		                // Matrix[53] += w_i*w_j*w_k* (lambda*N1_y*N4_x + miu*N1_x*N4_y) * det_J;
		                // Matrix[54] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_x*N4_x + miu*N1_y*N4_y + miu*N1_z*N4_z) * det_J;
		
		                // // L11
		                // Matrix[55] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_y*N4_y + miu*N4_x*N4_x + miu*N4_z*N4_z) * det_J;
		                // Matrix[56] += w_i*w_j*w_k* (lambda*N4_x*N4_y + miu*N4_y*N4_x) * det_J;
		                // Matrix[57] += w_i*w_j*w_k* (lambda*N3_z*N4_y + miu*N3_y*N4_z) * det_J;
		                // Matrix[58] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_y*N4_y + miu*N3_x*N4_x + miu*N3_z*N4_z) * det_J;
		                // Matrix[59] += w_i*w_j*w_k* (lambda*N3_x*N4_y + miu*N3_y*N4_x) * det_J;
		                // Matrix[60] += w_i*w_j*w_k* (lambda*N2_z*N4_y + miu*N2_y*N4_z) * det_J;
		                // Matrix[61] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_y*N4_y + miu*N2_x*N4_x + miu*N2_z*N4_z) * det_J;
		                // Matrix[62] += w_i*w_j*w_k* (lambda*N2_x*N4_y + miu*N2_y*N4_x) * det_J;
		                // Matrix[63] += w_i*w_j*w_k* (lambda*N1_z*N4_y + miu*N1_y*N4_z) * det_J;
		                // Matrix[64] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_y*N4_y + miu*N1_x*N4_x + miu*N1_z*N4_z) * det_J;
		                // Matrix[65] += w_i*w_j*w_k* (lambda*N1_x*N4_y + miu*N1_y*N4_x) * det_J;
		
		                // // L12
		                // Matrix[66] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_z*N4_z + miu*N4_x*N4_x + miu*N4_y*N4_y) * det_J;
		                // Matrix[67] += w_i*w_j*w_k* (lambda*N4_y*N4_z + miu*N4_z*N4_y) * det_J;
		                // Matrix[68] += w_i*w_j*w_k* (lambda*N4_x*N4_z + miu*N4_z*N4_x) * det_J;
		                // Matrix[69] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_z*N4_z + miu*N3_x*N4_x + miu*N3_y*N4_y) * det_J;
		                // Matrix[70] += w_i*w_j*w_k* (lambda*N3_y*N4_z + miu*N3_z*N4_y) * det_J;
		                // Matrix[71] += w_i*w_j*w_k* (lambda*N3_x*N4_z + miu*N3_z*N4_x) * det_J;
		                // Matrix[72] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_z*N4_z + miu*N2_x*N4_x + miu*N2_y*N4_y) * det_J;
		                // Matrix[73] += w_i*w_j*w_k* (lambda*N2_y*N4_z + miu*N2_z*N4_y) * det_J;
		                // Matrix[74] += w_i*w_j*w_k* (lambda*N2_x*N4_z + miu*N2_z*N4_x) * det_J;
		                // Matrix[75] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_z*N4_z + miu*N1_x*N4_x + miu*N1_y*N4_y) * det_J;
		                // Matrix[76] += w_i*w_j*w_k* (lambda*N1_y*N4_z + miu*N1_z*N4_y) * det_J;
		                // Matrix[77] += w_i*w_j*w_k* (lambda*N1_x*N4_z + miu*N1_z*N4_x) * det_J;
		
		                // // L13
		                // Matrix[78] += w_i*w_j*w_k* ((lambda + 2*miu)*N5_x*N5_x + miu*N5_y*N5_y + miu*N5_z*N5_z) * det_J;
		                // Matrix[79] += w_i*w_j*w_k* (lambda*N4_z*N5_x + miu*N4_x*N5_z) * det_J;
		                // Matrix[80] += w_i*w_j*w_k* (lambda*N4_y*N5_x + miu*N4_x*N5_y) * det_J;
		                // Matrix[81] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_x*N5_x + miu*N4_y*N5_y + miu*N4_z*N5_z) * det_J;
		                // Matrix[82] += w_i*w_j*w_k* (lambda*N3_z*N5_x + miu*N3_x*N5_z) * det_J;
		                // Matrix[83] += w_i*w_j*w_k* (lambda*N3_y*N5_x + miu*N3_x*N5_y) * det_J;
		                // Matrix[84] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_x*N5_x + miu*N3_y*N5_y + miu*N3_z*N5_z) * det_J;
		                // Matrix[85] += w_i*w_j*w_k* (lambda*N2_z*N5_x + miu*N2_x*N5_z) * det_J;
		                // Matrix[86] += w_i*w_j*w_k* (lambda*N2_y*N5_x + miu*N2_x*N5_y) * det_J;
		                // Matrix[87] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_x*N5_x + miu*N2_y*N5_y + miu*N2_z*N5_z) * det_J;
		                // Matrix[88] += w_i*w_j*w_k* (lambda*N1_z*N5_x + miu*N1_x*N5_z) * det_J;
		                // Matrix[89] += w_i*w_j*w_k* (lambda*N1_y*N5_x + miu*N1_x*N5_y) * det_J;
		                // Matrix[90] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_x*N5_x + miu*N1_y*N5_y + miu*N1_z*N5_z) * det_J;
		
		                // // L14
		                // Matrix[91] += w_i*w_j*w_k* ((lambda + 2*miu)*N5_y*N5_y + miu*N5_x*N5_x + miu*N5_z*N5_z) * det_J;
		                // Matrix[92] += w_i*w_j*w_k* (lambda*N5_x*N5_y + miu*N5_y*N5_x) * det_J;
		                // Matrix[93] += w_i*w_j*w_k* (lambda*N4_z*N5_y + miu*N4_y*N5_z) * det_J;
		                // Matrix[94] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_y*N5_y + miu*N4_x*N5_x + miu*N4_z*N5_z) * det_J;
		                // Matrix[95] += w_i*w_j*w_k* (lambda*N4_x*N5_y + miu*N4_y*N5_x) * det_J;
		                // Matrix[96] += w_i*w_j*w_k* (lambda*N3_z*N5_y + miu*N3_y*N5_z) * det_J;
		                // Matrix[97] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_y*N5_y + miu*N3_x*N5_x + miu*N3_z*N5_z) * det_J;
		                // Matrix[98] += w_i*w_j*w_k* (lambda*N3_x*N5_y + miu*N3_y*N5_x) * det_J;
		                // Matrix[99] += w_i*w_j*w_k* (lambda*N2_z*N5_y + miu*N2_y*N5_z) * det_J;
		                // Matrix[100] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_y*N5_y + miu*N2_x*N5_x + miu*N2_z*N5_z) * det_J;
		                // Matrix[101] += w_i*w_j*w_k* (lambda*N2_x*N5_y + miu*N2_y*N5_x) * det_J;
		                // Matrix[102] += w_i*w_j*w_k* (lambda*N1_z*N5_y + miu*N1_y*N5_z) * det_J;
		                // Matrix[103] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_y*N5_y + miu*N1_x*N5_x + miu*N1_z*N5_z) * det_J;
		                // Matrix[104] += w_i*w_j*w_k* (lambda*N1_x*N5_y + miu*N1_y*N5_x) * det_J;
		
		                // // L15
		                // Matrix[105] += w_i*w_j*w_k* ((lambda + 2*miu)*N5_z*N5_z + miu*N5_x*N5_x + miu*N5_y*N5_y) * det_J;
		                // Matrix[106] += w_i*w_j*w_k* (lambda*N5_y*N5_z + miu*N5_z*N5_y) * det_J;
		                // Matrix[107] += w_i*w_j*w_k* (lambda*N5_x*N5_z + miu*N5_z*N5_x) * det_J;
		                // Matrix[108] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_z*N5_z + miu*N4_x*N5_x + miu*N4_y*N5_y) * det_J;
		                // Matrix[109] += w_i*w_j*w_k* (lambda*N4_y*N5_z + miu*N4_z*N5_y) * det_J;
		                // Matrix[110] += w_i*w_j*w_k* (lambda*N4_x*N5_z + miu*N4_z*N5_x) * det_J;
		                // Matrix[111] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_z*N5_z + miu*N3_x*N5_x + miu*N3_y*N5_y) * det_J;
		                // Matrix[112] += w_i*w_j*w_k* (lambda*N3_y*N5_z + miu*N3_z*N5_y) * det_J;
		                // Matrix[113] += w_i*w_j*w_k* (lambda*N3_x*N5_z + miu*N3_z*N5_x) * det_J;
		                // Matrix[114] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_z*N5_z + miu*N2_x*N5_x + miu*N2_y*N5_y) * det_J;
		                // Matrix[115] += w_i*w_j*w_k* (lambda*N2_y*N5_z + miu*N2_z*N5_y) * det_J;
		                // Matrix[116] += w_i*w_j*w_k* (lambda*N2_x*N5_z + miu*N2_z*N5_x) * det_J;
		                // Matrix[117] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_z*N5_z + miu*N1_x*N5_x + miu*N1_y*N5_y) * det_J;
		                // Matrix[118] += w_i*w_j*w_k* (lambda*N1_y*N5_z + miu*N1_z*N5_y) * det_J;
		                // Matrix[119] += w_i*w_j*w_k* (lambda*N1_x*N5_z + miu*N1_z*N5_x) * det_J;
		
		                // // L16
		                // Matrix[120] += w_i*w_j*w_k* ((lambda + 2*miu)*N6_x*N6_x + miu*N6_y*N6_y + miu*N6_z*N6_z) * det_J;
		                // Matrix[121] += w_i*w_j*w_k* (lambda*N5_z*N6_x + miu*N5_x*N6_z) * det_J;
		                // Matrix[122] += w_i*w_j*w_k* (lambda*N5_y*N6_x + miu*N5_x*N6_y) * det_J;
		                // Matrix[123] += w_i*w_j*w_k* ((lambda + 2*miu)*N5_x*N6_x + miu*N5_y*N6_y + miu*N5_z*N6_z) * det_J;
		                // Matrix[124] += w_i*w_j*w_k* (lambda*N4_z*N6_x + miu*N4_x*N6_z) * det_J;
		                // Matrix[125] += w_i*w_j*w_k* (lambda*N4_y*N6_x + miu*N4_x*N6_y) * det_J;
		                // Matrix[126] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_x*N6_x + miu*N4_y*N6_y + miu*N4_z*N6_z) * det_J;
		                // Matrix[127] += w_i*w_j*w_k* (lambda*N3_z*N6_x + miu*N3_x*N6_z) * det_J;
		                // Matrix[128] += w_i*w_j*w_k* (lambda*N3_y*N6_x + miu*N3_x*N6_y) * det_J;
		                // Matrix[129] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_x*N6_x + miu*N3_y*N6_y + miu*N3_z*N6_z) * det_J;
		                // Matrix[130] += w_i*w_j*w_k* (lambda*N2_z*N6_x + miu*N2_x*N6_z) * det_J;
		                // Matrix[131] += w_i*w_j*w_k* (lambda*N2_y*N6_x + miu*N2_x*N6_y) * det_J;
		                // Matrix[132] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_x*N6_x + miu*N2_y*N6_y + miu*N2_z*N6_z) * det_J;
		                // Matrix[133] += w_i*w_j*w_k* (lambda*N1_z*N6_x + miu*N1_x*N6_z) * det_J;
		                // Matrix[134] += w_i*w_j*w_k* (lambda*N1_y*N6_x + miu*N1_x*N6_y) * det_J;
		                // Matrix[135] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_x*N6_x + miu*N1_y*N6_y + miu*N1_z*N6_z) * det_J;
		
		                // // L17
		                // Matrix[136] += w_i*w_j*w_k* ((lambda + 2*miu)*N6_y*N6_y + miu*N6_x*N6_x + miu*N6_z*N6_z) * det_J;
		                // Matrix[137] += w_i*w_j*w_k* (lambda*N6_x*N6_y + miu*N6_y*N6_x) * det_J;
		                // Matrix[138] += w_i*w_j*w_k* (lambda*N5_z*N6_y + miu*N5_y*N6_z) * det_J;
		                // Matrix[139] += w_i*w_j*w_k* ((lambda + 2*miu)*N5_y*N6_y + miu*N5_x*N6_x + miu*N5_z*N6_z) * det_J;
		                // Matrix[140] += w_i*w_j*w_k* (lambda*N5_x*N6_y + miu*N5_y*N6_x) * det_J;
		                // Matrix[141] += w_i*w_j*w_k* (lambda*N4_z*N6_y + miu*N4_y*N6_z) * det_J;
		                // Matrix[142] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_y*N6_y + miu*N4_x*N6_x + miu*N4_z*N6_z) * det_J;
		                // Matrix[143] += w_i*w_j*w_k* (lambda*N4_x*N6_y + miu*N4_y*N6_x) * det_J;
		                // Matrix[144] += w_i*w_j*w_k* (lambda*N3_z*N6_y + miu*N3_y*N6_z) * det_J;
		                // Matrix[145] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_y*N6_y + miu*N3_x*N6_x + miu*N3_z*N6_z) * det_J;
		                // Matrix[146] += w_i*w_j*w_k* (lambda*N3_x*N6_y + miu*N3_y*N6_x) * det_J;
		                // Matrix[147] += w_i*w_j*w_k* (lambda*N2_z*N6_y + miu*N2_y*N6_z) * det_J;
		                // Matrix[148] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_y*N6_y + miu*N2_x*N6_x + miu*N2_z*N6_z) * det_J;
		                // Matrix[149] += w_i*w_j*w_k* (lambda*N2_x*N6_y + miu*N2_y*N6_x) * det_J;
		                // Matrix[150] += w_i*w_j*w_k* (lambda*N1_z*N6_y + miu*N1_y*N6_z) * det_J;
		                // Matrix[151] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_y*N6_y + miu*N1_x*N6_x + miu*N1_z*N6_z) * det_J;
		                // Matrix[152] += w_i*w_j*w_k* (lambda*N1_x*N6_y + miu*N1_y*N6_x) * det_J;
		
		                // // L18
		                // Matrix[153] += w_i*w_j*w_k* ((lambda + 2*miu)*N6_z*N6_z + miu*N6_x*N6_x + miu*N6_y*N6_y) * det_J;
		                // Matrix[154] += w_i*w_j*w_k* (lambda*N6_y*N6_z + miu*N6_z*N6_y) * det_J;
		                // Matrix[155] += w_i*w_j*w_k* (lambda*N6_x*N6_z + miu*N6_z*N6_x) * det_J;
		                // Matrix[156] += w_i*w_j*w_k* ((lambda + 2*miu)*N5_z*N6_z + miu*N5_x*N6_x + miu*N5_y*N6_y) * det_J;
		                // Matrix[157] += w_i*w_j*w_k* (lambda*N5_y*N6_z + miu*N5_z*N6_y) * det_J;
		                // Matrix[158] += w_i*w_j*w_k* (lambda*N5_x*N6_z + miu*N5_z*N6_x) * det_J;
		                // Matrix[159] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_z*N6_z + miu*N4_x*N6_x + miu*N4_y*N6_y) * det_J;
		                // Matrix[160] += w_i*w_j*w_k* (lambda*N4_y*N6_z + miu*N4_z*N6_y) * det_J;
		                // Matrix[161] += w_i*w_j*w_k* (lambda*N4_x*N6_z + miu*N4_z*N6_x) * det_J;
		                // Matrix[162] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_z*N6_z + miu*N3_x*N6_x + miu*N3_y*N6_y) * det_J;
		                // Matrix[163] += w_i*w_j*w_k* (lambda*N3_y*N6_z + miu*N3_z*N6_y) * det_J;
		                // Matrix[164] += w_i*w_j*w_k* (lambda*N3_x*N6_z + miu*N3_z*N6_x) * det_J;
		                // Matrix[165] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_z*N6_z + miu*N2_x*N6_x + miu*N2_y*N6_y) * det_J;
		                // Matrix[166] += w_i*w_j*w_k* (lambda*N2_y*N6_z + miu*N2_z*N6_y) * det_J;
		                // Matrix[167] += w_i*w_j*w_k* (lambda*N2_x*N6_z + miu*N2_z*N6_x) * det_J;
		                // Matrix[168] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_z*N6_z + miu*N1_x*N6_x + miu*N1_y*N6_y) * det_J;
		                // Matrix[169] += w_i*w_j*w_k* (lambda*N1_y*N6_z + miu*N1_z*N6_y) * det_J;
		                // Matrix[170] += w_i*w_j*w_k* (lambda*N1_x*N6_z + miu*N1_z*N6_x) * det_J;
		
		                // // L19
		                // Matrix[171] += w_i*w_j*w_k* ((lambda + 2*miu)*N7_x*N7_x + miu*N7_y*N7_y + miu*N7_z*N7_z) * det_J;
		                // Matrix[172] += w_i*w_j*w_k* (lambda*N6_z*N7_x + miu*N6_x*N7_z) * det_J;
		                // Matrix[173] += w_i*w_j*w_k* (lambda*N6_y*N7_x + miu*N6_x*N7_y) * det_J;
		                // Matrix[174] += w_i*w_j*w_k* ((lambda + 2*miu)*N6_x*N7_x + miu*N6_y*N7_y + miu*N6_z*N7_z) * det_J;
		                // Matrix[175] += w_i*w_j*w_k* (lambda*N5_z*N7_x + miu*N5_x*N7_z) * det_J;
		                // Matrix[176] += w_i*w_j*w_k* (lambda*N5_y*N7_x + miu*N5_x*N7_y) * det_J;
		                // Matrix[177] += w_i*w_j*w_k* ((lambda + 2*miu)*N5_x*N7_x + miu*N5_y*N7_y + miu*N5_z*N7_z) * det_J;
		                // Matrix[178] += w_i*w_j*w_k* (lambda*N4_z*N7_x + miu*N4_x*N7_z) * det_J;
		                // Matrix[179] += w_i*w_j*w_k* (lambda*N4_y*N7_x + miu*N4_x*N7_y) * det_J;
		                // Matrix[180] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_x*N7_x + miu*N4_y*N7_y + miu*N4_z*N7_z) * det_J;
		                // Matrix[181] += w_i*w_j*w_k* (lambda*N3_z*N7_x + miu*N3_x*N7_z) * det_J;
		                // Matrix[182] += w_i*w_j*w_k* (lambda*N3_y*N7_x + miu*N3_x*N7_y) * det_J;
		                // Matrix[183] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_x*N7_x + miu*N3_y*N7_y + miu*N3_z*N7_z) * det_J;
		                // Matrix[184] += w_i*w_j*w_k* (lambda*N2_z*N7_x + miu*N2_x*N7_z) * det_J;
		                // Matrix[185] += w_i*w_j*w_k* (lambda*N2_y*N7_x + miu*N2_x*N7_y) * det_J;
		                // Matrix[186] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_x*N7_x + miu*N2_y*N7_y + miu*N2_z*N7_z) * det_J;
		                // Matrix[187] += w_i*w_j*w_k* (lambda*N1_z*N7_x + miu*N1_x*N7_z) * det_J;
		                // Matrix[188] += w_i*w_j*w_k* (lambda*N1_y*N7_x + miu*N1_x*N7_y) * det_J;
		                // Matrix[189] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_x*N7_x + miu*N1_y*N7_y + miu*N1_z*N7_z) * det_J;
		
		                // // L20
		                // Matrix[190] += w_i*w_j*w_k* ((lambda + 2*miu)*N7_y*N7_y + miu*N7_x*N7_x + miu*N7_z*N7_z) * det_J;
		                // Matrix[191] += w_i*w_j*w_k* (lambda*N7_x*N7_y + miu*N7_y*N7_x) * det_J;
		                // Matrix[192] += w_i*w_j*w_k* (lambda*N6_z*N7_y + miu*N6_y*N7_z) * det_J;
		                // Matrix[193] += w_i*w_j*w_k* ((lambda + 2*miu)*N6_y*N7_y + miu*N6_x*N7_x + miu*N6_z*N7_z) * det_J;
		                // Matrix[194] += w_i*w_j*w_k* (lambda*N6_x*N7_y + miu*N6_y*N7_x) * det_J;
		                // Matrix[195] += w_i*w_j*w_k* (lambda*N5_z*N7_y + miu*N5_y*N7_z) * det_J;
		                // Matrix[196] += w_i*w_j*w_k* ((lambda + 2*miu)*N5_y*N7_y + miu*N5_x*N7_x + miu*N5_z*N7_z) * det_J;
		                // Matrix[197] += w_i*w_j*w_k* (lambda*N5_x*N7_y + miu*N5_y*N7_x) * det_J;
		                // Matrix[198] += w_i*w_j*w_k* (lambda*N4_z*N7_y + miu*N4_y*N7_z) * det_J;
		                // Matrix[199] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_y*N7_y + miu*N4_x*N7_x + miu*N4_z*N7_z) * det_J;
		                // Matrix[200] += w_i*w_j*w_k* (lambda*N4_x*N7_y + miu*N4_y*N7_x) * det_J;
		                // Matrix[201] += w_i*w_j*w_k* (lambda*N3_z*N7_y + miu*N3_y*N7_z) * det_J;
		                // Matrix[202] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_y*N7_y + miu*N3_x*N7_x + miu*N3_z*N7_z) * det_J;
		                // Matrix[203] += w_i*w_j*w_k* (lambda*N3_x*N7_y + miu*N3_y*N7_x) * det_J;
		                // Matrix[204] += w_i*w_j*w_k* (lambda*N2_z*N7_y + miu*N2_y*N7_z) * det_J;
		                // Matrix[205] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_y*N7_y + miu*N2_x*N7_x + miu*N2_z*N7_z) * det_J;
		                // Matrix[206] += w_i*w_j*w_k* (lambda*N2_x*N7_y + miu*N2_y*N7_x) * det_J;
		                // Matrix[207] += w_i*w_j*w_k* (lambda*N1_z*N7_y + miu*N1_y*N7_z) * det_J;
		                // Matrix[208] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_y*N7_y + miu*N1_x*N7_x + miu*N1_z*N7_z) * det_J;
		                // Matrix[209] += w_i*w_j*w_k* (lambda*N1_x*N7_y + miu*N1_y*N7_x) * det_J;
		
		                // // L21
		                // Matrix[210] += w_i*w_j*w_k* ((lambda + 2*miu)*N7_z*N7_z + miu*N7_x*N7_x + miu*N7_y*N7_y) * det_J;
		                // Matrix[211] += w_i*w_j*w_k* (lambda*N7_y*N7_z + miu*N7_z*N7_y) * det_J;
		                // Matrix[212] += w_i*w_j*w_k* (lambda*N7_x*N7_z + miu*N7_z*N7_x) * det_J;
		                // Matrix[213] += w_i*w_j*w_k* ((lambda + 2*miu)*N6_z*N7_z + miu*N6_x*N7_x + miu*N6_y*N7_y) * det_J;
		                // Matrix[214] += w_i*w_j*w_k* (lambda*N6_y*N7_z + miu*N6_z*N7_y) * det_J;
		                // Matrix[215] += w_i*w_j*w_k* (lambda*N6_x*N7_z + miu*N6_z*N7_x) * det_J;
		                // Matrix[216] += w_i*w_j*w_k* ((lambda + 2*miu)*N5_z*N7_z + miu*N5_x*N7_x + miu*N5_y*N7_y) * det_J;
		                // Matrix[217] += w_i*w_j*w_k* (lambda*N5_y*N7_z + miu*N5_z*N7_y) * det_J;
		                // Matrix[218] += w_i*w_j*w_k* (lambda*N5_x*N7_z + miu*N5_z*N7_x) * det_J;
		                // Matrix[219] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_z*N7_z + miu*N4_x*N7_x + miu*N4_y*N7_y) * det_J;
		                // Matrix[220] += w_i*w_j*w_k* (lambda*N4_y*N7_z + miu*N4_z*N7_y) * det_J;
		                // Matrix[221] += w_i*w_j*w_k* (lambda*N4_x*N7_z + miu*N4_z*N7_x) * det_J;
		                // Matrix[222] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_z*N7_z + miu*N3_x*N7_x + miu*N3_y*N7_y) * det_J;
		                // Matrix[223] += w_i*w_j*w_k* (lambda*N3_y*N7_z + miu*N3_z*N7_y) * det_J;
		                // Matrix[224] += w_i*w_j*w_k* (lambda*N3_x*N7_z + miu*N3_z*N7_x) * det_J;
		                // Matrix[225] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_z*N7_z + miu*N2_x*N7_x + miu*N2_y*N7_y) * det_J;
		                // Matrix[226] += w_i*w_j*w_k* (lambda*N2_y*N7_z + miu*N2_z*N7_y) * det_J;
		                // Matrix[227] += w_i*w_j*w_k* (lambda*N2_x*N7_z + miu*N2_z*N7_x) * det_J;
		                // Matrix[228] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_z*N7_z + miu*N1_x*N7_x + miu*N1_y*N7_y) * det_J;
		                // Matrix[229] += w_i*w_j*w_k* (lambda*N1_y*N7_z + miu*N1_z*N7_y) * det_J;
		                // Matrix[230] += w_i*w_j*w_k* (lambda*N1_x*N7_z + miu*N1_z*N7_x) * det_J;
		
		                // // L22
		                // Matrix[231] += w_i*w_j*w_k* ((lambda + 2*miu)*N8_x*N8_x + miu*N8_y*N8_y + miu*N8_z*N8_z) * det_J;
		                // Matrix[232] += w_i*w_j*w_k* (lambda*N7_z*N8_x + miu*N7_x*N8_z) * det_J;
		                // Matrix[233] += w_i*w_j*w_k* (lambda*N7_y*N8_x + miu*N7_x*N8_y) * det_J;
		                // Matrix[234] += w_i*w_j*w_k* ((lambda + 2*miu)*N7_x*N8_x + miu*N7_y*N8_y + miu*N7_z*N8_z) * det_J;
		                // Matrix[235] += w_i*w_j*w_k* (lambda*N6_z*N8_x + miu*N6_x*N8_z) * det_J;
		                // Matrix[236] += w_i*w_j*w_k* (lambda*N6_y*N8_x + miu*N6_x*N8_y) * det_J;
		                // Matrix[237] += w_i*w_j*w_k* ((lambda + 2*miu)*N6_x*N8_x + miu*N6_y*N8_y + miu*N6_z*N8_z) * det_J;
		                // Matrix[238] += w_i*w_j*w_k* (lambda*N5_z*N8_x + miu*N5_x*N8_z) * det_J;
		                // Matrix[239] += w_i*w_j*w_k* (lambda*N5_y*N8_x + miu*N5_x*N8_y) * det_J;
		                // Matrix[240] += w_i*w_j*w_k* ((lambda + 2*miu)*N5_x*N8_x + miu*N5_y*N8_y + miu*N5_z*N8_z) * det_J;
		                // Matrix[241] += w_i*w_j*w_k* (lambda*N4_z*N8_x + miu*N4_x*N8_z) * det_J;
		                // Matrix[242] += w_i*w_j*w_k* (lambda*N4_y*N8_x + miu*N4_x*N8_y) * det_J;
		                // Matrix[243] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_x*N8_x + miu*N4_y*N8_y + miu*N4_z*N8_z) * det_J;
		                // Matrix[244] += w_i*w_j*w_k* (lambda*N3_z*N8_x + miu*N3_x*N8_z) * det_J;
		                // Matrix[245] += w_i*w_j*w_k* (lambda*N3_y*N8_x + miu*N3_x*N8_y) * det_J;
		                // Matrix[246] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_x*N8_x + miu*N3_y*N8_y + miu*N3_z*N8_z) * det_J;
		                // Matrix[247] += w_i*w_j*w_k* (lambda*N2_z*N8_x + miu*N2_x*N8_z) * det_J;
		                // Matrix[248] += w_i*w_j*w_k* (lambda*N2_y*N8_x + miu*N2_x*N8_y) * det_J;
		                // Matrix[249] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_x*N8_x + miu*N2_y*N8_y + miu*N2_z*N8_z) * det_J;
		                // Matrix[250] += w_i*w_j*w_k* (lambda*N1_z*N8_x + miu*N1_x*N8_z) * det_J;
		                // Matrix[251] += w_i*w_j*w_k* (lambda*N1_y*N8_x + miu*N1_x*N8_y) * det_J;
		                // Matrix[252] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_x*N8_x + miu*N1_y*N8_y + miu*N1_z*N8_z) * det_J;
		
		                // // L23
		                // Matrix[253] += w_i*w_j*w_k* ((lambda + 2*miu)*N8_y*N8_y + miu*N8_x*N8_x + miu*N8_z*N8_z) * det_J;
		                // Matrix[254] += w_i*w_j*w_k* (lambda*N8_x*N8_y + miu*N8_y*N8_x) * det_J;
		                // Matrix[255] += w_i*w_j*w_k* (lambda*N7_z*N8_y + miu*N7_y*N8_z) * det_J;
		                // Matrix[256] += w_i*w_j*w_k* ((lambda + 2*miu)*N7_y*N8_y + miu*N7_x*N8_x + miu*N7_z*N8_z) * det_J;
		                // Matrix[257] += w_i*w_j*w_k* (lambda*N7_x*N8_y + miu*N7_y*N8_x) * det_J;
		                // Matrix[258] += w_i*w_j*w_k* (lambda*N6_z*N8_y + miu*N6_y*N8_z) * det_J;
		                // Matrix[259] += w_i*w_j*w_k* ((lambda + 2*miu)*N6_y*N8_y + miu*N6_x*N8_x + miu*N6_z*N8_z) * det_J;
		                // Matrix[260] += w_i*w_j*w_k* (lambda*N6_x*N8_y + miu*N6_y*N8_x) * det_J;
		                // Matrix[261] += w_i*w_j*w_k* (lambda*N5_z*N8_y + miu*N5_y*N8_z) * det_J;
		                // Matrix[262] += w_i*w_j*w_k* ((lambda + 2*miu)*N5_y*N8_y + miu*N5_x*N8_x + miu*N5_z*N8_z) * det_J;
		                // Matrix[263] += w_i*w_j*w_k* (lambda*N5_x*N8_y + miu*N5_y*N8_x) * det_J;
		                // Matrix[264] += w_i*w_j*w_k* (lambda*N4_z*N8_y + miu*N4_y*N8_z) * det_J;
		                // Matrix[265] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_y*N8_y + miu*N4_x*N8_x + miu*N4_z*N8_z) * det_J;
		                // Matrix[266] += w_i*w_j*w_k* (lambda*N4_x*N8_y + miu*N4_y*N8_x) * det_J;
		                // Matrix[267] += w_i*w_j*w_k* (lambda*N3_z*N8_y + miu*N3_y*N8_z) * det_J;
		                // Matrix[268] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_y*N8_y + miu*N3_x*N8_x + miu*N3_z*N8_z) * det_J;
		                // Matrix[269] += w_i*w_j*w_k* (lambda*N3_x*N8_y + miu*N3_y*N8_x) * det_J;
		                // Matrix[270] += w_i*w_j*w_k* (lambda*N2_z*N8_y + miu*N2_y*N8_z) * det_J;
		                // Matrix[271] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_y*N8_y + miu*N2_x*N8_x + miu*N2_z*N8_z) * det_J;
		                // Matrix[272] += w_i*w_j*w_k* (lambda*N2_x*N8_y + miu*N2_y*N8_x) * det_J;
		                // Matrix[273] += w_i*w_j*w_k* (lambda*N1_z*N8_y + miu*N1_y*N8_z) * det_J;
		                // Matrix[274] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_y*N8_y + miu*N1_x*N8_x + miu*N1_z*N8_z) * det_J;
		                // Matrix[275] += w_i*w_j*w_k* (lambda*N1_x*N8_y + miu*N1_y*N8_x) * det_J;
		
		                // // L24
		                // Matrix[276] += w_i*w_j*w_k* ((lambda + 2*miu)*N8_z*N8_z + miu*N8_x*N8_x + miu*N8_y*N8_y) * det_J;
		                // Matrix[277] += w_i*w_j*w_k* (lambda*N8_y*N8_z + miu*N8_z*N8_y) * det_J;
		                // Matrix[278] += w_i*w_j*w_k* (lambda*N8_x*N8_z + miu*N8_z*N8_x) * det_J;
		                // Matrix[279] += w_i*w_j*w_k* ((lambda + 2*miu)*N7_z*N8_z + miu*N7_x*N8_x + miu*N7_y*N8_y) * det_J;
		                // Matrix[280] += w_i*w_j*w_k* (lambda*N7_y*N8_z + miu*N7_z*N8_y) * det_J;
		                // Matrix[281] += w_i*w_j*w_k* (lambda*N7_x*N8_z + miu*N7_z*N8_x) * det_J;
		                // Matrix[282] += w_i*w_j*w_k* ((lambda + 2*miu)*N6_z*N8_z + miu*N6_x*N8_x + miu*N6_y*N8_y) * det_J;
		                // Matrix[283] += w_i*w_j*w_k* (lambda*N6_y*N8_z + miu*N6_z*N8_y) * det_J;
		                // Matrix[284] += w_i*w_j*w_k* (lambda*N6_x*N8_z + miu*N6_z*N8_x) * det_J;
		                // Matrix[285] += w_i*w_j*w_k* ((lambda + 2*miu)*N5_z*N8_z + miu*N5_x*N8_x + miu*N5_y*N8_y) * det_J;
		                // Matrix[286] += w_i*w_j*w_k* (lambda*N5_y*N8_z + miu*N5_z*N8_y) * det_J;
		                // Matrix[287] += w_i*w_j*w_k* (lambda*N5_x*N8_z + miu*N5_z*N8_x) * det_J;
		                // Matrix[288] += w_i*w_j*w_k* ((lambda + 2*miu)*N4_z*N8_z + miu*N4_x*N8_x + miu*N4_y*N8_y) * det_J;
		                // Matrix[289] += w_i*w_j*w_k* (lambda*N4_y*N8_z + miu*N4_z*N8_y) * det_J;
		                // Matrix[290] += w_i*w_j*w_k* (lambda*N4_x*N8_z + miu*N4_z*N8_x) * det_J;
		                // Matrix[291] += w_i*w_j*w_k* ((lambda + 2*miu)*N3_z*N8_z + miu*N3_x*N8_x + miu*N3_y*N8_y) * det_J;
		                // Matrix[292] += w_i*w_j*w_k* (lambda*N3_y*N8_z + miu*N3_z*N8_y) * det_J;
		                // Matrix[293] += w_i*w_j*w_k* (lambda*N3_x*N8_z + miu*N3_z*N8_x) * det_J;
		                // Matrix[294] += w_i*w_j*w_k* ((lambda + 2*miu)*N2_z*N8_z + miu*N2_x*N8_x + miu*N2_y*N8_y) * det_J;
		                // Matrix[295] += w_i*w_j*w_k* (lambda*N2_y*N8_z + miu*N2_z*N8_y) * det_J;
		                // Matrix[296] += w_i*w_j*w_k* (lambda*N2_x*N8_z + miu*N2_z*N8_x) * det_J;
		                // Matrix[297] += w_i*w_j*w_k* ((lambda + 2*miu)*N1_z*N8_z + miu*N1_x*N8_x + miu*N1_y*N8_y) * det_J;
		                // Matrix[298] += w_i*w_j*w_k* (lambda*N1_y*N8_z + miu*N1_z*N8_y) * det_J;
		                // Matrix[299] += w_i*w_j*w_k* (lambda*N1_x*N8_z + miu*N1_z*N8_x) * det_J;
            		}
		}
	}
	unsigned s = 0;
    	for (unsigned int ii = 0; ii < 24; ++ii) {
        	for (int jj = ii; jj > -1; --jj) {
            	Matrix[s] = K[ii][jj];
            	++s;
        	}
    	}
}

//  Calculate element stress
void CH8::ElementStress(double* stress, double* Displacement)
{
	C3DMaterial* material_ = dynamic_cast<C3DMaterial*>(ElementMaterial_);

	double E = material_->E;
	double mu = material_->mu;
    double lambda = E * mu/((1 + mu) * (1 - 2 * mu));
    double miu = E/(2*(1 + mu));

	for (unsigned int i = 0; i < 48; i++) {
		stress[i] = 0.0;
	}

    // Rebuild local displacement
    double displacement[24];


    for (unsigned int i = 0; i < 24; i++)
    {
        if (LocationMatrix_[i])
        {
            displacement[i] = Displacement[LocationMatrix_[i] - 1];
        }
        else
        {
            displacement[i] = 0;
        }
    }
    
    double w[2] = {0};
    double gp[2] = {0};
    GaussianQuadrature H8gauss;
    unsigned int ngp = H8gauss.ngp;

    H8gauss.Gauss_Calculate();
    for(unsigned int i = 0; i<ngp; i++)
    {
        w[i] = H8gauss.g_w[i];
        gp[i] = H8gauss.g_p[i];
    }

    
	for (unsigned i = 0; i < 2; i++)
	{
		for (unsigned j = 0; j < 2; j++)
		{
		            for (unsigned k = 0; k < 2; k++)
		            {
		                // Coordinates of the stress superconvergence point
		 
		                double p_i = gp[i];
		                double p_j = gp[j];
		                double p_k = gp[k];
		
		                // Jacobian matrix and det
		                double J_11 = 1.0/8 * ((p_j-1)*(1-p_k)*(nodes_[0]->XYZ[0]) + (1-p_j)*(1-p_k)*(nodes_[1]->XYZ[0]) + (1+p_j)*(1-p_k)*(nodes_[2]->XYZ[0]) + (-p_j-1)*(1-p_k)*(nodes_[3]->XYZ[0]) + (p_j-1)*(1+p_k)*(nodes_[4]->XYZ[0]) + (1-p_j)*(1+p_k)*(nodes_[5]->XYZ[0]) + (1+p_j)*(1+p_k)*(nodes_[6]->XYZ[0]) + (-1-p_j)*(1+p_k)*(nodes_[7]->XYZ[0]));
		                double J_12 = 1.0/8 * ((p_i-1)*(1-p_k)*(nodes_[0]->XYZ[0]) + (1+p_i)*(p_k-1)*(nodes_[1]->XYZ[0]) + (1+p_i)*(1-p_k)*(nodes_[2]->XYZ[0]) + (1-p_i)*(1-p_k)*(nodes_[3]->XYZ[0]) + (p_i-1)*(1+p_k)*(nodes_[4]->XYZ[0]) + (-1-p_i)*(1+p_k)*(nodes_[5]->XYZ[0]) + (1+p_i)*(1+p_k)*(nodes_[6]->XYZ[0]) + (1-p_i)*(1+p_k)*(nodes_[7]->XYZ[0]));
		                double J_13 = 1.0/8 * ((p_i-1)*(1-p_j)*(nodes_[0]->XYZ[0]) + (1+p_i)*(p_j-1)*(nodes_[1]->XYZ[0]) + (-1-p_i)*(1+p_j)*(nodes_[2]->XYZ[0]) + (p_i-1)*(1+p_j)*(nodes_[3]->XYZ[0]) + (1-p_i)*(1-p_j)*(nodes_[4]->XYZ[0]) + (1+p_i)*(1-p_j)*(nodes_[5]->XYZ[0]) + (1+p_i)*(1+p_j)*(nodes_[6]->XYZ[0]) + (1-p_i)*(1+p_j)*(nodes_[7]->XYZ[0]));
		                double J_21 = 1.0/8 * ((p_j-1)*(1-p_k)*(nodes_[0]->XYZ[1]) + (1-p_j)*(1-p_k)*(nodes_[1]->XYZ[1]) + (1+p_j)*(1-p_k)*(nodes_[2]->XYZ[1]) + (-p_j-1)*(1-p_k)*(nodes_[3]->XYZ[1]) + (p_j-1)*(1+p_k)*(nodes_[4]->XYZ[1]) + (1-p_j)*(1+p_k)*(nodes_[5]->XYZ[1]) + (1+p_j)*(1+p_k)*(nodes_[6]->XYZ[1]) + (-1-p_j)*(1+p_k)*(nodes_[7]->XYZ[1]));
		                double J_22 = 1.0/8 * ((p_i-1)*(1-p_k)*(nodes_[0]->XYZ[1]) + (1+p_i)*(p_k-1)*(nodes_[1]->XYZ[1]) + (1+p_i)*(1-p_k)*(nodes_[2]->XYZ[1]) + (1-p_i)*(1-p_k)*(nodes_[3]->XYZ[1]) + (p_i-1)*(1+p_k)*(nodes_[4]->XYZ[1]) + (-1-p_i)*(1+p_k)*(nodes_[5]->XYZ[1]) + (1+p_i)*(1+p_k)*(nodes_[6]->XYZ[1]) + (1-p_i)*(1+p_k)*(nodes_[7]->XYZ[1]));
		                double J_23 = 1.0/8 * ((p_i-1)*(1-p_j)*(nodes_[0]->XYZ[1]) + (1+p_i)*(p_j-1)*(nodes_[1]->XYZ[1]) + (-1-p_i)*(1+p_j)*(nodes_[2]->XYZ[1]) + (p_i-1)*(1+p_j)*(nodes_[3]->XYZ[1]) + (1-p_i)*(1-p_j)*(nodes_[4]->XYZ[1]) + (1+p_i)*(1-p_j)*(nodes_[5]->XYZ[1]) + (1+p_i)*(1+p_j)*(nodes_[6]->XYZ[1]) + (1-p_i)*(1+p_j)*(nodes_[7]->XYZ[1]));
		                double J_31 = 1.0/8 * ((p_j-1)*(1-p_k)*(nodes_[0]->XYZ[2]) + (1-p_j)*(1-p_k)*(nodes_[1]->XYZ[2]) + (1+p_j)*(1-p_k)*(nodes_[2]->XYZ[2]) + (-p_j-1)*(1-p_k)*(nodes_[3]->XYZ[2]) + (p_j-1)*(1+p_k)*(nodes_[4]->XYZ[2]) + (1-p_j)*(1+p_k)*(nodes_[5]->XYZ[2]) + (1+p_j)*(1+p_k)*(nodes_[6]->XYZ[2]) + (-1-p_j)*(1+p_k)*(nodes_[7]->XYZ[2]));
		                double J_32 = 1.0/8 * ((p_i-1)*(1-p_k)*(nodes_[0]->XYZ[2]) + (1+p_i)*(p_k-1)*(nodes_[1]->XYZ[2]) + (1+p_i)*(1-p_k)*(nodes_[2]->XYZ[2]) + (1-p_i)*(1-p_k)*(nodes_[3]->XYZ[2]) + (p_i-1)*(1+p_k)*(nodes_[4]->XYZ[2]) + (-1-p_i)*(1+p_k)*(nodes_[5]->XYZ[2]) + (1+p_i)*(1+p_k)*(nodes_[6]->XYZ[2]) + (1-p_i)*(1+p_k)*(nodes_[7]->XYZ[2]));
		                double J_33 = 1.0/8 * ((p_i-1)*(1-p_j)*(nodes_[0]->XYZ[2]) + (1+p_i)*(p_j-1)*(nodes_[1]->XYZ[2]) + (-1-p_i)*(1+p_j)*(nodes_[2]->XYZ[2]) + (p_i-1)*(1+p_j)*(nodes_[3]->XYZ[2]) + (1-p_i)*(1-p_j)*(nodes_[4]->XYZ[2]) + (1+p_i)*(1-p_j)*(nodes_[5]->XYZ[2]) + (1+p_i)*(1+p_j)*(nodes_[6]->XYZ[2]) + (1-p_i)*(1+p_j)*(nodes_[7]->XYZ[2]));
		                double det_J = J_11*(J_22*J_33-J_23*J_32)-J_12*(J_21*J_33-J_23*J_31)+J_13*(J_21*J_32-J_31*J_22);
		
		                // Inv of J
		                double J_inv_11 = (J_22*J_33-J_23*J_32) / det_J;
		                double J_inv_12 = - (J_21*J_33-J_23*J_31) / det_J;
		                double J_inv_13 = (J_21*J_32-J_31*J_22) / det_J;
		                double J_inv_21 = - (J_12*J_33-J_32*J_13) / det_J;
		                double J_inv_22 = (J_11*J_33-J_13*J_31) / det_J;
		                double J_inv_23 = - (J_11*J_32-J_31*J_12) / det_J;
		                double J_inv_31 = (J_12*J_23-J_22*J_13) / det_J;
		                double J_inv_32 = - (J_11*J_23-J_21*J_13) / det_J;
		                double J_inv_33 = (J_11*J_22-J_21*J_12) / det_J;
		
		                // Grad of shape function
		                double N1_x = 1.0/8 * (J_inv_11 * (p_j-1)*(1-p_k) + J_inv_12 * (p_i-1)*(1-p_k) + J_inv_13 * (p_i-1)*(1-p_j));
		                double N2_x = 1.0/8 * (J_inv_11 * (1-p_j)*(1-p_k) + J_inv_12 * (1+p_i)*(p_k-1) + J_inv_13 * (1+p_i)*(p_j-1));
		                double N3_x = 1.0/8 * (J_inv_11 * (1+p_j)*(1-p_k) + J_inv_12 * (1+p_i)*(1-p_k) + J_inv_13 * (-1-p_i)*(1+p_j));
		                double N4_x = 1.0/8 * (J_inv_11 * (1+p_j)*(p_k-1) + J_inv_12 * (1-p_i)*(1-p_k) + J_inv_13 * (p_i-1)*(1+p_j));
		                double N5_x = 1.0/8 * (J_inv_11 * (p_j-1)*(1+p_k) + J_inv_12 * (p_i-1)*(1+p_k) + J_inv_13 * (1-p_i)*(1-p_j));
		                double N6_x = 1.0/8 * (J_inv_11 * (1-p_j)*(1+p_k) + J_inv_12 * (-1-p_i)*(1+p_k) + J_inv_13 * (1+p_i)*(1-p_j));
		                double N7_x = 1.0/8 * (J_inv_11 * (1+p_j)*(1+p_k) + J_inv_12 * (1+p_i)*(1+p_k) + J_inv_13 * (1+p_i)*(1+p_j));
		                double N8_x = 1.0/8 * (J_inv_11 * (-1-p_j)*(1+p_k) + J_inv_12 * (1-p_i)*(1+p_k) + J_inv_13 * (1-p_i)*(1+p_j));
		                double N1_y = 1.0/8 * (J_inv_21 * (p_j-1)*(1-p_k) + J_inv_22 * (p_i-1)*(1-p_k) + J_inv_23 * (p_i-1)*(1-p_j));
		                double N2_y = 1.0/8 * (J_inv_21 * (1-p_j)*(1-p_k) + J_inv_22 * (1+p_i)*(p_k-1) + J_inv_23 * (1+p_i)*(p_j-1));
		                double N3_y = 1.0/8 * (J_inv_21 * (1+p_j)*(1-p_k) + J_inv_22 * (1+p_i)*(1-p_k) + J_inv_23 * (-1-p_i)*(1+p_j));
		                double N4_y = 1.0/8 * (J_inv_21 * (1+p_j)*(p_k-1) + J_inv_22 * (1-p_i)*(1-p_k) + J_inv_23 * (p_i-1)*(1+p_j));
		                double N5_y = 1.0/8 * (J_inv_21 * (p_j-1)*(1+p_k) + J_inv_22 * (p_i-1)*(1+p_k) + J_inv_23 * (1-p_i)*(1-p_j));
		                double N6_y = 1.0/8 * (J_inv_21 * (1-p_j)*(1+p_k) + J_inv_22 * (-1-p_i)*(1+p_k) + J_inv_23 * (1+p_i)*(1-p_j));
		                double N7_y = 1.0/8 * (J_inv_21 * (1+p_j)*(1+p_k) + J_inv_22 * (1+p_i)*(1+p_k) + J_inv_23 * (1+p_i)*(1+p_j));
		                double N8_y = 1.0/8 * (J_inv_21 * (-1-p_j)*(1+p_k) + J_inv_22 * (1-p_i)*(1+p_k) + J_inv_23 * (1-p_i)*(1+p_j));
		                double N1_z = 1.0/8 * (J_inv_31 * (p_j-1)*(1-p_k) + J_inv_32 * (p_i-1)*(1-p_k) + J_inv_33 * (p_i-1)*(1-p_j));
		                double N2_z = 1.0/8 * (J_inv_31 * (1-p_j)*(1-p_k) + J_inv_32 * (1+p_i)*(p_k-1) + J_inv_33 * (1+p_i)*(p_j-1));
		                double N3_z = 1.0/8 * (J_inv_31 * (1+p_j)*(1-p_k) + J_inv_32 * (1+p_i)*(1-p_k) + J_inv_33 * (-1-p_i)*(1+p_j));
		                double N4_z = 1.0/8 * (J_inv_31 * (1+p_j)*(p_k-1) + J_inv_32 * (1-p_i)*(1-p_k) + J_inv_33 * (p_i-1)*(1+p_j));
		                double N5_z = 1.0/8 * (J_inv_31 * (p_j-1)*(1+p_k) + J_inv_32 * (p_i-1)*(1+p_k) + J_inv_33 * (1-p_i)*(1-p_j));
		                double N6_z = 1.0/8 * (J_inv_31 * (1-p_j)*(1+p_k) + J_inv_32 * (-1-p_i)*(1+p_k) + J_inv_33 * (1+p_i)*(1-p_j));
		                double N7_z = 1.0/8 * (J_inv_31 * (1+p_j)*(1+p_k) + J_inv_32 * (1+p_i)*(1+p_k) + J_inv_33 * (1+p_i)*(1+p_j));
		                double N8_z = 1.0/8 * (J_inv_31 * (-1-p_j)*(1+p_k) + J_inv_32 * (1-p_i)*(1+p_k) + J_inv_33 * (1-p_i)*(1+p_j));
		
		                double nabla_N[8][3] = {{N1_x,N1_y,N1_z},{N2_x,N2_y,N2_z},{N3_x,N3_y,N3_z},{N4_x,N4_y,N4_z},{N5_x,N5_y,N5_z},{N6_x,N6_y,N6_z},{N7_x,N7_y,N7_z},{N8_x,N8_y,N8_z}};
		
		                for (unsigned int m = 0; m < 8; m++)
		                {
		                    stress[24*i+12*j+6*k+0] += (lambda + 2*miu)*nabla_N[m][0]*displacement[3*m+0] + lambda*nabla_N[m][1]*displacement[3*m+1] + lambda*nabla_N[m][2]*displacement[m*i+2];
		                    stress[24*i+12*j+6*k+1] += lambda*nabla_N[m][0]*displacement[3*m+0] + (lambda + 2*miu)*nabla_N[m][1]*displacement[3*m+1] + lambda*nabla_N[m][2]*displacement[3*m+2];
		                    stress[24*i+12*j+6*k+2] += lambda*nabla_N[m][0]*displacement[3*m+0] + lambda*nabla_N[m][1]*displacement[3*m+1] + (lambda + 2*miu)*nabla_N[m][2]*displacement[3*m+2];
		                    stress[24*i+12*j+6*k+3] += miu*nabla_N[m][1]*displacement[3*m+0] + miu*nabla_N[m][0]*displacement[3*m+1];
		                    stress[24*i+12*j+6*k+4] += miu*nabla_N[m][2]*displacement[3*m+1] + miu*nabla_N[m][1]*displacement[3*m+2];
		                    stress[24*i+12*j+6*k+5] += miu*nabla_N[m][2]*displacement[3*m+0] + miu*nabla_N[m][0]*displacement[3*m+2];
		
		                } 
		            }
			
		}
        
	}
}
