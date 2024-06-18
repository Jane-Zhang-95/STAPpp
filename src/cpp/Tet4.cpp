/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/


#include "Tet4.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CTet4::CTet4()
{
	NEN_ = 4;	// Each element has 4 nodes 
	nodes_ = new CNode*[NEN_]; 
    
    ND_ = 12;
    LocationMatrix_ = new unsigned int[ND_]; //n of Matrix 

	ElementMaterial_ = nullptr;
}

CTet4::~CTet4(){
	if(Bmat)
		delete [] Bmat;
}



//	Read element data from stream Input
bool CTet4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{ 
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4;	// Left node number and right node number 

	Input >> N1 >> N2 >> N3 >> N4 >> MSet; 
    ElementMaterial_ = dynamic_cast<C3DMaterial*>(MaterialSets) + MSet - 1; 
	nodes_[0] = &NodeList[N1 - 1]; 
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//	Write element data to stream 
void CTet4::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber 
		   << setw(9) << nodes_[2]->NodeNumber 
		   << setw(9) << nodes_[3]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

void CTet4::D_Elastic(C3DMaterial* material_)
{
	double Elas[36] = {0};
	double k_stress = (material_->E) / ((1.0 + material_->mu)*(1.0 - 2*(material_->mu)) );

	Elas[0] = k_stress * (1.0-(material_->mu));
	Elas[1] = k_stress * material_->mu;
	Elas[2] = k_stress * material_->mu;

	Elas[6] = k_stress * material_->mu;
	Elas[7] = k_stress * (1.0-(material_->mu));
	Elas[8] = k_stress * material_->mu;

	Elas[12] = k_stress * material_->mu;
	Elas[13] = k_stress * material_->mu;
	Elas[14] = k_stress * (1.0-(material_->mu));

	Elas[21] = k_stress * (1.0-2*(material_->mu))/2;
	Elas[28] = k_stress * (1.0-2*(material_->mu))/2;
	Elas[35] = k_stress * (1.0-2*(material_->mu))/2;

	//std::cout<<Elas[0]<<" "<<Elas[1]<<" "<<Elas[2]<<" "<<Elas[6]<<" "<<
	//Elas[7]<<" "<<Elas[8]<<" "<<Elas[12]<<" "<<Elas[13]<<" "<<Elas[14]<<" "<<Elas[21]<<" "<<
	//Elas[28]<<" "<<Elas[35]<<" "<<endl;

	for(unsigned int i = 0;i<36;i++){
		(material_->D)[i] = Elas[i];
	}
}

void CTet4::BmatElast(C3DMaterial* material_)
{
	double node0[3]={0,0,0};  
	double node1[3]={0,0,0};
	double node2[3]={0,0,0};
	double node3[3]={0,0,0};		
	double Tet4_b[4]={0,0,0,0};	// const abcd
	double Tet4_c[4]={0,0,0,0};	
	double Tet4_d[4]={0,0,0,0};
	double Tet4_a[4]={0,0,0,0};
	
	// calculate const abc
	for (unsigned int i = 0; i < 3; i++){
		node0[i] = nodes_[0]->XYZ[i];
		node1[i] = nodes_[1]->XYZ[i];
		node2[i] = nodes_[2]->XYZ[i];
		node3[i] = nodes_[3]->XYZ[i];
	}

	Tet4_a[0]= Calculate_Determinant(node1[0],node2[0],node3[0],node1[1],node2[1],node3[1],node1[2],node2[2],node3[2]);
	Tet4_b[0]= -Calculate_Determinant(1,1,1,node1[1],node2[1],node3[1],node1[2],node2[2],node3[2]);
	Tet4_c[0]= Calculate_Determinant(1,1,1,node1[0],node2[0],node3[0],node1[2],node2[2],node3[2]);
	Tet4_d[0]= -Calculate_Determinant(1,1,1,node1[0],node2[0],node3[0],node1[1],node2[1],node3[1]);
	Tet4_a[1]= -Calculate_Determinant(node2[0],node3[0],node0[0],node2[1],node3[1],node0[1],node2[2],node3[2],node0[2]);
	Tet4_b[1]= Calculate_Determinant(1,1,1,node2[1],node3[1],node0[1],node2[2],node3[2],node0[2]);
	Tet4_c[1]= -Calculate_Determinant(1,1,1,node2[0],node3[0],node0[0],node2[2],node3[2],node0[2]);
	Tet4_d[1]= Calculate_Determinant(1,1,1,node2[0],node3[0],node0[0],node2[1],node3[1],node0[1]);
	Tet4_a[2]= Calculate_Determinant(node3[0],node0[0],node1[0],node3[1],node0[1],node1[1],node3[2],node0[2],node1[2]);
	Tet4_b[2]= -Calculate_Determinant(1,1,1,node3[1],node0[1],node1[1],node3[2],node0[2],node1[2]);
	Tet4_c[2]= Calculate_Determinant(1,1,1,node3[0],node0[0],node1[0],node3[2],node0[2],node1[2]);
	Tet4_d[2]= -Calculate_Determinant(1,1,1,node3[0],node0[0],node1[0],node3[1],node0[1],node1[1]);
	Tet4_a[3]= -Calculate_Determinant(node0[0],node1[0],node2[0],node0[1],node1[1],node2[1],node0[2],node1[2],node2[2]);
	Tet4_b[3]= Calculate_Determinant(1,1,1,node0[1],node1[1],node2[1],node0[2],node1[2],node2[2]);
	Tet4_c[3]= -Calculate_Determinant(1,1,1,node0[0],node1[0],node2[0],node0[2],node1[2],node2[2]);
	Tet4_d[3]= Calculate_Determinant(1,1,1,node0[0],node1[0],node2[0],node0[1],node1[1],node2[1]);

	//std:: cout<<" "<<Tet4_b[0]<<" "<<Tet4_b[1]<<" "<<Tet4_b[2]<<" "<<Tet4_b[3]<<" "<<
	//Tet4_c[0]<<" "<<Tet4_c[1]<<" "<<Tet4_c[2]<<" "<<Tet4_c[3]<<" "<<
	//Tet4_d[0]<<" "<<Tet4_d[1]<<" "<<Tet4_d[2]<<" "<<Tet4_d[3]<<" "<<endl;
	
	//element V
	double V = (Tet4_a[0]+Tet4_a[1]+Tet4_a[2]+Tet4_a[3])/6.0;
	//std::cout<<Tet4_a[0]<<Tet4_a[1]<<Tet4_a[2]<<Tet4_a[3]<<V<<endl;
	//std::cout<<V<<endl;
	material_->Area = V;
	double k = 1 / (6*V);
	double B[72] = {0};
	for(unsigned int i = 0; i<4 ; i++){
		B[3*i] = k* Tet4_b[i] ;
		B[3*i+1] = 0 ;
		B[3*i+2] = 0 ;
		B[12+3*i] = 0 ;
		B[12+3*i+1] = k* Tet4_c[i] ;
		B[12+3*i+2] = 0 ;
		B[24+3*i] = 0 ;
		B[24+3*i+1] = 0 ;
		B[24+3*i+2] = k* Tet4_d[i] ;
		B[36+3*i] = k* Tet4_c[i] ;
		B[36+3*i+1] = k* Tet4_b[i] ;
		B[36+3*i+2] = 0 ;
		B[48+3*i] = k* Tet4_d[i] ;
		B[48+3*i+1] = 0 ;
		B[48+3*i+2] = k* Tet4_b[i] ;
		B[60+3*i] = 0 ;
		B[60+3*i+1] = k* Tet4_d[i] ;
		B[60+3*i+2] = k* Tet4_c[i] ;
	}

	for(unsigned int i = 0;i<72;i++){
		(*this).Bmat[i] = B[i];
	}
}


//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CTet4::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix()); 

	//	Calculate element stiffness matrix

	C3DMaterial* material_ = dynamic_cast<C3DMaterial*>(ElementMaterial_);	// Pointer to material of the element
	material_->D = new double[36];
	(*this).Bmat = new double[72];
	D_Elastic(material_);
	BmatElast(material_);
	double* Dmat = material_->D;
	//AtBDB
	double k = material_->Area;
	//std::cout<<k<<endl;
	double M[144] = {0};
	//std::cout<<k<<Dmat[0]<<endl;
	//std::cout<<Bmat[0]<<" "<<Bmat[1]<<" "<<Bmat[2]<<" "<<Bmat[3]<<endl;

	for(unsigned int i=0; i<12; i++){//column
		for(unsigned int j=0; j<12; j++){
			for(unsigned int i2 = 0; i2 < 6; i2++){
				for(unsigned int j2 = 0; j2 < 6; j2++){
					M[12*i+j] += k * Bmat[12*i2+i] * Dmat[6*i2+j2] * Bmat[12*j2+j];
				}
			}
		}
	}
	//std::cout<<M[0]<<M[1]<<M[2]<<M[3]<<M[4]<<M[5]<<M[6]<<endl;
	//std::cout<<M[7]<<M[8]<<M[9]<<M[10]<<M[11]<<endl;
	//std::cout<<M[13*11]<<M[130]<<M[13*9]<<endl;

	unsigned int num = 0;
	for(unsigned int j = 0; j < 12; j++){//column
		for(unsigned int i = j+1; i > 0; i--){
			Matrix[num] = M[12*(i-1)+j];
			num=num+1;
		}
	}
	//std::cout<<Matrix[45]<<Matrix[46]<<Matrix[47]<<Matrix[48]<<Matrix[49]<<Matrix[50]<<Matrix[51]<<Matrix[52]
	//<<Matrix[53]<<Matrix[54]<<Matrix[55]<<Matrix[56]<<Matrix[57]<<" "<<Matrix[66]<<Matrix[67]<<Matrix[68]<<Matrix[69]<<Matrix[70]
	//<<Matrix[18]<<Matrix[19]<<Matrix[20]<<Matrix[21]<<Matrix[22]<<Matrix[23]<<Matrix[24]<<" "<<Matrix[25]<<Matrix[26]<<Matrix[27]
	//<<endl;
	//std::cout<<Matrix[0]<<Matrix[1]<<Matrix[3]<<Matrix[6]<<Matrix[10]<<Matrix[15]<<Matrix[21]<<
	//" "<<Matrix[28]<<" " <<Matrix[36]<<" "<<Matrix[45]<<" "<<Matrix[66]<<" "<<Matrix[55]<<endl;
	delete[] Dmat;
}

// assemble body force
void CTet4::ElementForce(double* Force)
{
	clear(Force, 12);

	double bx1 = nodes_[0]->BODY[0];
	double bx2 = nodes_[1]->BODY[0];
	double bx3 = nodes_[2]->BODY[0];
	double bx4 = nodes_[3]->BODY[0];
	double by1 = nodes_[0]->BODY[1];
	double by2 = nodes_[1]->BODY[1];
	double by3 = nodes_[2]->BODY[1];
	double by4 = nodes_[3]->BODY[1];
	double bz1 = nodes_[0]->BODY[2];
	double bz2 = nodes_[1]->BODY[2];
	double bz3 = nodes_[2]->BODY[2];
	double bz4 = nodes_[3]->BODY[2];

	Force[0] = 2*bx1 + bx2 + bx3 + bx4;
	Force[1] = 2*by1 + by2 + by3 + by4;
	Force[2] = 2*bz1 + bz2 + bz3 + bz4;
	Force[3] = 2*bx2 + bx1 + bx3 + bx4;
	Force[4] = 2*by2 + by1 + by3 + by4;
	Force[5] = 2*bz2 + bz1 + bz3 + bz4;
	Force[6] = 2*bx3 + bx2 + bx1 + bx4;
	Force[7] = 2*by3 + by2 + by1 + by4;
	Force[8] = 2*bz3 + bz2 + bz1 + bz4;
	Force[9] = 2*bx4 + bx2 + bx3 + bx1;
	Force[10] = 2*by4 + by2 + by3 + by1;
	Force[11] = 2*bz4 + bz2 + bz3 + bz1;
}


//	Calculate element stress
void CTet4::ElementStress(double* stress, double* Displacement)
{ //stress and dis
	clear(stress,6);
	C3DMaterial* material_ = dynamic_cast<C3DMaterial*>(ElementMaterial_);	// Pointer to material of the element
	Bmat = new double[72];	
	BmatElast(material_);
	D_Elastic(material_);
	double B[72] = {0};
	for(unsigned int i = 0;i<72;i++){
		B[i] = ((*this).Bmat)[i];}
	//std::cout<<B[0]<<B[1]<<B[2]<<endl;
	double D[36] ={0};
	//std::cout<<D[0]<<endl;
	for(unsigned int i = 0;i<36;i++){
		D[i]= (material_->D)[i];
	} 
	//std::cout<<D[0]<<D[1]<<endl;
	double s[72] = {0}; //DB

	for (unsigned int i = 0; i < 6; i++){
		for (unsigned int j = 0; j < 12; j++){
			for (unsigned int k = 0; k < 6; k++){
				s[12*i+j] += D[6*i+k] * B[12*k+j];
			}
		}
	} //DB
	//std::cout<<s[0]<<s[1]<<s[2]<<s[3]<<s[4]<<s[5]<<s[6]<<endl;
	//stress = new double[3];
	//std::cout<<news[0]<<news[1]<<news[2]<<news[3]<<news[4]<<news[5]<<endl;
	for (unsigned int i = 0; i < 6; i++)
	{
		for( int j=0; j<12; j++){
			if (LocationMatrix_[j]){
			stress[i] += s[12*i + j] * Displacement[LocationMatrix_[j]-1];
			//std::cout<<news[3]<<news[4]<<news[5]<<news[6]<<news[7]<<news[8]<<endl;
			//std::cout<<Displacement[LocationMatrix_[j]-1]<<endl;
			}
		}
	} //DBd
}

double CTet4::Calculate_Determinant(double x1, double x2, double x3, double y1, double y2, 
	double y3, double z1, double z2, double z3){

	double det = x1*y2*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1 - y3*z2*x1 - z3*y1*x2;
	return det;
}