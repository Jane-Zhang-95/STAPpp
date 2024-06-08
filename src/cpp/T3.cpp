/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/


#include "T3.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CT3::CT3()
{
	NEN_ = 3;	// Each element has 3 nodes 
	nodes_ = new CNode*[NEN_]; 
    
    ND_ = 9;
    LocationMatrix_ = new unsigned int[ND_]; //n of Matrix 

	ElementMaterial_ = nullptr;
}

CT3::~CT3(){
	if(Bmat)
		delete [] Bmat;
}



//	Read element data from stream Input
bool CT3::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{ 
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3;	// Left node number and right node number 

	Input >> N1 >> N2 >> N3 >> MSet; 
    ElementMaterial_ = dynamic_cast<CT3Material*>(MaterialSets) + MSet - 1; 
	nodes_[0] = &NodeList[N1 - 1]; 
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];

	return true;
}

//	Write element data to stream 
void CT3::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(9) << nodes_[2]->NodeNumber 
		   << setw(12) << ElementMaterial_->nset << endl;
}

void CT3::D_Elastic(CT3Material* material_)
{
	unsigned int Mode = material_->Stress_Mode;
	double Elas[9] = {0};
	if( Mode == 0 ) //stress
	{
		double k_stress = (material_->E) / (1.0 - (material_->nu)*(material_->nu) );
		Elas[0] = k_stress * 1.0;
		Elas[1] = k_stress * material_->nu;
		Elas[2] = 0;
		Elas[3] = k_stress * material_->nu;
		Elas[4] = k_stress * 1.0;
		Elas[5] = 0;
		Elas[6] = 0;
		Elas[7] = 0;
		Elas[8] = k_stress * 0.5 * (1.0 - material_->nu);
	}
	else if( Mode == 1 ) //strain
	{
		double k_strain = (material_->E) / (1.0 + material_->nu)*(1.0 - 2.0 * material_->nu) ;
		Elas[0] = k_strain * (1.0 - material_->nu);
		Elas[1] = k_strain * material_->nu;
		Elas[2] = 0;
		Elas[3] = k_strain * material_->nu;
		Elas[4] = k_strain * (1.0 - material_->nu);
		Elas[5] = 0;
		Elas[6] = 0;
		Elas[7] = 0;
		Elas[8] = k_strain * 0.5 * (1.0 - 2.0 * material_->nu);
	}
	for(unsigned int i = 0;i<9;i++){
		(material_->D)[i] = Elas[i];
	}
}

void CT3::BmatElast2D(CT3Material* material_)
{
	double node0[2]={0,0};  
	double node1[2]={0,0};
	double node2[2]={0,0};		
	double T3_b[3]={0,0,0};	// const abc
	double T3_c[3]={0,0,0};	
	
	// calculate const abc
	for (unsigned int i = 0; i < 2; i++){
		node0[i] = nodes_[0]->XYZ[i];
		node1[i] = nodes_[1]->XYZ[i];
		node2[i] = nodes_[2]->XYZ[i];
	}
	
	T3_b[0]=node1[1]-node2[1];
	T3_b[1]=node2[1]-node0[1];
	T3_b[2]=node0[1]-node1[1];
	T3_c[0]=node2[0]-node1[0];
	T3_c[1]=node0[0]-node2[0];
	T3_c[2]=node1[0]-node0[0];
	//element A
	double Area = 0.5 * fabs((node0[0]-node1[0])*(node0[1]-node2[1])-(node0[0]-node2[0])*(node0[1]-node1[1]));
	material_->Area = Area;
	double k = 0.5 / Area;
	double B[18] = {0};
	B[0] = k * (T3_b[0]);
	B[1] = 0;
	B[2] = k * (T3_b[1]);
	B[3] = 0;
	B[4] = k * (T3_b[2]);
	B[5] = 0;
	B[6] = 0;
	B[7] = k * (T3_c[0]);
	B[8] = 0;
	B[9] = k * (T3_c[1]);
	B[10] = 0;
	B[11] = k * (T3_c[2]);
	B[12] = k * (T3_c[0]);
	B[13] = k * (T3_b[0]);
	B[14] = k * (T3_c[1]);
	B[15] = k * (T3_b[1]);
	B[16] = k * (T3_c[2]);
	B[17] = k * (T3_b[2]);

	for(unsigned int i = 0;i<18;i++){
		(*this).Bmat[i] = B[i];
	}
}

void CT3::NmatElast2D(double* Nmat, CT3Material* material_)
{
	double node0[2];  
	double node1[2];
	double node2[2];
	double T3_a[3];
	for (unsigned int i = 0; i < 2; i++){
		node0[i] = nodes_[0]->XYZ[i];
		node1[i] = nodes_[1]->XYZ[i];
		node2[i] = nodes_[2]->XYZ[i];
	}
	T3_a[0]=node1[0]*node2[1]-node2[0]*node1[1];
	T3_a[1]=node2[0]*node0[1]-node0[0]*node2[1];
	T3_a[2]=node0[0]*node1[1]-node1[0]*node0[1];
	
	double Area = material_->Area;
	Nmat[0] = (T3_a[0] + Bmat[12]*node0[1] + Bmat[13]*node0[0])/(2*Area);
	Nmat[1] = (T3_a[1] + Bmat[14]*node1[1] + Bmat[15]*node1[0])/(2*Area);
	Nmat[2] = (T3_a[2] + Bmat[16]*node2[1] + Bmat[17]*node2[0])/(2*Area);

}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CT3::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix()); 

	//	Calculate element stiffness matrix

	CT3Material* material_ = dynamic_cast<CT3Material*>(ElementMaterial_);	// Pointer to material of the element
	material_->D = new double[9];
	(*this).Bmat = new double[18];
	D_Elastic(material_);
	BmatElast2D(material_);
	double* Dmat = material_->D;
	//AtBDB
	double k = (material_->Area) * (material_->h_T3) ;
	double M[36];

	for(unsigned int i=0 ; i<36; i++){
		M[i]=0;}

	for(unsigned int i=0; i<6; i++){//column
		for(unsigned int j=0; j<6; j++){
			for(unsigned int i2 = 0; i2 < 3; i2++){
				for(unsigned int j2 = 0; j2 < 3; j2++){
					M[6*i+j] += k * Bmat[6*i2+i] * Dmat[3*i2+j2] * Bmat[6*j2+j];
				}
			}
		}
	}
	//std::cout<<M[0]<<M[1]<<M[2]<<M[3]<<M[4]<<M[5]<<M[6]<<endl;
	double newM[81];
	for(unsigned int i=0 ; i<81; i++){
		newM[i]=0;}

    unsigned int rowMap[] = {0, 1, 3, 4, 6, 7}; 

    for (unsigned int i = 0; i < 6; i++) {
		for (unsigned int j = 0; j < 6; j++) {
			newM[rowMap[i] * 9 + rowMap[j]] += M[i * 6 + j];
		}
    }
	unsigned int num = 0;
	for(unsigned int j = 0; j < 9; j++){//column
		for(unsigned int i = j+1; i > 0; i--){
			Matrix[num] = newM[9*(i-1)+j];
			num=num+1;
		}
	}
	//std::cout<<Matrix[0]<<Matrix[1]<<Matrix[2]<<Matrix[3]<<Matrix[4]<<Matrix[5]<<Matrix[6]<<Matrix[7]
	//<<Matrix[8]<<Matrix[9]<<Matrix[10]<<Matrix[11]<<Matrix[12]<<Matrix[13]<<Matrix[14]<<Matrix[15]<<Matrix[16]
	//<<endl;
}

// assemble body force
void CT3::ElementForce(double* Force)
{
	clear(Force, 9);

	double bx1 = nodes_[0]->BODY[0];
	double bx2 = nodes_[1]->BODY[0];
	double bx3 = nodes_[2]->BODY[0];
	double by1 = nodes_[0]->BODY[1];
	double by2 = nodes_[1]->BODY[1];
	double by3 = nodes_[2]->BODY[1];

	Force[0] = 2*bx1 + bx2 + bx3;
	Force[1] = 2*by1 + by2 + by3;
	Force[3] = 2*bx2 + bx1 + bx3;
	Force[4] = 2*by2 + by1 + by3;
	Force[6] = 2*bx3 + bx2 + bx1;
	Force[7] = 2*by3 + by2 + by1;
}

//	Calculate element stress
void CT3::ElementStress(double* stress, double* Displacement)
{ //stress and dis
	clear(stress,3);
	CT3Material* material_ = dynamic_cast<CT3Material*>(ElementMaterial_);	// Pointer to material of the element
	Bmat = new double[18];	
	BmatElast2D(material_);
	double B[18] = {0};
	for(unsigned int i = 0;i<18;i++){
		B[i] = ((*this).Bmat)[i];}
	//std::cout<<B[0]<<B[1]<<B[2]<<endl;
	double D[9] ={0};
	for(unsigned int i = 0;i<9;i++){
		D[i]= (material_->D)[i];
	} 
	//std::cout<<D[0]<<D[1]<<endl;
	double s[18] = {0}; //DB
	double news[27] = {0};

	for (unsigned int i = 0; i < 3; i++){
		for (unsigned int j = 0; j < 6; j++){
			for (unsigned int k = 0; k < 3; k++){
				s[6*i+j] += D[3*i+k] * B[6*k+j];
			}
		}
	} //DB
	//std::cout<<s[0]<<s[1]<<s[2]<<s[3]<<s[4]<<s[5]<<s[6]<<endl;
	//stress = new double[3];
	unsigned int rowMap[] = {0, 1, 3, 4, 6, 7}; 

    for (unsigned int i = 0; i < 3; i++) {
		for (unsigned int j = 0; j < 6; j++) {
			news[9*i + rowMap[j]] = s[i * 6 + j];
		}
    }
	//std::cout<<news[0]<<news[1]<<news[2]<<news[3]<<news[4]<<news[5]<<endl;
	for (unsigned int i = 0; i < 3; i++)
	{
		for( int j=0; j<9; j++){
			if (LocationMatrix_[j]){
			stress[i] += news[9*i + j] * Displacement[LocationMatrix_[j]-1];
			//std::cout<<news[3]<<news[4]<<news[5]<<news[6]<<news[7]<<news[8]<<endl;
			//std::cout<<Displacement[LocationMatrix_[j]-1]<<endl;
			}
		}
	} //DBd
}



