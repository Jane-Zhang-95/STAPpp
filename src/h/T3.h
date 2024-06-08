/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! Bar element class
class CT3 : public CElement
{
public:

	double* Bmat;

public:

//!	Constructor
	CT3();

//!	Desconstructor
	~CT3();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Calculate D
	virtual void D_Elastic(CT3Material* material_);

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//! B matrix
	virtual void BmatElast2D(CT3Material* material_);

//! N matrix
	virtual void NmatElast2D(double* Nmat, CT3Material* material_);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

//body
	virtual void ElementForce(double* Force);

//!	Calculate element stress
	//virtual double* T3_GAUSS(int ngp);
};
