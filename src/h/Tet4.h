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
class CTet4 : public CElement
{
public:

	double* Bmat;

public:

//!	Constructor
	CTet4();

//!	Desconstructor
	~CTet4();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Calculate D
	virtual void D_Elastic(C3DMaterial* material_);

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//! B matrix
	virtual void BmatElast(C3DMaterial* material_);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

//body
	virtual void ElementForce(double* Force);

	virtual double Calculate_Determinant(double x1, double x2, double x3, double y1, double y2, 
	double y3, double z1, double z2, double z3);

//!	Calculate element stress
	//virtual double* T3_GAUSS(int ngp);
};
