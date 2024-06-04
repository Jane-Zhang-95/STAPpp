#include "Q4.h"
#include "Gauss.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//  Constructor
CQ4::CQ4()
{
	NEN_ = 4;
	nodes_ = new CNode * [NEN_];

	ND_ = 12;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//  Desconstructor
CQ4::~CQ4()
{
}

//  Read element data from stream Input
bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int Mset;
	unsigned int N1, N2, N3, N4;

	Input >> N1 >> N2 >> N3 >> N4 >> Mset;
	ElementMaterial_ = dynamic_cast<C2DMaterial*>(MaterialSets) + Mset - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//  Write element data to stream
void CQ4::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(9) << nodes_[2]->NodeNumber
		   << setw(9) << nodes_[3]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

//  Calculate element stiffness matrix
//  Upper triangular matrix stored as an array column by colum starting from the diagonal element
void CQ4::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	C2DMaterial* material_ = dynamic_cast<C2DMaterial*>(ElementMaterial_);

	double E = material_->E;
	double v = material_->mu;
	bool plane_stress = material_->plane_stress;

	if (!plane_stress)	// plane strain condition
	{
		E = E / (1 - v * v);
		v = v / (1 - v);
	}

	for (unsigned int i = 0; i < 2; i++) {
		for (unsigned int j = 0; j < 2; j++) {
			// Gaussian quadrature weights and points
			double W_i = GaussianQuadrature::weights_2p[i];
			double W_j = GaussianQuadrature::weights_2p[j];
			double P_i = GaussianQuadrature::points_2p[i];
			double P_j = GaussianQuadrature::points_2p[j];

			// Jacobian matrix and det
			double J_11 = 1.0 / 4 * ((P_j - 1) * (nodes_[0]->XYZ[0]) + (1 - P_j) * (nodes_[1]->XYZ[0]) + (1 + P_j) * (nodes_[2]->XYZ[0]) + (-P_j - 1) * (nodes_[3]->XYZ[0]));
			double J_12 = 1.0 / 4 * ((P_j - 1) * (nodes_[0]->XYZ[1]) + (1 - P_j) * (nodes_[1]->XYZ[1]) + (1 + P_j) * (nodes_[2]->XYZ[1]) + (-P_j - 1) * (nodes_[3]->XYZ[1]));
			double J_21 = 1.0 / 4 * ((P_i - 1) * (nodes_[0]->XYZ[0]) + (-P_i - 1) * (nodes_[1]->XYZ[0]) + (1 + P_i) * (nodes_[2]->XYZ[0]) + (1 - P_i) * (nodes_[3]->XYZ[0]));
			double J_22 = 1.0 / 4 * ((P_i - 1) * (nodes_[0]->XYZ[1]) + (-P_i - 1) * (nodes_[1]->XYZ[1]) + (1 + P_i) * (nodes_[2]->XYZ[1]) + (1 - P_i) * (nodes_[3]->XYZ[1]));
			double det_J = J_11 * J_22 - J_12 * J_21;

			// Inv of J
			double J_inv_11 = J_22 / det_J;
			double J_inv_12 = -J_12 / det_J;
			double J_inv_21 = -J_21 / det_J;
			double J_inv_22 = J_11 / det_J;

			// Grad of shape function
			double N1_x = 1.0 / 4 * (J_inv_11 * (P_j - 1) + J_inv_12 * (P_i - 1));
			double N2_x = 1.0 / 4 * (J_inv_11 * (1 - P_j) + J_inv_12 * (-P_i - 1));
			double N3_x = 1.0 / 4 * (J_inv_11 * (1 + P_j) + J_inv_12 * (1 + P_i));
			double N4_x = 1.0 / 4 * (J_inv_11 * (-P_j - 1) + J_inv_12 * (1 - P_i));
			double N1_y = 1.0 / 4 * (J_inv_21 * (P_j - 1) + J_inv_22 * (P_i - 1));
			double N2_y = 1.0 / 4 * (J_inv_21 * (1 - P_j) + J_inv_22 * (-P_i - 1));
			double N3_y = 1.0 / 4 * (J_inv_21 * (1 + P_j) + J_inv_22 * (1 + P_i));
			double N4_y = 1.0 / 4 * (J_inv_21 * (-P_j - 1) + J_inv_22 * (1 - P_i));

			Matrix[0] += W_i * W_j * (N1_x * N1_x * E / (1 - v * v) + N1_y * N1_y * E / 2 / (1 + v)) * det_J;

			Matrix[1] += W_i * W_j * (N1_y * N1_y * E / (1 - v * v) + N1_x * N1_x * E / 2 / (1 + v)) * det_J;
			Matrix[2] += W_i * W_j * (N1_x * N1_y * E / 2 / (1 - v)) * det_J;

			Matrix[6] += W_i * W_j * (N2_x * N2_x * E / (1 - v * v) + N2_y * N2_y * E / 2 / (1 + v)) * det_J;
			Matrix[8] += W_i * W_j * (N1_y * N2_x * E * v / (1 - v * v) + N1_x * N2_y * E / 2 / (1 + v)) * det_J;
			Matrix[9] += W_i * W_j * (N1_x * N2_x * E / (1 - v * v) + N1_y * N2_y * E / 2 / (1 + v)) * det_J;

			Matrix[10] += W_i * W_j * (N2_y * N2_y * E / (1 - v * v) + N2_x * N2_x * E / 2 / (1 + v)) * det_J;
			Matrix[11] += W_i * W_j * (N2_x * N2_y * E / 2 / (1 - v)) * det_J;
			Matrix[13] += W_i * W_j * (N1_y * N2_y * E / (1 - v * v) + N1_x * N2_x * E / 2 / (1 + v)) * det_J;
			Matrix[14] += W_i * W_j * (N1_x * N2_y * E * v / (1 - v * v) + N1_y * N2_x * E / 2 / (1 + v)) *det_J;

			Matrix[21] += W_i * W_j * (N3_x * N3_x * E / (1 - v * v) + N3_y * N3_y * E / 2 / (1 + v)) * det_J;
			Matrix[23] += W_i * W_j * (N2_y * N3_x * E * v / (1 - v * v) + N2_x * N3_y * E / 2 / (1 + v)) * det_J;
			Matrix[24] += W_i * W_j * (N2_x * N3_x * E / (1 - v * v) + N2_y * N3_y * E / 2 / (1 + v)) * det_J;
			Matrix[26] += W_i * W_j * (N1_y * N3_x * E * v / (1 - v * v) + N1_x * N3_y * E / 2 / (1 + v)) * det_J;
			Matrix[27] += W_i * W_j * (N1_x * N3_x * E / (1 - v * v) + N1_y * N3_y * E / 2 / (1 + v)) * det_J;
			
			Matrix[28] += W_i * W_j * (N3_y * N3_y * E / (1 - v * v) + N3_x * N3_x * E / 2 / (1 + v)) * det_J;
			Matrix[29] += W_i * W_j * (N3_x * N3_y * E / 2 / (1 - v)) * det_J;
			Matrix[31] += W_i * W_j * (N2_y * N3_y * E / (1 - v * v) + N2_x * N3_x * E / 2 / (1 + v)) * det_J;
			Matrix[32] += W_i * W_j * (N2_x * N3_y * E * v / (1 - v * v) + N2_y * N3_x * E / 2 / (1 + v)) * det_J;
			Matrix[34] += W_i * W_j * (N1_y * N3_y * E / (1 - v * v) + N1_x * N3_x * E / 2 / (1 + v)) * det_J;
			Matrix[35] += W_i * W_j * (N1_x * N3_y * E * v / (1 - v * v) + N1_y * N3_x * E / 2 / (1 + v)) * det_J;
			
			Matrix[45] += W_i * W_j * (N4_x * N4_x * E / (1 - v * v) + N4_y * N4_y * E / 2 / (1 + v)) * det_J;
			Matrix[47] += W_i * W_j * (N3_y * N4_x * E * v / (1 - v * v) + N3_x * N4_y * E / 2 / (1 + v)) * det_J;
			Matrix[48] += W_i * W_j * (N3_x * N4_x * E / (1 - v * v) + N3_y * N4_y * E / 2 / (1 + v)) * det_J;
			Matrix[50] += W_i * W_j * (N2_y * N4_x * E * v / (1 - v * v) + N2_x * N4_y * E / 2 / (1 + v)) * det_J;
			Matrix[51] += W_i * W_j * (N2_x * N4_x * E / (1 - v * v) + N2_y * N4_y * E / 2 / (1 + v)) * det_J;
			Matrix[53] += W_i * W_j * (N1_y * N4_x * E * v / (1 - v * v) + N1_x * N4_y * E / 2 / (1 + v)) * det_J;
			Matrix[54] += W_i * W_j * (N1_x * N4_x * E / (1 - v * v) + N1_y * N4_y * E / 2 / (1 + v)) * det_J;

			Matrix[55] += W_i * W_j * (N4_y * N4_y * E / (1 - v * v) + N4_x * N4_x * E / 2 / (1 + v)) * det_J;
			Matrix[56] += W_i * W_j * (N4_x * N4_y * E / 2 / (1 - v)) * det_J;
			Matrix[58] += W_i * W_j * (N3_y * N4_y * E / (1 - v * v) + N3_x * N4_x * E / 2 / (1 + v)) * det_J;
			Matrix[59] += W_i * W_j * (N3_x * N4_y * E * v / (1 - v * v) + N3_y * N4_x * E / 2 / (1 + v)) * det_J;
			Matrix[61] += W_i * W_j * (N2_y * N4_y * E / (1 - v * v) + N2_x * N4_x * E / 2 / (1 + v)) * det_J;
			Matrix[62] += W_i * W_j * (N2_x * N4_y * E * v / (1 - v * v) + N2_y * N4_x * E / 2 / (1 + v)) * det_J;
			Matrix[64] += W_i * W_j * (N1_y * N4_y * E / (1 - v * v) + N1_x * N4_x * E / 2 / (1 + v)) * det_J;
			Matrix[65] += W_i * W_j * (N1_x * N4_y * E * v / (1 - v * v) + N1_y * N4_x * E / 2 / (1 + v)) * det_J;
		}
	}
	
}

//  Calculate element stress
void CQ4::ElementStress(double* stress, double* Displacement)
{
	C2DMaterial* material_ = dynamic_cast<C2DMaterial*>(ElementMaterial_);

	double E = material_->E;
	double v = material_->mu;
	bool plane_stress = material_->plane_stress;

	if (!plane_stress)	// plane strain condition
	{
		E = E / (1 - v * v);
		v = v / (1 - v);
	}

	for (unsigned int i = 0; i < 24; i++) {
		stress[i] = 0.0;
	}

	for (unsigned i = 0; i < 2; i++)
	{
		for (unsigned j = 0; j < 2; j++)
		{
			// Coordinates of the stress superconvergence point
			double P_i = GaussianQuadrature::points_2p[i];
			double P_j = GaussianQuadrature::points_2p[j];

			// Jacobian matrix and det
			double J_11 = 1.0 / 4 * ((P_j - 1) * (nodes_[0]->XYZ[0]) + (1 - P_j) * (nodes_[1]->XYZ[0]) + (1 + P_j) * (nodes_[2]->XYZ[0]) + (-P_j - 1) * (nodes_[3]->XYZ[0]));
			double J_12 = 1.0 / 4 * ((P_j - 1) * (nodes_[0]->XYZ[1]) + (1 - P_j) * (nodes_[1]->XYZ[1]) + (1 + P_j) * (nodes_[2]->XYZ[1]) + (-P_j - 1) * (nodes_[3]->XYZ[1]));
			double J_21 = 1.0 / 4 * ((P_i - 1) * (nodes_[0]->XYZ[0]) + (-P_i - 1) * (nodes_[1]->XYZ[0]) + (1 + P_i) * (nodes_[2]->XYZ[0]) + (1 - P_i) * (nodes_[3]->XYZ[0]));
			double J_22 = 1.0 / 4 * ((P_i - 1) * (nodes_[0]->XYZ[1]) + (-P_i - 1) * (nodes_[1]->XYZ[1]) + (1 + P_i) * (nodes_[2]->XYZ[1]) + (1 - P_i) * (nodes_[3]->XYZ[1]));
			double det_J = J_11 * J_22 - J_12 * J_21;

			// Inv of J
			double J_inv_11 = J_22 / det_J;
			double J_inv_12 = -J_12 / det_J;
			double J_inv_21 = -J_21 / det_J;
			double J_inv_22 = J_11 / det_J;

			// Grad of shape function
			double N1_x = 1.0 / 4 * (J_inv_11 * (P_j - 1) + J_inv_12 * (P_i - 1));
			double N2_x = 1.0 / 4 * (J_inv_11 * (1 - P_j) + J_inv_12 * (-P_i - 1));
			double N3_x = 1.0 / 4 * (J_inv_11 * (1 + P_j) + J_inv_12 * (1 + P_i));
			double N4_x = 1.0 / 4 * (J_inv_11 * (-P_j - 1) + J_inv_12 * (1 - P_i));
			double N1_y = 1.0 / 4 * (J_inv_21 * (P_j - 1) + J_inv_22 * (P_i - 1));
			double N2_y = 1.0 / 4 * (J_inv_21 * (1 - P_j) + J_inv_22 * (-P_i - 1));
			double N3_y = 1.0 / 4 * (J_inv_21 * (1 + P_j) + J_inv_22 * (1 + P_i));
			double N4_y = 1.0 / 4 * (J_inv_21 * (-P_j - 1) + J_inv_22 * (1 - P_i));

			// Rebuild local displacement
			double displacement[12];

			for (unsigned int i = 0; i < 12; i++)
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

			stress[(i * 2 + j) * 6 + 0] = E / (1 - v * v) *
				(N1_x * displacement[0] + N2_x * displacement[3] + N3_x * displacement[6] + N4_x * displacement[9])
				+ E * v / (1 - v * v) *
				(N1_y * displacement[1] + N2_y * displacement[4] + N3_y * displacement[7] + N4_y * displacement[10]);

			stress[(i * 2 + j) * 6 + 1] = E * v / (1 - v * v) *
				(N1_x * displacement[0] + N2_x * displacement[3] + N3_x * displacement[6] + N4_x * displacement[9])
				+ E / (1 - v * v) *
				(N1_y * displacement[1] + N2_y * displacement[4] + N3_y * displacement[7] + N4_y * displacement[10]);

			if (!plane_stress)
			{
				stress[(i * 2 + j) * 6 + 2] = v / (1 + v) * (stress[0] + stress[1]);
			}
			else
			{
				stress[(i * 2 + j) * 6 + 2] = 0;
			}

			stress[(i * 2 + j) * 6 + 3] = E / 2 / (1 + v) *
				(N1_y * displacement[0] + N2_y * displacement[3] + N3_y * displacement[6] + N4_y * displacement[9]
					+ N1_x * displacement[1] + N2_x * displacement[4] + N3_x * displacement[7] + N4_x * displacement[10]);
		}
	}
}