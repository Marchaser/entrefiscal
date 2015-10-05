#include "mex.h"
#include <cstring>
#include <cmath>
// #include "mkl.h"
#ifdef USE_OPENMP
#include "omp.h"
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#define _OccPolicy prhs[0]
#define _ApManagerCell prhs[1]
#define _ApManagerLeftShare prhs[2]
#define _EpsilonTrans prhs[3]
#define _ZetaTrans prhs[4]

	double* occPolicy = mxGetPr(_OccPolicy);
	int* apManagerCell = (int*)mxGetData(_ApManagerCell);
	double* apManagerLeftShare = mxGetPr(_ApManagerLeftShare);
	double* epsilonTrans = mxGetPr(_EpsilonTrans);
	double* zetaTrans = mxGetPr(_ZetaTrans);

	const int* distDim = mxGetDimensions(_OccPolicy);
	int epsilonPts = distDim[0];
	int zetaPts = distDim[1];
	int aPts = distDim[2];
	plhs[0] = mxCreateNumericArray(3, distDim, mxDOUBLE_CLASS, mxREAL);
	double* survival = mxGetPr(plhs[0]);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for (int ie = 0; ie < epsilonPts; ++ie)
	{
		for (int izeta = 0; izeta < zetaPts; ++izeta)
		{
			for (int ia = 0; ia < aPts; ++ia)
			{
#define Mat(data,ie,izeta,ia) data[(ia)*zetaPts*epsilonPts + (izeta)*epsilonPts + (ie)]
				int iManagerCell = Mat(apManagerCell, ie, izeta, ia);
				double iManagerLeftShare = Mat(apManagerLeftShare, ie, izeta, ia);
#define EpsilonTransMat(ie,ieprime) epsilonTrans[(ie)+(ieprime)*epsilonPts]
#define ZetaTransMat(izeta,izetaprime) zetaTrans[(izeta)+(izetaprime)*zetaPts]
				for (int ieprime = 0; ieprime < epsilonPts; ++ieprime)
				{
					for (int izetaprime = 0; izetaprime < zetaPts; ++izetaprime)
					{
						Mat(survival, ie, izeta, ia) +=
							(Mat(occPolicy, ieprime, izetaprime, iManagerCell)*iManagerLeftShare
							+ Mat(occPolicy, ieprime, izetaprime, iManagerCell + 1)*(1 - iManagerLeftShare))
							* EpsilonTransMat(ie, ieprime)*ZetaTransMat(izeta, izetaprime);
					}
				}
			}
		}
	}
}

