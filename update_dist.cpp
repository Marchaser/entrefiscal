#include "mex.h"
#include <cstring>
#include <cmath>
// #include "mkl.h"
#ifdef USE_OPENMP
#include "omp.h"
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#define _Dist prhs[0]
#define _Occ prhs[1]
#define _ApManagerCell prhs[2]
#define _ApWorkerCell prhs[3]
#define _ApManagerLeftShare prhs[4]
#define _ApWorkerLeftShare prhs[5]
#define _EpsilonTrans prhs[6]
#define _ZetaTrans prhs[7]

	double* dist = mxGetPr(_Dist);
	double* occ = mxGetPr(_Occ);
	int* apManagerCell = (int*)mxGetData(_ApManagerCell);
	int* apWorkerCell = (int*)mxGetData(_ApWorkerCell);
	double* apManagerLeftShare = (double*)mxGetPr(_ApManagerLeftShare);
	double* apWorkerLeftShare = (double*)mxGetPr(_ApWorkerLeftShare);
	double* epsilonTrans = (double*)mxGetPr(_EpsilonTrans);
	double* zetaTrans = (double*)mxGetPr(_ZetaTrans);

	const int* distDim = mxGetDimensions(_Dist);
	int epsilonPts = distDim[0];
	int zetaPts = distDim[1];
	int aPts = distDim[2];
	plhs[0] = mxCreateNumericArray(3, distDim, mxDOUBLE_CLASS, mxREAL);
	double* tDist = mxGetPr(plhs[0]);
	memset(tDist, 0, sizeof(double)*epsilonPts*zetaPts*aPts);
	double* tDistDistributed = (double*)malloc(sizeof(double)*epsilonPts*epsilonPts*zetaPts*aPts);
	memset(tDistDistributed, 0, sizeof(double)*epsilonPts*epsilonPts*zetaPts*aPts);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for (int ie = 0; ie < epsilonPts; ++ie)
	{
		double* iTDistDistributed = tDistDistributed + ie*epsilonPts*zetaPts*aPts;
		for (int izeta = 0; izeta < zetaPts; ++izeta)
		{
			for (int ia = 0; ia < aPts; ++ia)
			{
#define Mat(data,ie,izeta,ia) data[(ia)*zetaPts*epsilonPts + (izeta)*epsilonPts + (ie)]
				double iDist = Mat(dist, ie, izeta, ia);
				int iManagerCell = Mat(apManagerCell, ie, izeta, ia);
				int iWorkerCell = Mat(apWorkerCell, ie, izeta, ia);
				double iManagerLeftShare = Mat(apManagerLeftShare, ie, izeta, ia);
				double iWorkerLeftShare = Mat(apWorkerLeftShare, ie, izeta, ia);
				double iOcc = Mat(occ, ie, izeta, ia);
#define EpsilonTransMat(ie,ieprime) epsilonTrans[(ie)+(ieprime)*epsilonPts]
#define ZetaTransMat(izeta,izetaprime) zetaTrans[(izeta)+(izetaprime)*zetaPts]
				for (int ieprime = 0; ieprime < epsilonPts; ++ieprime)
				{
					for (int izetaprime = 0; izetaprime < zetaPts; ++izetaprime)
					{
						Mat(iTDistDistributed, ieprime, izetaprime, iManagerCell) += iDist*iManagerLeftShare*iOcc*EpsilonTransMat(ie, ieprime)*ZetaTransMat(izeta, izetaprime);
						Mat(iTDistDistributed, ieprime, izetaprime, iManagerCell + 1) += iDist*(1 - iManagerLeftShare)*iOcc*EpsilonTransMat(ie, ieprime)*ZetaTransMat(izeta, izetaprime);
						Mat(iTDistDistributed, ieprime, izetaprime, iWorkerCell) += iDist*iWorkerLeftShare*(1-iOcc)*EpsilonTransMat(ie, ieprime)*ZetaTransMat(izeta, izetaprime);
						Mat(iTDistDistributed, ieprime, izetaprime, iWorkerCell + 1) += iDist*(1 - iWorkerLeftShare)*(1-iOcc)*EpsilonTransMat(ie, ieprime)*ZetaTransMat(izeta, izetaprime);
					}
				}
			}
		}
	}
	for (int ie = 0; ie < epsilonPts; ++ie)
	{
		double* iTDistDistributed = tDistDistributed + ie*epsilonPts*zetaPts*aPts;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
		for (int i = 0; i < epsilonPts*zetaPts*aPts; ++i)
			tDist[i] += iTDistDistributed[i];
	}
	free(tDistDistributed);
}

