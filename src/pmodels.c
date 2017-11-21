#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "defs.h"
#include "mem.h"
#include "common.h"
#include "pmodels.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CMWeight *CreateWeightModel(uint32_t size){
  uint32_t n;
  CMWeight *CMW    = (CMWeight *) Calloc(1, sizeof(CMWeight));
  CMW->totModels   = size;
  CMW->totalWeight = 0;
  CMW->weight      = (double *) Calloc(CMW->totModels, sizeof(double));
  for(n = 0 ; n < CMW->totModels ; ++n)
    CMW->weight[n] = 1.0 / CMW->totModels;
  return CMW;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ResetWeightModel(CMWeight *CMW){
  uint32_t n;
  double fraction = 1.0 / CMW->totModels;
  CMW->totalWeight = 0;
  for(n = 0 ; n < CMW->totModels ; ++n)
    CMW->weight[n] = fraction;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void RenormalizeWeights(CMWeight *CMW){
  uint32_t n;
  for(n = 0 ; n < CMW->totModels ; ++n) 
    CMW->weight[n] /= CMW->totalWeight;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void CalcDecayment(CMWeight *CMW, PModel **PM, uint8_t sym, double gamma){
  uint32_t n;
  CMW->totalWeight = 0;
  for(n = 0 ; n < CMW->totModels ; ++n){
    CMW->weight[n] = Power(CMW->weight[n], gamma) * (double) PM[n]->freqs[sym] 
    / PM[n]->sum;
    CMW->totalWeight += CMW->weight[n];
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void RemoveWeightModel(CMWeight *CMW){
  Free(CMW->weight);
  Free(CMW);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ComputeMXProbs(FloatPModel *PT, PModel *MX, uint32_t nSym){
  uint32_t x;
  MX->sum = 0;
  for(x = 0 ; x < nSym ; ++x)
    MX->sum += MX->freqs[x] = 1 + (unsigned) (PT->freqs[x] * MX_PMODEL);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PModel *CreatePModel(U32 n){
  PModel *PM = (PModel *) Calloc(1, sizeof(PModel));
  PM->freqs  = (U32    *) Calloc(n, sizeof(U32));
  return PM;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void RemovePModel(PModel *PM){
  Free(PM->freqs);
  Free(PM);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FloatPModel *CreateFloatPModel(U32 n){
  FloatPModel *F = (FloatPModel *) Calloc(1, sizeof(FloatPModel));
  F->freqs = (double *) Calloc(n, sizeof(double));
  return F;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void RemoveFPModel(FloatPModel *FM){
  Free(FM->freqs);
  Free(FM);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ComputeWeightedFreqs(double w, PModel *P, FloatPModel *PT, uint32_t nSym){
  uint32_t x;
  double f = w / P->sum;
  for(x = 0 ; x < nSym ; ++x)
    PT->freqs[x] += (double) P->freqs[x] * f;
  }


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

