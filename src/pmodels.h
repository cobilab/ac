#ifndef PMODELS_H_INCLUDED
#define PMODELS_H_INCLUDED

#include "defs.h"

#define MX_PMODEL      65535
#define DEFAULT_GAMMA  0.90

typedef struct{
  U32       *freqs;
  U32       sum;
  }
PModel;

typedef struct{
  double    *freqs;
  }
FloatPModel;

typedef struct{
  uint32_t  totModels;
  double    *weight;
  double    totalWeight;
  }
CMWeight;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PModel      *CreatePModel        (U32);
FloatPModel *CreateFloatPModel   (U32);
void        RemovePModel         (PModel *);
void        RemoveFPModel        (FloatPModel *);
void        ComputeMXProbs       (FloatPModel *, PModel *, uint32_t);
void        ComputeWeightedFreqs (double, PModel *, FloatPModel *, uint32_t);
CMWeight    *CreateWeightModel   (uint32_t);
void        ResetWeightModel     (CMWeight *);
void        RenormalizeWeights   (CMWeight *);
void        CalcDecayment        (CMWeight *, PModel **, uint8_t, double);
void        RemoveWeightModel    (CMWeight *);
double      PModelSymbolNats     (PModel *, uint32_t);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
