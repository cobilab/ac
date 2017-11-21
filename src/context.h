#ifndef CONTEXT_H_INCLUDED
#define CONTEXT_H_INCLUDED

#include "defs.h"
#include "buffer.h"

#define WINDOW_SIZE           6    // IT WILL ACCEPT X SUBSTITUTIONS IN 6
#define ARRAY_MODE            0
#define HASH_TABLE_MODE       1
#define HASH_TABLE_BEGIN_CTX  15
//#define HASH_SIZE             33554471
#define HASH_SIZE             16777259

#if defined(PREC32B)
  #define MAX_HASH_CTX        28 
#elif defined(PREC16B)
  #define MAX_HASH_CTX        20 
#else
  #define MAX_HASH_CTX        16
#endif

#define  TMM_ON               1
#define  TMM_OFF              0

typedef U16  ACC;                  // Size of context counters for arrays
typedef U8   HCC;                // Size of context counters for hash tables
typedef U8   ENTMAX;                // Entry size (nKeys for each hIndex)

typedef struct{
  uint64_t  key;            //The key stored in this entry 
  HCC       *counters;      //The context counters 
  }
Entry;

typedef struct{
  uint32_t  size;           //Size of the hash table
  uint8_t   *entrySize;     //Number of keys in this entry
  Entry     **entries;      //The heads of the hash table lists
  HCC       **zeroCounters;  
  }
HashTable;

typedef struct{
  ACC       *counters;
  }
Array;

typedef struct{
  uint32_t  in;
  CBUF      *seq;           // BUFFER FOR EDITED SEQUENCE
  uint8_t   *mask;          // BUFFER FOR FAILS & HITS
  uint64_t  idx;            // INDEX FOR UPDATE
  uint64_t  idx2;           // AUXILIAR INDEX FOR UPDATE
  uint32_t  threshold;      // DISCARD ABOVE THIS VALUE
  uint32_t  eDen;           // ALPHA DENOMINATOR FOR THIS MODEL
  }
Correct;

typedef struct{
  U32        ctx;           // Current depth of context template
  U64        nPModels;      // Maximum number of probability models
  U32        alphaDen;      // Denominator of alpha
  U32        maxCount;      // Counters /= 2 if one counter >= maxCount
  U64        multiplier;
  U8         ref;
  U32        mode;
  HashTable  hTable;
  Array      array;
  U64        pModelIdx;     // IDX
  U32        edits;         // EDITS
  U32        nSym;          // EDITS
  Correct    SUBS;          // EDITS
  }
CModel;

typedef struct{
  U32        *freqs;
  U32        sum;
  }
PModel;

typedef struct{
  double     *freqs;
  }
FloatPModel;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t         BestId               (uint32_t *, uint32_t, uint32_t);
void            HitSUBS              (CModel *);
void            FailSUBS             (CModel *);
void            FreeCModel           (CModel *);
void            GetPModelIdx         (U8 *, CModel *);
uint64_t        GetPModelIdxCorr     (U8 *, CModel *, uint64_t);
void            CorrectCModelSUBS    (CModel *, PModel *, uint8_t);
PModel          *CreatePModel        (U32);
FloatPModel     *CreateFloatPModel   (U32);
void            ResetCModelIdx       (CModel *);
void            UpdateCModelCounter  (CModel *, U32, U64);
CModel          *CreateCModel        (U32, U32, U8, U32, U32, U32);
void            ComputePModel        (CModel *, PModel *, uint64_t, uint32_t);
void            ComputeWeightedFreqs (double, PModel *, FloatPModel *, uint32_t);
double          PModelSymbolNats     (PModel *, U32);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
