#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "defs.h"
#include "mem.h"
#include "common.h"
#include "context.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static uint64_t XHASH(uint64_t x){
  return (x * 786433 + 196613) % 68719476735;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static uint64_t YHASH(uint64_t y){
  return (y * 786491 + 216617) % 66719476787;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static uint64_t ZHASH(uint64_t z){
  z = (~z) + (z << 21);
  z = z    ^ (z >> 24);
  z = (z   + (z << 3)) + (z << 8);
  z = z    ^ (z >> 14);
  z = (z   + (z << 2)) + (z << 4);
  z = z    ^ (z >> 28);
  z = z    + (z << 31);
  return z;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void InitHashTable(CModel *M){ 
  M->hTable.entries      = (Entry   **) Calloc(M->hTable.size, sizeof(Entry *));
  M->hTable.entrySize    = (uint8_t  *) Calloc(M->hTable.size, sizeof(uint8_t));
  M->hTable.zeroCounters = (HCC     **) Calloc(M->nSym, sizeof(HCC *));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void FreeCModel(CModel *M){
  U32 k, n;
  if(M->mode == HASH_TABLE_MODE){
    for(k = 0 ; k < M->hTable.size ; ++k){
      if(M->hTable.entrySize[k] > 0){
        for(n = 0 ; n < M->hTable.entrySize[k] ; ++n){
          Free(M->hTable.entries[k][n].counters);
          }
        Free(M->hTable.entries[k]);
        }
      }
    Free(M->hTable.entries);
    Free(M->hTable.entrySize);
    }
  else // TABLE_MODE
    Free(M->array.counters);
  Free(M);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void InitArray(CModel *M){
  M->array.counters = (ACC *) Calloc(M->nPModels*M->nSym, sizeof(ACC));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void InsertKey(HashTable *H, U32 hi, U64 key, U8 s){
  H->entries[hi] = (Entry *) Realloc(H->entries[hi], (H->entrySize[hi] + 1) * 
  sizeof(Entry), sizeof(Entry));
  H->entries[hi][H->entrySize[hi]].counters = (HCC *) Calloc(s, sizeof(HCC));
  H->entries[hi][H->entrySize[hi]].key = key;
  H->entrySize[hi]++;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static HCC *GetHCCounters(HashTable *H, U64 key){

 key = ZHASH(key);

 uint64_t n, hi = key % H->size;              //The hash index

 for(n = 0 ; n < H->entrySize[hi] ; ++n)
   if(H->entries[hi][n].key == key)     // If key found
     return H->entries[hi][n].counters;

  return NULL;
  } 

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PModel *CreatePModel(U32 n){
  PModel *P = (PModel *) Calloc(1, sizeof(PModel));
  P->freqs  = (U32    *) Calloc(n, sizeof(U32));
  return P;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FloatPModel *CreateFloatPModel(U32 n){
  FloatPModel *F = (FloatPModel *) Calloc(1, sizeof(FloatPModel));
  F->freqs = (double *) Calloc(n, sizeof(double));
  return F;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateCModelCounter(CModel *M, U32 sym, U64 im){
  U32 n;
  ACC *AC;
  U64 idx = ZHASH(im);

  if(M->mode == HASH_TABLE_MODE){
    uint64_t s;
    uint64_t hIndex = idx % M->hTable.size;

    for(n = 0 ; n < M->hTable.entrySize[hIndex] ; ++n)
      if(M->hTable.entries[hIndex][n].key == idx){
        if(++M->hTable.entries[hIndex][n].counters[sym] == 255) 
          for(s = 0 ; s < M->nSym ; ++s)
            M->hTable.entries[hIndex][n].counters[s] /= 2;
        return;
        }

    InsertKey(&M->hTable, hIndex, idx, M->nSym);
    M->hTable.entries[hIndex][n].counters[sym]++;
    }
  else{
    AC = &M->array.counters[im * M->nSym];
    if(++AC[sym] == M->maxCount){    
      for(n = 0 ; n < M->nSym ; ++n)
        AC[n] /= 2;
      }
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CModel *CreateCModel(U32 ctx, U32 aDen, U8 ref, U32 edits, U32 eDen, U32 nSym){
  CModel *M = (CModel *) Calloc(1, sizeof(CModel));
  U64    prod = 1, *mult;
  U32    n;

  if(ctx > MAX_HASH_CTX){
    fprintf(stderr, "Error: context size cannot be greater than %d\n", 
    MAX_HASH_CTX);
    exit(1);
    }

  M->nSym        = nSym;
  mult           = (U64 *) Calloc(ctx, sizeof(U64));
  M->nPModels    = (U64) pow(M->nSym, ctx);
  M->ctx         = ctx;
  M->alphaDen    = aDen;
  M->edits       = edits;
  M->pModelIdx   = 0;
  M->ref         = ref == 0 ? 0 : 1;

  if((ULL)(M->nPModels) * M->nSym * sizeof(ACC) >> 20 > MAX_ARRAY_MEMORY){
    M->mode     = HASH_TABLE_MODE;
    M->maxCount = DEFAULT_MAX_COUNT >> 8;
    M->hTable.size = HASH_SIZE;
    InitHashTable(M);
    }
  else{
    M->mode     = ARRAY_MODE;
    M->maxCount = DEFAULT_MAX_COUNT;
    InitArray(M);
    }

  for(n = 0 ; n < M->ctx ; ++n){
    mult[n] = prod;
    prod *= M->nSym;
    }

  M->multiplier = mult[M->ctx-1];

  if(edits != 0){
    M->SUBS.seq       = CreateCBuffer(BUFFER_SIZE, BGUARD);
    M->SUBS.in        = 0;
    M->SUBS.idx       = 0;
    M->SUBS.mask      = (uint8_t *) Calloc(BGUARD, sizeof(uint8_t));
    M->SUBS.threshold = edits;
    M->SUBS.eDen      = eDen;
    }

  Free(mult);
  return M;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t BestId(uint32_t *f, uint32_t sum, uint32_t nSym){
  if(sum == nSym) return -2; // ZERO COUNTERS

  uint32_t x, best = 0, max = f[0];
  for(x = 1 ; x < nSym ; ++x)
    if(f[x] > max){
      max = f[x];
      best = x;
      }

  for(x = 0 ; x < nSym ; ++x) if(best != x && max == f[x]) return -1;

  return best;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ResetCModelIdx(CModel *M){
  M->pModelIdx   = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void GetPModelIdx(U8 *p, CModel *M){
  M->pModelIdx = ((M->pModelIdx-*(p-M->ctx)*M->multiplier)*M->nSym)+*p;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint64_t GetPModelIdxCorr(U8 *p, CModel *M, uint64_t idx){
  return (((idx-*(p-M->ctx)*M->multiplier)*M->nSym)+*p);
  }
 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void FailSUBS(CModel *M){
  uint32_t x, fails = 0;
  for(x = 0 ; x < M->ctx ; ++x)
    if(M->SUBS.mask[x] != 0)
      ++fails;
  if(fails <= M->SUBS.threshold)
    ShiftBuffer(M->SUBS.mask, M->ctx, 1);
  else 
    M->SUBS.in = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void HitSUBS(CModel *M){
  ShiftBuffer(M->SUBS.mask, M->ctx, 0);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void CorrectCModelSUBS(CModel *M, PModel *P, uint8_t sym){
  int32_t best = BestId(P->freqs, P->sum, M->nSym);
  switch(best){
    case -2:  // IT IS A ZERO COUNTER [NOT SEEN BEFORE]
      if(M->SUBS.in != 0)
        HitSUBS(M); // HIT WORKS BETTER WITH AMINO ACIDS
    break;
    case -1:  // IT HAS AT LEAST TWO MAXIMUM FREQS [CONFUSION FREQS]
      if(M->SUBS.in != 0)
        FailSUBS(M);
    break;
    default:  // IT HAS ONE MAXIMUM FREQ
      if(M->SUBS.in == 0){ // IF IS OUT
        M->SUBS.in = 1;
        memset(M->SUBS.mask, 0, M->ctx);
        }
      else{ // IF IS IN
        if(best == sym) HitSUBS(M);
        else{
          FailSUBS(M);
          M->SUBS.seq->buf[M->SUBS.seq->idx] = best; 
          } // UPDATE BUFFER WITH NEW SYMBOL
        }
    }
  UpdateCBuffer(M->SUBS.seq);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ComputePModel(CModel *M, PModel *P, uint64_t idx, uint32_t aDen){
  ACC *ac;
  HCC *hc;
  uint32_t x;
  switch(M->mode){

    case HASH_TABLE_MODE:
     if(!(hc = GetHCCounters(&M->hTable, idx)))
       hc = (HCC *) M->hTable.zeroCounters;
     P->sum = 0;
     for(x = 0 ; x < M->nSym ; ++x){
       P->freqs[x] = 1 + aDen * hc[x];
       P->sum += P->freqs[x];
       }
    break;

    case ARRAY_MODE:
      ac = &M->array.counters[idx*M->nSym];
      P->sum = 0;
      for(x = 0 ; x < M->nSym ; ++x){
        P->freqs[x] = 1 + aDen * ac[x];
        P->sum += P->freqs[x];
        }
    break;

    default:
    fprintf(stderr, "Error: not implemented!\n");
    exit(1);
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ComputeWeightedFreqs(double w, PModel *P, FloatPModel *PT, uint32_t nSym){
  uint32_t x;
  double f = w / P->sum;
  for(x = 0 ; x < nSym ; ++x)
    PT->freqs[x] += (double) P->freqs[x] * f;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double PModelSymbolNats(PModel *P, U32 s){
  return log((double) P->sum / P->freqs[s]);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
