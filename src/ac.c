#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include "mem.h"
#include "defs.h"
#include "msg.h"
#include "buffer.h"
#include "alphabet.h"
#include "levels.h"
#include "common.h"
#include "pmodels.h"
#include "tolerant.h"
#include "context.h"
#include "bitio.h"
#include "arith.h"
#include "arith_aux.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S O R - - - - - - - - - - - - - -

int Compress(Parameters *P, CModel **cModels, uint8_t id, uint32_t 
refNModels, INF *I){
  FILE        *Reader  = Fopen(P->tar[id], "r");
  char        *name    = concatenate(P->tar[id], ".co");
  FILE        *Writter = Fopen(name, "w");
  uint32_t    n, k, x, cModel, totModels, idxPos;
  int32_t     idx = 0;
  uint64_t    i, size = 0;
  uint8_t     *readerBuffer, sym, irSym, *pos, type = 0, 
              header = 1, line = 0, dna = 0;
  double      se = 0;
  PModel      **pModel, *MX;
  FloatPModel *PT;
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  CMWeight    *WM;

  if(P->verbose)
    fprintf(stderr, "Analyzing data and creating models ...\n");

  #ifdef ESTIMATE
  FILE *IAE = NULL;
  char *IAEName = NULL;
  if(P->estim == 1){
    IAEName = concatenate(P->tar[id], ".iae");
    IAE = Fopen(IAEName, "w");
    }
  #endif
  
  _bytes_output = 0;
  size = NBytesInFile(Reader);

  // BUILD ALPHABET
  ALPHABET *AL = CreateAlphabet(P->low);
  LoadAlphabet(AL, Reader);
  PrintAlphabet(AL);

  // ADAPT ALPHABET FOR NON FREQUENT SYMBOLS
  AdaptAlphabetNonFrequent(AL, Reader);
  
  // EXTRA MODELS DERIVED FROM EDITS
  totModels = P->nModels;
  for(n = 0 ; n < P->nModels ; ++n) 
    if(P->model[n].edits != 0){
      totModels++;
      }

  fprintf(stderr, "Using %u probabilistic models\n", totModels);

  pModel        = (PModel  **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n]   = CreatePModel(AL->cardinality);
  MX            = CreatePModel(AL->cardinality);
  PT            = CreateFloatPModel(AL->cardinality);
  WM            = CreateWeightModel(totModels);
  readerBuffer  = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));

  for(n = 0 ; n < P->nModels ; ++n){
    if(P->model[n].type == TARGET){
      cModels[n] = CreateCModel(P->model[n].ctx, P->model[n].den, TARGET, 
      P->model[n].edits, P->model[n].eDen, AL->cardinality, P->model[n].gamma,
      P->model[n].eGamma);
      }
    }

  if(P->verbose){
    fprintf(stderr, "Done!\n");
    fprintf(stderr, "Compressing target sequence %d [%"PRIu64" symbols] ...\n", 
    id + 1, size);
    }

  startoutputtingbits();
  start_encode();

  WriteNBits(WATERMARK,                32, Writter);
  WriteNBits(P->checksum,              46, Writter);
  WriteNBits(size,                     46, Writter);

  // PRE HEADER : NON FREQUENT SYMBOLS
  WriteNBits(P->low,                   32, Writter);
  WriteNBits(AL->nLow,                 32, Writter);
  for(x = 0 ; x < AL->nLow ; ++x){
    WriteNBits(AL->lowAlpha[x],         8, Writter);
    }

  // REMAP ALPHABET
  //// ResetAlphabet(AL);
  // PrintAlphabet(AL);

  WriteNBits(AL->cardinality,          16, Writter);
  for(x = 0 ; x < AL->cardinality ; ++x)
    WriteNBits(AL->toChars[x],          8, Writter);
  WriteNBits(P->nModels,               16, Writter);
  for(n = 0 ; n < P->nModels ; ++n){
    WriteNBits(cModels[n]->ctx,         5, Writter);
    WriteNBits(cModels[n]->alphaDen,   11, Writter);
    WriteNBits((int)(cModels[n]->gamma * 65536), 17, Writter);
    WriteNBits(cModels[n]->edits,       7, Writter);
    if(cModels[n]->edits != 0){
      WriteNBits((int)(cModels[n]->eGamma * 65536), 17, Writter);
      WriteNBits(cModels[n]->TM->den,     9, Writter);
      }
    WriteNBits(P->model[n].type,        1, Writter);
    }

  I[id].header = _bytes_output;

  // GIVE SPECIFIC GAMMA:
  int pIdx = 0;
  for(n = 0 ; n < P->nModels ; ++n){
    WM->gamma[pIdx++] = cModels[n]->gamma;
    if(P->model[n].edits != 0){
      WM->gamma[pIdx++] = cModels[n]->eGamma;
      }
    }

  i = 0;
  while((k = fread(readerBuffer, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){

      CalcProgress(size, ++i);
      
      symBuf->buf[symBuf->idx] = sym = AL->revMap[ readerBuffer[idxPos] ];
      memset((void *)PT->freqs, 0, AL->cardinality * sizeof(double));

      n = 0;
      pos = &symBuf->buf[symBuf->idx-1];
      for(cModel = 0 ; cModel < P->nModels ; ++cModel){
        CModel *CM = cModels[cModel];
        GetPModelIdx(pos, CM);
        ComputePModel(CM, pModel[n], CM->pModelIdx, CM->alphaDen);
        ComputeWeightedFreqs(WM->weight[n], pModel[n], PT, CM->nSym);
        if(CM->edits != 0){
          ++n;
          CM->TM->seq->buf[CM->TM->seq->idx] = sym;
          CM->TM->idx = GetPModelIdxCorr(CM->TM->seq->buf+
          CM->TM->seq->idx-1, CM, CM->TM->idx);
          ComputePModel(CM, pModel[n], CM->TM->idx, CM->TM->den);
          ComputeWeightedFreqs(WM->weight[n], pModel[n], PT, CM->nSym);
          }
        ++n;
        }

      ComputeMXProbs(PT, MX, AL->cardinality);

      AESym(sym, (int *)(MX->freqs), (int) MX->sum, Writter);
      #ifdef ESTIMATE
      if(P->estim != 0)
        fprintf(IAE, "%.3g\n", PModelSymbolNats(MX, sym) / M_LN2);
      #endif

      CalcDecayment(WM, pModel, sym);

      for(n = 0 ; n < P->nModels ; ++n)
        if(cModels[n]->ref == TARGET)
          UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);

      RenormalizeWeights(WM);

      n = 0;
      for(cModel = 0 ; cModel < P->nModels ; ++cModel){
        if(cModels[cModel]->edits != 0)
          UpdateTolerantModel(cModels[cModel]->TM, pModel[++n], sym);
        ++n;
        }

      UpdateCBuffer(symBuf);
      }

  finish_encode(Writter);
  doneoutputtingbits(Writter);
  fclose(Writter);

  se = PrintSE(AL);

  #ifdef ESTIMATE
  if(P->estim == 1){
    fclose(IAE);
    Free(IAEName);
    }
  #endif

  Free(MX);
  Free(name);
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == REFERENCE)
      ResetCModelIdx(cModels[n]);
    else
      RemoveCModel(cModels[n]);
  for(n = 0 ; n < totModels ; ++n){
    RemovePModel(pModel[n]);
    }
  Free(pModel);

  Free(PT);
  Free(readerBuffer);
  RemoveCBuffer(symBuf);
  RemoveAlphabet(AL);
  RemoveWeightModel(WM);
  int card = AL->cardinality;
  fclose(Reader);

  if(P->verbose == 1)
    fprintf(stderr, "Done!                          \n");  // SPACES ARE VALID 

  I[id].bytes = _bytes_output;
  I[id].size  = i;
  I[id].se    = se;
  return card;
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - R E F E R E N C E - - - - - - - - - - - - -

CModel **LoadReference(Parameters *P){
  FILE      *Reader = Fopen(P->ref, "r");
  uint32_t  n, k, idxPos;
  uint64_t  nSymbols = 0;
  int32_t   idx = 0;
  uint8_t   *readerBuffer, sym;
  CBUF      *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  CModel    **cModels;
  
  #ifdef PROGRESS
  uint64_t  i = 0;
  #endif

  if(P->verbose == 1)
    fprintf(stdout, "Building reference model ...\n");

  // BUILD ALPHABET
  ALPHABET *AL = CreateAlphabet(P->low);
  LoadAlphabet(AL, Reader);
  PrintAlphabet(AL);

  readerBuffer  = (uint8_t *) Calloc(BUFFER_SIZE + 1, sizeof(uint8_t));
  cModels       = (CModel **) Malloc(P->nModels * sizeof(CModel *)); 
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == REFERENCE)
      cModels[n] = CreateCModel(P->model[n].ctx, P->model[n].den, REFERENCE, 
      P->model[n].edits, P->model[n].eDen, AL->cardinality, P->model[n].gamma,
      P->model[n].eGamma);

  nSymbols = NBytesInFile(Reader);

  P->checksum = 0;
  while((k = fread(readerBuffer, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){

      symBuf->buf[symBuf->idx] = sym = AL->revMap[ readerBuffer[idxPos] ];
      P->checksum = (P->checksum + (uint8_t) sym);

      for(n = 0 ; n < P->nModels ; ++n)
        if(P->model[n].type == REFERENCE){
          GetPModelIdx(symBuf->buf+symBuf->idx-1, cModels[n]);
          UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);
          }

      UpdateCBuffer(symBuf);

      #ifdef PROGRESS
      CalcProgress(nSymbols, ++i);
      #endif
      }
 
  P->checksum %= CHECKSUMGF; 
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == REFERENCE)
      ResetCModelIdx(cModels[n]);
  Free(readerBuffer);
  RemoveCBuffer(symBuf);
  RemoveAlphabet(AL);
  fclose(Reader);

  if(P->verbose == 1)
    fprintf(stdout, "Done!                          \n");  // SPACES ARE VALID  
  else
    fprintf(stdout, "                               \n");  // SPACES ARE VALID

  return cModels;
  }
  
//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t main(int argc, char *argv[]){
  char        **p = *&argv, **xargv, *xpl = NULL;
  CModel      **refModels;
  int32_t     xargc = 0, cardinality = 1;
  uint32_t    n, k, refNModels;
  uint64_t    totalBytes, headerBytes, totalSize;
  clock_t     stop = 0, start = clock();
  double      se_average;
  
  Parameters  *P;
  INF         *I;

  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  if((P->help = ArgsState(DEFAULT_HELP, p, argc, "-h")) == 1 || argc < 2){
    PrintMenu();
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_VERSION, p, argc, "-V")){
    PrintVersion();
    return EXIT_SUCCESS;
    }

  if(ArgsState(0, p, argc, "-s")){
    PrintLevels(); 
    return EXIT_SUCCESS;
    }

  P->verbose = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v" );
  P->force   = ArgsState  (DEFAULT_FORCE,   p, argc, "-f" );
  P->estim   = ArgsState  (0,               p, argc, "-e" );
  P->level   = ArgsNum    (0,   p, argc, "-l", MIN_LEVEL, MAX_LEVEL);
  P->low     = ArgsNum    (200, p, argc, "-t", MIN_THRESHOLD, MAX_THRESHOLD);

  P->nModels = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-rm") == 0 || strcmp(argv[n], "-tm") == 0)
      P->nModels += 1;

  if(P->nModels == 0 && P->level == 0)
    P->level = DEFAULT_LEVEL;
  
  if(P->level != 0){
    xpl = GetLevels(P->level);
    xargc = StrToArgv(xpl, &xargv);
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-rm") == 0 || strcmp(xargv[n], "-tm") == 0)
        P->nModels += 1;
    }

  if(P->nModels == 0){
    fprintf(stderr, "Error: at least you need to use a context model!\n");
    return 1;
    }

  P->model = (ModelPar *) Calloc(P->nModels, sizeof(ModelPar));

  k = 0;
  refNModels = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-rm") == 0){
      P->model[k++] = ArgsUniqModel(argv[n+1], 1);
      ++refNModels;
      }
  if(P->level != 0){
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-rm") == 0){
        P->model[k++] = ArgsUniqModel(xargv[n+1], 1);
        ++refNModels;
        }
    }

  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-tm") == 0)
      P->model[k++] = ArgsUniqModel(argv[n+1], 0);
  if(P->level != 0){
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-tm") == 0)
        P->model[k++] = ArgsUniqModel(xargv[n+1], 0);
    }

  P->ref      = ArgsString (NULL, p, argc, "-r");
  P->nTar     = ReadFNames (P, argv[argc-1]);
  P->checksum = 0;
  if(P->verbose) 
    PrintArgs(P);

  if(refNModels == 0)
    refModels = (CModel **) Malloc(P->nModels * sizeof(CModel *));
  else{
    if(P->ref == NULL){
      fprintf(stderr, "Error: using reference model(s) in nonexistent "
      "reference sequence!\n");
      exit(1);
      }
    refModels = LoadReference(P);
    if(P->verbose)
      fprintf(stderr, "Checksum: %"PRIu64"\n", P->checksum);
    }

  I = (INF *) Calloc(P->nTar, sizeof(INF));

  se_average  = 0;
  totalSize   = 0;
  totalBytes  = 0;
  headerBytes = 0;
  cardinality = 1;
  for(n = 0 ; n < P->nTar ; ++n){
    cardinality  = Compress(P, refModels, n, refNModels, I);
    totalSize   += I[n].size;
    totalBytes  += I[n].bytes;
    headerBytes += I[n].header;
    se_average  += I[n].se;
    }
  se_average /= P->nTar;

  if(P->nTar > 1)
    for(n = 0 ; n < P->nTar ; ++n){
      fprintf(stdout, "File %d compressed bytes: %"PRIu64" (", n+1, (uint64_t) 
      I[n].bytes);
      PrintHRBytes(I[n].bytes);
      fprintf(stdout, ") , Normalized Dissimilarity Rate: %.6g\n", 
      (8.0*I[n].bytes)/(log2(cardinality)*I[n].size));

      fprintf(stdout, "Shannon entropy: %.6g\n", I[n].se);
      }

  fprintf(stdout, "Total bytes: %"PRIu64" (", totalBytes);
  PrintHRBytes(totalBytes);
  fprintf(stdout, "), %.5g bps, %.5g bps w/ no header\n",
  ((8.0*totalBytes)/totalSize), ((8.0*(totalBytes-headerBytes))/totalSize)); 

  fprintf(stdout, "Normalized Dissimilarity Rate: %.6g\n", (8.0*totalBytes)/
  (log2(cardinality)*totalSize));  

  fprintf(stdout, "Average Shannon entropy: %.6g\n", se_average);

  stop = clock();
  fprintf(stdout, "Spent %g sec.\n", ((double)(stop-start))/CLOCKS_PER_SEC);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
