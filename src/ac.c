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
  double      *cModelWeight, cModelTotalWeight = 0;
  uint8_t     *readerBuffer, sym, irSym, *pos, type = 0, 
              header = 1, line = 0, dna = 0;
  PModel      **pModel, *MX;
  FloatPModel *PT;
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);

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
      #ifdef RUN_SUBS
      totModels++;
      #endif
      }

  pModel        = (PModel  **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n]   = CreatePModel(AL->cardinality);
  MX            = CreatePModel(AL->cardinality);
  PT            = CreateFloatPModel(AL->cardinality);
  readerBuffer  = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  cModelWeight  = (double   *) Calloc(totModels, sizeof(double));
  
  for(n = 0 ; n < totModels ; ++n)
    cModelWeight[n] = 1.0 / totModels;

  for(n = 0 ; n < P->nModels ; ++n){
    if(P->model[n].type == TARGET){
      cModels[n] = CreateCModel(P->model[n].ctx, P->model[n].den, TARGET, 
      P->model[n].edits, P->model[n].eDen, AL->cardinality);
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
  WriteNBits((int) (P->gamma * 65536), 32, Writter);
  WriteNBits(P->nModels,               16, Writter);
  for(n = 0 ; n < P->nModels ; ++n){
    WriteNBits(cModels[n]->ctx,        16, Writter);
    WriteNBits(cModels[n]->alphaDen,   16, Writter);
    WriteNBits(cModels[n]->edits,       8, Writter);
    WriteNBits(cModels[n]->SUBS.eDen,  32, Writter);
    WriteNBits(P->model[n].type,        1, Writter);
    }

  I[id].header = _bytes_output;

  i = 0;
  while((k = fread(readerBuffer, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){

      CalcProgress(size, ++i);
      
      /*
      if(IsLowChar(AL, readerBuffer[idxPos]) == 1){
        #ifdef ESTIMATE
        if(P->estim != 0)
          fprintf(IAE, "%.3g\n", log2(AL->cardinality));
        #endif
        continue;
        }
      */

      symBuf->buf[symBuf->idx] = sym = AL->revMap[ readerBuffer[idxPos] ];
      memset((void *)PT->freqs, 0, AL->cardinality * sizeof(double));

      n = 0;
      pos = &symBuf->buf[symBuf->idx-1];
      for(cModel = 0 ; cModel < P->nModels ; ++cModel){
        CModel *CM = cModels[cModel];
        GetPModelIdx(pos, CM);
        ComputePModel(CM, pModel[n], CM->pModelIdx, CM->alphaDen);
        ComputeWeightedFreqs(cModelWeight[n], pModel[n], PT, CM->nSym);
        if(CM->edits != 0){
          ++n;
          CM->SUBS.seq->buf[CM->SUBS.seq->idx] = sym;
          CM->SUBS.idx = GetPModelIdxCorr(CM->SUBS.seq->buf+
          CM->SUBS.seq->idx-1, CM, CM->SUBS.idx);
          ComputePModel(CM, pModel[n], CM->SUBS.idx, CM->SUBS.eDen);
          ComputeWeightedFreqs(cModelWeight[n], pModel[n], PT, CM->nSym);
          }
        ++n;
        }

      MX->sum = 0;
      for(x = 0 ; x < AL->cardinality ; ++x){
        MX->sum += MX->freqs[x] = 1 + (unsigned) (PT->freqs[x] * MX_PMODEL);
        }

      AESym(sym, (int *)(MX->freqs), (int) MX->sum, Writter);
      #ifdef ESTIMATE
      if(P->estim != 0)
        fprintf(IAE, "%.3g\n", PModelSymbolNats(MX, sym) / M_LN2);
      #endif

      cModelTotalWeight = 0;
      for(n = 0 ; n < totModels ; ++n){
        cModelWeight[n] = Power(cModelWeight[n], P->gamma) * (double) 
        pModel[n]->freqs[sym] / pModel[n]->sum;
        cModelTotalWeight += cModelWeight[n];
        }

      for(n = 0 ; n < P->nModels ; ++n)
        if(cModels[n]->ref == TARGET)
          UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);

      for(n = 0 ; n < totModels ; ++n)
        cModelWeight[n] /= cModelTotalWeight; // RENORMALIZE THE WEIGHTS

      n = 0;
      for(cModel = 0 ; cModel < P->nModels ; ++cModel){
        if(cModels[cModel]->edits != 0)
          CorrectCModelSUBS(cModels[cModel], pModel[++n], sym);
        ++n;
        }

      UpdateCBuffer(symBuf);
      }

  finish_encode(Writter);
  doneoutputtingbits(Writter);
  fclose(Writter);


  #ifdef ESTIMATE
  if(P->estim == 1){
    fclose(IAE);
    Free(IAEName);
    }
  #endif

  Free(MX);
  Free(name);
  Free(cModelWeight);
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == REFERENCE)
      ResetCModelIdx(cModels[n]);
    else
      FreeCModel(cModels[n]);
  for(n = 0 ; n < totModels ; ++n){
    Free(pModel[n]->freqs);
    Free(pModel[n]);
    }
  Free(pModel);
  Free(PT);
  Free(readerBuffer);
  RemoveCBuffer(symBuf);
  RemoveAlphabet(AL);
  int card = AL->cardinality;
  fclose(Reader);

  if(P->verbose == 1)
    fprintf(stderr, "Done!                          \n");  // SPACES ARE VALID 

  I[id].bytes = _bytes_output;
  I[id].size  = i;
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
      P->model[n].edits, P->model[n].eDen, AL->cardinality);

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
  double      gamma;
  
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

  gamma = DEFAULT_GAMMA;
  for(n = 1 ; n < xargc ; ++n) 
    if(strcmp(xargv[n], "-g") == 0) 
      gamma = atof(xargv[n+1]);

  P->gamma    = ArgsDouble (gamma, p, argc, "-g");
  P->gamma    = ((int)(P->gamma * 65536)) / 65536.0;
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

  totalSize   = 0;
  totalBytes  = 0;
  headerBytes = 0;
  cardinality = 1;
  for(n = 0 ; n < P->nTar ; ++n){
    cardinality  = Compress(P, refModels, n, refNModels, I);
    totalSize   += I[n].size;
    totalBytes  += I[n].bytes;
    headerBytes += I[n].header;
    }

  if(P->nTar > 1)
    for(n = 0 ; n < P->nTar ; ++n){
      fprintf(stdout, "File %d compressed bytes: %"PRIu64" (", n+1, (uint64_t) 
      I[n].bytes);
      PrintHRBytes(I[n].bytes);

      fprintf(stdout, ") , Normalized Dissimilarity Rate: %.6g\n", 
      (8.0*I[n].bytes)/(log2(cardinality)*I[n].size));
      }


  fprintf(stdout, "Total bytes: %"PRIu64" (", totalBytes);
  PrintHRBytes(totalBytes);
  fprintf(stdout, "), %.4g bpb, %.4g bps w/ no header, Normalized Dissimilarity" 
  " Rate: %.6g\n", ((8.0*totalBytes)/totalSize), ((8.0*(totalBytes-headerBytes))
  /totalSize), (8.0*totalBytes)/(log2(cardinality)*totalSize));  
  stop = clock();
  fprintf(stdout, "Spent %g sec.\n", ((double)(stop-start))/CLOCKS_PER_SEC);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
