#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include "mem.h"
#include "msg.h"
#include "defs.h"
#include "buffer.h"
#include "alphabet.h"
#include "common.h"
#include "context.h"
#include "bitio.h"
#include "arith.h"
#include "arith_aux.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - D E C O M P R E S S O R - - - - - - - - - - - -

void Decompress(Parameters *P, CModel **cModels, uint8_t id){
  FILE        *Reader  = Fopen(P->tar[id], "r");
  char        *name    = ReplaceSubStr(P->tar[id], ".co", ".de"); 
  FILE        *Writter = Fopen(name, "w");
  uint64_t    nSymbols = 0;
  uint32_t    n, k, x, cModel, totModels;
  int32_t     idx = 0, idxOut = 0;
  uint8_t     *outBuffer, sym = 0, *pos;
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  PModel      **pModel, *MX;
  FloatPModel *PT;
  CMWeight    *WM;
  uint64_t    i = 0;

  if(P->verbose)
    fprintf(stderr, "Decompressing %"PRIu64" symbols of target %d ...\n", 
    P[id].size, id + 1);

  startinputtingbits();
  start_decode(Reader);

  P[id].watermark        = ReadNBits(32, Reader);
  garbage                = ReadNBits(46, Reader);
  P[id].size             = ReadNBits(46, Reader);
  P[id].low              = ReadNBits(32, Reader);
  ALPHABET *AL = CreateAlphabet(P[id].low);
  AL->length   = P[id].size;
  AL->nLow               = ReadNBits(32, Reader);
  for(x = 0 ; x < AL->nLow ; ++x)
    AL->lowAlpha[x]      = ReadNBits( 8, Reader);
  AL->cardinality        = ReadNBits(16, Reader);
  for(x = 0 ; x < 256 ; ++x)
    AL->revMap[x] = INVALID_SYMBOL;
  for(x = 0 ; x < AL->cardinality ; ++x){
    AL->toChars[x]       = ReadNBits( 8, Reader);
    AL->revMap[(uint8_t) AL->toChars[x]] = x;
    }
  P[id].nModels          = ReadNBits(16, Reader);
  for(k = 0 ; k < P[id].nModels ; ++k){
    P[id].model[k].ctx   = ReadNBits( 5, Reader);
    P[id].model[k].den   = ReadNBits(11, Reader);
    P[id].model[k].gamma = ReadNBits(17, Reader) / 65536.0;
    P[id].model[k].edits = ReadNBits( 7, Reader);
    P[id].model[k].eDen  = ReadNBits( 9, Reader);
    P[id].model[k].type  = ReadNBits( 1, Reader);
    }

  PrintAlphabet(AL);

  // EXTRA MODELS DERIVED FROM EDITS
  totModels = P[id].nModels;
  for(n = 0 ; n < P[id].nModels ; ++n)
    if(P[id].model[n].edits != 0)
      ++totModels; // TOLERANT

  nSymbols      = P[id].size;
  pModel        = (PModel  **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n]   = CreatePModel(AL->cardinality);
  MX            = CreatePModel(AL->cardinality);
  PT            = CreateFloatPModel(AL->cardinality);
  outBuffer     = (uint8_t  *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  WM            = CreateWeightModel(totModels);

  for(n = 0 ; n < P[id].nModels ; ++n){
    if(P[id].model[n].type == TARGET)
      cModels[n] = CreateCModel(P[id].model[n].ctx , P[id].model[n].den, 
      TARGET, P[id].model[n].edits, P[id].model[n].eDen, AL->cardinality,
      P[id].model[n].gamma);
    }

  // GIVE SPECIFIC GAMMA:
  int pIdx = 0;
  for(n = 0 ; n < P[id].nModels ; ++n){
    WM->gamma[pIdx++] = cModels[n]->gamma;
    if(P[id].model[n].edits != 0){
      WM->gamma[pIdx++] = cModels[n]->gamma;
      }
    }

  i = 0;
  while(nSymbols--){
    CalcProgress(P[id].size, ++i);

    memset((void *)PT->freqs, 0, AL->cardinality * sizeof(double));

    n = 0;
    pos = &symBuf->buf[symBuf->idx-1];
    for(cModel = 0 ; cModel < P[id].nModels ; ++cModel){
      GetPModelIdx(pos, cModels[cModel]);
      ComputePModel(cModels[cModel], pModel[n], cModels[cModel]->pModelIdx,
      cModels[cModel]->alphaDen);
      ComputeWeightedFreqs(WM->weight[n], pModel[n], PT, AL->cardinality);
      if(cModels[cModel]->edits != 0){ // SUBSTITUTIONAL HANDLING
        ++n;
        cModels[cModel]->SUBS.idx = GetPModelIdxCorr(cModels[cModel]->SUBS.seq
        ->buf+cModels[cModel]->SUBS.seq->idx-1, cModels[cModel], cModels[cModel]
        ->SUBS.idx);
        ComputePModel(cModels[cModel], pModel[n], cModels[cModel]->SUBS.idx,
        cModels[cModel]->SUBS.eDen);
        ComputeWeightedFreqs(WM->weight[n], pModel[n], PT, AL->cardinality);
        }
      ++n;
      }

    ComputeMXProbs(PT, MX, AL->cardinality);

    symBuf->buf[symBuf->idx] = sym = ArithDecodeSymbol(AL->cardinality, (int *) 
    MX->freqs, (int) MX->sum, Reader);
    outBuffer[idxOut] = AL->toChars[sym];

    for(n = 0 ; n < P[id].nModels ; ++n)
      if(cModels[n]->edits != 0){
        cModels[n]->SUBS.seq->buf[cModels[n]->SUBS.seq->idx] = sym;
        }         

    CalcDecayment(WM, pModel, sym);

    for(n = 0 ; n < P[id].nModels ; ++n){
      if(P[id].model[n].type == TARGET)
        UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);
      }

    RenormalizeWeights(WM);

    n = 0;
    for(cModel = 0 ; cModel < P[id].nModels ; ++cModel){
      if(cModels[cModel]->edits != 0){      // CORRECT CMODEL CONTEXTS
        CorrectCModelSUBS(cModels[cModel], pModel[++n], sym);
        }
      ++n;
      }

    if(++idxOut == BUFFER_SIZE){
      fwrite(outBuffer, 1, idxOut, Writter);
      idxOut = 0;
      }

    UpdateCBuffer(symBuf);
    }

  if(idxOut != 0) 
    fwrite(outBuffer, 1, idxOut, Writter);

  finish_decode();
  doneinputtingbits();

  fclose(Writter);
  Free(MX);
  Free(name);
  for(n = 0 ; n < P[id].nModels ; ++n){
    if(P[id].model[n].type == REFERENCE)
      ResetCModelIdx(cModels[n]);
    else
      RemoveCModel(cModels[n]);
    }

  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(pModel[n]);
  Free(pModel);

  RemoveWeightModel(WM);
  Free(PT);
  Free(outBuffer);
  RemoveCBuffer(symBuf);
  RemoveAlphabet(AL);
  fclose(Reader);

  if(P->verbose == 1)
    fprintf(stderr, "Done!                          \n");  // SPACES ARE VALID 
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

  uint64_t  i = 0;

  if(P->verbose == 1)
    fprintf(stderr, "Building reference model ...\n");

  // BUILD ALPHABET
  ALPHABET *AL = CreateAlphabet(P->low);
  LoadAlphabet(AL, Reader);
  PrintAlphabet(AL);

  readerBuffer  = (uint8_t *) Calloc(BUFFER_SIZE + 1, sizeof(uint8_t));
  cModels       = (CModel **) Malloc(P->nModels * sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == REFERENCE)
      cModels[n] = CreateCModel(P->model[n].ctx, P->model[n].den, REFERENCE,
      P->model[n].edits, P->model[n].eDen, AL->cardinality, P->model[n].gamma);

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

      CalcProgress(nSymbols, ++i);
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
    fprintf(stderr, "Done!                          \n");  // SPACES ARE VALID  
  else
    fprintf(stderr, "                               \n");  // SPACES ARE VALID

  return cModels;
  }

  
//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t main(int argc, char *argv[]){
  char        **p = *&argv;
  CModel      **refModels; 
  uint32_t    n, k, *checksum, refNModels = 0, garbage, cardinality;
  Parameters  *P;
  FILE        *Reader = NULL;
  uint8_t     help, verbose, force, nTar = 1;
  clock_t     stop = 0, start = clock();
  
  if((help = ArgsState(DEFAULT_HELP, p, argc, "-h")) == 1 || argc < 2){
    PrintMenuD();
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_VERSION, p, argc, "-V")){
    PrintVersion();  
    return EXIT_SUCCESS;
    }

  verbose  = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v");
  force    = ArgsState  (DEFAULT_FORCE,   p, argc, "-f");

  for(n = 0 ; n != strlen(argv[argc-1]) ; ++n)
    if(argv[argc-1][n] == ':')
      ++nTar;

  P        = (Parameters *) Malloc(nTar * sizeof(Parameters));
  checksum = (uint32_t   *) Calloc(nTar , sizeof(uint32_t));

  P[0].force   = force;
  P[0].verbose = verbose;
  P[0].nTar    = ReadFNames (P, argv[argc-1]);
  P[0].ref     = ArgsString (NULL, p, argc, "-r");
  for(n = 0 ; n < nTar ; ++n){
    Reader = Fopen(P[0].tar[n], "r");
    startinputtingbits();
    start_decode(Reader);

    refNModels = 0;
    P[n].watermark = ReadNBits(32, Reader);
    if(P[n].watermark != WATERMARK){
      fprintf(stderr, "Error: Invalid compressed file to decompress!\n");
      return 1;
      }
    checksum[n]    = ReadNBits(46, Reader);
    P[n].size      = ReadNBits(46, Reader);
    ///////////////////////////////////////////////////////////
    P[n].low       = ReadNBits(32, Reader);
    uint32_t nLow  = ReadNBits(32, Reader);
    uint32_t x;
    for(x = 0 ; x < nLow ; ++x)
      garbage      = ReadNBits( 8, Reader);
    ///////////////////////////////////////////////////////////
    cardinality    = ReadNBits(16, Reader);
    for(k = 0 ; k < cardinality ; ++k)
      garbage      = ReadNBits(8,  Reader);
    P[n].nModels   = ReadNBits(16, Reader);
    P[n].model     = (ModelPar *) Calloc(P[n].nModels, sizeof(ModelPar));
    for(k = 0 ; k < P[n].nModels ; ++k){
      P[n].model[k].ctx   = ReadNBits( 5, Reader); 
      P[n].model[k].den   = ReadNBits(11, Reader); 
      P[n].model[k].gamma = ReadNBits(17, Reader) / 65536.0; 
      P[n].model[k].edits = ReadNBits( 7, Reader); 
      P[n].model[k].eDen  = ReadNBits( 9, Reader); 
      P[n].model[k].type  = ReadNBits( 1, Reader);
      if(P[n].model[k].type == 1)
        ++refNModels;
      }

    finish_decode();
    doneinputtingbits();
    fclose(Reader);
    }

  if(P->verbose)
    PrintArgs(P);
 
  if(refNModels > 0 && P[0].ref == NULL){
    fprintf(stderr, "Error: using reference model(s) in nonexistent "
    "reference sequence!\n");
    exit(1);
    }

  if(refNModels != 0)
    refModels = LoadReference(P);
  else
    refModels = (CModel **) Malloc(P->nModels * sizeof(CModel *));

  if(P->verbose && refNModels != 0)
    fprintf(stderr, "Checksum: %"PRIu64"\n", P->checksum); 

  for(n = 0 ; n < nTar ; ++n){
    if(refNModels != 0){
      if(CmpCheckSum(checksum[n], P[0].checksum) == 0)
        Decompress(P, refModels, n);
      }
    else
      Decompress(P, refModels, n);
    }

  stop = clock();
  fprintf(stderr, "Spent %g sec.\n", ((double)(stop-start))/CLOCKS_PER_SEC);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
