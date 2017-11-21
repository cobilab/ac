#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "levels.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// GET PARAMETERS FROM LEVELS
//
char *GetLevels(uint8_t l){
  switch(l){
    case 1: return LEVEL_1;
    case 2: return LEVEL_2;
    case 3: return LEVEL_3;
    case 4: return LEVEL_4;
    case 5: return LEVEL_5;
    default: fprintf(stderr, "[x] Unknown level!\n");
    exit(1);
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// PRINT LEVELS FROM LEVELS
//
void PrintLevels(void){
  fprintf(stderr, "Level%s\n", LEVEL_1);
  fprintf(stderr, "Level%s\n", LEVEL_2);
  fprintf(stderr, "Level%s\n", LEVEL_3);
  fprintf(stderr, "Level%s\n", LEVEL_4);
  fprintf(stderr, "Level%s\n", LEVEL_5);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

