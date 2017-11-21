#ifndef LEVELS_H_INCLUDED
#define LEVELS_H_INCLUDED

#include "defs.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
#define LEVEL_1 " 1: -rm 7:1:0/0 "
#define LEVEL_2 " 2: -rm 10:500:1/100 -g 0.9 "
#define LEVEL_3 " 3: -rm 3:1:0/0 -rm 12:200:1/10 -g 0.9 "
#define LEVEL_4 " 4: -rm 1:1:0/0 -rm 3:1:0/0 -rm 12:200:1/10 -g 0.9 "
#define LEVEL_5 " 5: -rm 1:1:0/0 -rm 3:1:0/0 -rm 5:1:0/0 -rm 13:500:2/10 -g 0.9 "

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char    *GetLevels  (uint8_t);
void    PrintLevels (void);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

