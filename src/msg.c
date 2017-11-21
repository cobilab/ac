#include "msg.h"
#include <stdio.h>
#include <stdlib.h>

void ModelsExplanation(void){
  fprintf(stderr,
  "                                                                       \n"
  "  -rm <c>:<d>:<m/e>      reference context model (ex:-rm 13:100:0/0),  \n"
  "  -rm <c>:<d>:<m/e>      reference context model (ex:-rm 18:1000:1/50),\n"
  "  ...                                                                  \n"
  "  -tm <c>:<d>:<m/e>      target context model (ex:-tm 4:1:0/0),        \n"
  "  -tm <c>:<d>:<m/e>      target context model (ex:-tm 18:20:2/10),     \n"
  "  ...                                                                  \n"
  "                         target and reference templates use <c> for    \n"
  "                         context-order size, <d> for alpha (1/<d>),    \n"
  "                         <m> to the maximum sets the allowed mutations,\n"
  "                         on the context without being discarded (for   \n"
  "                         deep contexts), under the estimator <e>,      \n");
  } 

void PrintMenuD(void){
  fprintf(stderr,
  "Usage: AD [OPTION]... -r [FILE]  [FILE]:[...]                          \n"
  "Decompression of amino acid sequences (compressed by AC).              \n"
  "                                                                       \n"
  "Non-mandatory arguments:                                               \n"
  "                                                                       \n"
  "  -h                    give this help,                                \n"
  "  -v                    verbose mode (more information),               \n"
  "                                                                       \n"
  "  -r <FILE>             reference file,                                \n"
  "                                                                       \n"
  "Mandatory arguments:                                                   \n"
  "                                                                       \n"
  "  <FILE>                file to uncompress (last argument). For        \n"
  "                        more files use splitting \":\" characters.     \n"
  "                                                                       \n"
  "Report bugs to <pratas@ua.pt>.                                         \n");
  }

void PrintMenu(void){
  fprintf(stderr,
  "Usage: AC [OPTION]... -r [FILE]  [FILE]:[...]                          \n"
  "Compression of amino acid sequences.                                   \n"
  "                                                                       \n"
  "Non-mandatory arguments:                                               \n"
  "                                                                       \n"
  "  -h                     give this help,                               \n"
  "  -s                     show AC compression levels,                   \n"
  "  -v                     verbose mode (more information),              \n"
  "  -V                     display version number,                       \n"
  "  -f                     force overwrite of output,                    \n"
  "  -l <level>             level of compression [1;5] (lazy -tm setup),  \n"
  "  -t <threshold>         threshold frequency to discard from alphabet, \n"
  "  -g <gamma>             mixture decayment forgetting factor. It is    \n"
  "                         a real value in the interval [0;1),           \n");
  #ifdef ESTIMATE
  fprintf(stderr,
  "  -e                     it creates a file with the extension \".iae\" \n"
  "                         with the respective information content.      \n");
  #endif
  ModelsExplanation();
  fprintf(stderr,
  "                                                                       \n"
  "  -r <FILE>              reference file (\"-rm\" are loaded here),     \n"
  "                                                                       \n"
  "Mandatory arguments:                                                   \n"
  "                                                                       \n"
  "  <FILE>                 file to compress (last argument). For more    \n"
  "                         files use splitting \":\" characters.         \n"
  "                                                                       \n"
  "Report bugs to <pratas@ua.pt>.                                           \n");
  }


void PrintVersion(void){
  fprintf(stderr,
  "                                                                       \n"
  "                          ================                             \n"
  "                          |  AC  &  AD %u.%u |                        \n"
  "                          ================                             \n"
  "                                                                       \n"
  "Copyright (C) 2017-2018 University of Aveiro. This is a Free software. \n"
  "You may redistribute copies of it under the terms of the GNU - General \n"
  "Public License v3 <http://www.gnu.org/licenses/gpl.html>. There is NOT \n"
  "ANY WARRANTY, to the extent permitted by law. Developed and Written by \n"
  "Diogo Pratas.\n\n", VERSION, RELEASE);
  }

