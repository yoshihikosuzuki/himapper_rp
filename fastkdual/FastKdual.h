/********************************************************************************************
 *
 *  FastKdual: dual self/relative profile support for FastK.
 *
 *  This interface is deliberately separate from FastK.h so that the FastK executable and
 *  its source files retain their original behaviour.
 *
 ********************************************************************************************/

#ifndef _FASTK_DUAL
#define _FASTK_DUAL

#include "FastK.h"

extern char *DUAL_ROOT;

void FastKdual_Mark_Histogram_Touched(void);
void Sorting_Dual(char *path, char *root, char *relative_root);

#endif
