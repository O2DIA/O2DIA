#ifndef R_R_H
# include <R.h>
#endif

#ifndef R_EXT_DYNLOAD_H_
# include <R_ext/Rdynload.h>
#endif


#include <Rinternals.h>
#include <stdlib.h> // for NULL

/* register native routines ------------------------------------------------ */


/* NO .C calls */
/* NO .Call calls */
 
/* .Fortran calls */
void F77_NAME(inito2mod)  (void (* steadyparms)(int *, double *));
void F77_NAME(inito2forc) (void (* steadyforcs)(int *, double *));

void F77_NAME(o2mod)  (int *, double *, double *, double *, double *, int *);

R_FortranMethodDef FEntries[] = {
    {"inito2mod",      (DL_FUNC) &F77_SUB(inito2mod),      1},
    {"inito2forc",     (DL_FUNC) &F77_SUB(inito2forc),     1},
    {"o2mod",          (DL_FUNC) &F77_SUB(o2mod),          6},
    {NULL, NULL, 0}
};

/* Initialization ---------------------------------------------------------- */
void R_init_O2DIA(DllInfo *dll) {

  R_registerRoutines(dll, NULL, NULL, FEntries, NULL);

  // the following line protects against accidentially finding entry points

  R_useDynamicSymbols(dll, FALSE); // disable dynamic searching
}
