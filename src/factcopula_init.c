#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void gauleg(void *, void *, void *, void *, void *);

/* .Fortran calls */
extern void F77_NAME(bb11fact)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(bb1m1f1t)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(bifactnllk)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(gum1fact)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(latupdate)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(latupdate3)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(latupdatebifact3)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lfrk1fact)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pfactnllk)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strbb1frk2)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strbb1frk2nllk)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strbb1gum2)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strbb1gum2nllk)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strfrk1)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strfrk2)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strfrk2nllk)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strfrkbb1)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strfrkgum1)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strgum1)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strgum2)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strgum2nllk)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strgumbb1)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strgumfrk2)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strgumfrk2nllk)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strt1)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strt2)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strt2nllk)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(strtbb1)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(t1fact)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(tbifactnllk)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(tpfactnllk)(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"gauleg", (DL_FUNC) &gauleg, 5},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"bb11fact",         (DL_FUNC) &F77_NAME(bb11fact),         11},
    {"bb1m1f1t",         (DL_FUNC) &F77_NAME(bb1m1f1t),         14},
    {"bifactnllk",       (DL_FUNC) &F77_NAME(bifactnllk),        8},
    {"gum1fact",         (DL_FUNC) &F77_NAME(gum1fact),         11},
    {"latupdate",        (DL_FUNC) &F77_NAME(latupdate),         9},
    {"latupdate3",       (DL_FUNC) &F77_NAME(latupdate3),        9},
    {"latupdatebifact3", (DL_FUNC) &F77_NAME(latupdatebifact3), 13},
    {"lfrk1fact",        (DL_FUNC) &F77_NAME(lfrk1fact),        11},
    {"pfactnllk",        (DL_FUNC) &F77_NAME(pfactnllk),         7},
    {"strbb1frk2",       (DL_FUNC) &F77_NAME(strbb1frk2),       13},
    {"strbb1frk2nllk",   (DL_FUNC) &F77_NAME(strbb1frk2nllk),   11},
    {"strbb1gum2",       (DL_FUNC) &F77_NAME(strbb1gum2),       13},
    {"strbb1gum2nllk",   (DL_FUNC) &F77_NAME(strbb1gum2nllk),   11},
    {"strfrk1",          (DL_FUNC) &F77_NAME(strfrk1),          13},
    {"strfrk2",          (DL_FUNC) &F77_NAME(strfrk2),          13},
    {"strfrk2nllk",      (DL_FUNC) &F77_NAME(strfrk2nllk),      11},
    {"strfrkbb1",        (DL_FUNC) &F77_NAME(strfrkbb1),        13},
    {"strfrkgum1",       (DL_FUNC) &F77_NAME(strfrkgum1),       13},
    {"strgum1",          (DL_FUNC) &F77_NAME(strgum1),          13},
    {"strgum2",          (DL_FUNC) &F77_NAME(strgum2),          13},
    {"strgum2nllk",      (DL_FUNC) &F77_NAME(strgum2nllk),      11},
    {"strgumbb1",        (DL_FUNC) &F77_NAME(strgumbb1),        13},
    {"strgumfrk2",       (DL_FUNC) &F77_NAME(strgumfrk2),       13},
    {"strgumfrk2nllk",   (DL_FUNC) &F77_NAME(strgumfrk2nllk),   11},
    {"strt1",            (DL_FUNC) &F77_NAME(strt1),            14},
    {"strt2",            (DL_FUNC) &F77_NAME(strt2),            14},
    {"strt2nllk",        (DL_FUNC) &F77_NAME(strt2nllk),        12},
    {"strtbb1",          (DL_FUNC) &F77_NAME(strtbb1),          15},
    {"t1fact",           (DL_FUNC) &F77_NAME(t1fact),           12},
    {"tbifactnllk",      (DL_FUNC) &F77_NAME(tbifactnllk),       9},
    {"tpfactnllk",       (DL_FUNC) &F77_NAME(tpfactnllk),        8},
    {NULL, NULL, 0}
};

void R_init_FactorCopulaModel(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
