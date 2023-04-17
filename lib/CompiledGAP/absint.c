/*Written by Paul Smith*/

#include "src/compiled.h"  /* the GAP headers */
#include <stdlib.h>        /* for abs */


/***************** The new GAP kernel functions ***************/

/** 
GAP kernel C function to calculate ane return the absolute value of an integer. 
This is assumed to be a GAP small integer (i.e. less than 2^28-1 or 2^60-1 
depending on whether the machine is 32- or 64-bit). 
@param self The standard GAP first parameter 
@param n The standard GAP first parameter 
@return The maximum small integer (according to this kernel module).
**/
Obj FuncABSINT_HAP(Obj self, Obj n)
{
  Int Cn;
  Cn = INT_INTOBJ(n);     /* Convert the GAP object n into a C integer */
  Cn = abs(Cn);           /* Get the absolute value of this integer */
  return INTOBJ_INT(Cn);  /* Convert it back to a GAP object and return it */
}


/******************** The interface to GAP ***************/

/**
Details of the functions to make available to GAP. 
This is used in InitKernel() and InitLibrary()
*/
static StructGVarFunc GVarFuncs[] = 
{
  {"AbsIntt_HAP",      /* The function name in GAP */
   1,                 /* The number of parameters */
   "n",               /* The names of the parameters */
   FuncABSINT_HAP,    /* The C function to call */
   "absint.c:FuncABSINT_HAP"  /* A user-friendly description of where
                                         this function is */
  }, 

  { 0 } /* Finish with an empty entry */
};



/**
The first function to be called when the library is loaded by the kernel.
**/
static Int InitKernel(StructInitInfo* module)
{
  /* init filters and functions                                          */
  InitHdlrFuncsFromTable( GVarFuncs );
  
  /* return success                                                      */
  return 0;
}


/**
The second function to be called when the library is loaded by the kernel.
**/
static Int InitLibrary(StructInitInfo* module)
{
    /* init filters and functions                                          */
    InitGVarFuncsFromTable( GVarFuncs );

    /* return success                                                      */
    return 0;
}


/**
Information about this library, returned when the library is loaded, 
for example by Init__Dynamic(). This contains details of the library name,
and the further initialisation functions to call.
**/
static StructInitInfo module = {
#ifdef STATICMODULE
 /* type        = */ MODULE_STATIC,
#else
 /* type        = */ MODULE_DYNAMIC,
#endif
 /* name        = */ "absolute value of an integer",
 /* revision_c  = */ 0,
 /* revision_h  = */ 0,
 /* version     = */ 0,
 /* crc         = */ 0,
 /* initKernel  = */ InitKernel,
 /* initLibrary = */ InitLibrary,
 /* checkInit   = */ 0,
 /* preSave     = */ 0,
 /* postSave    = */ 0,
 /* postRestore = */ 0
};


#ifndef STATICGAP
/** 
Function called by GAP as soon as the library is dynamically loaded. 
This returns the StructInitInfo data for this library
**/
StructInitInfo * Init__Dynamic (void)
{
 return &module;
}
#endif
StructInitInfo * Init__linbox(void)
{
  return &module;
}

