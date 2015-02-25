#pragma once

#include <Windows.h>
#include <tchar.h>
#include <atlconv.h>
#include <time.h>

/* disable c runtime warnings for "unsafe" strings*/
//#define _CRT_SECURE_NO_WARNINGS

// Number of "worker" threads, usually equal to the number of cpu cores minus 1.
#define N_WORKERS 1

// numerical constants
#define _PI  3.1415926
#define _TINY  1E-6
#define _BIG  1E10
#define Viscos  1.1E-5

#define CBIG  1e8
#define CSMALL  1e-6

// constants for DW hydraulic equations
#define A1  .314159265359e04  /* 1000*PI */
#define A2  0.157079632679e04  /* 500*PI  */
#define A3  0.502654824574e02  /* 16*PI   */
#define A4  6.283185307        /* 2*PI    */
#define A8  4.61841319859      /* 5.74*(PI/4)^.9 */
#define A9  -8.685889638e-01   /* -2/ln(10)      */
#define AA  -1.5634601348      /* -2*.9*2/ln(10) */
#define AB  3.28895476345e-03  /* 5.74/(4000^.9) */
#define AC  -5.14214965799e-03 /* AA*AB */


#define MAX_NET_ID_LEN  256 ///> Network node/link identification string
#define MAX_COMP_TYPE_STR_SIZE 64


// 
#define SQR(x) ((x)*(x))
#ifndef SGN
	#define SGN(x) ((x)>0 ? 1 : ( (x)==0? 0 :-1))
#endif

/// max size of error txt string
#define MAX_ERR_STRING_SIZE 1024

/// max size of error txt string
#define MAX_FILE_NAME_SIZE 1024


/// max size of error txt prefix
#define MAX_ERR_PREFIX_SIZE 128

/// initialize global error/warning/info (EWI) reporting settings
// Enable reporting to a log file
void init_ewi(TCHAR const * log_filepath);

/// common Error reporting function, 
/**default behavior is sending output to stderr. No initialization required
 . More reporting targets can be provided by modifying this function*/
void ewi(TCHAR const* txt, int thread_no=0);

/// end ewi system
void end_ewi();


