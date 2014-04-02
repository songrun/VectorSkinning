#ifndef _STAM_EIGEN_READER_H_
#define _STAM_EIGEN_READER_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef
   struct
   {
      double * val;
      double * vecI;
      double ** Phi;
   } EVALSTRUCT;

EVALSTRUCT ** read_eval ( int * pNmax );
void print_eval_N ( EVALSTRUCT ** ev, int N );
void print_eval ( EVALSTRUCT ** ev, int Nmax );

#ifdef __cplusplus
}
#endif

#endif // _STAM_EIGEN_READER_H_