#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "stam_eigen_reader.h"

#define IS_BIG_ENDIAN (!*(unsigned char *)&(uint16_t){1})

#define IX(i,j,n) ((i)+(n)*(j))


EVALSTRUCT ** read_eval ( int * pNmax )
{
	EVALSTRUCT ** ev;
	FILE * f;
	uint32_t Nmax, i, N, K;

	if ( !(f = fopen ( IS_BIG_ENDIAN ? "data/ccdata50.dat" : "data/ccdata50NT.dat", "r" )) ) {
		printf( "cannot open the Stam data file.\n" );
		return ( NULL );
	}

	fread ( &Nmax, sizeof(int32_t), 1, f );

	ev = (EVALSTRUCT **) malloc ( (Nmax-2)*sizeof(EVALSTRUCT *) );

	for ( i=0 ; i<Nmax-2 ; i++ )
	{
		N = i+3;
		K = 2*N+8;

		ev[i] = (EVALSTRUCT *) malloc ( sizeof(EVALSTRUCT) );
		ev[i]->val = (double *) malloc ( K*sizeof(double) );
		ev[i]->vecI = (double *) malloc ( K*K*sizeof(double) );
		ev[i]->Phi = (double **) malloc ( 3*sizeof(double *) );
		ev[i]->Phi[0] = (double *) malloc ( K*16*sizeof(double) );
		ev[i]->Phi[1] = (double *) malloc ( K*16*sizeof(double) );
		ev[i]->Phi[2] = (double *) malloc ( K*16*sizeof(double) );

		fread ( ev[i]->val, sizeof(double), K, f );
		fread ( ev[i]->vecI, sizeof(double), K*K, f );
		fread ( ev[i]->Phi[0], sizeof(double), K*16, f );
		fread ( ev[i]->Phi[1], sizeof(double), K*16, f );
		fread ( ev[i]->Phi[2], sizeof(double), K*16, f );
	}

	fclose ( f );

	*pNmax = Nmax;

	return ( ev );
}

void print_eval_N ( EVALSTRUCT ** ev, int N ) 
{
	int i, j, k, l, K;
	K = 2*N+8;
	i = N-3; 
	
	printf ( "hello world" );
	printf ( "N=%d\n\n", N );

	printf ( "Eigenvalues:\n\n" );
	for ( j=0 ; j<K ; j++ ) printf ( "%6.3f ", ev[i]->val[j] );
	printf ( "\n\n" );

	printf ( "Inverse of Eigenvectors:\n\n" );

	for ( j=0 ; j<K ; j++ )
	{
		for ( k=0 ; k<K ; k++ ) printf ( "%6.3f ", ev[i]->vecI[IX(j,k,K)] );
		printf ( "\n" );
	}
	printf ( "\n\n" );

	printf ( "Coefficients of the Eigenbasis functions:\n" );

	for ( k=0 ; k<3 ; k++ )
	{
		printf ( "k=%d:\n", j );
		for ( j=0 ; j<K ; j++ )
		{
			for ( l=0 ; l<16 ; l++ ) printf("%6.3f ",ev[i]->Phi[k][IX(j,l,K)]);
				printf ( "\n" );
		}
		printf ( "\n" );
	}
}

void print_eval ( EVALSTRUCT ** ev, int Nmax )
{
	int i, N;
	printf ( "Nmax = %d\n\n", Nmax );

	for ( i=0 ; i<Nmax-2 ; i++ )
	{
		N = i+3;
		print_eval_N( ev, N );
	}

	return;
}