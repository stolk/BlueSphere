#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <immintrin.h>

#include "pseudorand.h"
#include "threadpool.h"

#define NUMCONCURRENTTASKS 16

#define DUMPVERTS


typedef float floatx16 __attribute__((vector_size(4*16)));


#define MAXHIT	(1<<21)		// Maximum number of samples we will generate.
#define MAXCAN	(MAXHIT>>12)	// Maximum number of candidates that we will consider for each sample.

static float* listx;
static float* listy;
static float* listz;

static float* candx;
static float* candy;
static float* candz;

static float* candv;

static int sz=0;
static int numcand=0;


static int   choseni[ NUMCONCURRENTTASKS ];
static float chosenv[ NUMCONCURRENTTASKS ];


static threadpool_t* threadpool = 0;


void mkrandom(int i)
{
	float x,y,z;
	float ls;
	do
	{
		x = -1.0f + 2.0f * pseudorand_float();
		y = -1.0f + 2.0f * pseudorand_float();
		z = -1.0f + 2.0f * pseudorand_float();
		ls = x*x + y*y + z*z;
	} while( ls >= 1.0f );

	const float scl = 1.0f / sqrtf( ls );
	candx[i] = x * scl;
	candy[i] = y * scl;
	candz[i] = z * scl;
}


void promote(int i)
{
	listx[sz] = candx[i];
	listy[sz] = candy[i];
	listz[sz] = candz[i];
	sz++;
}


int evaluate_candidate_range(int s0, int s1)
{
	assert( s0 != s1 );
	const int numsets = s1 - s0;

	floatx16* cax = (floatx16*) ( candx + s0 * 16 );
	floatx16* cay = (floatx16*) ( candy + s0 * 16 );
	floatx16* caz = (floatx16*) ( candz + s0 * 16 );
	floatx16* cav = (floatx16*) ( candv + s0 * 16 );

	// Initialize all worst case dotprods as -1.
	for ( int i=0; i<numsets; ++i )
		cav[i] = _mm512_set1_ps( -1.0f );

	// Iterate over every promoted vector that is in the list.
	for ( int i=0; i<sz; ++i )
	{
		const floatx16 lx = _mm512_set1_ps( listx[i] );
		const floatx16 ly = _mm512_set1_ps( listy[i] );
		const floatx16 lz = _mm512_set1_ps( listz[i] );

		floatx16* cx = cax;
		floatx16* cy = cay;
		floatx16* cz = caz;
		floatx16* cv = cav;

		// Compute dot product with all candidates for this vector.
		for (int s=0; s<numsets; ++s)
		{
			floatx16 d = *cx * lx + *cy * ly + *cz * lz;
			*cv = _mm512_max_ps( *cv, d );
			cx++;
			cy++;
			cz++;
			cv++;
		}
	}

	// Now find the candidate with the best (lowest) value.
	const int c0 = s0*16;
	const int c1 = s1*16;

	float lowv = FLT_MAX;
	int lowi = -1;

	for ( int i=c0; i<c1; ++i )
	{
		lowi = candv[i] < lowv ?        i : lowi;
		lowv = candv[i] < lowv ? candv[i] : lowv;
	}
	//fprintf( stderr,"Best in range %d..%d = %d at %f\n", c0, c1, lowi, lowv );
	return lowi;
}


static void evaluate_slice( argument_t* arg )
{
        const int rank = (int) (long) arg->arg;
	const int slicesz = numcand / 256;
	assert( slicesz >= 1 );

	const int s0 = (rank+0) * slicesz;
	const int s1 = (rank+1) * slicesz;

	const int best = evaluate_candidate_range( s0, s1 );
	assert( best >= 0 );
	choseni[ rank ] = best;
	chosenv[ rank ] = candv[ best ];
}


static int evaluate_all_slices( void )
{
	for ( int i=0; i<NUMCONCURRENTTASKS; ++i )
		choseni[i] = -1;

	for ( int i=0; i<NUMCONCURRENTTASKS; ++i )
	{
		argument_t arg = { (void*)(long) i, NO_DELETE, 0 };
		threadpool_add( threadpool, evaluate_slice, arg, MEDIUM );
	}
	threadpool_wait( threadpool );
	float lowv = FLT_MAX;
	int lowi = -1;
	for ( int i=0; i<NUMCONCURRENTTASKS; ++i )
	{
		lowi = chosenv[i] < lowv ? choseni[i] : lowi;
		lowv = chosenv[i] < lowv ? chosenv[i] : lowv;
	}
	//fprintf( stderr,"Chose %d with value %f\n", lowi, lowv );
	assert( lowi >= 0 && lowi < numcand );
	return lowi;
}


#if defined(DUMPVERTS)
void dump_vertices( const char* fname, int count, float* dirsx, float* dirsy, float* dirsz )
{
	FILE* f = fopen( fname, "wb" );
	assert(f);

	struct per_photon
	{
		float pos[3];
		uint32_t clr;
	};

	const size_t photonsz = sizeof( struct per_photon );
	const size_t filesz = count * photonsz;
	fprintf( stderr, "Dumping 0x%x vertices of size %zu each, for a total of %zu bytes.\n", count, photonsz, filesz );
	struct per_photon* photons = (struct per_photon*) malloc( filesz );
	assert( photons );

	for ( int i=0; i<count; ++i )
	{
		photons[i].pos[0] = dirsx[i];
		photons[i].pos[1] = dirsy[i];
		photons[i].pos[2] = dirsz[i];
		photons[i].clr    = 0xff004060;
	}
	const int written = fwrite( photons, filesz, 1, f );
	assert( written == 1 );
	fclose(f);
}
#endif


void dump_floats( const char* fname, int count, float* v )
{
	FILE* f = fopen( fname, "wb" );
	assert(f);

	fwrite( v, sizeof(float), count, f );
	fclose(f);
}



int main(int argc, char* argv[])
{
	fprintf(stderr,"MAXHIT %d\n", MAXHIT);
	fprintf(stderr,"MAXCAN %d\n", MAXCAN);

	// Promoted candidates, part of the list.
	listx = (float*) aligned_alloc( 64, MAXHIT*sizeof(float) );
	listy = (float*) aligned_alloc( 64, MAXHIT*sizeof(float) );
	listz = (float*) aligned_alloc( 64, MAXHIT*sizeof(float) );

	// Candidate x/y/z coordinates.
	candx = (float*) aligned_alloc( 64, MAXCAN*sizeof(float) );
	candy = (float*) aligned_alloc( 64, MAXCAN*sizeof(float) );
	candz = (float*) aligned_alloc( 64, MAXCAN*sizeof(float) );

	// Candidate max dot products with list vectors.
	candv = (float*) aligned_alloc( 64, MAXCAN*sizeof(float) );

	threadpool = threadpool_create( NUMCONCURRENTTASKS );

	// Initialize a list with 1 coordinate in it.
	mkrandom(0);
	promote (0);

	while (sz < MAXHIT)
	{
		//fprintf( stderr, "%d ", sz );
		numcand = sz>>12;
		assert( numcand < MAXCAN );
		numcand = numcand < 256 ? 256 : numcand;
		while ( numcand & 255 )
			numcand++;

		for ( int i=0; i<numcand; ++i )
			mkrandom(i);
		fprintf(stderr,"Created %d candidates for sample %d...\n", numcand, sz);

		int chosen = evaluate_all_slices();
		promote( chosen );
	}

	threadpool_free( threadpool );

	// Save the results

	fprintf( stderr,"sz = %d\n", sz );

	dump_floats( "out.x", sz, listx );
	dump_floats( "out.y", sz, listy );
	dump_floats( "out.z", sz, listz );

	dump_vertices( "out.pts", sz, listx, listy, listz );

	return 0;
}

