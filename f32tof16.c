#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>	// for strstr()
#include <math.h>	// for sqrtf()
#include <immintrin.h>	// for _mm256_cvtps_ph


void convert_to_half_floats
(
	int count,
	float* __restrict__ src,
	short* __restrict__ dst		// really fp16, not short ints.
)
{
	assert( ( count & 7 ) == 0 );
	const int numbatches = count/8;
	for ( int b=0; b<numbatches; ++b )
	{
		// Load up 8 floating points values.
		const __m256  v32 = _mm256_load_ps( src + 8*b );
		// Convert to 16 bit half floats.
		const __m128i v16 = _mm256_cvtps_ph( v32, 0 );
		// Write into destination memory.
		__m128i* tgt = (__m128i*) ( dst + 8*b );
		_mm_store_si128( tgt, v16 );
	}
}


void convert_to_full_floats
(
	int count,
	short* __restrict__ src,	// really fp16, not short ints.
	float* __restrict__ dst
)
{
	assert( ( count & 7 ) == 0 );
	const int numbatches = count/8;
	for ( int b=0; b<numbatches; ++b )
	{
		// Load up 8 __half floats.
		const __m128i v16 = _mm_load_si128( (__m128i*) (src + 8*b) );
		// Convert to 32bit floats.
		const __m256  v32 = _mm256_cvtph_ps( v16 );
		// Write to destination memory.
		_mm256_store_ps( dst + 8*b, v32 );
	}
}


void random_shuffle( int count, float* dx, float* dy, float* dz )
{
	for ( int i=0; i<count; ++i )
	{
		int j = rand() % count;
		float t0 = dx[ i ];
		float t1 = dy[ i ];
		float t2 = dz[ i ];
		dx[ i ] = dx[ j ];
		dy[ i ] = dy[ j ];
		dz[ i ] = dz[ j ];
		dx[ j ] = t0;
		dy[ j ] = t1;
		dz[ j ] = t2;
	}
	fprintf( stderr, "Dirs are randomly shuffled.\n" );
}




int main(int argc, char* argv[])
{
	if ( argc!=6 )
	{
		fprintf( stderr, "Usage: %s count coords.x coords.y coords.z outfile\n", argv[0] );
		return -1;
	}

	const int count = atoi( argv[1] );
	assert(count>1);
	const char* oname = argv[ 5 ];

	float* coords32[3] = {0,0,0};
	short* coords16[3] = {0,0,0};
	for ( int i=0; i<3; ++i )
	{
		coords32[i] = aligned_alloc( 64, count*sizeof(float) );
		coords16[i] = aligned_alloc( 64, count*sizeof(short) );
		const char* iname = argv[2+i];
		FILE* f = fopen( iname, "rb" );
		if ( !f )
		{
			fprintf( stderr, "Failed to open '%s' for reading.\n", iname );
			return -1;
		}
		const int numread = fread( coords32[i], sizeof(float), count, f );
		assert( numread == count );
		fclose(f);
		convert_to_half_floats( count, coords32[i], coords16[i] );
	}

	FILE* f = fopen( oname, "wb" );
	if ( !f )
	{
		fprintf( stderr, "Failed to open '%s' for writing.\n", oname );
		return -2;
	}

	for ( int i=0; i<3; ++i )
	{
		const int numwritten = fwrite( coords16[i], sizeof(short), count, f );
		assert( numwritten == count );
	}

	fclose(f);

	fprintf( stderr, "Wrote %d coordinates as 16bit floats in SoA form to: %s\n", count, oname );

	return 0;
}


