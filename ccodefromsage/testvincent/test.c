#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

void forward(double     c1,
             double     rho,
             int        iterations,
             int        samples,
             double     ll_x, 
             double     ll_y, 
             double     dx, 
             int        rows, 
             int        cols,
             uint64_t (*m)[cols]
             )
{
  int pixels = 0;
  int i, j, s, z;
  for( i = 0; i < rows; ++i )
    for( j = 0; j < cols; ++j )
    {
      if( ! (m[i][j] & 0xffffffff) )
        continue;

      pixels += 1;

      for( s = 0; s < samples; ++s )
      {
        double x = ll_x + j * dx + dx * rand() / RAND_MAX;
        double y = ll_y + i * dx + dx * rand() / RAND_MAX;

        for( z = 0; z < iterations; ++z )
        {
          double t = x;

          x = c1 - x*x + rho*y; 
          y = t;
        }

        int row = (y - ll_y + dx/2)/dx;
        int col = (x - ll_x + dx/2)/dx;

        if( 0 <= row && row < rows && 0 <= col && col < cols )
	  m[row][col] += UINT64_C(1) << 32;
      }
    }

  for( i = 0; i < rows; ++i )
    for( j = 0; j < rows; ++j )
      m[i][j] >>= 32;
  
  printf( "processed %u pixels\n", pixels );
}

