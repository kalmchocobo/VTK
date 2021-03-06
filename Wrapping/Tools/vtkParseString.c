/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkParseString.c

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*-------------------------------------------------------------------------
  Copyright (c) 2012 David Gobbi.

  Contributed to the VisualizationToolkit by the author in April 2012
  under the terms of the Visualization Toolkit 2008 copyright.
-------------------------------------------------------------------------*/

#include "vtkParseString.h"
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------
 * String utility methods
 *
 * Strings are centrally allocated and are const.  They should not
 * be freed until the parse is complete and all the data structures
 * generated by the parse have been freed.
 */

/* allocate a string of n+1 bytes */
void vtkParse_InitStringCache(StringCache *cache)
{
  cache->NumberOfChunks = 0;
  cache->Chunks = NULL;
  cache->ChunkSize = 0;
  cache->Position = 0;
}

/* allocate a string of n+1 bytes */
char *vtkParse_NewString(StringCache *cache, size_t n)
{
  size_t nextPosition;
  char *cp;

  if (cache->ChunkSize == 0)
    {
    cache->ChunkSize = 8176;
    }

  // align next start position on an 8-byte boundary
  nextPosition = (((cache->Position + n + 8) | 7 ) - 7);

  if (cache->NumberOfChunks == 0 || nextPosition > cache->ChunkSize)
    {
    if (n + 1 > cache->ChunkSize)
      {
      cache->ChunkSize = n + 1;
      }
    cp = (char *)malloc(cache->ChunkSize);

    /* if empty, alloc for the first time */
    if (cache->NumberOfChunks == 0)
      {
      cache->Chunks = (char **)malloc(sizeof(char *));
      }
    /* if count is power of two, reallocate with double size */
    else if ((cache->NumberOfChunks & (cache->NumberOfChunks-1)) == 0)
      {
      cache->Chunks = (char **)realloc(
        cache->Chunks, (2*cache->NumberOfChunks)*sizeof(char *));
      }

    cache->Chunks[cache->NumberOfChunks++] = cp;

    cache->Position = 0;
    nextPosition = (((n + 8) | 7) - 7);
    }

  cp = &cache->Chunks[cache->NumberOfChunks-1][cache->Position];
  cp[0] = '\0';

  cache->Position = nextPosition;

  return cp;
}

/* free all allocated strings */
void vtkParse_FreeStringCache(StringCache *cache)
{
  unsigned long i;

  for (i = 0; i < cache->NumberOfChunks; i++)
    {
    free(cache->Chunks[i]);
    }
  if (cache->Chunks)
    {
    free(cache->Chunks);
    }

  cache->Chunks = NULL;
  cache->NumberOfChunks = 0;
}

/* duplicate the first n bytes of a string and terminate it */
const char *vtkParse_CacheString(StringCache *cache, const char *in, size_t n)
{
  char *res = NULL;

  res = vtkParse_NewString(cache, n);
  strncpy(res, in, n);
  res[n] = '\0';

  return res;
}
