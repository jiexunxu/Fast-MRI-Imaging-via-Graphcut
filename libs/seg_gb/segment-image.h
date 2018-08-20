/*
Copyright (C) 2006 Pedro Felzenszwalb

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

#ifndef SEGMENT_IMAGE
#define SEGMENT_IMAGE

#include <cstdlib>
#include "image.h"
#include "misc.h"
#include "filter.h"
#include "segment-graph.h"

// random color
rgb random_rgb(){ 
  rgb c;
  
  c.r = (uchar)rand();
  c.g = (uchar)rand();
  c.b = (uchar)rand();

  return c;
}

// dissimilarity measure between pixels
static inline float diff(image<float> *r, image<float> *g, image<float> *b,
			 int x1, int y1, int x2, int y2) {
  return sqrt(square(imRef(r, x1, y1)-imRef(r, x2, y2)) +
	      square(imRef(g, x1, y1)-imRef(g, x2, y2)) +
	      square(imRef(b, x1, y1)-imRef(b, x2, y2)));
}

/*
 * Segment an image
 *
 * Returns a color image representing the segmentation.
 *
 * im: image to segment.
 * sigma: to smooth the image.
 * c: constant for treshold function.
 * min_size: minimum component size (enforced by post-processing stage).
 * num_ccs: number of connected components in the segmentation.
 */
void segment_image(image<rgb> *im, float sigma, float c, int min_size, int *num_ccs, uint32_t *output, int numRows) {
  int width = im->width();
  int height = im->height();

  image<float> *r = new image<float>(width, height);
  image<float> *g = new image<float>(width, height);
  image<float> *b = new image<float>(width, height);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      imRef(r, x, y) = imRef(im, x, y).r;
      imRef(g, x, y) = imRef(im, x, y).g;
      imRef(b, x, y) = imRef(im, x, y).b;
    }
  }

  // build graph
  edge *edges = new edge[width*height*3];
  int num = 0;
  int colTraversed=0;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      if (colTraversed==numRows)
        colTraversed=0;
      else
        colTraversed++;

      if (x < width-1) {
	edges[num].a = y * width + x;
	edges[num].b = y * width + (x+1);
	edges[num].w = diff(r, g, b, x, y, x+1, y);
	num++;
      }

      if ((y < height-1) && (colTraversed!=numRows)){
	edges[num].a = y * width + x;
	edges[num].b = (y+1) * width + x;
	edges[num].w = diff(r, g, b, x, y, x, y+1);
	num++;
      }

	 if(y<height-numRows){
     edges[num].a = y * width + x;
	edges[num].b = (y+numRows) * width + x;
	edges[num].w = diff(r, g, b, x, y, x, y+numRows);
	num++;
       }
    }
  }
  delete r;
  delete g;
  delete b;

  // segment
  universe *u = segment_graph(width*height, num, edges, c);
  
  // post process small components
  for (int i = 0; i < num; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
      u->join(a, b);
  }
  delete [] edges;
  *num_ccs = u->num_sets();
  
  // output the regions
  for (int y = 0; y < height-1; y++) {
    for (int x = 0; x < width-1; x++) {
      output[x * height + y] = (uint32_t)u->find(y * width + x);
    }
  }  

  delete u;

  return;
}

#endif