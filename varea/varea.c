/* file: varea.c */
/* code by Dr. C.S. Botzler */
/* if you use this code for your work, please cite: */
/*   Botzler, 2004, PhD Thesis, Ludwig-Maximilian University, Munich */
/*   and Shewchuk, 1996, in First Workshop on Applied Computational */
/*   Geometry, ACM (http://www.cs.cmu.edu/~quake/triangle.html) */

/* May 2007: Changed allocation of space for nodes to be dynamic */

/* varea - Find the nucleus and compute the area of voronoi cells
 *
 * Compile: make
 *
 * Uses input from triangle **Add some ref here**
 *
 * Syntax:
 * varea fieldname fieldsize nmc
 *    nmc: Number of Monte-Carlo iterations
 *
 * Needs 1 inputfile:
 *
 * fieldname.nuclei:
 *       Format: id x y z MUNICS_ID
 *
 * $Id: varea.c,v 1.9 2001/12/13 19:01:19 snigula Exp $
 */

//#define DEBUG 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>

/* create (double) random numbers within [0.0, 1.0] */
#define drand() ((double) rand() / (RAND_MAX + 1.0))

#define REAL double
#include "triangle.h"

//#define MAXNODES 200

struct nucleus 
{
  double x;
  double y;
  double z;
  int munid;
};

struct node 
{
  double x;
  double y;
  int numedges;
  int valid;
};
struct node * vnodes;
int n_vnodes;

struct edge
{
  int from;
  int to;
};
struct edge * vedges;
int n_vedges;

struct vcell
{
  int nucleus;
  int numnodes;
  int n_nodes;
  int * nodes;
  double * distances;
  double area;
  int on_edge; /* 0: in, 1: out */
};
struct vcell * vcells;

int maxnodes;
const int nodeblock = 10; // Number of nodes to be allocated

void add_node( index, node, distance )
{
  if( vcells[index].numnodes + 1 > vcells[index].n_nodes ) { // Reallocate more nodespace!
    vcells[index].n_nodes+=nodeblock;
    vcells[index].nodes = (int *)realloc( vcells[index].nodes, vcells[index].n_nodes*sizeof( int ) );
    vcells[index].distances = (double *)realloc( vcells[index].distances, vcells[index].n_nodes*sizeof( double ) );
  }

  vcells[index].nodes[vcells[index].numnodes]=node;
  vcells[index].distances[vcells[index].numnodes]=distance;
  vcells[index].numnodes++;

  if( vcells[index].n_nodes > maxnodes ) maxnodes = vcells[index].n_nodes;
  
}

int fsize;
double ** distances;
double ** tmpdistances;

void calc_vor (struct nucleus* nuclei, int n_nuclei, int verbose);
int line_counter (FILE * fp1);
void write_files( char * fname, int n_nuclei );

void usage (void) 
{
  fprintf (stderr, "varea fieldname fieldsize nmc\n");
  fprintf (stderr, "  nmc: Number of Monte-Carlo iteration\n");
}

int main (int argc, char* argv[])
{
  struct nucleus * nuclei;
  int n_nuclei;
  struct nucleus * mcnuclei;

  int nmc, n_areas;
  double tmp_area;
  double* mcdensities;
  double x,y, tmpx, tmpy;

  int new_obj;
   
  FILE* nuc;
  FILE* out;
  FILE* dbgout;

  int i,j,k,count;
  double zlow, zhigh;

  char fname[256];
  char froot[256];

  srand(time(NULL));

  if (argc != 4)
    {
      usage ();
      exit (EXIT_FAILURE);
    }
   
  fsize = atoi( argv[2] );
  nmc = atoi( argv[3] );

  /* Read in nuclei information */
   
  sprintf (fname, "%s.nuclei", argv[1]);

  nuc = fopen (fname, "r");
   
  n_nuclei = line_counter (nuc);
  rewind (nuc);
   
  nuclei = (struct nucleus *) malloc (n_nuclei * sizeof (struct nucleus));
  mcnuclei = (struct nucleus *) malloc (n_nuclei * sizeof (struct nucleus));
  vcells = (struct vcell *) malloc (n_nuclei * sizeof (struct vcell));

  memset( nuclei, 0, n_nuclei * sizeof (struct nucleus) );
  memset( mcnuclei, 0, n_nuclei * sizeof (struct nucleus) );
  memset( vcells, 0, n_nuclei * sizeof (struct vcell) );

  fprintf (stderr, "Reading %d nuclei...", n_nuclei);
   
  for (i=0; i<n_nuclei; i++)
    {
      // Initialize vcells

      vcells[i].numnodes=0;
      vcells[i].nodes= (int*) malloc( nodeblock * sizeof( int ) );
      vcells[i].distances= (double*) malloc( nodeblock * sizeof( double ) );
      vcells[i].n_nodes=nodeblock;

      if (fscanf (nuc, "%lf %lf %lf %d \n", 
                  & nuclei[i].x,
                  & nuclei[i].y,
                  & nuclei[i].z,
                  & nuclei[i].munid) == EOF )
        {
          fprintf (stderr, "\nFormat error in file %s! \n", fname);
          exit (-1);
        }       
    }
   
  fprintf (stderr, "done\n" );
  fclose( nuc );

  /* Monte-Carlo */
  
  fprintf (stderr, "Estimating background using Monte-Carlo simulations..." );
  mcdensities = (double*) malloc (nmc*sizeof (double));
  
  sprintf( fname, "%s.vor.areas", argv[1] );
  out = fopen (fname, "w");
  
  for (i=0; i<nmc; i++)
    {
      /* Generate random data */
      for (j=0; j<n_nuclei; j++)
        {          
          mcnuclei[j].x=drand()*fsize;
          mcnuclei[j].y=drand()*fsize;
        }
      
      mcdensities[i] = 0;
      n_areas = 0;

      calc_vor (mcnuclei, n_nuclei, 0);

      for (j=0; j<n_nuclei; j++)
        {
          if (!vcells[j].on_edge)
            {
              if (vcells[j].area < 0) {
                fprintf (stderr, "\nAy Caramba!! %f\n", vcells[j].area);
                sprintf( fname, "%s.mc.%d", argv[1], (int)(drand()*1000) );
                fprintf (stderr, "Saving to %s.*\n", fname);
                write_files( fname, n_nuclei );
                sprintf (fname, "%s.nuclei", fname );
                dbgout = fopen( fname, "w" );
                for (k=0; k<n_nuclei; k++)
                  {
                    fprintf (dbgout, "%f\t%f\n", mcnuclei[k].x, mcnuclei[k].y);
                  } 
                fflush (dbgout);
                fclose (dbgout);
              } else { 
                mcdensities[i] += 1./vcells[j].area;
                n_areas++;
              }
            }
          vcells[j].numnodes = 0;
        }
      
      if( n_areas )
	mcdensities[i] = mcdensities[i] / (double)n_areas;
      
      
      for (j=0; j<n_vnodes; j++)
        free(distances[j]);
      free(distances);
      
#ifdef DEBUG
      for (j=0; j<n_vnodes; j++)
        free(tmpdistances[j]);
      free(tmpdistances);
#endif
      
      free (vnodes);
      n_vnodes = 0;
      free (vedges);
      n_vedges = 0;

      if( mcdensities[i] > 1. ) {
	// We had a nan, discard run
	--i; // Do it again
      } else {
	fprintf (out, "%.15f\n", mcdensities[i]);
      }
    }
  
  fflush (out);
  fclose (out);
  
  fprintf (stderr, "done\n" );
  
  calc_vor( nuclei, n_nuclei, 1);
  
  /* Write results to file */
  
  write_files( argv[1], n_nuclei );

  return 0;
}

/* counts the lines of the source file */

int line_counter (FILE * fp1)
{
  int c;
  int counter = 0;

  while ((c = fgetc(fp1)) != EOF)
    {
      if (c == '\n') 
        counter ++;     
    }
  
  return counter;
}

void write_files( char * fnam, int n_nuclei )
{
  FILE* nuc;
  FILE* vnode_out;
  FILE* vedge_out;

  int i,j,k;

  char fname[256];

  fprintf (stderr, "Writing output files...");
  
  sprintf( fname, "%s.vor.edge", fnam );
  vedge_out = fopen (fname, "w");
  
  for (i=0; i<n_vedges; i++)
    {
      fprintf (vedge_out, "%d %d %d\n", i, vedges[i].to, vedges[i].from);
    }

  fflush (vedge_out);
  fclose (vedge_out);

  sprintf( fname, "%s.vor.node", fnam );
  vnode_out = fopen (fname, "w");
  
  for (i=0; i<n_vnodes; i++)
    {
      fprintf (vnode_out, "%d %f %f\n", i, vnodes[i].x, vnodes[i].y);
    }

  fflush (vnode_out);
  fclose (vnode_out);

#ifdef DEBUG
  sprintf( fname, "%s.vor.distances", fnam );
  vnode_out = fopen (fname, "w");
  
  fprintf (vnode_out, "\t\t");

  for (i=0; i<n_nuclei; i++)
    {
      fprintf (vnode_out, "%d\t\t", i);
    }
  
  fprintf (vnode_out, "\n");
      
  for (i=0; i<n_vnodes; i++)
    {
      fprintf (vnode_out, "%d\t\t", i);
      for (j=0; j<n_nuclei; j++)
        {
          fprintf (vnode_out, "%f\t", tmpdistances[i][j]);
        }
      fprintf (vnode_out, "\n");
    }

  fflush (vnode_out);
  fclose (vnode_out);
#endif

  sprintf( fname, "%s.vor.cells", fnam );
  vnode_out = fopen (fname, "w");
  
  for (i=0; i<n_nuclei; i++)
    {
      fprintf (vnode_out, "%d %d %d %f\t\t", 
               i, 
               vcells[i].numnodes, 
               vcells[i].on_edge,
               vcells[i].area );

#ifdef DEBUG
      for( j=0; j<vcells[i].numnodes; j++ )
        {
          fprintf (vnode_out, "%d\t\t", vcells[i].nodes[j]);
        }
      for( j=vcells[i].numnodes; j<vcells[i].n_nodes; j++ )
        {
          fprintf (vnode_out, "-1\t\t" );
        }
      for( j=0; j<vcells[i].numnodes; j++ )
        {
          fprintf (vnode_out, "%f\t\t", vcells[i].distances[j]);
        }
      for( j=vcells[i].numnodes; j<vcells[i].n_nodes; j++ )
        {
          fprintf (vnode_out, "-1\t\t" );
        }
#endif
      fprintf (vnode_out, "\n");
    }

  fflush (vnode_out);
  fclose (vnode_out);

  fprintf (stderr, "done\n" );

}



void calc_vor (struct nucleus* nuclei, int n_nuclei, int verbose)
{
  struct triangulateio tri_in;
  struct triangulateio tri_out;
  struct triangulateio tri_vor;

  double tmp1, tmp2;
  double mindist;
  int mem;
  int i, j, k;

  int test;
   
  int firstnode;
  int secondnode;
  int currnode;
  int nextnode;
  int cont;
  //  int tmpnodes[MAXNODES];
  int * tmpnodes;
  int nodes_visited;

  double xmin, xmax, x, y, tmp_xmin;
  int left, right;

  double convexa, convexb;
  
  int n_upper, n_lower;
/*   int upper[MAXNODES]; */
/*   int lower[MAXNODES]; */
  int * upper;
  int * lower;
  int lastnode;
  double a, b;
  double xlow, xhigh, ylow, yhigh, quot;

  if (verbose)
     fprintf (stderr, "Calculating voronoi cells..." );

  tri_in.pointlist = (double *) malloc (n_nuclei*2*sizeof(double));
  tri_in.pointattributelist = NULL;
  tri_in.numberofpoints = n_nuclei;
  tri_in.numberofpointattributes = 0;

  tri_out.pointlist=NULL;

  tri_vor.pointlist=NULL;
  tri_vor.edgelist=NULL;
  tri_vor.normlist=NULL;

  for (i=0, j=0; i<n_nuclei; i++, j++)
    {
      tri_in.pointlist[j] = nuclei[i].x;
      tri_in.pointlist[++j] = nuclei[i].y;
      tri_in.pointmarkerlist = NULL;
    }

  triangulate ("zQNEv", &tri_in, &tri_out, &tri_vor);

  if (verbose)
  fprintf (stderr, "done\n");

  n_vnodes = tri_vor.numberofpoints;
  n_vedges = tri_vor.numberofedges;

  vnodes = (struct node *) malloc (n_vnodes * sizeof (struct node));
  vedges = (struct edge *) malloc (n_vedges * sizeof (struct edge));
  
  for (i=0, j=0; i<n_vnodes; i++, j++)
    {
       vnodes[i].x = tri_vor.pointlist[j];
       vnodes[i].y = tri_vor.pointlist[++j];

       if (vnodes[i].x < 0 || vnodes[i].x > fsize ||
           vnodes[i].y < 0 || vnodes[i].y > fsize)
          vnodes[i].valid=0;
       else
          vnodes[i].valid=1;   

       vnodes[i].numedges=0;
    }
  
  for (i=0, j=0; i<n_vedges; i++, j++)
    {
       vedges[i].from = tri_vor.edgelist[j];
       vedges[i].to = tri_vor.edgelist[++j];

       if ( vedges[i].to != -1 )
	 if (!vnodes[vedges[i].to].valid)
	   vedges[i].to = -1;

       if ( vedges[i].from != -1 )
	 {
	   if (!vnodes[vedges[i].from].valid)
	     {
	       vedges[i].from = vedges[i].to;
	       vedges[i].to = -1;
	     }
	 }

       if ( vedges[i].to > -1 )
	 vnodes[vedges[i].to].numedges++;
       if ( vedges[i].from > -1 )
	 vnodes[vedges[i].from].numedges++;
    }

  /* Free unused space */

  free (tri_in.pointlist);
  free (tri_out.pointlist);
  free (tri_vor.pointlist);
  free (tri_vor.edgelist);
  free (tri_vor.normlist);
  
  /* Test node consistency */
  
  if (verbose)
  fprintf (stderr, "Testing node consistency...");

  test = 0;
  for (i=0; i<n_vnodes; i++) {
    for (j=i+1; j<n_vnodes; j++) {
      if( fabs(vnodes[i].x-vnodes[j].x)<DBL_EPSILON &&
          fabs(vnodes[i].y-vnodes[j].y)<DBL_EPSILON )
        test++;
    }
  }

  if (test) {
     if (verbose)
        fprintf (stderr, "failed for %d of %d nodes!\nExiting...\n", test, n_vnodes);
     else
        fprintf (stderr, "Consistency check failed for %d of %d nodes!\nExiting...\n", test, n_vnodes);
     exit (-1);
  } else {
     if (verbose)
        fprintf (stderr, "done\n");
  }

  /* Find the voronoi cells */

  if (verbose)
     fprintf (stderr, "Constructing voronoi cells...");

  distances = (double**) malloc (n_vnodes * sizeof (double*));
  for (i=0; i<n_vnodes; i++) {
     distances[i] = (double*) malloc (n_nuclei * sizeof (double));
  }

#ifdef DEBUG
  tmpdistances = (double**) malloc (n_vnodes * sizeof (double*));
  for (i=0; i<n_vnodes; i++) {
     tmpdistances[i] = (double*) malloc (n_nuclei * sizeof (double));
  }
#endif

  for (i=0; i<n_vnodes; i++) {
     for (j=0; j<n_nuclei; j++) {
        tmp1 = vnodes[i].x - nuclei[j].x;
        tmp2 = vnodes[i].y - nuclei[j].y;
        distances[i][j] = sqrt( tmp1*tmp1 + tmp2*tmp2 );
#ifdef DEBUG
        tmpdistances[i][j] = sqrt( tmp1*tmp1 + tmp2*tmp2 );
#endif
     }
  }

  for (i=0; i<n_vnodes; i++) 
    {
      for (j=1; j<=vnodes[i].numedges; j++)
        {
          mindist = DBL_MAX;
          mem = -1;
          for (k=0; k<n_nuclei; k++)
            {
              if (distances[i][k]<mindist) 
                {
                  mindist = distances[i][k];
                  mem = k;
                }
            }

	  if ( mem == -1 ) exit(-1);

/* 	  fprintf( stderr, "%d %d\n", mem, vcells[mem].numnodes ); */
/*           vcells[mem].distances[vcells[mem].numnodes] = distances[i][mem]; */
/*           vcells[mem].nodes[vcells[mem].numnodes] = i; */
/* 	  if ( mem == 721 ) fprintf( stderr, "%d\n", vcells[mem].numnodes ); */
/*           vcells[mem].numnodes++; */
	  add_node( mem, i, distances[i][mem] );
          vcells[mem].nucleus=mem;
          
          distances[i][mem]=DBL_MAX;          
        }
    }


  if (verbose)
     fprintf (stderr, "done\n" );

  // At this point we know the maximum number of nodes used, so we can create the arrays
  if (verbose)
    fprintf( stderr, "Using a maximum number of %d nodes\n", maxnodes );
  upper = (int*)malloc( maxnodes * sizeof(int) );
  lower = (int*)malloc( maxnodes * sizeof(int) );
  tmpnodes = (int*)malloc( maxnodes * sizeof(int) );

  /* Find edge nuclei */

  if (verbose)
     fprintf (stderr, "Removing edge nuclei...");

  for (i=0; i<n_nuclei; i++ )
    {

      memset( tmpnodes, 0, maxnodes*sizeof(int) );

      nodes_visited = 0;
      /* Copy the node information */
      cont = 0;
      for (j=0; j<vcells[i].numnodes; j++)
        {
          tmpnodes[j] = vcells[i].nodes[j];
        }

      for (j=vcells[i].numnodes; j<vcells[i].n_nodes; j++)
        {
          tmpnodes[j] = -2;
        }

      secondnode = -2;
      firstnode = tmpnodes[0];
      currnode = firstnode;

      while (1) 
        {
          nextnode = -2;
          for (j=0; j<n_vedges; j++)
            {
              if( vedges[j].from == currnode ) 
                {
                  for (k=0; k<vcells[i].numnodes; k++)
                    {
                      if( vedges[j].to == tmpnodes[k] )
                        {
                          nextnode = tmpnodes[k];
                          if( currnode == secondnode && nextnode == firstnode ) { 
                            nextnode = -2;
                            continue;
                          }
                          
                          tmpnodes[k] = -2;
                          nodes_visited++;
                          if( secondnode == -2 ) secondnode=nextnode;
                          break;
                        }
                    }
                  if( nextnode != -2 ) break;
                }
              if( vedges[j].to == currnode ) 
                {
                  for (k=0; k<vcells[i].numnodes; k++)
                    {
                      if( vedges[j].from == tmpnodes[k] )
                        {
                          nextnode = tmpnodes[k];
                          if( currnode == secondnode && nextnode == firstnode ) { 
                            nextnode = -2;
                            continue;
                          }             
                          
                          tmpnodes[k] = -2;
                          nodes_visited++;
                          if( secondnode == -2 ) secondnode=nextnode;
                          break;
                        }
                    }
                  if( nextnode != -2 ) break;
                }
            }
          
          if( nextnode == firstnode && nodes_visited != vcells[i].numnodes) {
            fprintf( stderr, "Something ugly happened down here!" );
            fprintf( stderr, "%d %d\n", nodes_visited, vcells[i].numnodes );
            vcells[i].on_edge = 1;
            break;
          }
          if( nextnode == firstnode && nodes_visited == vcells[i].numnodes) {
            vcells[i].on_edge = 0;        
            break;
          }
          
          if( nextnode == -2 && nodes_visited == vcells[i].numnodes ) {
            vcells[i].on_edge = 0;
            break;
          }
          
          if( nextnode == -2 ) {
            vcells[i].on_edge = 1;
            break;
          }
          
          currnode = nextnode;
        }
    }        
      
  
  if (verbose)
     fprintf (stderr, "done\n" );

  /* Compute area of voronoi cells */


  if (verbose)
     fprintf (stderr, "Computing area of %d voronoi cells...", n_nuclei);

  for (i=0; i<n_nuclei; i++) 
    {

      memset( upper, 0, maxnodes*sizeof(int) );
      memset( lower, 0, maxnodes*sizeof(int) );
      memset( tmpnodes, 0, maxnodes*sizeof(int) );

      if (vcells[i].on_edge) 
        {
          vcells[i].area = -1;
          continue;
        }

      xmin = fsize;
      xmax = 0;

      for (j=0; j<vcells[i].numnodes; j++)
        {
          x = vnodes[vcells[i].nodes[j]].x;
          y = vnodes[vcells[i].nodes[j]].y;

          if (x<xmin) 
            {
              xmin = x;
              left=vcells[i].nodes[j]; 
            }
          if (x>xmax)
            {
              xmax = x;
              right=vcells[i].nodes[j]; 
            }

          if (fabs(x-xmin)<DBL_EPSILON && y > vnodes[left].y) 
            {
              xmin = x; 
              left=vcells[i].nodes[j]; 
            }
          if (fabs(x-xmax)<DBL_EPSILON && y > vnodes[right].y) 
            {
              xmax = x; 
              right=vcells[i].nodes[j]; 
            }
        }
          
      convexa = (vnodes[right].y - vnodes[left].y)/(vnodes[right].x - vnodes[left].x);
      convexb = vnodes[left].y - convexa*vnodes[left].x;

      n_upper = 0;
      n_lower = 0;
      
      for (j=0; j<vcells[i].numnodes; j++)
        {
          tmpnodes[j] = vcells[i].nodes[j];
          if( tmpnodes[j] == left || tmpnodes[j] == right )
            tmpnodes[j] = -2;
        }
      for (j=vcells[i].numnodes; j<vcells[i].n_nodes; j++)
        {
          tmpnodes[j] = -2;
        }

      nodes_visited = 2;

      while( nodes_visited < vcells[i].numnodes )
        {
          tmp_xmin = DBL_MAX;
          mem = -2;
          for (j=0; j<vcells[i].numnodes; j++)
            {
              if ( tmpnodes[j] != -2) 
                { 
                  x = vnodes[tmpnodes[j]].x;
                  if (x<tmp_xmin) 
                    { 
                      tmp_xmin = x; 
                      mem = j;
                    }
                }
            }
          nextnode=tmpnodes[mem]; 
          tmpnodes[mem] = -2;
          nodes_visited++;
          
          if (vnodes[nextnode].y < (convexa*vnodes[nextnode].x + convexb))
            lower[n_lower++]=nextnode;
          else
            upper[n_upper++]=nextnode;

        }

#ifdef VERB_DEBUG
      fprintf( stderr, "%d  %d ", i, left );
      for( j=0; j<n_lower; j++ )
        fprintf( stderr, "%d ", lower[j] );
      fprintf( stderr, "   ---    " );
      for( j=0; j<n_upper; j++ )
        fprintf( stderr, "%d ", upper[j] );
      fprintf( stderr, "%d\n", right );      
#endif

      lower[n_lower++] = right;
      upper[n_upper++] = right;

      /* Upper integration */
      
      vcells[i].area = 0;
      lastnode = left;
      
      for( j=0; j<n_upper; j++ )
        {
#ifdef VERB_DEBUG
          fprintf( stderr, "up: From: %d to %d\n", lastnode, upper[j] );
#endif
          xlow = vnodes[lastnode].x;
          ylow = vnodes[lastnode].y;
          xhigh = vnodes[upper[j]].x;
          yhigh = vnodes[upper[j]].y;

          quot = (xhigh - xlow);
          if (fabs( quot ) > DBL_EPSILON )
            {
	      a = (yhigh - ylow)/(xhigh - xlow);
	      b = yhigh - a*xhigh;

	      vcells[i].area += 0.5*a*(xhigh*xhigh - xlow*xlow) + b*(xhigh - xlow);
	    }
          lastnode = upper[j];
        }

      /* Lower integration */
      
      lastnode = left;
      
      for( j=0; j<n_lower; j++ )
        {
#ifdef VERB_DEBUG
          fprintf( stderr, "low: From: %d to %d\n", lastnode, lower[j] );
#endif
          xlow = vnodes[lastnode].x;
          ylow = vnodes[lastnode].y;
          xhigh = vnodes[lower[j]].x;
          yhigh = vnodes[lower[j]].y;

          quot = (xhigh - xlow);
          if (fabs( quot ) > DBL_EPSILON )
            {
              a = (yhigh - ylow)/(xhigh - xlow);
              b = yhigh - a*xhigh;

              vcells[i].area -= 0.5*a*(xhigh*xhigh - xlow*xlow) + b*(xhigh - xlow);
            }
          lastnode = lower[j];
        }

      if (vcells[i].area < 0) {
        if (vcells[i].numnodes!=0) {
          fprintf (stderr, "\n" );
          fprintf (stderr, "Ugly cell out of %d!\n", n_nuclei );
          fprintf (stderr, "vcells[%d].nucleus=[%d,(%f,%f)]\n", 
                   i, vcells[i].nucleus, nuclei[vcells[i].nucleus].x, 
                 nuclei[vcells[i].nucleus].y );
          fprintf (stderr, "vcells[%d].numnodes=%d\n", i, vcells[i].numnodes );
          fprintf (stderr, "vcells[%d].area=%d\n", i, vcells[i].area );
          fprintf (stderr, "vcells[%d].on_edge=%d\n", i, vcells[i].on_edge );
          fprintf (stderr, "vcells[%d].nodes=\n", i);
          for (j=0; j<vcells[i].numnodes; j++) 
            fprintf( stderr, "[%f,%f,%d,%d]\n", vnodes[vcells[i].nodes[j]].x,
                     vnodes[vcells[i].nodes[j]].y,
                     vnodes[vcells[i].nodes[j]].numedges,
                     vnodes[vcells[i].nodes[j]].valid);          
          fprintf( stderr, "\n" );
        } else {
          vcells[i].area = -1;
          vcells[i].on_edge = 1;
        }
      }
    }

  if (verbose)
     fprintf (stderr, "done\n" );
}
