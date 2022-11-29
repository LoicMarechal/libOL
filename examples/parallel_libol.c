

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                      LIB OCTREE LOCALISATION V1.81                         */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Description:         Parallel localization on a surface mesh               */
/* Author:              Loic MARECHAL                                         */
/* Creation date:       oct 02 2020                                           */
/* Last modification:   mar 11 2022                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <libmeshb7.h>
#include <lplib3.h>
#include <libol1.h>


/*----------------------------------------------------------------------------*/
/* Local defines                                                              */
/*----------------------------------------------------------------------------*/

#define INX 100
#define INY 100
#define INZ 100
#define MUL 1
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define POW(a) ((a)*(a))


/*----------------------------------------------------------------------------*/
/* Global variables                                                           */
/*----------------------------------------------------------------------------*/

int64_t OctIdx;
double MinCrd[3], IncCrd[3], AvgDis[ MaxPth ];


/*----------------------------------------------------------------------------*/
/* Prototypes of local procedures                                             */
/*----------------------------------------------------------------------------*/

void ParallelGridSearch(int, int, int, void *);


/*----------------------------------------------------------------------------*/
/* Load a mesh, build an octree & perform some localisations                  */
/*----------------------------------------------------------------------------*/

int main()
{
   int i, j, NmbVer, NmbTri, ver, dim, (*TriTab)[3], ref, NmbThr, GrdTyp;
   int64_t MshIdx, ParIdx;
   double t, mid, siz, avg = 0., (*VerTab)[3], MaxCrd[3];


   /*---------------------------------------*/
   /* Open, allocate and read the mesh file */
   /*---------------------------------------*/

   t = GetWallClock();
   printf("\nMesh reading: ");

   if(!(MshIdx = GmfOpenMesh("../sample_meshes/test.meshb", GmfRead, &ver, &dim)))
   {
     puts("Cannot open test.meshb.");
     return(1);
   }

   NmbVer = (int)GmfStatKwd(MshIdx, GmfVertices);
   NmbTri = (int)GmfStatKwd(MshIdx, GmfTriangles);
 
   if( !NmbVer || !NmbTri || (dim != 3) )
   {
      puts("Incompatible mesh.");
      return(1);
   }

   VerTab = malloc((NmbVer+1) * 3 * sizeof(double));
   assert(VerTab);
   GmfGetBlock(MshIdx, GmfVertices, 1, NmbVer, 0, NULL, NULL,
               GmfDoubleVec, 3, VerTab[1], VerTab[ NmbVer ],
               GmfInt, &ref, &ref );

   TriTab = malloc((NmbTri+1) * 3 * sizeof(int));
   assert(TriTab);
   GmfGetBlock(MshIdx, GmfTriangles, 1, NmbTri, 0, NULL, NULL,
               GmfIntVec, 3, TriTab[1], TriTab[ NmbTri ],
               GmfInt, &ref, &ref );

   GmfCloseMesh(MshIdx);
   printf(" %g s\n", GetWallClock() - t);
   printf("vertices = %d, triangles = %d, dimension = %d, version = %d\n",
            NmbVer, NmbTri, dim, ver);


   /*---------------------------------------------------------*/
   /* Build an octree from this mesh and perform some queries */
   /*---------------------------------------------------------*/

   NmbThr = GetNumberOfCores();

   if(!(ParIdx = InitParallel(NmbThr)))
   {
      puts("Error initializing the LPlib.");
      exit(1);
   }

   SetExtendedAttributes(ParIdx, SetInterleavingFactor, 16);

   if(!(GrdTyp = NewType(ParIdx, INX)))
   {
      puts("Error while creating the grid data type.");
      exit(1);
   }

   printf("\nOctree building: ");
   t = clock();
   OctIdx = LolNewOctree(  NmbVer, VerTab[1], VerTab[2],
                           0, NULL, NULL,
                           NmbTri, TriTab[1], TriTab[2],
                           0, NULL, NULL,
                           0, NULL, NULL,
                           0, NULL, NULL,
                           0, NULL, NULL,
                           0, NULL, NULL ,
                           1, NmbThr);

   if(!OctIdx)
   {
      puts("The octree building failled.");
      exit(1);
   }

   printf(" %g s\n", (double)(clock() - t) / CLOCKS_PER_SEC);

   // Compute the bounding box and the testing grid steps
   for(i=0;i<3;i++)
      MaxCrd[i] = MinCrd[i] = VerTab[1][i];

   for(i=1;i<=NmbVer;i++)
      for(j=0;j<3;j++)
         if(VerTab[i][j] < MinCrd[j])
            MinCrd[j] = VerTab[i][j];
      else if(VerTab[i][j] > MaxCrd[j])
         MaxCrd[j] = VerTab[i][j];

   // Scale the grid
   for(i=0;i<3;i++)
   {
      mid = (MaxCrd[i] + MinCrd[i]) / 2.;
      siz = (MaxCrd[i] - MinCrd[i]) / 2.;
      MinCrd[i] = mid - siz * MUL;
      MaxCrd[i] = mid + siz * MUL;
   }

   IncCrd[0] = (MaxCrd[0] - MinCrd[0]) / INX;
   IncCrd[1] = (MaxCrd[1] - MinCrd[1]) / INY;
   IncCrd[2] = (MaxCrd[2] - MinCrd[2]) / INZ;

   for(i=0;i<NmbThr;i++)
      AvgDis[i] = 0.;

   // Call parallel GetNearest
   t = GetWallClock();
   LaunchParallel(ParIdx, GrdTyp, 0, ParallelGridSearch, NULL);
   t = GetWallClock() - t;

   for(i=0;i<NmbThr;i++)
      avg += AvgDis[i];

   printf("Search time: %g s on %d threads\n", t, NmbThr);
   printf("nb samples = %d, mean dist = %g\n", INX*INY*INZ, avg / (INX*INY*INZ));


   /*------------------*/
   /* Cleanup memories */ 
   /*------------------*/

   printf("\nMemory used = %zd bytes\n", LolFreeOctree(OctIdx));
   free(VerTab);
   free(TriTab);

   return(0);
}


/*----------------------------------------------------------------------------*/
/* Parallel call to GetNearest from the grid points to the mesh triangles     */
/*----------------------------------------------------------------------------*/

void ParallelGridSearch(int BegIdx, int EndIdx, int PthIdx, void *nil)
{
   int i, j, k;
   double dis, TstCrd[3], avg = 0.;

   for(i=BegIdx; i<=EndIdx; i++)
   {
      TstCrd[0] = MinCrd[0] + (i-1) * IncCrd[0];

      for(j=0;j<INY;j++)
      {
         TstCrd[1] = MinCrd[1] + j * IncCrd[1];

         for(k=0;k<INZ;k++)
         {
            TstCrd[2] = MinCrd[2] + k * IncCrd[2];
            LolGetNearest(OctIdx, LolTypTri, TstCrd, &dis, 0, NULL, NULL, PthIdx);
            avg += sqrt(dis);
         }
      }
   }

   AvgDis[ PthIdx ] += avg;
}
