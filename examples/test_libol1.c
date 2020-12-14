

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                      LIB OCTREE LOCALISATION V1.64                         */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Description:         Basic localization test on a surface mesh             */
/* Author:              Loic MARECHAL                                         */
/* Creation date:       mar 16 2012                                           */
/* Last modification:   dec 14 2020                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <libmeshb7.h>
#include <libol1.h>


/*----------------------------------------------------------------------------*/
/* Local defines                                                              */
/*----------------------------------------------------------------------------*/

#define BufSiz 1000000
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define POW(a) ((a)*(a))
#define CUB(a) ((a)*(a)*(a))
#define INC 100


/*----------------------------------------------------------------------------*/
/* Load a mesh, build an octree & perform some localisations                  */
/*----------------------------------------------------------------------------*/

int main()
{
   int i, j, k, cpt=0, NmbVer, NmbEdg, NmbTri, NmbQad, NmbTet, NmbItm, ver, dim;
   int (*EdgTab)[2], (*TriTab)[3], (*QadTab)[4], (*TetTab)[4], buf[ BufSiz ];
   int ref, idx;
   int64_t MshIdx, OctIdx;
   double crd1[3] = {0,0,-10 };
   double crd2[3] = {0.002928, 0.079575, 0.006978};
   double crd3[3] = {-.0054, .1488, -.0067};
   double dis, (*VerTab)[3], MinCrd[3], MaxCrd[3], IncCrd[3], AvgDis=0;
   time_t t, t2, MinTim = INT_MAX, MaxTim = 0;


   /*---------------------------------------*/
   /* Open, allocate and read the mesh file */
   /*---------------------------------------*/

   t = clock();
   puts("\nRead mesh");

   if(!(MshIdx = GmfOpenMesh("../sample_meshes/test.meshb", GmfRead, &ver, &dim)))
   //if(!(MshIdx = GmfOpenMesh("../sample_meshes/hexes.meshb", GmfRead, &ver, &dim)))
   {
     puts("Cannot open test.meshb.");
     return(1);
   }

   NmbVer = GmfStatKwd(MshIdx, GmfVertices);
   NmbEdg = GmfStatKwd(MshIdx, GmfEdges);
   NmbTri = GmfStatKwd(MshIdx, GmfTriangles);
   NmbQad = GmfStatKwd(MshIdx, GmfQuadrilaterals);
   NmbTet = GmfStatKwd(MshIdx, GmfTetrahedra);
   printf("vertices = %d, edges = %d, triangles = %d, quads = %d, tets = %d, dimension = %d, version = %d\n",
            NmbVer, NmbEdg, NmbTri, NmbQad, NmbTet, dim, ver);

   if( !NmbVer || (dim != 3) )
   {
      puts("Incompatible mesh.");
      return(1);
   }

   VerTab = malloc((NmbVer+1) * 3 * sizeof(double));
   GmfGetBlock(MshIdx, GmfVertices, 1, NmbVer, 0, NULL, NULL,
               GmfDoubleVec, 3, VerTab[1], VerTab[ NmbVer ],
               GmfInt, &ref, &ref );

   if(NmbEdg)
   {
      EdgTab = malloc((NmbEdg+1) * 2 * sizeof(int));
      GmfGetBlock(MshIdx, GmfEdges, 1, NmbEdg, 0, NULL, NULL,
                  GmfIntVec, 2, EdgTab[1], EdgTab[ NmbEdg ],
                  GmfInt, &ref, &ref);
   }

   if(NmbTri)
   {
      TriTab = malloc((NmbTri+1) * 3 * sizeof(int));
      GmfGetBlock(MshIdx, GmfTriangles, 1, NmbTri, 0, NULL, NULL,
                  GmfIntVec, 3, TriTab[1], TriTab[ NmbTri ],
                  GmfInt, &ref, &ref );
   }

   if(NmbQad)
   {
      QadTab = malloc((NmbQad+1) * 4 * sizeof(int));
      GmfGetBlock(MshIdx, GmfQuadrilaterals, 1, NmbQad, 0, NULL, NULL,
                  GmfIntVec, 4, QadTab[1], QadTab[ NmbQad ],
                  GmfInt, &ref, &ref );
   }

   if(NmbTet)
   {
      TetTab = malloc((NmbTet+1) * 4 * sizeof(int));
      GmfGetBlock(MshIdx, GmfTetrahedra, 1, NmbTet, 0, NULL, NULL,
                  GmfIntVec, 4, TetTab[1], TetTab[ NmbTet ],
                  GmfInt, &ref, &ref );
   }

   GmfCloseMesh(MshIdx);
   printf(" %g s\n", (double)(clock() - t) / CLOCKS_PER_SEC);


   /*---------------------------------------------------------*/
   /* Build an octree from this mesh and perform some queries */
   /*---------------------------------------------------------*/

   puts("\nBuild the octree : ");
   t = clock();
   OctIdx = LolNewOctree(  NmbVer, VerTab[1], VerTab[2],
                           NmbEdg, EdgTab[1], EdgTab[2],
                           NmbTri, TriTab[1], TriTab[2],
                           NmbQad, QadTab[1], QadTab[2],
                           NmbTet, TetTab[1], TetTab[2],
                           0, NULL, NULL,
                           0, NULL, NULL,
                           0, NULL, NULL ,
                           1, 1);

   printf(" %g s\n", (double)(clock() - t) / CLOCKS_PER_SEC);

   for(i=0;i<3;i++)
   {
      MinCrd[i] = crd2[i] - .0001;
      MaxCrd[i] = crd2[i] + .0001;
   }

   puts("\nSearch for vertices in a bounding box :");
   NmbItm = LolGetBoundingBox(OctIdx, LolTypVer, BufSiz, buf, MinCrd, MaxCrd, 0);

   for(i=0;i<NmbItm;i++)
      printf(" vertex : %d\n", buf[i]);

   if(NmbEdg)
   {
      puts("\nSearch for the closest edge from a given point :");
      idx = LolGetNearest(OctIdx, LolTypEdg, crd1, &dis, 0., NULL, NULL, 0);
      printf(" closest edge = %d (%d - %d), distance = %g\n",
               idx, EdgTab[ idx ][0], EdgTab[ idx ][1], dis);
      LolProjectVertex(OctIdx, crd1, LolTypEdg, idx, MinCrd, 0);
      printf(" projection on the closest edge: %g %g %g\n",
               MinCrd[0], MinCrd[1], MinCrd[2] );
      printf(" distance from image = %g\n",\
               sqrt(POW(MinCrd[0] - crd1[0])
                  + POW(MinCrd[1] - crd1[1])
                  + POW(MinCrd[2] - crd1[2]) ));
   }

   if(NmbTri)
   {
      puts("\nSearch for triangles in a bounding box :");

      NmbItm = LolGetBoundingBox(OctIdx, LolTypTri, BufSiz, buf, MinCrd, MaxCrd, 0);

      for(i=0;i<NmbItm;i++)
         printf(" triangle : %d\n", buf[i]);

      puts("\nSearch for the closest triangle from a given point :");
      idx = LolGetNearest(OctIdx, LolTypTri, crd1, &dis, 0., NULL, NULL, 0);
      printf(" closest triangle = %d, distance = %g\n", idx, dis);
      LolProjectVertex(OctIdx, crd1, LolTypTri, idx, MinCrd, 0);
      printf(" projection on the closest triangle: %g %g %g\n",
               MinCrd[0], MinCrd[1], MinCrd[2] );
      printf(" distance from image = %g\n",\
               sqrt(POW(MinCrd[0] - crd1[0])
                  + POW(MinCrd[1] - crd1[1])
                  + POW(MinCrd[2] - crd1[2]) ));

      for(i=1;i<=NmbVer;i++)
         for(j=0;j<3;j++)
            if(VerTab[i][j] < MinCrd[j])
               MinCrd[j] = VerTab[i][j];
         else if(VerTab[i][j] > MaxCrd[j])
            MaxCrd[j] = VerTab[i][j];

      for(i=0;i<3;i++)
         IncCrd[i] = (MaxCrd[i] - MinCrd[i]) / INC;

      crd1[0] = MinCrd[0];

      t = clock();

      for(i=0;i<INC;i++)
      {
         crd1[1] = MinCrd[1];

         for(j=0;j<INC;j++)
         {
            crd1[2] = MinCrd[2];

            for(k=0;k<INC;k++)
            {
               t2 = clock();
               idx = LolGetNearest(OctIdx, LolTypTri, crd1, &dis, 0, NULL, NULL, 0);
               t2 = clock() - t2;
               MinTim = MIN(MinTim, t2);
               MaxTim = MAX(MaxTim, t2);
               AvgDis += sqrt(dis);
               crd1[2] += IncCrd[2];
            }

            crd1[1] += IncCrd[1];
         }

         crd1[0] += IncCrd[0];
      }

      printf("nb samples = %d, mean dist = %g, total time = %g s, min time = %g s, max time = %g s\n",
            CUB(INC), AvgDis / CUB(INC), (double)(clock() - t) / CLOCKS_PER_SEC,
            (double)MinTim / CLOCKS_PER_SEC, (double)MaxTim / CLOCKS_PER_SEC);
   }

   if(NmbQad)
   {
      puts("\nSearch for quads in a bounding box :");

      NmbItm = LolGetBoundingBox(OctIdx, LolTypQad, BufSiz, buf, MinCrd, MaxCrd, 0);

      for(i=0;i<NmbItm;i++)
         printf(" quad : %d\n", buf[i]);

      puts("\nSearch for the closest quad from a given point :");
      idx = LolGetNearest(OctIdx, LolTypQad, crd1, &dis, 0., NULL, NULL, 0);
      printf(" closest quad = %d, distance = %g\n", idx, dis);
      LolProjectVertex(OctIdx, crd1, LolTypQad, idx, MinCrd, 0);
      printf(" projection on the closest quad: %g %g %g\n",
               MinCrd[0], MinCrd[1], MinCrd[2] );
      printf(" distance from image = %g\n",\
               sqrt(POW(MinCrd[0] - crd1[0])
                  + POW(MinCrd[1] - crd1[1])
                  + POW(MinCrd[2] - crd1[2]) ));
   }

   if(NmbTet)
   {
      crd1[0] = 0.5;
      crd1[1] = 0.;
      crd1[2] = 0.;
      idx = LolGetNearest(OctIdx, LolTypTet, crd1, &dis, 0., NULL, NULL, 0);
      printf("tetra %d, dis = %g\n", idx, dis);

      MinCrd[0] = -2;
      MaxCrd[0] = 2;
      MinCrd[1] = -2;
      MaxCrd[1] = 2;
      MinCrd[2] = -.6;
      MaxCrd[2] = -.5;

      puts("\nSearch for vertices in a bounding box :");
      NmbItm = LolGetBoundingBox(OctIdx, LolTypTet, BufSiz, buf, MinCrd, MaxCrd, 0);
      printf("NmbTet = %d\n", NmbItm);

      if(!(MshIdx = GmfOpenMesh("../sample_meshes/cut.meshb", GmfWrite, 2, 3)))
      {
        puts("Cannot open test.meshb.");
        return(1);
      }

      GmfSetKwd(MshIdx, GmfVertices, NmbVer);
      GmfSetBlock(MshIdx, GmfVertices, 1, NmbVer, 0, NULL, NULL,
                  GmfDoubleVec, 3, VerTab[1], VerTab[ NmbVer ],
                  GmfInt, &ref, &ref );

      GmfSetKwd(MshIdx, GmfTetrahedra, NmbItm);
      for(i=0;i<NmbItm;i++)
         GmfSetLin(  MshIdx, GmfTetrahedra,
                     TetTab[ buf[i] ][0], TetTab[ buf[i] ][1],
                     TetTab[ buf[i] ][2], TetTab[ buf[i] ][3], 2);

      GmfCloseMesh(MshIdx);
   }


   /*------------------*/
   /* Cleanup memories */ 
   /*------------------*/

   printf("\nFree octree %lld\n", OctIdx);
   printf(" memory used = %zd bytes\n", LolFreeOctree(OctIdx));
   free(VerTab);

   if(NmbEdg)
      free(EdgTab);

   if(NmbTri)
      free(TriTab);

   if(NmbTet)
      free(TetTab);

   return(0);
}
