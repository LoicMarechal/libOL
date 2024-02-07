

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                      LIB OCTREE LOCALISATION V1.82                         */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Description:         Searching for the nearest triangle                    */
/*                      belonging to a given reference                        */
/* Author:              Loic MARECHAL                                         */
/* Creation date:       feb 07 2024                                           */
/* Last modification:   feb 07 2024                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <libmeshb7.h>
#include <libol1.h>


/*----------------------------------------------------------------------------*/
/* Mesh structure                                                             */
/*----------------------------------------------------------------------------*/

typedef struct
{
   int NmbVer, NmbTri, NmbQad, TstRef, (*TriTab)[4], (*QadTab)[5];
   double (*VerCrd)[3];
}MshSct;


/*----------------------------------------------------------------------------*/
/* Filter out triangles that do not belong to the requested reference         */
/* The filtering procedure had only two input parameters:                     */
/* one pointing to a global structure that contains all the user's data       */
/* and the 2nd one that is the triangle's index being processed by the libOL  */
/* With these two arguments, it is up to the user's filter to decide whether  */
/* to reject (return 0) or accept (return 1) this triangle                    */
/*----------------------------------------------------------------------------*/

int FltTri(void *UsrDat, int TriIdx)
{
   MshSct *TriMsh = (MshSct *)UsrDat;

   if(TriMsh->TstRef == TriMsh->TriTab[ TriIdx ][3])
      return(1);
   else
      return(0);
}


/*----------------------------------------------------------------------------*/
/* Load a triangulated mesh, build an octree from it, load a quad mesh        */
/* and project its nodes a the triangulated mesh with the selected referrence */
/*----------------------------------------------------------------------------*/

int main()
{
   int i, j, ref, idx, ver, dim;
   int64_t MshIdx, OctIdx;
   double dis, PrjCrd[3];
   MshSct TriMsh, QadMsh;
   char *PrjNam[2] = {"../sample_meshes/projection1.meshb",
                      "../sample_meshes/projection2.meshb"};


   /*-------------------------------------------------------*/
   /* Open, allocate and read the triangulated surface mesh */
   /*-------------------------------------------------------*/

   puts("Reading the triangulated mesh.");

   // Try to open the triangulated mesh file for reading
   if(!(MshIdx = GmfOpenMesh("../sample_meshes/spheres_tri.meshb", GmfRead, &ver, &dim)))
   {
     puts("Cannot open spheres_tri.meshb.");
     return(1);
   }

   // Get the mesh size
   TriMsh.NmbVer = (int)GmfStatKwd(MshIdx, GmfVertices);
   TriMsh.NmbTri = (int)GmfStatKwd(MshIdx, GmfTriangles);

   if( !TriMsh.NmbVer || !TriMsh.NmbTri || (dim != 3) )
   {
      puts("Incompatible mesh.");
      return(1);
   }

   // Allocate and read the vertices coordinates
   TriMsh.VerCrd = malloc( (TriMsh.NmbVer + 1) * 3 * sizeof(double) );
   assert(TriMsh.VerCrd);
   GmfGetBlock(MshIdx, GmfVertices, 1, TriMsh.NmbVer, 0, NULL, NULL,
               GmfDoubleVec, 3, TriMsh.VerCrd[1], TriMsh.VerCrd[ TriMsh.NmbVer ],
               GmfInt, &ref, &ref );

   // Allocate and read the triangles' nodes and reference
   TriMsh.TriTab = malloc( (TriMsh.NmbTri + 1) * 4 * sizeof(int) );
   assert(TriMsh.TriTab);
   GmfGetBlock(MshIdx, GmfTriangles, 1, TriMsh.NmbTri, 0, NULL, NULL,
               GmfIntVec, 4, TriMsh.TriTab[1], TriMsh.TriTab[ TriMsh.NmbTri ]);

   GmfCloseMesh(MshIdx);
   printf("vertices = %d, triangles = %d\n", TriMsh.NmbVer, TriMsh.NmbTri);


   /*---------------------------------------------------------*/
   /* Build an octree from the vertices and triangles         */
   /*---------------------------------------------------------*/

   puts("Building the octree.");
   OctIdx = LolNewOctree(  TriMsh.NmbVer, TriMsh.VerCrd[1], TriMsh.VerCrd[2],
                           0, NULL, NULL,
                           TriMsh.NmbTri, TriMsh.TriTab[1], TriMsh.TriTab[2],
                           0, NULL, NULL,
                           0, NULL, NULL,
                           0, NULL, NULL,
                           0, NULL, NULL,
                           0, NULL, NULL ,
                           1, 1);

   if(!OctIdx)
   {
      puts("The octree building failled.");
      exit(1);
   }


   /*--------------------------------------------------------*/
   /* Open, allocate and read the quadrilateral surface mesh */
   /*--------------------------------------------------------*/

   puts("Reading the quadrilateral mesh");

   // Try to open the quadrilateral mesh file for reading
   if(!(MshIdx = GmfOpenMesh("../sample_meshes/sphere_quad.meshb", GmfRead, &ver, &dim)))
   {
     puts("Cannot open sphere_quad.meshb.");
     return(1);
   }

   // Get the mesh size
   QadMsh.NmbVer = GmfStatKwd(MshIdx, GmfVertices);
   QadMsh.NmbQad = GmfStatKwd(MshIdx, GmfQuadrilaterals);

   if( !QadMsh.NmbVer || !QadMsh.NmbQad || (dim != 3) )
   {
      puts("Incompatible mesh.");
      return(1);
   }

   // Allocate and read the vertices coordinates
   QadMsh.VerCrd = malloc( (QadMsh.NmbVer + 1) * 3 * sizeof(double) );
   assert(QadMsh.VerCrd);
   GmfGetBlock(MshIdx, GmfVertices, 1, QadMsh.NmbVer, 0, NULL, NULL,
               GmfDoubleVec, 3, QadMsh.VerCrd[1], QadMsh.VerCrd[ QadMsh.NmbVer ],
               GmfInt, &ref, &ref );

   // Allocate and read the quads' nodes and reference
   QadMsh.QadTab = malloc( (QadMsh.NmbQad + 1) * 5 * sizeof(int) );
   assert(QadMsh.QadTab);
   GmfGetBlock(MshIdx, GmfQuadrilaterals, 1, QadMsh.NmbQad, 0, NULL, NULL,
               GmfIntVec, 5, QadMsh.QadTab[1], QadMsh.QadTab[ QadMsh.NmbQad ]);

   GmfCloseMesh(MshIdx);
   printf("vertices = %d, quads = %d\n", QadMsh.NmbVer, QadMsh.NmbQad);


   /*-----------------------------------------------------------------*/
   /* Project the quad mesh vertices on the closest ref 1/2 triangles */
   /*-----------------------------------------------------------------*/

   for(ref = 1; ref <= 2; ref++)
   {
      printf("Projecting the quads on the triangulated surface with ref %d :\n", ref);

      TriMsh.TstRef = ref;

      for(i=1;i<=QadMsh.NmbVer;i++)
      {
         // For each quad mesh's vertex, the closest triangle with
         // the requested reference is searched from the octree
         idx = LolGetNearest(OctIdx, LolTypTri, QadMsh.VerCrd[i], &dis, 0., FltTri, &TriMsh, 0);

         if(!idx)
            continue;

         // Compute the projection of vertex on the closest triangle
         LolProjectVertex(OctIdx, QadMsh.VerCrd[i], LolTypTri, idx, PrjCrd, 0);

         for(j=0;j<3;j++)
            QadMsh.VerCrd[i][j] = PrjCrd[j];
      }

      printf("Saving the projected mesh %s\n", PrjNam[ ref - 1 ]);

      // Try to create a quad mesh file
      if(!(MshIdx = GmfOpenMesh(PrjNam[ ref - 1 ], GmfWrite, 2, 3)))
      {
        printf("Cannot create file %s\n", PrjNam[ ref - 1 ]);
        return(1);
      }

      // Save the vertices' coordinates
      GmfSetKwd(MshIdx, GmfVertices, QadMsh.NmbVer);
      GmfSetBlock(MshIdx, GmfVertices, 1, QadMsh.NmbVer, 0, NULL, NULL,
                  GmfDoubleVec, 3, QadMsh.VerCrd[1], QadMsh.VerCrd[ QadMsh.NmbVer ],
                  GmfInt, &ref, &ref );

      // Save the quads' nodes and reference
      GmfSetKwd(MshIdx, GmfQuadrilaterals, QadMsh.NmbQad);
      GmfSetBlock(MshIdx, GmfQuadrilaterals, 1, QadMsh.NmbQad, 0, NULL, NULL,
                  GmfIntVec, 4, QadMsh.QadTab[1], QadMsh.QadTab[ QadMsh.NmbQad ],
                  GmfInt, &ref, &ref );

      GmfCloseMesh(MshIdx);
   }


   /*------------------*/
   /* Cleanup memories */ 
   /*------------------*/

   printf("Memory used by the octree = %zd bytes\n", LolFreeOctree(OctIdx));
   free(TriMsh.VerCrd);
   free(QadMsh.VerCrd);
   free(TriMsh.TriTab);
   free(QadMsh.QadTab);

   return(0);
}
