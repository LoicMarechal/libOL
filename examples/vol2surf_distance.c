

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                      LIB OCTREE LOCALISATION V1.84                         */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Description:         For each volume vertex, get the distance from         */
/*                      the surface mesh and store the results in a sizemap   */
/* Author:              Loic MARECHAL                                         */
/* Creation date:       mar 31 2025                                           */
/* Last modification:   mar 12 2026                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <libmeshb8.h>
#include <lplib4.h>
#include <libol1.h>


/*----------------------------------------------------------------------------*/
/* Mesh structure                                                             */
/*----------------------------------------------------------------------------*/

typedef struct
{
   int      ver, dim, NmbVer, NmbTri, (*TriTab)[3];
   double   (*VerTab)[3], *SizMap;
   char     *FilNam;
   int64_t  FilIdx;
}MshSct;

typedef struct
{
   int      NmbCpu;
   int64_t  LplIdx, OctIdx, VerTyp;
   MshSct   *VolMsh, *SrfMsh, *SolMsh;
}ParSct;


/*----------------------------------------------------------------------------*/
/* Prototypes of local procedures                                             */
/*----------------------------------------------------------------------------*/

void     ParallelSearch(int, int, int, ParSct *);
MshSct  *ReadMesh(char *);
int      WriteSol(MshSct *, char *);
void     CpyVec(double *, double *);
double   GetTetVol(double [4][3]);


/*----------------------------------------------------------------------------*/
/* Load a mesh, build an octree & perform some localisations                  */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   int i, j;
   double t;
   char  *VolNam, *SrfNam;
   ParSct par;
   MshSct *VolMsh, *SrfMsh, SolMsh;


   /*---------------------------------------*/
   /* Parse the command line and print help */
   /*---------------------------------------*/

   if(ArgCnt != 4)
   {
      puts("\nvol2surf_distance  VolumeMesh  SurfaceMesh  OutputSizeMap\n\n");
      exit(0);
   }
   else
   {
      VolNam = *++ArgVec;
      SrfNam = *++ArgVec;
      SolMsh.FilNam = *++ArgVec;
   }


   /*----------------------------------*/
   /* Allocate and read the mesh files */
   /*----------------------------------*/

   t = GetWallClock();
   puts("");
   printf("Reading the meshes        : ");
   VolMsh = ReadMesh(VolNam);
   SrfMsh = ReadMesh(SrfNam);
   printf(" %g s\n", GetWallClock() - t);

   SolMsh.NmbVer = VolMsh->NmbVer;
   SolMsh.SizMap = malloc( (VolMsh->NmbVer + 1) * sizeof(double) );
   par.VolMsh = VolMsh;
   par.SrfMsh = SrfMsh;
   par.SolMsh = &SolMsh;


   /*-------------------------------------------*/
   /* Initialize the parallelism with the LPlib */
   /*-------------------------------------------*/

   par.NmbCpu = 4 * GetNumberOfCores();

   if(!(par.LplIdx = InitParallel(par.NmbCpu)))
   {
      puts("Error initializing the LPlib.");
      exit(1);
   }

   if(!(par.VerTyp = NewType(par.LplIdx, VolMsh->NmbVer)))
   {
      puts("Error while creating the vertex data type.");
      exit(1);
   }


   /*---------------------------------------------------------*/
   /* Build an octree from this mesh and perform some queries */
   /*---------------------------------------------------------*/

   printf("Building the octree       : ");
   t = GetWallClock();
   par.OctIdx = LolNewOctree( SrfMsh->NmbVer, SrfMsh->VerTab[1], SrfMsh->VerTab[2],
                              0, NULL, NULL,
                              SrfMsh->NmbTri, SrfMsh->TriTab[1], SrfMsh->TriTab[2],
                              0, NULL, NULL,
                              0, NULL, NULL,
                              0, NULL, NULL,
                              0, NULL, NULL,
                              0, NULL, NULL ,
                              1, par.NmbCpu);

   if(!par.OctIdx)
   {
      puts("Failled to build the octree.");
      exit(1);
   }

   printf(" %g s\n", GetWallClock() - t);

   // Call parallel GetNearest
   printf("Computing the distances   : ");
   t = GetWallClock();
   LaunchParallel(par.LplIdx, par.VerTyp, 0, ParallelSearch, &par);
   printf(" %g s\n", GetWallClock() - t);


   /*-------------------------------------*/
   /* Write the distances in a .solb file */
   /*-------------------------------------*/

   printf("Writing the solution      : ");
   t = GetWallClock();
   WriteSol(&SolMsh, SolMsh.FilNam);
   printf(" %g s\n", GetWallClock() - t);


   /*------------------*/
   /* Cleanup memories */ 
   /*------------------*/

   printf("Memory used by the octree : %zd bytes\n\n", LolFreeOctree(par.OctIdx));
   StopParallel(par.LplIdx);
   

   return(0);
}


/*----------------------------------------------------------------------------*/
/* Read a volume or a surface mesh and allocate all needed memories           */
/*----------------------------------------------------------------------------*/

MshSct *ReadMesh(char *FilNam)
{
   int ref;
   MshSct *msh;

   msh = calloc(1, sizeof(MshSct));
   assert(msh);
   msh->FilNam = FilNam;

   if(!(msh->FilIdx = GmfOpenMesh(FilNam, GmfRead, &msh->ver, &msh->dim)))
   {
      printf("Cannot open mesh file %s.", FilNam);
      exit(1);
   }

   msh->NmbVer = GmfStatKwd(msh->FilIdx, GmfVertices);
   msh->NmbTri = GmfStatKwd(msh->FilIdx, GmfTriangles);
 
   if( !msh->NmbVer || (msh->dim != 3) )
   {
      puts("Incompatible mesh.");
      exit(2);
   }

   msh->VerTab = malloc( (msh->NmbVer + 1) * 3 * sizeof(double) );
   assert(msh->VerTab);

   GmfGetBlock(msh->FilIdx, GmfVertices, 1, msh->NmbVer, 0, NULL, NULL,
               GmfDoubleVec, 3, msh->VerTab[1], msh->VerTab[ msh->NmbVer ],
               GmfInt, &ref, &ref );

   if(msh->NmbTri)
   {
      msh->TriTab = malloc( (msh->NmbTri + 1) * 3 * sizeof(int) );
      assert(msh->TriTab);

      GmfGetBlock(msh->FilIdx, GmfTriangles, 1, msh->NmbTri, 0, NULL, NULL,
                  GmfIntVec, 3, msh->TriTab[1], msh->TriTab[ msh->NmbTri ],
                  GmfInt, &ref, &ref );
   }

   GmfCloseMesh(msh->FilIdx);

   return(msh);
}

/*----------------------------------------------------------------------------*/
/* Write the size map                                                         */
/*----------------------------------------------------------------------------*/

int WriteSol(MshSct *msh, char *FilNam)
{
   int TypTab[1];

   if(!(msh->FilIdx = GmfOpenMesh(FilNam, GmfWrite, 3, 3)))
      exit(1);

   // Write the distances
   TypTab[0] = GmfSca;
   GmfSetKwd(msh->FilIdx, GmfSolAtVertices, msh->NmbVer, 1, TypTab);
   GmfSetBlock(msh->FilIdx, GmfSolAtVertices, 1, msh->NmbVer, 0, NULL, NULL,
               GmfDouble, &msh->SizMap[1], &msh->SizMap[ msh->NmbVer ]);

   GmfCloseMesh(msh->FilIdx);

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Parallel call to GetNearest from the volume points to the mesh triangles   */
/*----------------------------------------------------------------------------*/

void ParallelSearch(int BegIdx, int EndIdx, int PthIdx, ParSct *par)
{
   int i, j, TriIdx;
   double TetCrd[4][3];

   for(i=BegIdx;i<=EndIdx;i++)
   {
      // Search for the closest surface triangle from each volume vertex
      TriIdx = LolGetNearest( par->OctIdx, LolTypTri, par->VolMsh->VerTab[i],
                              &par->SolMsh->SizMap[i], 0, NULL, NULL, PthIdx );

      if(!TriIdx)
      {
         printf("Could not find the nearest triangle from vertex %d\n", i);
         continue;
      }

      CpyVec(par->VolMsh->VerTab[i], TetCrd[3]);

      for(j=0;j<3;j++)
         CpyVec(par->SrfMsh->VerTab[ par->SrfMsh->TriTab[ TriIdx ][j] ], TetCrd[j]);

      if(GetTetVol(TetCrd) < 0.)
         par->SolMsh->SizMap[i] = -par->SolMsh->SizMap[i];
   }
}


/*----------------------------------------------------------------------------*/
/* Vector functions                                                           */
/*----------------------------------------------------------------------------*/

void CpyVec(double *u, double *v)
{
   v[0] = u[0];
   v[1] = u[1];
   v[2] = u[2];
}

double GetTetVol(double TetCrd[4][3])
{
   double mat[9];

   // Create the linear system, each column is filled with a vector
   mat[0] = TetCrd[1][0] - TetCrd[0][0];
   mat[1] = TetCrd[2][0] - TetCrd[0][0];
   mat[2] = TetCrd[3][0] - TetCrd[0][0];

   mat[3] = TetCrd[1][1] - TetCrd[0][1];
   mat[4] = TetCrd[2][1] - TetCrd[0][1];
   mat[5] = TetCrd[3][1] - TetCrd[0][1];

   mat[6] = TetCrd[1][2] - TetCrd[0][2];
   mat[7] = TetCrd[2][2] - TetCrd[0][2];
   mat[8] = TetCrd[3][2] - TetCrd[0][2];

   // Return the determinant of the matrix
   return(  mat[0] * (mat[4]*mat[8] - mat[5]*mat[7])
         +  mat[1] * (mat[5]*mat[6] - mat[3]*mat[8])
         +  mat[2] * (mat[3]*mat[7] - mat[4]*mat[6]) );
}
