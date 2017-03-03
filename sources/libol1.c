

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                      LIB OCTREE LOCALISATION V1.40                         */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    Description:         Octree for mesh localization                       */
/*    Author:              Loic MARECHAL                                      */
/*    Creation date:       mar 16 2012                                        */
/*    Last modification:   mar 03 2017                                        */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* ANSI C headers                                                             */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <fcntl.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "libol1.h"


/*----------------------------------------------------------------------------*/
/* Defines and macros                                                         */
/*----------------------------------------------------------------------------*/

#define MaxItmOct 20
#define MaxOctLvl 10
#define MinGrdLvl 3
#define ItmPerBuc 100
#define MemBlkSiz 100000
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define POW(a) ((a)*(a))
#define CUB(a) ((a)*(a)*(a))


/*----------------------------------------------------------------------------*/
/* Local structures                                                           */
/*----------------------------------------------------------------------------*/

typedef struct
{
   LolInt idx;
   double crd[3];
}VerSct;

typedef struct
{
   VerSct *ver[2];
   double tng[3], siz;
   LolInt idx;
}EdgSct;

typedef struct
{
   EdgSct edg[3];
   VerSct *ver[3];
   double nrm[3];
   LolInt idx;
}TriSct;

typedef struct
{
   VerSct *ver[4];
   EdgSct edg[4];
   TriSct tri[2];
   LolInt idx;
}QadSct;

typedef struct
{
   VerSct *ver[4];
   EdgSct edg[6];
   TriSct tri[4];
   LolInt idx;
}TetSct;

typedef struct
{
   VerSct *ver[5];
   EdgSct edg[8];
   TriSct tri[5][2];
   LolInt idx;
}PyrSct;

typedef struct
{
   VerSct *ver[6];
   EdgSct edg[9];
   TriSct tri[5][2];
   LolInt idx;
}PriSct;

typedef struct
{
   VerSct *ver[8];
   EdgSct edg[12];
   TriSct tri[6][2];
   LolInt idx;
}HexSct;

typedef struct LnkSctPtr
{
   LolInt typ, idx;
   struct LnkSctPtr *nex;
}LnkSct;

typedef struct MemSctPtr
{
   LolInt siz;
   void *adr;
   struct MemSctPtr *nex;
}MemSct;

typedef struct
{
   VerSct ver[3];
   EdgSct edg;
   TriSct tri;
   QadSct qad;
   TetSct tet;
   PyrSct pyr;
   PriSct pri;
   HexSct hex;
   LolInt NmbVer, NmbEdg, NmbTri, NmbQad, NmbTet;
   LolInt NmbPyr, NmbPri, NmbHex, NmbItm;
   char *FlgTab, *VerCrd, VerSiz, *EdgIdx, EdgSiz, *TriIdx, TriSiz;
   char *QadIdx, QadSiz, TetIdx, TetSiz, *PyrIdx, PyrSiz;
   char *PriIdx, PriSiz, *HexIdx, HexSiz;
}MshSct;

typedef struct OctSctPtr
{
   union
   {
      LnkSct *lnk;
      struct OctSctPtr *son;
   };
   unsigned char cpt[ LolNmbTyp ], lvl, sub;
}OctSct;

typedef struct
{
   OctSct *oct;
   unsigned int tag;
   char pos[3];
}BucSct;

typedef struct
{
   VerSct ver[8];
   HexSct hex;
   unsigned int tag;
   LolInt MaxItm, MaxLvl, NmbFreOct, NmbOct, GrdLvl, NmbBuc;
   LolInt MemUse;
   double eps, MaxSiz, MinSiz, BucSiz, bnd[2][3];
   OctSct oct, *CurOctBlk;
   BucSct *grd, **stk;
   MshSct *msh;
   LnkSct *NexFreLnk;
   MemSct *NexMem;
}OctMshSct;


/*----------------------------------------------------------------------------*/
/* Prototypes of octree procedures                                            */
/*----------------------------------------------------------------------------*/

static void SetMshBox(OctMshSct *, MshSct *);
static void AddVer(MshSct *, OctMshSct *, OctSct *, double [3], double [3]);
static void AddEdg(MshSct *, OctMshSct *, OctSct *, double [3], double [3]);
static void AddTri(MshSct *, OctMshSct *, OctSct *, double [3], double [3]);
static void SubOct(MshSct *, OctMshSct *, OctSct *, double [3], double [3]);
static void LnkItm(OctMshSct *, OctSct *, LolInt, LolInt);
static OctSct *GetCrd(OctSct *, int, double [3], double [3], double [3]);
static void GetBox(  OctSct *, LolInt, LolInt *, LolInt, LolInt *, char *, \
                     double [2][3], double, double [3], double [3] );
static LolInt BoxIntBox(double [2][3], double [2][3], double);
static void SetEdg(MshSct *, LolInt);
static void SetTri(MshSct *, LolInt);
static void SetSonCrd(LolInt, double [3], double [3], double [3], double [3]);
static void GetOctLnk(  MshSct *, LolInt, double [3], LolInt *, double *, \
                        OctSct *, double [3], double [3] );
static void GetBucBox(OctMshSct *, BucSct *, double [3], double [3]);
static BucSct *GetBucNgb(OctMshSct *, BucSct *, LolInt);
static double DisVerOct(double [3], double [3], double [3]);
static LolInt VerInsOct(double [3], double [3], double [3]);
static double *GetPtrCrd(MshSct *, LolInt);
static LolInt *GetPtrEdg(MshSct *, LolInt);
static LolInt *GetPtrTri(MshSct *, LolInt);


/*----------------------------------------------------------------------------*/
/* Prototypes of meshing procedures                                           */
/*----------------------------------------------------------------------------*/

static LolInt EdgIntEdg(EdgSct *, EdgSct *, VerSct *, double);
static double DisVerTri(MshSct *, double [3], LolInt);
static double GetTriSrf(TriSct *);
static double DisVerEdg(double [3], EdgSct *);
static void GetTriVec(TriSct *, double [3]);
static void SetTriNrm(TriSct *);
static void SetTmpHex(HexSct *, double [3], double [3]);
static LolInt  VerInsHex(VerSct *, HexSct *);
static LolInt EdgIntHex(EdgSct *, HexSct *, double);
static LolInt TriIntHex(TriSct *, HexSct *, double);
static LolInt EdgIntTri(TriSct *, EdgSct *, VerSct *, double);
static LolInt VerInsTri(TriSct *, VerSct *, double);
static LolInt VerInsEdg(EdgSct *, VerSct *, double);
static void SetEdgTng(EdgSct *);


/*----------------------------------------------------------------------------*/
/* Prototypes of geometric procedures                                         */
/*----------------------------------------------------------------------------*/

static void PrjVerLin(double [3], double [3], double [3], double [3]);
static double PrjVerPla(double [3], double [3], double [3], double [3]);
static void LinCmbVec3(double, double [3], double, double [3], double [3]);
static void CpyVec(double [3], double [3]);
static void AddVec2(double [3], double [3]);
static void SubVec3(double [3], double [3], double [3]);
static void AddScaVec1(double, double [3]);
static void AddScaVec2(double, double [3], double [3]);
static void MulVec1(double, double [3]);
static void MulVec2(double, double [3], double [3]);
static void NrmVec(double [3]);
static void CrsPrd(double [3], double [3], double [3]);
static double DotPrd(double [3], double [3]);
static double dis(double [3], double [3]);
static double DisPow(double [3], double [3]);
static double DisVerPla(double [3], double [3], double [3]);
static double GetNrmVec(double [3]);


/*----------------------------------------------------------------------------*/
/* Prototypes of memory handling procedures                                   */
/*----------------------------------------------------------------------------*/

static void *NewMem(OctMshSct *, LolInt);
static void FreAllMem(OctMshSct *);


/*----------------------------------------------------------------------------*/
/* Allocate and build a new octree from user's data                           */
/*----------------------------------------------------------------------------*/

int64_t LolNewOctree(   LolInt NmbVer, double *PtrCrd1, double *PtrCrd2, \
                        LolInt NmbEdg, LolInt *PtrEdg1, LolInt *PtrEdg2, \
                        LolInt NmbTri, LolInt *PtrTri1, LolInt *PtrTri2, \
                        LolInt NmbQad, LolInt *PtrQad1, LolInt *PtrQad2, \
                        LolInt NmbTet, LolInt *PtrTet1, LolInt *PtrTet2, \
                        LolInt NmbPyr, LolInt *PtrPyr1, LolInt *PtrPyr2, \
                        LolInt NmbPri, LolInt *PtrPri1, LolInt *PtrPri2, \
                        LolInt NmbHex, LolInt *PtrHex1, LolInt *PtrHex2 )
{
   LolInt i, j, k;
   const int tvpe[12][2] = {  {3,2}, {0,1}, {4,5}, {7,6}, {3,7}, {2,6}, \
                              {1,5}, {0,4}, {3,0}, {7,4}, {6,5}, {2,1} };
   const int tvpf[6][4] = {   {3,0,4,7}, {5,1,2,6}, {3,2,1,0}, \
                              {5,6,7,4},{3,7,6,2}, {5,4,0,1} };
   double crd[3];
   BucSct *buc;
   OctMshSct *OctMsh = NULL;
   MshSct *msh;

   // Setup a single octant octree
   assert((OctMsh = calloc(1, sizeof(OctMshSct))));

   // Setup the mesh structure
   msh = NewMem(OctMsh, sizeof(MshSct));

   msh->NmbVer = NmbVer;
   msh->NmbItm = NmbVer;
   msh->VerCrd = (char *)PtrCrd1;
   msh->VerSiz = (char *)PtrCrd2 - (char *)PtrCrd1;

   msh->NmbEdg = NmbEdg;
   msh->NmbItm = MAX(msh->NmbItm, NmbEdg);
   msh->EdgIdx = (char *)PtrEdg1;
   msh->EdgSiz = (char *)PtrEdg2 - (char *)PtrEdg1;

   msh->NmbTri = NmbTri;
   msh->NmbItm = MAX(msh->NmbItm, NmbTri);
   msh->TriIdx = (char *)PtrTri1;
   msh->TriSiz = (char *)PtrTri2 - (char *)PtrTri1;

   msh->FlgTab = NewMem(OctMsh, (msh->NmbItm + 1) * sizeof(char) );
   memset(msh->FlgTab, 0, (msh->NmbItm + 1) * sizeof(char));

   for(i=0;i<3;i++)
   {
      msh->tri.ver[i] = &msh->ver[i];
      msh->tri.edg[i].ver[0] = &msh->ver[ (i+1)%3 ];
      msh->tri.edg[i].ver[1] = &msh->ver[ (i+2)%3 ];
   }

   // Setup the octree structure
   SetMshBox(OctMsh, msh);
   OctMsh->msh = msh;
   OctMsh->MaxItm = MaxItmOct;

   // Set the grid size depending on the number of entities in the mesh
   OctMsh->GrdLvl = MAX(MinGrdLvl, log((msh->NmbTri + msh->NmbEdg + msh->NmbVer) \
                                       / ItmPerBuc) / (3 * log(2)));

   OctMsh->NmbBuc = 1<<OctMsh->GrdLvl;

   OctMsh->grd = NewMem(OctMsh, CUB(OctMsh->NmbBuc) * sizeof(BucSct));
   OctMsh->stk = NewMem(OctMsh, CUB(OctMsh->NmbBuc) * sizeof(void *));

   // Setup the temporary edge for local geometric calculations
   for(i=0;i<2;i++)
      msh->edg.ver[i] = &msh->ver[i];

   // Setup the temporary triangle for local geometric calculations
   for(i=0;i<3;i++)
   {
      msh->tri.ver[i] = &msh->ver[i];
      msh->tri.edg[i].ver[0] = &msh->ver[ (i+1)%3 ];
      msh->tri.edg[i].ver[1] = &msh->ver[ (i+2)%3 ];
   }

   // Setup the temporary hex for local geometric calculations
   for(i=0;i<8;i++)
      OctMsh->hex.ver[i] = &OctMsh->ver[i];

   for(i=0;i<12;i++)
   {
      OctMsh->hex.edg[i].ver[0] =  &OctMsh->ver[ tvpe[i][0] ];
      OctMsh->hex.edg[i].ver[1] =  &OctMsh->ver[ tvpe[i][0] ];
   }

   for(i=0;i<6;i++)
      for(j=0;j<3;j++)
      {
         OctMsh->hex.tri[i][0].ver[j] =  &OctMsh->ver[ tvpf[i][j] ];
         OctMsh->hex.tri[i][1].ver[j] =  &OctMsh->ver[ tvpf[i][ (j+2)%4 ] ];
      }

   for(i=0;i<4;i++)
   {
      OctMsh->hex.edg[i].tng[0] = 1;
      OctMsh->hex.edg[i].tng[1] = 0;
      OctMsh->hex.edg[i].tng[2] = 0;

      OctMsh->hex.edg[i+4].tng[0] = 0;
      OctMsh->hex.edg[i+4].tng[1] = 1;
      OctMsh->hex.edg[i+4].tng[2] = 0;

      OctMsh->hex.edg[i+8].tng[0] = 0;
      OctMsh->hex.edg[i+8].tng[1] = 0;
      OctMsh->hex.edg[i+8].tng[2] = 1;
   }

   for(i=0;i<2;i++)
   {
      OctMsh->hex.tri[0][i].nrm[0] = 1;
      OctMsh->hex.tri[0][i].nrm[1] = 0;
      OctMsh->hex.tri[0][i].nrm[2] = 0;

      OctMsh->hex.tri[1][i].nrm[0] = -1;
      OctMsh->hex.tri[1][i].nrm[1] = 0;
      OctMsh->hex.tri[1][i].nrm[2] = 0;

      OctMsh->hex.tri[2][i].nrm[0] = 0;
      OctMsh->hex.tri[2][i].nrm[1] = 1;
      OctMsh->hex.tri[2][i].nrm[2] = 0;

      OctMsh->hex.tri[3][i].nrm[0] = 0;
      OctMsh->hex.tri[3][i].nrm[1] = -1;
      OctMsh->hex.tri[3][i].nrm[2] = 0;

      OctMsh->hex.tri[4][i].nrm[0] = 0;
      OctMsh->hex.tri[4][i].nrm[1] = 0;
      OctMsh->hex.tri[4][i].nrm[2] = 1;

      OctMsh->hex.tri[5][i].nrm[0] = 0;
      OctMsh->hex.tri[5][i].nrm[1] = 0;
      OctMsh->hex.tri[5][i].nrm[2] = -1;
   }

   // Insert each vertices in the octree
   for(i=1;i<=msh->NmbVer;i++)
   {
      msh->ver[0].idx = i;
      CpyVec(GetPtrCrd(msh, i), msh->ver[0].crd);
      AddVer(msh, OctMsh, &OctMsh->oct, OctMsh->bnd[0], OctMsh->bnd[1]);
   }

   // Insert each edges in the octree
   for(i=1;i<=msh->NmbEdg;i++)
   {
      SetEdg(msh, i);
      AddEdg(msh, OctMsh, &OctMsh->oct, OctMsh->bnd[0], OctMsh->bnd[1]);
   }

   // Insert each triangles in the octree
   for(i=1;i<=msh->NmbTri;i++)
   {
      SetTri(msh, i);
      AddTri(msh, OctMsh, &OctMsh->oct, OctMsh->bnd[0], OctMsh->bnd[1]);
   }

   // Setup an acceleration grid whose buckets point to an octant
   OctMsh->BucSiz = (OctMsh->bnd[1][0] - OctMsh->bnd[0][0]) / (double)OctMsh->NmbBuc;
   crd[0] = OctMsh->bnd[0][0] + OctMsh->BucSiz / 2.;

   for(i=0;i<OctMsh->NmbBuc;i++)
   {
      crd[1] = OctMsh->bnd[0][1] + OctMsh->BucSiz / 2.;

      for(j=0;j<OctMsh->NmbBuc;j++)
      {
         crd[2] = OctMsh->bnd[0][2] + OctMsh->BucSiz / 2.;

         for(k=0;k<OctMsh->NmbBuc;k++)
         {
            buc = &OctMsh->grd[ i * POW(OctMsh->NmbBuc) + j * OctMsh->NmbBuc + k ];
            buc->oct = GetCrd(&OctMsh->oct, OctMsh->GrdLvl, \
                              crd, OctMsh->bnd[0], OctMsh->bnd[1]);
            buc->pos[0] = i;
            buc->pos[1] = j;
            buc->pos[2] = k;
            buc->tag = 0;
            crd[2] += OctMsh->BucSiz;
         }

         crd[1] += OctMsh->BucSiz;
      }

      crd[0] += OctMsh->BucSiz;
   }

   // Return the octree's unique ID
   return((int64_t)OctMsh);
}


/*----------------------------------------------------------------------------*/
/* Free the octants and the links                                             */
/*----------------------------------------------------------------------------*/

LolInt LolFreeOctree(int64_t OctIdx)
{
   LolInt i;
   OctMshSct *OctMsh = (OctMshSct *)OctIdx;
   LolInt MemUse = OctMsh->MemUse;

   FreAllMem(OctMsh);
   memset(OctMsh, 0, sizeof(OctMshSct));

   return(MemUse);
}


/*----------------------------------------------------------------------------*/
/* Search the octree for triangles included in this box                       */
/*----------------------------------------------------------------------------*/

LolInt LolGetBoundingBox(  int64_t OctIdx, LolInt typ, LolInt MaxItm, \
                           LolInt *ItmTab, double MinCrd[3], double MaxCrd[3] )
{
   LolInt i, NmbItm = 0;
   double box[2][3] = { {MinCrd[0], MinCrd[1], MinCrd[2]}, \
                        {MaxCrd[0], MaxCrd[1], MaxCrd[2]} };
   OctMshSct *OctMsh = (OctMshSct *)OctIdx;

   GetBox(  &OctMsh->oct, typ, &NmbItm, MaxItm, ItmTab, OctMsh->msh->FlgTab, \
            box, OctMsh->eps, OctMsh->bnd[0], OctMsh->bnd[1] );

   for(i=0;i<NmbItm;i++)
      OctMsh->msh->FlgTab[ ItmTab[i] ] = 0;

   return(NmbItm);
}


/*----------------------------------------------------------------------------*/
/* Recusrsive coordinates search                                              */
/*----------------------------------------------------------------------------*/

static OctSct *GetCrd(  OctSct *oct, int MaxLvl, double VerCrd[3], \
                        double MinCrd[3], double MaxCrd[3] )
{
   LolInt SonIdx;
   double MidCrd[3], OctMin[3], OctMax[3], SonMin[3], SonMax[3];

   CpyVec(MinCrd, OctMin);
   CpyVec(MaxCrd, OctMax);

   // Search for the smallest octant containing
   // this vertex but stop at the grid level

   while(oct->sub && (oct->lvl < MaxLvl))
   {
      LinCmbVec3(.5, OctMin, .5, OctMax, MidCrd);

      SonIdx =    ((VerCrd[0] < MidCrd[0]) ? 0 : 1) \
               |  ((VerCrd[1] < MidCrd[1]) ? 0 : 2) \
               |  ((VerCrd[2] < MidCrd[2]) ? 0 : 4);

      SetSonCrd(SonIdx, SonMin, SonMax, OctMin, OctMax);
      CpyVec(SonMin, OctMin);
      CpyVec(SonMax, OctMax);
      oct = oct->son + SonIdx;
   }

   return(oct);
}


/*----------------------------------------------------------------------------*/
/* Recusrsive box search                                                      */
/*----------------------------------------------------------------------------*/

static void GetBox(  OctSct *oct, LolInt typ, LolInt *NmbItm, LolInt MaxItm, \
                     LolInt *ItmTab, char *FlgTab, double box[2][3], \
                     double eps, double MinCrd[3], double MaxCrd[3] )
{
   LolInt i;
   LnkSct *lnk;
   double xmid = (MinCrd[0] + MaxCrd[0])/2.;
   double ymid = (MinCrd[1] + MaxCrd[1])/2.;
   double zmid = (MinCrd[2] + MaxCrd[2])/2.;
   double son[8][2][3] = { \
      { {MinCrd[0], MinCrd[1], MinCrd[2]}, {xmid, ymid, zmid} }, \
      { {xmid, MinCrd[1], MinCrd[2]}, {MaxCrd[0], ymid, zmid} }, \
      { {MinCrd[0], ymid, MinCrd[2]}, {xmid, MaxCrd[1], zmid} }, \
      { {xmid, ymid, MinCrd[2]}, {MaxCrd[0], MaxCrd[1], zmid} }, \
      { {MinCrd[0], MinCrd[1], zmid}, {xmid, ymid, MaxCrd[2]} }, \
      { {xmid, MinCrd[1], zmid}, {MaxCrd[0], ymid, MaxCrd[2]} }, \
      { {MinCrd[0], ymid, zmid}, {xmid, MaxCrd[1], MaxCrd[2]} }, \
      { {xmid, ymid, zmid}, {MaxCrd[0], MaxCrd[1], MaxCrd[2]} } };

   if(oct->sub)
   {
      // Recursively intersect the box with the octree
      for(i=0;i<8;i++)
         if(BoxIntBox(box, son[i], eps))
            GetBox(  oct->son+i, typ, NmbItm, MaxItm, ItmTab, FlgTab, \
                     box, eps, son[i][0], son[i][1] );
   }
   else if((lnk = oct->lnk) && (*NmbItm < MaxItm) )
   {
      // When a terminal octant is reached, add its linked entities
      // to the table and flag them to avoid duplicates
      do
      {
         if(lnk->typ != typ)
            continue;

         if(!FlgTab[ lnk->idx ])
         {
            ItmTab[ (*NmbItm)++ ] = lnk->idx;
            FlgTab[ lnk->idx ] = 1;
         }
      }while( (lnk = lnk->nex) && (*NmbItm < MaxItm) );
   }
}


/*----------------------------------------------------------------------------*/
/* Search for the nearest item from a vertex in the octree                    */
/*----------------------------------------------------------------------------*/

LolInt LolGetNearest(int64_t OctIdx, LolInt typ, double VerCrd[3], \
                     double *MinDis, double MaxDis)
{
   OctMshSct *OctMsh = (OctMshSct *)OctIdx;
   LolInt i, ins=0, out=0, MinItm = 0, ini[3];
   double MinCrd[3], MaxCrd[3];
   MshSct *msh = OctMsh->msh;
   BucSct *IniBuc, *buc, *ngb;
   OctMsh->tag++;

   if(MaxDis > 0.)
      *MinDis = MaxDis * MaxDis;
   else
      *MinDis = DBL_MAX;

   // Get the vertex's integer coordinates in the grid
   // and clip it if it stands outside the bounding box

   for(i=0;i<3;i++)
   {
      ini[i] = (VerCrd[i] - OctMsh->bnd[0][i]) / OctMsh->BucSiz;
      ini[i] = MAX(0, ini[i]);
      ini[i] = MIN(OctMsh->NmbBuc-1, ini[i]);
   }

   IniBuc = &OctMsh->grd[ ini[0] * POW(OctMsh->NmbBuc) \
                        + ini[1] * OctMsh->NmbBuc + ini[2] ];

   // Push the octant containing the starting point on the lifo stack
   OctMsh->stk[ ins++ ] = IniBuc;
   IniBuc->tag = OctMsh->tag;

   // Flood fill processing of the grid : 
   // check octant's contents distance against the closest item
   while(ins > out)
   {
      buc = OctMsh->stk[ out++ ];
      GetBucBox(OctMsh, buc, MinCrd, MaxCrd);
      GetOctLnk(msh, typ, VerCrd, &MinItm, MinDis, buc->oct, MinCrd, MaxCrd);

      // Push unprocessed neighbours on the stack as long as they are not too far
      for(i=0;i<6;i++)
      {
         if( !(ngb = GetBucNgb(OctMsh, buc, i)) || (ngb->tag == OctMsh->tag) )
            continue;

         GetBucBox(OctMsh, buc, MinCrd, MaxCrd);

         if(DisVerOct(VerCrd, MinCrd, MaxCrd) < *MinDis)
         {
            OctMsh->stk[ ins++ ] = ngb;
            ngb->tag = OctMsh->tag;
         }
      }
   }

   *MinDis = sqrt(*MinDis);

   return(MinItm);
}


/*----------------------------------------------------------------------------*/
/* Project a vertex on a given enitity: vertex, edge or triangle              */
/*----------------------------------------------------------------------------*/

LolInt LolProjectVertex(int64_t OctIdx, double *VerCrd, \
                        LolInt typ, LolInt MinItm, double *MinCrd)
{
   OctMshSct *OctMsh = (OctMshSct *)OctIdx;
   MshSct *msh = OctMsh->msh;
   VerSct TmpVer;
   LolInt i, EdgFlg = 0;
   double CurDis, MinDis = DBL_MAX;

   if(typ == LolTypVer)
   {
      // Vertex case, there is only one possible projection:
      // the target vertex itself
      CpyVec(GetPtrCrd(msh, MinItm), MinCrd);
      return(1);
   }
   else if(typ == LolTypEdg)
   {
      // Edge case, the closest position may be on the edge itself
      SetEdg(msh, MinItm);
      PrjVerLin(VerCrd, msh->edg.ver[0]->crd, msh->edg.tng, TmpVer.crd);

      if(VerInsEdg(&msh->edg, &TmpVer, OctMsh->eps))
      {
         CpyVec(TmpVer.crd, MinCrd);
         return(2);
      }

      // Or one of its two vertices
      if(dis(VerCrd, msh->edg.ver[0]->crd) < dis(VerCrd, msh->edg.ver[1]->crd))
         CpyVec(msh->edg.ver[0]->crd, MinCrd);
      else
         CpyVec(msh->edg.ver[1]->crd, MinCrd);

      return(1);
   }
   else if(typ == LolTypTri)
   {
      // Triangle case, the closest position may be on the triangle itself
      SetTri(msh, MinItm);
      PrjVerPla(VerCrd, msh->tri.ver[0]->crd, msh->tri.nrm, TmpVer.crd);

      if(VerInsTri(&msh->tri, &TmpVer, OctMsh->eps))
      {
         CpyVec(TmpVer.crd, MinCrd);
         return(3);
      }

      // Or fall inside one of its three edges
      for(i=0;i<3;i++)
      {
         PrjVerLin(VerCrd, msh->tri.edg[i].ver[0]->crd, msh->tri.edg->tng, TmpVer.crd);

         if(VerInsEdg(msh->tri.edg, &TmpVer, OctMsh->eps) \
         && (dis(VerCrd, TmpVer.crd) < MinDis) )
         {
            CpyVec(TmpVer.crd, MinCrd);
            EdgFlg = 2;
         }
      }

      if(EdgFlg)
         return(2);

      // Or one of the three vertices
      for(i=0;i<3;i++)
      {
         CurDis = dis(VerCrd, msh->tri.ver[i]->crd);

         if(CurDis < MinDis)
         {
            MinDis = CurDis;
            CpyVec(msh->tri.ver[i]->crd, MinCrd);
         }
      }

      return(1);
   }
   else
      return(0);
}


/*----------------------------------------------------------------------------*/
/* Extract the bounding box from a grid's bucket                              */
/*----------------------------------------------------------------------------*/

static void GetBucBox(  OctMshSct *OctMsh, BucSct *buc, \
                        double MinCrd[3], double MaxCrd[3] )
{
   int i;

   for(i=0;i<3;i++)
   {
      MinCrd[i] = OctMsh->bnd[0][i] + buc->pos[i] * OctMsh->BucSiz;
      MaxCrd[i] = OctMsh->bnd[0][i] + (buc->pos[i]+1) * OctMsh->BucSiz;
   }

}

/*----------------------------------------------------------------------------*/
/* Get a bucket's neighbour from the grid                                     */
/*----------------------------------------------------------------------------*/

static BucSct *GetBucNgb(OctMshSct *OctMsh, BucSct *buc, LolInt dir)
{
   if( (dir == 0) && (buc->pos[0] > 0) )
      return(&OctMsh->grd[ (buc->pos[0]-1) * POW(OctMsh->NmbBuc) \
                           + buc->pos[1] * OctMsh->NmbBuc + buc->pos[2] ]);

   if( (dir == 1) && (buc->pos[0] < OctMsh->NmbBuc-1) )
      return(&OctMsh->grd[ (buc->pos[0]+1) * POW(OctMsh->NmbBuc) \
                           + buc->pos[1] * OctMsh->NmbBuc + buc->pos[2] ]);

   if( (dir == 2) && (buc->pos[1] > 0) )
      return(&OctMsh->grd[ buc->pos[0] * POW(OctMsh->NmbBuc) \
                           + (buc->pos[1]-1) * OctMsh->NmbBuc + buc->pos[2] ]);

   if( (dir == 3) && (buc->pos[1] < OctMsh->NmbBuc-1) )
      return(&OctMsh->grd[ buc->pos[0] * POW(OctMsh->NmbBuc) \
                           + (buc->pos[1]+1) * OctMsh->NmbBuc + buc->pos[2] ]);

   if( (dir == 4) && (buc->pos[2] > 0) )
      return(&OctMsh->grd[ buc->pos[0] * POW(OctMsh->NmbBuc) \
                           + buc->pos[1] * OctMsh->NmbBuc + buc->pos[2]-1 ]);

   if( (dir == 5) && (buc->pos[2] < OctMsh->NmbBuc-1) )
      return(&OctMsh->grd[ buc->pos[0] * POW(OctMsh->NmbBuc) \
                           + buc->pos[1] * OctMsh->NmbBuc + buc->pos[2]+1 ]);

   return(NULL);
}


/*----------------------------------------------------------------------------*/
/* Compute the distance between a point and an octant                         */
/*----------------------------------------------------------------------------*/

static double DisVerOct(double VerCrd[3], double MinCrd[3], double MaxCrd[3])
{
   LolInt i;
   double ClpCrd[3];

   // Project the vertex on the octant's surface
   for(i=0;i<3;i++)
   {
      ClpCrd[i] = MAX(VerCrd[i], MinCrd[i]);
      ClpCrd[i] = MIN(ClpCrd[i], MaxCrd[i]);
   }

   // The distance between the projection and the vertex is the shortest
   // if this latter stands OUTSIDE the octant
   return(DisPow(ClpCrd, VerCrd));
}


/*----------------------------------------------------------------------------*/
/* Search for the nearest item from a vertex from an octant                   */
/*----------------------------------------------------------------------------*/

static void GetOctLnk(  MshSct *msh, LolInt typ, double VerCrd[3], \
                        LolInt *MinItm, double *MinDis, OctSct *oct, \
                        double MinCrd[3], double MaxCrd[3] )
{
   LolInt i, *IdxTab;
   double CurDis, SonMin[3], SonMax[3];
   LnkSct *lnk;

   if(oct->sub)
   {
      // Check each sons recursively as long as the minimum distance
      // between the vertex and the octant is lower than
      // the current distance from the closest entity

      for(i=0;i<8;i++)
      {
         SetSonCrd(i, SonMin, SonMax, MinCrd, MaxCrd);

         if(DisVerOct(VerCrd, SonMin, SonMax) <= *MinDis)
            GetOctLnk(msh, typ, VerCrd, MinItm, MinDis, oct->son+i, SonMin, SonMax);
      }
   }
   else if((lnk = oct->lnk))
   {
      // When a leaf octant is reached, compute the distance
      // between its linked enities and the vertex

      do
      {
         if(lnk->typ != typ)
            continue;

         if(lnk->typ == LolTypVer)
            CurDis = DisPow(VerCrd, GetPtrCrd(msh, lnk->idx));
         else if(lnk->typ == LolTypEdg)
         {
            SetEdg(msh, lnk->idx);
            CurDis = DisVerEdg(VerCrd, &msh->edg);
            IdxTab = GetPtrEdg(msh, lnk->idx);
         }
         else if(lnk->typ == LolTypTri)
            CurDis = DisVerTri(msh, VerCrd, lnk->idx);

         if(CurDis < *MinDis)
         {
            *MinItm = lnk->idx;
            *MinDis = CurDis;
         }
      }while((lnk = lnk->nex));
   }
}


/*----------------------------------------------------------------------------*/
/* Setup a edge structure from its index                                      */
/*----------------------------------------------------------------------------*/

static void SetEdg(MshSct *msh, LolInt EdgIdx)
{
   LolInt j, *IdxTab;

   /* Setup the temporary edge structure with this edge's ID */

   msh->edg.idx = EdgIdx;
   IdxTab = GetPtrEdg(msh, EdgIdx);

   for(j=0;j<2;j++)
      CpyVec(GetPtrCrd(msh, IdxTab[j]), msh->edg.ver[j]->crd);

   SetEdgTng(&msh->edg);
}


/*----------------------------------------------------------------------------*/
/* Setup a triangle structure from its index                                  */
/*----------------------------------------------------------------------------*/

static void SetTri(MshSct *msh, LolInt TriIdx)
{
   LolInt j, *IdxTab;

   // Setup the temporary triangle structure with this triangle's ID
   msh->tri.idx = TriIdx;
   IdxTab = GetPtrTri(msh, TriIdx);

   for(j=0;j<3;j++)
      CpyVec(GetPtrCrd(msh, IdxTab[j]), msh->tri.ver[j]->crd);

   SetTriNrm(&msh->tri);

   for(j=0;j<3;j++)
      SetEdgTng(&msh->tri.edg[j]);
}


/*----------------------------------------------------------------------------*/
/* Get a vertex's coordinates table                                           */
/*----------------------------------------------------------------------------*/

static double *GetPtrCrd(MshSct *msh, LolInt idx)
{
   char *adr = &msh->VerCrd[ (idx-1) * msh->VerSiz ];
   return((double *)adr);
}


/*----------------------------------------------------------------------------*/
/* Get a edge's vertices index table                                          */
/*----------------------------------------------------------------------------*/

static LolInt *GetPtrEdg(MshSct *msh, LolInt idx)
{
   char *adr = &msh->EdgIdx[ (idx-1) * msh->EdgSiz ];
   return((LolInt *)adr);
}


/*----------------------------------------------------------------------------*/
/* Get a triangle's vertices index table                                      */
/*----------------------------------------------------------------------------*/

static LolInt *GetPtrTri(MshSct *msh, LolInt idx)
{
   char *adr = &msh->TriIdx[ (idx-1) * msh->TriSiz ];
   return((LolInt *)adr);
}


/*----------------------------------------------------------------------------*/
/* Create a mesh with one hexa enclosing the whole surface                    */
/*----------------------------------------------------------------------------*/

static void SetMshBox(OctMshSct *box, MshSct *msh)
{
   LolInt i, j;
   double MinCrd[3], MaxCrd[3], MidCrd[3], *CrdTab, siz;

   // Compute the bounding box (rectangular)
   CpyVec(GetPtrCrd(msh, 1), MinCrd);
   CpyVec(GetPtrCrd(msh, 1), MaxCrd);

   for(i=2;i<=msh->NmbVer;i++)
   {
      CrdTab = GetPtrCrd(msh, i);

      for(j=0;j<3;j++)
      {
         MinCrd[j] = MIN(MinCrd[j], CrdTab[j]);
         MaxCrd[j] = MAX(MaxCrd[j], CrdTab[j]);
      }
   }

   // Compute the bounding box
   siz = MAX(MaxCrd[0] - MinCrd[0], MaxCrd[1] - MinCrd[1]);
   siz = MAX(siz, MaxCrd[2] - MinCrd[2]);
   box->eps = siz * FLT_EPSILON;
   box->MinSiz = box->MaxSiz = siz;

   // Move the center 1/1000th away
   LinCmbVec3(.5, MinCrd, .5, MaxCrd, MidCrd);
   AddScaVec1(siz * FLT_EPSILON, MidCrd);
   AddScaVec2(-1.02 * siz / 2., MidCrd, box->bnd[0]);
   AddScaVec2( 1.02 * siz / 2., MidCrd, box->bnd[1]);
}


/*----------------------------------------------------------------------------*/
/* Add a vertex to leaf octants                                               */
/*----------------------------------------------------------------------------*/

static void AddVer(  MshSct *msh, OctMshSct *OctMsh, OctSct *oct, \
                     double MinCrd[3], double MaxCrd[3] )
{
   LolInt i;
   double SonMin[3], SonMax[3];

   if(oct->sub)
   {
      for(i=0;i<8;i++)
      {
         SetSonCrd(i, SonMin, SonMax, MinCrd, MaxCrd);

         if(VerInsOct(msh->ver[0].crd, SonMin, SonMax))
            AddVer(msh, OctMsh, oct->son+i, SonMin, SonMax);
      }
   }
   else
   {
      LnkItm(OctMsh, oct, LolTypVer, msh->ver[0].idx);

      if( (oct->lvl < OctMsh->GrdLvl) || ((oct->cpt[ LolTypVer ] \
         >= OctMsh->MaxItm) && (oct->lvl < MaxOctLvl)) )
      {
         SubOct(msh, OctMsh, oct, MinCrd, MaxCrd);
      }
   }
}


/*----------------------------------------------------------------------------*/
/* Add an edge to leaf octants                                                */
/*----------------------------------------------------------------------------*/

static void AddEdg(  MshSct *msh, OctMshSct *OctMsh, OctSct *oct, \
                     double MinCrd[3], double MaxCrd[3] )
{
   LolInt i;
   double SonMin[3], SonMax[3];

   if(oct->sub)
   {
      for(i=0;i<8;i++)
      {
         SetSonCrd(i, SonMin, SonMax, MinCrd, MaxCrd);
         SetTmpHex(&OctMsh->hex, SonMin, SonMax);

         if(EdgIntHex(&msh->edg, &OctMsh->hex, OctMsh->eps))
            AddEdg(msh, OctMsh, oct->son+i, SonMin, SonMax);
      }
   }
   else
   {
      LnkItm(OctMsh, oct, LolTypEdg, msh->edg.idx);

      if( (oct->lvl < OctMsh->GrdLvl) \
      || ((oct->cpt[ LolTypEdg ] >= OctMsh->MaxItm) && (oct->lvl < MaxOctLvl)) )
      {
         SubOct(msh, OctMsh, oct, MinCrd, MaxCrd);
      }
   }
}


/*----------------------------------------------------------------------------*/
/* Add a triangle to leaf octants                                             */
/*----------------------------------------------------------------------------*/

static void AddTri(  MshSct *msh, OctMshSct *OctMsh, OctSct *oct, \
                     double MinCrd[3], double MaxCrd[3] )
{
   LolInt i;
   double SonMin[3], SonMax[3];

   if(oct->sub)
   {
      for(i=0;i<8;i++)
      {
         SetSonCrd(i, SonMin, SonMax, MinCrd, MaxCrd);
         SetTmpHex(&OctMsh->hex, SonMin, SonMax);

         if(TriIntHex(&msh->tri, &OctMsh->hex, OctMsh->eps))
            AddTri(msh, OctMsh, oct->son+i, SonMin, SonMax);
      }
   }
   else
   {
      LnkItm(OctMsh, oct, LolTypTri, msh->tri.idx);

      if( (oct->lvl < OctMsh->GrdLvl) \
      || ((oct->cpt[ LolTypTri ] >= OctMsh->MaxItm) && (oct->lvl < MaxOctLvl)) )
      {
         SubOct(msh, OctMsh, oct, MinCrd, MaxCrd);
      }
   }
}


/*----------------------------------------------------------------------------*/
/* Subdivide an octant and its content to its sons                            */
/*----------------------------------------------------------------------------*/

static void SubOct(  MshSct *msh, OctMshSct *OctMsh, OctSct *oct, \
                     double MinCrd[3], double MaxCrd[3] )
{
   LolInt i, *IdxTab;
   double SonMin[3], SonMax[3];
   LnkSct *lnk , *OctLnk = oct->lnk;
   OctSct *son;

   if(!OctMsh->NmbFreOct)
   {
      OctMsh->CurOctBlk = NewMem(OctMsh, MemBlkSiz * 8 * sizeof(OctSct));
      OctMsh->NmbFreOct = MemBlkSiz;
   }

   oct->son = &OctMsh->CurOctBlk[ (MemBlkSiz - OctMsh->NmbFreOct--) * 8 ];
   oct->sub = 1;
   OctMsh->NmbOct+=8;

   for(i=0;i<8;i++)
   {
      son = oct->son+i;
      son->lnk = NULL;
      son->son = NULL;
      son->cpt[0] = son->cpt[1] = son->sub = 0;
      son->lvl = oct->lvl + 1;
   }

   OctMsh->MinSiz = MIN(OctMsh->MinSiz, (MaxCrd[0] - MinCrd[0])/2.);
   OctMsh->MaxLvl = MAX(OctMsh->MaxLvl, oct->lvl+1);

   while((lnk = OctLnk))
   {
      if(lnk->typ == LolTypVer)
      {
         msh->ver[0].idx = lnk->idx;
         CpyVec(GetPtrCrd(msh, lnk->idx), msh->ver[0].crd);

         for(i=0;i<8;i++)
         {
            SetSonCrd(i, SonMin, SonMax, MinCrd, MaxCrd);

            if(VerInsOct(msh->ver[0].crd, SonMin, SonMax))
               LnkItm(OctMsh, oct->son+i, LolTypVer, lnk->idx);
         }
      }
      else if(lnk->typ == LolTypEdg)
      {
         msh->edg.idx = lnk->idx;
         IdxTab = GetPtrEdg(msh, lnk->idx);

         for(i=0;i<2;i++)
            CpyVec(GetPtrCrd(msh, IdxTab[i]), msh->edg.ver[i]->crd);

         SetEdgTng(&msh->edg);

         for(i=0;i<8;i++)
         {
            SetSonCrd(i, SonMin, SonMax, MinCrd, MaxCrd);
            SetTmpHex(&OctMsh->hex, SonMin, SonMax);

            if(EdgIntHex(&msh->edg, &OctMsh->hex, OctMsh->eps))
               LnkItm(OctMsh, oct->son+i, LolTypEdg, lnk->idx);
         }
      }
      else if(lnk->typ == LolTypTri)
      {
         msh->tri.idx = lnk->idx;
         IdxTab = GetPtrTri(msh, lnk->idx);

         for(i=0;i<3;i++)
            CpyVec(GetPtrCrd(msh, IdxTab[i]), msh->tri.ver[i]->crd);

         SetTriNrm(&msh->tri);

         for(i=0;i<3;i++)
         {
            CpyVec(msh->tri.ver[ (i+1)%3 ]->crd, msh->tri.edg[i].ver[0]->crd);
            CpyVec(msh->tri.ver[ (i+2)%3 ]->crd, msh->tri.edg[i].ver[1]->crd);
            SetEdgTng(&msh->tri.edg[i]);
         }

         for(i=0;i<8;i++)
         {
            SetSonCrd(i, SonMin, SonMax, MinCrd, MaxCrd);
            SetTmpHex(&OctMsh->hex, SonMin, SonMax);

            if(TriIntHex(&msh->tri, &OctMsh->hex, OctMsh->eps))
               LnkItm(OctMsh, oct->son+i, LolTypTri, lnk->idx);
         }
      }

      OctLnk = lnk->nex;
      lnk->nex = OctMsh->NexFreLnk;
      OctMsh->NexFreLnk = lnk;
   }
}


/*----------------------------------------------------------------------------*/
/* Add an entity to an octant linked list                                     */
/*----------------------------------------------------------------------------*/

static void LnkItm(OctMshSct *OctMsh, OctSct *oct, LolInt typ, LolInt idx)
{
   LolInt i;
   LnkSct *lnk;

   if(!OctMsh->NexFreLnk)
   {
      OctMsh->NexFreLnk = NewMem(OctMsh, MemBlkSiz * sizeof(LnkSct));

      for(i=0;i<MemBlkSiz;i++)
         OctMsh->NexFreLnk[i].nex = &OctMsh->NexFreLnk[ i+1 ];

      OctMsh->NexFreLnk[ MemBlkSiz - 1 ].nex = NULL;
   }

   lnk = OctMsh->NexFreLnk;
   OctMsh->NexFreLnk = lnk->nex;
   lnk->typ = typ;
   lnk->idx = idx;
   lnk->nex = oct->lnk;
   oct->lnk = lnk;
   oct->cpt[ typ ]++;
}


/*----------------------------------------------------------------------------*/
/* Build an octant strucutre from the two corner points                       */
/*----------------------------------------------------------------------------*/

static void SetSonCrd(  LolInt SonIdx, double SonMin[3], double SonMax[3], \
                        double MinCrd[3], double MaxCrd[3] )
{
   double MidCrd[3];

   LinCmbVec3(.5, MinCrd, .5, MaxCrd, MidCrd);

   switch(SonIdx)
   {
      case 0 : {
         SonMin[0] = MinCrd[0];
         SonMin[1] = MinCrd[1];
         SonMin[2] = MinCrd[2];
         SonMax[0] = MidCrd[0];
         SonMax[1] = MidCrd[1];
         SonMax[2] = MidCrd[2];
      }return;

      case 1 : {
         SonMin[0] = MidCrd[0];
         SonMin[1] = MinCrd[1];
         SonMin[2] = MinCrd[2];
         SonMax[0] = MaxCrd[0];
         SonMax[1] = MidCrd[1];
         SonMax[2] = MidCrd[2];
      }return;

      case 2 : {
         SonMin[0] = MinCrd[0];
         SonMin[1] = MidCrd[1];
         SonMin[2] = MinCrd[2];
         SonMax[0] = MidCrd[0];
         SonMax[1] = MaxCrd[1];
         SonMax[2] = MidCrd[2];
      }return;

      case 3 : {
         SonMin[0] = MidCrd[0];
         SonMin[1] = MidCrd[1];
         SonMin[2] = MinCrd[2];
         SonMax[0] = MaxCrd[0];
         SonMax[1] = MaxCrd[1];
         SonMax[2] = MidCrd[2];
      }return;

      case 4 : {
         SonMin[0] = MinCrd[0];
         SonMin[1] = MinCrd[1];
         SonMin[2] = MidCrd[2];
         SonMax[0] = MidCrd[0];
         SonMax[1] = MidCrd[1];
         SonMax[2] = MaxCrd[2];
      }return;

      case 5 : {
         SonMin[0] = MidCrd[0];
         SonMin[1] = MinCrd[1];
         SonMin[2] = MidCrd[2];
         SonMax[0] = MaxCrd[0];
         SonMax[1] = MidCrd[1];
         SonMax[2] = MaxCrd[2];
      }return;

      case 6 : {
         SonMin[0] = MinCrd[0];
         SonMin[1] = MidCrd[1];
         SonMin[2] = MidCrd[2];
         SonMax[0] = MidCrd[0];
         SonMax[1] = MaxCrd[1];
         SonMax[2] = MaxCrd[2];
      }return;

      case 7 : {
         SonMin[0] = MidCrd[0];
         SonMin[1] = MidCrd[1];
         SonMin[2] = MidCrd[2];
         SonMax[0] = MaxCrd[0];
         SonMax[1] = MaxCrd[1];
         SonMax[2] = MaxCrd[2];
      }return;
   }
}


/*----------------------------------------------------------------------------*/
/* Test if a vertex is inside an octant                                       */
/*----------------------------------------------------------------------------*/

static LolInt VerInsOct(double VerCrd[3], double MinCrd[3], double MaxCrd[3])
{
   LolInt i;

   for(i=0;i<3;i++)
      if( (VerCrd[i] > MaxCrd[i]) || (VerCrd[i] < MinCrd[i]) )
         return(0);

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Setup a temporary test hex from the bounding box limits                    */
/*----------------------------------------------------------------------------*/

static void SetTmpHex(HexSct *hex, double MinCrd[3], double MaxCrd[3])
{
   hex->ver[0]->crd[0] = MinCrd[0];
   hex->ver[0]->crd[1] = MinCrd[1];
   hex->ver[0]->crd[2] = MaxCrd[2];
   hex->ver[1]->crd[0] = MaxCrd[0];
   hex->ver[1]->crd[1] = MinCrd[1];
   hex->ver[1]->crd[2] = MaxCrd[2];
   hex->ver[2]->crd[0] = MaxCrd[0];
   hex->ver[2]->crd[1] = MinCrd[1];
   hex->ver[2]->crd[2] = MinCrd[2];
   hex->ver[3]->crd[0] = MinCrd[0];
   hex->ver[3]->crd[1] = MinCrd[1];
   hex->ver[3]->crd[2] = MinCrd[2];
   hex->ver[4]->crd[0] = MinCrd[0];
   hex->ver[4]->crd[1] = MaxCrd[1];
   hex->ver[4]->crd[2] = MaxCrd[2];
   hex->ver[5]->crd[0] = MaxCrd[0];
   hex->ver[5]->crd[1] = MaxCrd[1];
   hex->ver[5]->crd[2] = MaxCrd[2];
   hex->ver[6]->crd[0] = MaxCrd[0];
   hex->ver[6]->crd[1] = MaxCrd[1];
   hex->ver[6]->crd[2] = MinCrd[2];
   hex->ver[7]->crd[0] = MinCrd[0];
   hex->ver[7]->crd[1] = MaxCrd[1];
   hex->ver[7]->crd[2] = MinCrd[2];
}


/*----------------------------------------------------------------------------*/
/* Compute/test the intersection between an edge and a hex                    */
/*----------------------------------------------------------------------------*/

static LolInt EdgIntHex(EdgSct *edg, HexSct *hex, double eps)
{
   LolInt i;
   VerSct IntVer;

   // Test if an edge's vertex is included in the octant
   if(VerInsHex(edg->ver[0], hex) || VerInsHex(edg->ver[1], hex))
      return(1);

   // Compute the intersections between the edge and the hex's faces
   for(i=0;i<6;i++)
      if( EdgIntTri(&hex->tri[i][0], edg, &IntVer, eps) \
      ||  EdgIntTri(&hex->tri[i][1], edg, &IntVer, eps) )
      {
         return(1);
      }

   return(0);
}


/*----------------------------------------------------------------------------*/
/* Test if an octant is intersected by a triangle                             */
/*----------------------------------------------------------------------------*/

static LolInt TriIntHex(TriSct *tri, HexSct *hex, double eps)
{
   LolInt i, j, pos, neg;
   double CurDis;
   VerSct IntVer;

   // If there is no intersection between the bounding box
   // of the triangle and the octant it is no use to test,
   // the triangle doesn't intersect the octant

   if((  (tri->ver[0]->crd[0] < hex->ver[3]->crd[0]) \
      && (tri->ver[1]->crd[0] < hex->ver[3]->crd[0]) \
      && (tri->ver[2]->crd[0] < hex->ver[3]->crd[0]) ) \
   || (  (tri->ver[0]->crd[0] > hex->ver[5]->crd[0]) \
      && (tri->ver[1]->crd[0] > hex->ver[5]->crd[0]) \
      && (tri->ver[2]->crd[0] > hex->ver[5]->crd[0]) ) \
   || (  (tri->ver[0]->crd[1] < hex->ver[3]->crd[1]) \
      && (tri->ver[1]->crd[1] < hex->ver[3]->crd[1]) \
      && (tri->ver[2]->crd[1] < hex->ver[3]->crd[1]) ) \
   || (  (tri->ver[0]->crd[1] > hex->ver[5]->crd[1]) \
      && (tri->ver[1]->crd[1] > hex->ver[5]->crd[1]) \
      && (tri->ver[2]->crd[1] > hex->ver[5]->crd[1]) ) \
   || (  (tri->ver[0]->crd[2] < hex->ver[3]->crd[2]) \
      && (tri->ver[1]->crd[2] < hex->ver[3]->crd[2]) \
      && (tri->ver[2]->crd[2] < hex->ver[3]->crd[2]) ) \
   || (  (tri->ver[0]->crd[2] > hex->ver[5]->crd[2]) \
      && (tri->ver[1]->crd[2] > hex->ver[5]->crd[2]) \
      && (tri->ver[2]->crd[2] > hex->ver[5]->crd[2]) ) )
   {
      return(0);
   }

   // Test if a triangle's vertex is included in the octant
   for(i=0;i<3;i++)
      if(VerInsHex(tri->ver[i], hex))
         return(1);

   // Check whether the triangle plane intersects the octant
   pos = neg = 0;
    
   for(i=0;i<8;i++)
   {
      CurDis = DisVerPla(hex->ver[i]->crd, tri->ver[0]->crd, tri->nrm);
    
      if(CurDis < -eps)
         neg = 1;
      else if(CurDis > eps)
         pos = 1;
      else
         pos = neg = 1;
   }
    
   if(!pos || !neg)
      return(0);

   // Compute the intersections between the triangle and the hex edges
   for(i=0;i<12;i++)
      if(EdgIntTri(tri, &hex->edg[i], &IntVer, eps))
         return(1);

   // Compute the intersections between the triangle edges and the hex faces
   for(i=0;i<6;i++)
      for(j=0;j<3;j++)
         if(EdgIntTri(&hex->tri[i][0], &tri->edg[j], &IntVer, eps) \
         || EdgIntTri(&hex->tri[i][1], &tri->edg[j], &IntVer, eps))
         {
            return(1);
         }

   return(0);
}


/*----------------------------------------------------------------------------*/
/* Test if a vertex is inside a hex                                           */
/*----------------------------------------------------------------------------*/

static LolInt VerInsHex(VerSct *ver, HexSct *hex)
{
   LolInt i;

   for(i=0;i<3;i++)
      if( (ver->crd[i] > hex->ver[5]->crd[i]) \
      ||  (ver->crd[i] < hex->ver[3]->crd[i]) )
      {
         return(0);
      }

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Test if an edge intersects a triangle                                      */
/*----------------------------------------------------------------------------*/

static LolInt EdgIntTri(TriSct *tri, EdgSct *edg, VerSct *IntVer, double eps)
{
   LolInt i, NmbVer = 0;
   double sgn[2];
   VerSct *ver=NULL;
   EdgSct edg2;

   // Compute the distance between the edge's vertices and the triangle's plane
   for(i=0;i<2;i++)
   {
      sgn[i] = DisVerPla(edg->ver[i]->crd, tri->ver[0]->crd, tri->nrm);

      if(fabs(sgn[i]) < eps)
      {
         ver = edg->ver[i];
         NmbVer++;
      }
   }

   /* This leads to 3 possible cases */

   switch(NmbVer)
   {
      // Both vertices are far from the plane : 
      // if they stand in the same side of the plane, there is no intersection,
      // otherwise we compute the intersection between the edge and the plane
      // and test wether it falls inside the triangle or not

      case 0 :
      {
         // Test if both vertices stand on opposite sides of the plane
         if(sgn[0] * sgn[1] < 0)
         {
            LinCmbVec3( fabs(sgn[0]) / (fabs(sgn[0]) \
                     +  fabs(sgn[1])), edg->ver[1]->crd, \
                        fabs(sgn[1]) / (fabs(sgn[0]) \
                     +  fabs(sgn[1])), edg->ver[0]->crd, IntVer->crd);

            return(VerInsTri(tri, IntVer, eps));
         }
      }break;

      case 1 :
      {
         if(VerInsTri(tri, ver, eps))
         {
            CpyVec(ver->crd, IntVer->crd);
            return(1);
         }
      }break;

      case 2 :
      {
         for(i=0;i<3;i++)
         {
            edg2.ver[0] = tri->ver[i];
            edg2.ver[1] = tri->ver[ (i+1)%3 ];
            SetEdgTng(&edg2);

            if(EdgIntEdg(edg, &edg2, IntVer, eps))
               return(1);
         }
      }break;
   }

   return(0);
}


/*----------------------------------------------------------------------------*/
/* Test if a vertex is included in an inflated triangle                       */
/*----------------------------------------------------------------------------*/

static LolInt VerInsTri(TriSct *tri, VerSct *ver, double eps)
{
   LolInt i, ins = 1;
   double vec[3][3], nrm[3];
   VerSct img;
   EdgSct edg;

   // Project the vertex on the triangle's plane and check the distance
   if(PrjVerPla(ver->crd, tri->ver[0]->crd, tri->nrm, img.crd) > eps)
      return(0);

   // Compare the normals of the three sub triangles against
   // the original triangle's normal. If one of them is opposite,
   // the vertex does not belong to the triangle

   for(i=0;i<3;i++)
      SubVec3(tri->ver[i]->crd, img.crd, vec[i]);

   for(i=0;i<3;i++)
   {
      CrsPrd(vec[ (i+1)%3 ], vec[i], nrm);

      if(DotPrd(nrm, tri->nrm) <= 0)
      {
         ins = 0;
         break;
      }
   }

   if(ins)
      return(1);

   // If the vertex is not inside the triangle,
   // it may be close to one of its edges

   for(i=0;i<3;i++)
   {
      edg.ver[0] = tri->ver[i];
      edg.ver[1] = tri->ver[ (i+1)%3 ];
      SetEdgTng(&edg);

      if(VerInsEdg(&edg, &img, eps))
         return(1);
   }

   return(0);
}


/*----------------------------------------------------------------------------*/
/* Compute the intersection between coplanar edges                            */
/*----------------------------------------------------------------------------*/

static LolInt EdgIntEdg(EdgSct *edg1, EdgSct *edg2, VerSct *IntVer, double eps)
{
   LolInt i, NmbVer = 0;
   double siz[2];
   VerSct img, *ver=NULL;

   // Compute the distance between the edge1's vertices
   // and the edge2 support line
   for(i=0;i<2;i++)
   {
      PrjVerLin(edg2->ver[i]->crd, edg1->ver[0]->crd, edg1->tng, img.crd);
      siz[i] = dis(edg2->ver[i]->crd, img.crd);

      if(siz[i] < eps)
         NmbVer++;
   }

   // This leads to 3 possible cases
   if(NmbVer < 2)
   {
      LinCmbVec3( siz[0]/(siz[0]+siz[1]), edg2->ver[1]->crd, \
                  siz[1]/(siz[0]+siz[1]), edg2->ver[0]->crd, IntVer->crd);

      return(VerInsEdg(edg1, IntVer, eps));
   }
   else
   {
      ver = NULL;

      if(VerInsEdg(edg2, edg1->ver[0], eps))
         ver = edg1->ver[0];
      else if(VerInsEdg(edg2, edg1->ver[1], eps))
         ver = edg1->ver[1];
      else if(VerInsEdg(edg1, edg2->ver[0], eps))
         ver = edg2->ver[0];
      else if(VerInsEdg(edg1, edg2->ver[1], eps))
         ver = edg2->ver[1];

      if(ver)
      {
         CpyVec(ver->crd, IntVer->crd);
         return(1);
      }
   }

   return(0);
}


/*----------------------------------------------------------------------------*/
/* Test if a vertex belongs to an edge                                        */
/*----------------------------------------------------------------------------*/

static LolInt VerInsEdg(EdgSct *edg, VerSct *ver, double eps)
{
   LolInt i;
   double u[3], v[3];
   VerSct img;

   // Project the vertex on the edge's support line
   PrjVerLin(ver->crd, edg->ver[0]->crd, edg->tng, img.crd);

   if(DisPow(ver->crd, img.crd) > POW(eps))
      return(0);

   // Check if the image belongs to the edge
   SubVec3(img.crd, edg->ver[0]->crd, u);
   SubVec3(img.crd, edg->ver[1]->crd, v);

   if(DotPrd(u, v) < 0)
      return(1);

   // If the vertex does not belong to the edge,
   // it may be close to one of its vertices
   for(i=0;i<2;i++)
      if(DisPow(img.crd, edg->ver[i]->crd) < POW(eps))
         return(1);

   return(0);
}


/*----------------------------------------------------------------------------*/
/* Compute the distance between a vertex and a triangle                       */
/*----------------------------------------------------------------------------*/

static double DisVerTri(MshSct *msh, double VerCrd[3], LolInt TriIdx)
{
   LolInt i, *IdxTab, cod = 0, inc = 1;
   double dis1, dis2, TriSrf, SubSrf, TotSrf=0.;
   VerSct img, TriVer[3];
   EdgSct edg;
   TriSct SubTri, tri;

   // Fetch the triangle and its vertices data
   IdxTab = GetPtrTri(msh, TriIdx);

   for(i=0;i<3;i++)
   {
      tri.ver[i] = &TriVer[i];
      CpyVec(GetPtrCrd(msh, IdxTab[i]), TriVer[i].crd);
   }

   // Compute the triangle's normal and surface
   // and project the coordinates on its plane
   GetTriVec(&tri, tri.nrm);

   if((TriSrf = GetNrmVec(tri.nrm)))
      MulVec1(1./TriSrf, tri.nrm);

   SubTri.ver[2] = &img;

   dis1 = PrjVerPla(VerCrd, tri.ver[0]->crd, tri.nrm, img.crd);

   // Compute the barycentric coordinates and check the projection's position
   for(i=0;i<3;i++)
   {
      SubTri.ver[0] = tri.ver[ (i+1)%3 ];
      SubTri.ver[1] = tri.ver[ (i+2)%3 ];

      GetTriVec(&SubTri, SubTri.nrm);
      SubSrf = GetNrmVec(SubTri.nrm);
      TotSrf += SubSrf;

      if(DotPrd(SubTri.nrm, tri.nrm) < 0.)
         cod |= inc;

      inc = inc << 1;
   }

   // If the sum of the sub triangles surfaces is equal
   // to the main triangle's one, the projection lies inside the triangle
   if( (TotSrf - TriSrf) < .00001 * (TotSrf + TriSrf))
      return(POW(dis1));

   // Otherwise, compute the distance between an edge
   // or a vertex depending on the position code
   switch(cod)
   {
      // Facing edge 0 (1-2) or the barycentric coordinates are degenerate
      case 0 : case 1 :
      {
         edg.ver[0] = tri.ver[1];
         edg.ver[1] = tri.ver[2];
         SetEdgTng(&edg);
         return(DisVerEdg(VerCrd, &edg));
      }

      // Facing edge 1 (2-0)
      case 2 :
      {
         edg.ver[0] = tri.ver[2];
         edg.ver[1] = tri.ver[0];
         SetEdgTng(&edg);
         return(DisVerEdg(VerCrd, &edg));
      }

      // Facing vertex 2
      case 3 : return(DisPow(VerCrd, tri.ver[2]->crd));

      // Facing edge 2 (0-1)
      case 4 :
      {
         edg.ver[0] = tri.ver[0];
         edg.ver[1] = tri.ver[1];
         SetEdgTng(&edg);
         return(DisVerEdg(VerCrd, &edg));
      }

      // Facing vertex 1
      case 5 : return(DisPow(VerCrd, tri.ver[1]->crd));

      // Facing vertex 0
      default : case 6 : return(DisPow(VerCrd, tri.ver[0]->crd));
   }
}


/*----------------------------------------------------------------------------*/
/* Compute the triangle's surface                                             */
/*----------------------------------------------------------------------------*/

static double GetTriSrf(TriSct *tri)
{
   double nrm[3];

   // Compute the cross-product vector and get its size
   GetTriVec(tri, nrm);
   return(GetNrmVec(nrm) / 2.);
}


/*----------------------------------------------------------------------------*/
/* Compute the distance between a vertex and an edge                          */
/*----------------------------------------------------------------------------*/

static double DisVerEdg(double VerCrd[3], EdgSct *edg)
{
   double dis0, dis1, TotSiz = 0.;
   VerSct img;

   PrjVerLin(VerCrd, edg->ver[0]->crd, edg->tng, img.crd);
   TotSiz = dis(edg->ver[0]->crd, img.crd) + dis(edg->ver[1]->crd, img.crd);

   if( (TotSiz - edg->siz) < .00001 * (TotSiz + edg->siz))
      return(DisPow(VerCrd, img.crd));

   dis0 = DisPow(VerCrd, edg->ver[0]->crd);
   dis1 = DisPow(VerCrd, edg->ver[1]->crd);

   return(MIN(dis0, dis1));
}


/*----------------------------------------------------------------------------*/
/* Compute the triangle's normal vector                                       */
/*----------------------------------------------------------------------------*/

static void GetTriVec(TriSct *tri, double w[3])
{
   double u[3], v[3];

   // Compute the vector product
   SubVec3(tri->ver[1]->crd, tri->ver[0]->crd, u);
   SubVec3(tri->ver[2]->crd, tri->ver[0]->crd, v);
   CrsPrd(v, u, w);
}


/*----------------------------------------------------------------------------*/
/* Compute and set the triangle normal vector                                 */
/*----------------------------------------------------------------------------*/

static void SetTriNrm(TriSct *tri)
{
   // Compute the cross-product vector and normalize it
   GetTriVec(tri, tri->nrm);
   NrmVec(tri->nrm);
}


/*----------------------------------------------------------------------------*/
/* Compute and store the edge's unit tangent                                  */
/*----------------------------------------------------------------------------*/

static void SetEdgTng(EdgSct *edg)
{
   SubVec3(edg->ver[1]->crd, edg->ver[0]->crd, edg->tng);
   edg->siz = GetNrmVec(edg->tng);

   if(edg->siz)
      MulVec1(1./edg->siz, edg->tng);
}


/*----------------------------------------------------------------------------*/
/* Compute the normal projection of a point on a line                         */
/*----------------------------------------------------------------------------*/

static void PrjVerLin(  double VerCrd[3], double LinCrd[3], \
                        double LinTng[3], double ImgCrd[3] )
{
   double dp, u[3];

   // Compute the scalar product of the projected vector and the edge's vector
   SubVec3(VerCrd, LinCrd, u);
   dp = DotPrd(u, LinTng);
   LinCmbVec3(1., LinCrd, dp, LinTng, ImgCrd);
}


/*----------------------------------------------------------------------------*/
/* Compute the normal projection of a vertex on a plane                       */
/*----------------------------------------------------------------------------*/

static double PrjVerPla(double VerCrd[3], double PlaCrd[3], \
                        double PlaNrm[3], double ImgCrd[3] )
{
   double DisPla, u[3];

   // Compute the scalar product between the unit normal N and a vector V 
   // defined by the vertex to project and the base plane vertex

   SubVec3(PlaCrd, VerCrd, u);
   DisPla = DotPrd(u, PlaNrm);
   MulVec2(DisPla, PlaNrm, ImgCrd);
   AddVec2(VerCrd, ImgCrd);

   // Return the absolute ditance on the fly
   return(fabs(DisPla));
}


/*----------------------------------------------------------------------------*/
/* Compute the distance between a vertex and a plane                          */
/*----------------------------------------------------------------------------*/

static double DisVerPla(double VerCrd[3], double PlaCrd[3], double PlaNrm[3])
{
   double vec[3];

   SubVec3(VerCrd, PlaCrd, vec);
   return(DotPrd(vec, PlaNrm));
}


/*----------------------------------------------------------------------------*/
/* Test the intersection of two bounding boxex + epsilon                      */
/*----------------------------------------------------------------------------*/

static LolInt BoxIntBox(double box1[2][3], double box2[2][3], double eps)
{
   if((  ((box1[0][0] > box2[0][0] - eps) && (box1[0][0] < box2[1][0] + eps)) \
      || ((box1[1][0] > box2[0][0] - eps) && (box1[1][0] < box2[1][0] + eps)) \
      || ((box1[0][0] < box2[0][0]      ) && (box1[1][0] > box2[1][0]      ))) \
   && (  ((box1[0][1] > box2[0][1] - eps) && (box1[0][1] < box2[1][1] + eps)) \
      || ((box1[1][1] > box2[0][1] - eps) && (box1[1][1] < box2[1][1] + eps)) \
      || ((box1[0][1] < box2[0][1]      ) && (box1[1][1] > box2[1][1]      ))) \
   && (  ((box1[0][2] > box2[0][2] - eps) && (box1[0][2] < box2[1][2] + eps)) \
      || ((box1[1][2] > box2[0][2] - eps) && (box1[1][2] < box2[1][2] + eps)) \
      || ((box1[0][2] < box2[0][2]      ) && (box1[1][2] > box2[1][2]      ))) ) \
   {
      return(1);
   }

   return(0);
}


/*----------------------------------------------------------------------------*/
/* Various basic operations on vectors                                        */
/*----------------------------------------------------------------------------*/

// Euclidian distance
static double dis(double a[3], double b[3])
{
   int i;
   double siz = 0;

   for(i=0;i<3;i++)
      siz += POW(b[i] - a[i]);

   return(sqrt(siz));
}

// Euclidian distance to the power of two
static double DisPow(double a[3], double b[3])
{
   int i;
   double siz = 0;

   for(i=0;i<3;i++)
      siz += POW(b[i] - a[i]);

   return(siz);
}

// W = U - V
static void SubVec3(double u[3], double v[3], double w[3])
{
   int i;

   for(i=0;i<3;i++)
      w[i] = u[i] - v[i];
}

// Euclidian norm
static void NrmVec(double u[3])
{
   int i;
   double dp = 0;

   for(i=0;i<3;i++)
      dp += u[i] * u[i];

   if(dp < DBL_MIN)
      return;

   dp = 1. / sqrt(dp);

   for(i=0;i<3;i++)
      u[i] *= dp;
}

// Dot Product
static double DotPrd(double u[3], double v[3])
{
   int i;
   double dp = 0;

   for(i=0;i<3;i++)
      dp += u[i] * v[i];

   return(dp);
}

// Cross product
static void CrsPrd(double u[3], double v[3], double w[3])
{
   w[0] = u[1] * v[2] - u[2] * v[1];
   w[1] = u[2] * v[0] - u[0] * v[2];
   w[2] = u[0] * v[1] - u[1] * v[0];
}

// Linear combinaison: W = a*U + b*V
static void LinCmbVec3(double w1, double v1[3], double w2, double v2[3], double v3[3])
{
   int i;

   for(i=0;i<3;i++)
      v3[i] = w1 * v1[i] + w2 * v2[i];
}

// V = U
static void CpyVec(double u[3], double v[3])
{
   int i;

   for(i=0;i<3;i++)
      v[i] = u[i];
}

// V = V + U
static void AddVec2(double u[3], double v[3])
{
   int i;

   for(i=0;i<3;i++)
      v[i] += u[i];
}

// U = U + [s;s;s]
static void AddScaVec1(double s, double u[3])
{
   int i;

   for(i=0;i<3;i++)
      u[i] += s;
}

// V = U + [s;s;s]
static void AddScaVec2(double s, double u[3], double v[3])
{
   int i;

   for(i=0;i<3;i++)
      v[i] = u[i] + s;
}

// U = w*U
static void MulVec1(const double w, double u[3])
{
   u[0] *= w;
   u[1] *= w;
   u[2] *= w;
}

// V = w*U
static void MulVec2(double w, double u[3], double v[3])
{
   int i;

   for(i=0;i<3;i++)
      v[i] = w * u[i];
}

// Return Euclidian norm
static double GetNrmVec(double u[3])
{
   return(sqrt(POW(u[0]) + POW(u[1]) + POW(u[2])));
}


/*----------------------------------------------------------------------------*/
/* Allocate a chunk of memory and link it to memory allocator's list          */
/*----------------------------------------------------------------------------*/

static void *NewMem(OctMshSct *OctMsh, LolInt siz)
{
   MemSct *mem;
   assert((mem = malloc(sizeof(MemSct))));
   assert((mem->adr = malloc(siz)));
   mem->siz = siz;
   mem->nex = OctMsh->NexMem;
   OctMsh->NexMem = mem;
   OctMsh->MemUse += siz;

   return(mem->adr);
}


/*----------------------------------------------------------------------------*/
/* Cylct through the linked list of allocated block and free them all         */
/*----------------------------------------------------------------------------*/

static void FreAllMem(OctMshSct *OctMsh)
{
   MemSct *mem = OctMsh->NexMem, *nex;

   while(mem)
   {
      OctMsh->MemUse -= mem->siz;
      nex = mem->nex;
      free(mem->adr);
      free(mem);
      mem = nex;
   }
}


/*----------------------------------------------------------------------------*/
/* Fortran 77 API                                                             */
/*----------------------------------------------------------------------------*/

int64_t call(lolnewoctreef77)(LolInt *NmbVer, double *VerTab1, double *VerTab2, \
                              LolInt *NmbEdg, LolInt *EdgTab1, LolInt *EdgTab2, \
                              LolInt *NmbTri, LolInt *TriTab1, LolInt *TriTab2, \
                              LolInt *NmbQad, LolInt *QadTab1, LolInt *QadTab2, \
                              LolInt *NmbTet, LolInt *TetTab1, LolInt *TetTab2, \
                              LolInt *NmbPyr, LolInt *PyrTab1, LolInt *PyrTab2, \
                              LolInt *NmbPri, LolInt *PriTab1, LolInt *PriTab2, \
                              LolInt *NmbHex, LolInt *HexTab1, LolInt *HexTab2 )
{
   return(LolNewOctree(*NmbVer, VerTab1, VerTab2, *NmbEdg, EdgTab1, EdgTab2, \
                        *NmbTri, TriTab1, TriTab2, *NmbQad, QadTab1, QadTab2, \
                        *NmbTet, TetTab1, TetTab2, *NmbPyr, PyrTab1, PyrTab2, \
                        *NmbPri, PriTab1, PriTab2, *NmbHex, HexTab1, HexTab2 ));
}

LolInt call(lolfreeoctreef77)(int64_t *OctIdx)
{
   return(LolFreeOctree(*OctIdx));
}   

LolInt call(lolgetboundingboxf77)(int64_t *OctIdx, LolInt *typ, \
            LolInt *MaxItm,LolInt *ItmTab, double *MinCrd, double *MaxCrd)
{
   return(LolGetBoundingBox(*OctIdx, *typ, *MaxItm, ItmTab, MinCrd, MaxCrd));
}

LolInt call(lolgetnearestf77)(int64_t *OctIdx, LolInt *typ, double *MinCrd, \
            double *MinDis, double *MaxDis)
{
   return(LolGetNearest(*OctIdx, *typ, MinCrd, MinDis, *MaxDis));
}
