

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                      LIB OCTREE LOCALISATION V1.60                         */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    Description:         Octree for mesh localization                       */
/*    Author:              Loic MARECHAL                                      */
/*    Creation date:       mar 16 2012                                        */
/*    Last modification:   jun 19 2019                                        */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Public defines : type of elements and 32/64 bits int                       */
/*----------------------------------------------------------------------------*/

enum TypTag {  LolTypVer, LolTypEdg, LolTypTri, LolTypQad,
               LolTypTet, LolTypPyr, LolTypPri, LolTypHex, LolNmbTyp };


/*----------------------------------------------------------------------------*/
/* Global defines                                                             */
/*----------------------------------------------------------------------------*/

#ifdef INT64
#define itg int64_t
#define utg uint64_t
#else
#define itg int32_t
#define utg uint32_t
#endif


/*----------------------------------------------------------------------------*/
/* Public procedures                                                          */
/*----------------------------------------------------------------------------*/

int64_t LolNewOctree       (itg, double *, double *, itg, itg *, itg *,
                            itg, itg *, itg *, itg, itg *, itg *,
                            itg, itg *, itg *, itg, itg *, itg *,
                            itg, itg *, itg *, itg, itg *, itg *);
size_t  LolFreeOctree      (int64_t);
itg     LolGetBoundingBox  (int64_t , itg , itg, itg *, double [3], double [3]);
itg     LolGetNearest      (int64_t, itg, double *, double *,
                            double, itg (void *, itg), void * );
itg     LolIntersectSurface(int64_t, double *, double *, double *,
                            double, itg (void *, itg), void * );
itg     LolProjectVertex   (int64_t, double *, itg, itg, double *);


/*----------------------------------------------------------------------------*/
/* Fortran 77 API                                                             */
/*----------------------------------------------------------------------------*/

#if defined(F77_NO_UNDER_SCORE)
#define call(x) x
#else
#define call(x) x ## _
#endif
