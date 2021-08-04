

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                      LIB OCTREE LOCALISATION V1.79                         */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    Description:         Octree for mesh localization                       */
/*    Author:              Loic MARECHAL                                      */
/*    Creation date:       mar 16 2012                                        */
/*    Last modification:   jul 16 2021                                        */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* ANSI C headers                                                             */
/*----------------------------------------------------------------------------*/

#include <stdint.h>
#include <stddef.h>


/*----------------------------------------------------------------------------*/
/* Public defines : type of elements and 32/64 bits int and floating points   */
/*----------------------------------------------------------------------------*/

enum TypTag {  LolTypVer=1, LolTypEdg, LolTypTri, LolTypQad,
               LolTypTet,   LolTypPyr, LolTypPri, LolTypHex, LolNmbTyp };


/*----------------------------------------------------------------------------*/
/* Set integer and floating point sizes:                                      */
/* default values are 32-bit integers and 64-bit floating points              */
/* define INT64 or REAL32 to override default settings                        */
/*----------------------------------------------------------------------------*/

#ifdef INT64
#define itg int64_t
#define utg uint64_t
#else
#define itg int32_t
#define utg uint32_t
#endif

#ifdef REAL32
#define fpn float
#else
#define fpn double
#endif


/*----------------------------------------------------------------------------*/
/* Public procedures                                                          */
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

int64_t  LolNewOctree         (itg, const fpn *, const fpn *,
                               itg, const itg *, const itg *,
                               itg, const itg *, const itg *,
                               itg, const itg *, const itg *,
                               itg, const itg *, const itg *,
                               itg, const itg *, const itg *,
                               itg, const itg *, const itg *,
                               itg, const itg *, const itg *,
                               itg, itg);
int64_t  LolNewOctreeFromSTL  (itg, const fpn *, const fpn *, itg, itg);
size_t   LolFreeOctree        (int64_t);
itg      LolGetBoundingBox    (int64_t , itg , itg, itg *, fpn *, fpn *, itg);
itg      LolGetNearest        (int64_t, itg, fpn *, fpn *, fpn, itg (void *, itg),
                               void * , itg);
itg      LolIntersectSurface  (int64_t, fpn *, fpn *, fpn *,
                               fpn, itg (void *, itg), void *, itg);
itg      LolIsInside          (int64_t, fpn *, fpn *, itg);
itg      LolProjectVertex     (int64_t, fpn *, itg, itg, fpn *, itg);
itg      LolCheckIntersections(int64_t, itg, itg *);

#ifdef __cplusplus
}
#endif


/*----------------------------------------------------------------------------*/
/* Fortran 77 API                                                             */
/*----------------------------------------------------------------------------*/

#if defined(F77_NO_UNDER_SCORE)
#define call(x) x
#else
#define call(x) x ## _
#endif
