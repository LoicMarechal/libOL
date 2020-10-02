

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                      LIB OCTREE LOCALISATION V1.63                         */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    Description:         Octree for mesh localization                       */
/*    Author:              Loic MARECHAL                                      */
/*    Creation date:       mar 16 2012                                        */
/*    Last modification:   oct 02 2020                                        */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Public defines : type of elements and 32/64 bits int and floating points   */
/*----------------------------------------------------------------------------*/

enum TypTag {  LolTypVer, LolTypEdg, LolTypTri, LolTypQad,
               LolTypTet, LolTypPyr, LolTypPri, LolTypHex, LolNmbTyp };


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

int64_t LolNewOctree       (itg, fpn *, fpn *, itg, itg *, itg *,
                            itg, itg *, itg *, itg, itg *, itg *,
                            itg, itg *, itg *, itg, itg *, itg *,
                            itg, itg *, itg *, itg, itg *, itg *, itg, itg);
size_t  LolFreeOctree      (int64_t);
itg     LolGetBoundingBox  (int64_t , itg , itg, itg *, fpn [3], fpn [3]);
itg     LolGetNearest      (int64_t, itg, fpn *, fpn *, fpn, itg (void *, itg),
                            void * , itg);
itg     LolIntersectSurface(int64_t, fpn *, fpn *, fpn *,
                            fpn, itg (void *, itg), void * );
itg     LolProjectVertex   (int64_t, fpn *, itg, itg, fpn *);


/*----------------------------------------------------------------------------*/
/* Fortran 77 API                                                             */
/*----------------------------------------------------------------------------*/

#if defined(F77_NO_UNDER_SCORE)
#define call(x) x
#else
#define call(x) x ## _
#endif
