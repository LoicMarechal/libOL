

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                      LIB OCTREE LOCALISATION V1.40                         */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    Description:         Octree for mesh localization                       */
/*    Author:              Loic MARECHAL                                      */
/*    Creation date:       mar 16 2012                                        */
/*    Last modification:   jan 30 2017                                        */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Public defines : type of elements and 32/64 bits int                       */
/*----------------------------------------------------------------------------*/

enum TypTag {  LolTypVer, LolTypEdg, LolTypTri, LolTypQad, \
               LolTypTet, LolTypPyr, LolTypPri, LolTypHex, LolNmbTyp };

#ifdef i8
#define LolInt int64_t
#else
#define LolInt int32_t
#endif


/*----------------------------------------------------------------------------*/
/* Public procedures                                                          */
/*----------------------------------------------------------------------------*/

int64_t LolNewOctree(LolInt, double *, double *, LolInt, LolInt *, LolInt *, \
                     LolInt, LolInt *, LolInt *, LolInt, LolInt *, LolInt *, \
                     LolInt, LolInt *, LolInt *, LolInt, LolInt *, LolInt *, \
                     LolInt, LolInt *, LolInt *, LolInt, LolInt *, LolInt * );

LolInt LolFreeOctree(int64_t);
LolInt LolGetBoundingBox(  int64_t , LolInt , LolInt, LolInt *, \
                           double [3], double [3] );

LolInt LolGetNearest(int64_t, LolInt, double [3], double *, double);
LolInt LolProjectVertex(int64_t, double *, LolInt, LolInt, double *);


/*----------------------------------------------------------------------------*/
/* Fortran 77 API                                                             */
/*----------------------------------------------------------------------------*/

#if defined(F77_NO_UNDER_SCORE)
#define call(x) x
#else
#define call(x) x ## _
#endif
