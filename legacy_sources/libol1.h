

/*----------------------------------------------------------*/
/*															*/
/*				LIB OCTREE LOCALISATION V1.30				*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		Octree for mesh localization		*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		mar 16 2012							*/
/*	Last modification:	mar 29 2016							*/
/*															*/
/*----------------------------------------------------------*/


/*----------------------------------------------------------*/
/* Public defines : type of elements and 32/64 bits int		*/
/*----------------------------------------------------------*/

#define TypVer 0
#define TypTri 1

#ifdef i8
#define lng long long
#else
#define lng int
#endif


/*----------------------------------------------------------*/
/* Public procedures										*/
/*----------------------------------------------------------*/

long long NewOctree(lng, double *, double *, lng, lng *, lng *);
lng FreeOctree(long long);
lng GetBoundingBox(long long , lng , lng, lng *, double [3], double [3]);
lng GetNearest(long long, lng, double [3], double *, double);


/*----------------------------------------------------------*/
/* Fortran 77 API											*/
/*----------------------------------------------------------*/

#if defined(F77_NO_UNDER_SCORE)
#define call(x) x
#else
#define call(x) x ## _
#endif
