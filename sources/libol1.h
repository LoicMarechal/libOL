

/*----------------------------------------------------------*/
/*															*/
/*				LIB OCTREE LOCALISATION V1.40				*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		Octree for mesh localization		*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		mar 16 2012							*/
/*	Last modification:	jan 29 2016							*/
/*															*/
/*----------------------------------------------------------*/


/*----------------------------------------------------------*/
/* Public defines : type of elements and 32/64 bits int		*/
/*----------------------------------------------------------*/

enum TypTag {LolTypVer, LolTypEdg, LolTypTri, LolTypQad, LolTypTet, LolTypPyr, LolTypPri, LolTypHex, LolNmbTyp};

#ifdef i8
#define LolInt long long
#else
#define LolInt int
#endif


/*----------------------------------------------------------*/
/* Public procedures										*/
/*----------------------------------------------------------*/

long long LolNewOctree( LolInt, double *, double *, LolInt, LolInt *, LolInt *, \
                        LolInt, LolInt *, LolInt *, LolInt, LolInt *, LolInt *, \
                        LolInt, LolInt *, LolInt *, LolInt, LolInt *, LolInt *, \
                        LolInt, LolInt *, LolInt *, LolInt, LolInt *, LolInt * );

LolInt LolFreeOctree(long long);
LolInt LolGetBoundingBox(long long , LolInt , LolInt, LolInt *, double [3], double [3]);
LolInt LolGetNearest(long long, LolInt, double [3], double *, double);


/*----------------------------------------------------------*/
/* Fortran 77 API											*/
/*----------------------------------------------------------*/

#if defined(F77_NO_UNDER_SCORE)
#define call(x) x
#else
#define call(x) x ## _
#endif
