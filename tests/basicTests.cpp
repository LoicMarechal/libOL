#include <libol1.h>
#include <iostream>
#include <assert.h>
#include <math.h>

void check(bool b)
{
  if (not b){
    assert(false); // in debug, we have the backtrace
    exit(1); // in release, the code stops and return error code
  }
}


int testOnTriangles()
{
  const int nb_vert = 4;
  const int nb_tri = 2;
  double coords[ 13 ] = {0,  0,0,0,  1,0,0,  1,1,0,  0,1,0  };
  int def_tri[ 7 ] = {0, 1,2,4, 2,3,4 };

  auto octreeIdx = LolNewOctree(nb_vert, coords+1, coords+4,  //itg NmbVer, const fpn *PtrCrd1, const fpn *PtrCrd2,
				0, nullptr, nullptr,          //      itg NmbEdg, const itg *PtrEdg1, const itg *PtrEdg2,
				nb_tri, def_tri+1, def_tri+4,  //      itg NmbTri, const itg *PtrTri1, const itg *PtrTri2,
				0, nullptr, nullptr,          //      itg NmbQad, const itg *PtrQad1, const itg *PtrQad2,
				0, nullptr, nullptr,          //      itg NmbTet, const itg *PtrTet1, const itg *PtrTet2,
				0, nullptr, nullptr,          //      itg NmbPyr, const itg *PtrPyr1, const itg *PtrPyr2,
				0, nullptr, nullptr,          //      itg NmbPri, const itg *PtrPri1, const itg *PtrPri2,
				0, nullptr, nullptr,          //      itg NmbHex, const itg *PtrHex1, const itg *PtrHex2,
				1, 1//itg BasIdx, itg NmbThr
				);

  double dist;
  double proj[3];
  int infoProj;
  double pt[3];
  double projExact[3];
  int eltId;

  auto normDiff = [](double p1[3], double p2[3])->double {
    return std::sqrt( (p1[0]-p2[0])*(p1[0]-p2[0])
                     +(p1[1]-p2[1])*(p1[1]-p2[1])
                     +(p1[2]-p2[2])*(p1[2]-p2[2]));
  };

  pt[0]=1.5; pt[1]=-0.5;  pt  [2]=0;
  eltId = LolGetNearest(octreeIdx, LolTypTri, pt, &dist, -1, nullptr, nullptr, 0);
  LolProjectVertex(octreeIdx, pt, LolTypTri, eltId, proj, &infoProj, 0);
  check( eltId == 2);
  check( infoProj==LolVertex0 ) ;
  projExact[0] = 1; projExact[1] = 0; projExact[2] = 0;
  check( normDiff(proj, projExact) < 1e-12 );


  pt[0]=1.1; pt[1]=1.1;  pt[2]=0;
  eltId = LolGetNearest(octreeIdx, LolTypTri, pt, &dist, -1, nullptr, nullptr, 0);
  LolProjectVertex(octreeIdx, pt, LolTypTri, eltId, proj, &infoProj, 0);
  projExact[0] = 1; projExact[1] = 1; projExact[2] = 0;
  check( eltId == 2);
  check(infoProj==LolVertex1);
  check( normDiff(proj, projExact) < 1e-12 );


  pt[0]=-0.1; pt[1]=1.1;  pt[2]=0;
  eltId = LolGetNearest(octreeIdx, LolTypTri, pt, &dist, -1, nullptr, nullptr,0);
  LolProjectVertex(octreeIdx, pt, LolTypTri, eltId, proj, &infoProj, 0);
  projExact[0] = 0; projExact[1] = 1; projExact[2] = 0;
  check( eltId == 2);
  check(infoProj==LolVertex2);
  check( normDiff(proj, projExact) < 1e-12 );


  pt[0]=0; pt[1]=1.1;  pt[2]=0;
  eltId = LolGetNearest(octreeIdx, LolTypTri, pt, &dist, -1, nullptr, nullptr,0);
  LolProjectVertex(octreeIdx, pt, LolTypTri, eltId, proj, &infoProj, 0);
  projExact[0] = 0; projExact[1] = 1; projExact[2] = 0;
  check( eltId == 2);
  check(infoProj==LolVertex2);
  check( normDiff(proj, projExact) < 1e-12 );


  pt[0]=1.1; pt[1]=0.5;  pt[2]=0.0;
  eltId = LolGetNearest(octreeIdx, LolTypTri, pt, &dist, -1, nullptr, nullptr,0);
  LolProjectVertex(octreeIdx, pt, LolTypTri, eltId, proj, &infoProj, 0);
  projExact[0] = 1; projExact[1] = 0.5; projExact[2] = 0;
  check( eltId == 2);
  check(infoProj==LolEdge2);
  check( normDiff(proj, projExact) < 1e-12 );


  pt[0]=0.5; pt[1]=1.1;  pt[2]=0;
  eltId = LolGetNearest(octreeIdx, LolTypTri, pt, &dist, -1, nullptr, nullptr,0);
  LolProjectVertex(octreeIdx, pt, LolTypTri, eltId, proj, &infoProj, 0);
  projExact[0] = 0.5; projExact[1] = 1; projExact[2] = 0;
  check( eltId == 2);
  check(infoProj==LolEdge0);

  return 0;
}




int main()
{
  std::cout << "testOnTriangles ...";
  testOnTriangles();
  std::cout << " DONE\n";

  return 0;
}
