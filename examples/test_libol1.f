      include 'libmeshb7.ins'
      include 'libol1.ins'

      integer i, j, NmbVer, NmbTri, NmbItm, ver, dim, ref, idx
      integer res, EdgTab(2,10000000), TriTab(3,10000000), buf(10000)
      integer*8 mem, MshIdx, OctIdx
      real*8 VerTab(3,10000000), crd1(3), box(3,2), dis

c     Open the mesh file
      MshIdx = gmfopenmesh('../sample_meshes/test.meshb'
     +,GmfRead,ver,dim)
      print*, 'MshIdx = ', MshIdx, ' version = ', ver, ' dim = ', dim

      if(ver.lt.2) STOP ' version < 2'
      if(dim.ne.3) STOP ' dimension <> 3'

c     Check memory bounds
      NmbVer = gmfstatkwd(MshIdx, GmfVertices)
      print*, 'NmbVer = ', NmbVer
      if(NmbVer.le.0) STOP ' NmbVer <= 0'
      if(NmbVer.gt.10000000) STOP ' NmbVer > 10000000'

      NmbEdg = gmfstatkwd(MshIdx, GmfEdges)
      print*, 'NmbEdg = ', NmbEdg
      if(NmbEdg.gt.10000000) STOP ' NmbEdg > 10000000'

      NmbTri = gmfstatkwd(MshIdx, GmfTriangles)
      print*, 'NmbTri = ', NmbTri
      if(NmbTri.gt.10000000) STOP ' NmbTri > 10000000'

c     Read the vertices
      res = gmfgotokwd(MshIdx, GmfVertices)

      do i = 1, NmbVer
          res = gmfgetlin(MshIdx, GmfVertices
     +,VerTab(1,i), VerTab(2,i), VerTab(3,i), ref)
      end do

c     Read the edges
      res = gmfgotokwd(MshIdx, GmfEdges)

      do i = 1, NmbEdg
          res = gmfgetlin(MshIdx GmfEdges
     +, EdgTab(1,i), EdgTab(2,i), ref)
      end do

c     Read the triangles
      res = gmfgotokwd(MshIdx, GmfTriangles)

      do i = 1, NmbTri
          res = gmfgetlin(MshIdx, GmfTriangles
     +, TriTab(1,i), TriTab(2,i), TriTab(3,i), ref)
      end do

c     Close the mesh file
      res = gmfclosemesh(MshIdx)


c     Build the octree
      OctIdx = lolnewoctree(
     +  NmbVer, VerTab(1,1), VerTab(1,2)
     +, NmbEdg, EdgTab(1,1), EdgTab(1,2)
     +, NmbTri, TriTab(1,1), TriTab(1,2)
     +, 0, %val(0), %val(0)
     +, 0, %val(0), %val(0)
     +, 0, %val(0), %val(0)
     +, 0, %val(0), %val(0)
     +, 0, %val(0), %val(0))
      print*, 'OctIdx = ', OctIdx

c     Find the closest vertex and triangle from a given set of coordinates
      crd1(1) = 0
      crd1(2) = 0
      crd1(3) = 0

      idx = lolgetnearest(OctIdx, LolTypVer, crd1, dis, 0D0)
      print*, 'vertex closest from ', crd1, 'is ', idx

      idx = lolgetnearest(OctIdx, LolTypTri, crd1, dis, 0D0)
      print*, 'triangle closest from ', crd1, 'is ', idx

c     Find the triangles included in a bounding box
      box(1,1) = -0.000893 - 0.0001
      box(2,1) =  0.206075 - 0.0001
      box(3,1) =  0.012927 - 0.0001
      box(1,2) = -0.000893 + 0.0001
      box(2,2) =  0.206075 + 0.0001
      box(3,2) =  0.012927 + 0.0001

      NmbItm = lolgetboundingbox(OctIdx, LolTypTri, 10000, buf,
     +box(1,1), box(1,2))
      print*, NmbItm, ' triangles in box ', box(1:3,1), box(1:3,2)

      do i=1,NmbItm
         print*, buf(i)
      enddo

      mem = lolfreeoctree(OctIdx)
      print*, 'memory used = ', mem

      end
