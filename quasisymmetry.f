C quasisymmetry.f
C 2019 David S. Goodsell, the Scripps Research Institute
C
C Licensed under the Apache License, Version 2.0 (the "License");
C you may not use this file except in compliance with the License.
C You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
C Unless required by applicable law or agreed to in writing, software
C distributed under the License is distributed on an "AS IS" BASIS,
C WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
C See the License for the specific language governing permissions and
C
C This work was supported in part by the US National Institutes of Health R01-GM120604
C and P50-GM103368 to the HIVE Center.
C
C Creates a quasisymmetrical lattice with user-defined H and K and edge length=redge
C Projects points onto a sphere, excludes points with large deformations (at vertices)
C Selects ngag/6 points to represent gag protein hexamers
C Writes a standard PDB-format coordinate file
C
C reads from unit 5:
C quasisymmetry.pdb      # output file name
C 7,1,2000,7.5           # H,K,ngag,redge
C 1234                   # random seed
C
C requires two input data files:
C icos.mat -- transformation matrices for icosahedral symmetry
C Tlattice.pdb -- coordinates for quasisymmetrical tiling of triangular faces
C
	character*80 filename
	character*30 stuffin(1000),stuff(1000),stuffout(1000)
	real*4 xin(1000),yin(1000),zin(1000)
	real*4 hin(1000),kin(1000)
	real*4 h(1000),k(1000)
	real*4 hout(100000),kout(100000)
	real*4 x(1000),y(1000),z(1000)
	real*4 xout(100000),yout(100000),zout(100000)
	integer*4 edge(10000,6)
	real*4 rdist(100000),rcount(100000),rmin(100000),rmax(100000)
	integer*4 occ(100000,100),neighbor(100000,6)

c	initialize edge array
	do i=1,10000
	do j=1,6
	  edge(i,j)=0
	enddo
	enddo

C file "Tlattice.pdb" includes a triangular mesh that will be translated
C and clipped to give the final structure, indexed by H and K

 	i=1
	open(1,file="Tlattice.pdb")
 1	read(1,2,end=9) stuffin(i),xin(i),yin(i),zin(i),hin(i),kin(i)
 	  i=i+1
 	  goto 1
 2	format(a30,3f8.3,2f6.2)
 9	nin=i-1

C read H,K, and number of gag polyproteins
C typical values: 7,1,2000
	read(5,*) filename
	open(2,file=filename)
	read(5,*) rh,rk,ngag
	read(5,*) iseed

	do i=1,iseed
	r=rand()
	enddo

C **H along y!
	do i=1,nin
	if ((hin(i).eq.rh).and.(kin(i).eq.rk)) then
	  rht=yin(i)
	  rkt=xin(i)
	endif
	enddo

	rang=-atan(rkt/rht)
	do i=1,nin
	  xnew=xin(i)*cos(rang)+yin(i)*sin(rang)
	  ynew=-xin(i)*sin(rang)+yin(i)*cos(rang)
	  xin(i)=xnew
	  yin(i)=ynew
	enddo

	redge=sqrt(rht*rht+rkt*rkt)
	ricos=redge*0.951057
	rdihedral=20.905*3.14159/180.

C **clip to triangular face**
	n=0
	do i=1,nin
	  if (yin(i).gt.redge) goto 10
	  if (xin(i).lt.0.) goto 10
	  if (yin(i).gt.redge/2.) then
	    rxtri=(redge-yin(i))*0.866/0.5
	    if (xin(i).gt.rxtri) goto 10
	  endif
	  if (yin(i).le.redge/2.) then
	    rxtri=yin(i)*0.866/0.5
	    if (xin(i).gt.rxtri) goto 10
	  endif
	  n=n+1
	  x(n)=xin(i)
	  y(n)=yin(i)
	  z(n)=0.
	  h(n)=hin(i)
	  k(n)=kin(i)
	  stuff(n)=stuffin(i)
 10	continue
	enddo


C **apply icosahedral dihedral to face**
	do i=1,n
	  y(i)=y(i)-redge/2.
	  z(i)=redge*0.80902-x(i)*sin(rdihedral)
	  x(i)=x(i)*cos(rdihedral)
	enddo	

		
C now that we have one face, apply transformation for the whole icosahedron

	open(1,file="icos.mat",form="formatted")

			
	is=0
	iout=0
 111	read(1,200,end=99) rx1,rx2,rx3,tx
	read(1,200) ry1,ry2,ry3,ty
	read(1,200) rz1,rz2,rz3,tz
	if (is.eq.0) then
	  is=1
	else
	  is=0
	endif
		
 	do i=1,n
 	  xn=x(i)*rx1+y(i)*rx2+z(i)*rx3+tx
 	  yn=x(i)*ry1+y(i)*ry2+z(i)*ry3+ty
 	  zn=x(i)*rz1+y(i)*rz2+z(i)*rz3+tz
	do j=1,iout
	  rd=(xn-xout(j))**2+(yn-yout(j))**2+(zn-zout(j))**2
	  if (rd.lt.0.3) then
	    if (h(i)+k(i).lt.hout(j)+kout(j)) then
	       xout(j)=xn
	       yout(j)=yn
	       zout(j)=zn
	       hout(j)=h(i)
	       kout(j)=k(i)
	       stuffout(j)=stuff(i)
	     endif
	     goto 222
	   endif
	enddo
	  iout=iout+1
	  xout(iout)=xn
	  yout(iout)=yn
	  zout(iout)=zn
	  hout(iout)=h(i)
	  kout(iout)=k(i)
	  stuffout(iout)=stuff(i)
 222	continue
 	enddo

	goto 111


 200	format(23x,3f10.6,f15.5)
  99	continue	

C calculate edge array before transforming to a sphere (while all
C points are still spaced by 1 A)

	do i=1,iout-1
	 iedge=0
	do j=i+1,iout
	  r=sqrt((xout(i)-xout(j))**2+
     &           (yout(i)-yout(j))**2+
     &           (zout(i)-zout(j))**2)
	  if ((r.gt.0.9).and.(r.lt.1.1)) then
	    iedge=iedge+1
	    edge(i,iedge)=j
	  endif
	enddo
	enddo
	
C transform to the surface of a sphere
	rtot=0
	do i=1,iout
	  r=sqrt(xout(i)**2+yout(i)**2+zout(i)**2)
	  rtot=rtot+r
	enddo
	rave=rtot/float(iout)
	do i=1,iout
	  r=sqrt(xout(i)**2+yout(i)**2+zout(i)**2)
	  xout(i)=xout(i)*rave/r
	  yout(i)=yout(i)*rave/r
	  zout(i)=zout(i)*rave/r
	enddo

	do i=1,iout
	  rcount(i)=0.
	  rsum=0.
	  rmin(i)=9999.
	  rmax(i)=-9999.
	do j=1,iout
	  if (i.eq.j) goto 77
	  r=sqrt((xout(i)-xout(j))**2+
     &           (yout(i)-yout(j))**2+
     &           (zout(i)-zout(j))**2)
	  if (r.lt.1.2) then
	    rcount(i)=rcount(i)+1.
	    rsum=rsum+r
	    rmin(i)=min(r,rmin(i))
	    rmax(i)=max(r,rmax(i))
	    neighbor(i,int(rcount(i)))=j
	  endif
 77	continue
	enddo
	  rdist(i)=rsum/rcount(i)
	enddo

C this section selects a subset of points, based on the average distance
C and a total number of points--this is an attempt to create sets that look
C like the Briggs EM

	do i=1,iout
	  r=rand()
	do j=1,50
	  occ(i,j)=0
	  if (rdist(i).lt.0.8) occ(i,j)=9
	  if ((rdist(i).lt.0.9).and.(r.lt.0.9)) occ(i,j)=9
	enddo
	enddo

	occ(10,1)=1
	
	rout=float(iout)

	do iadd=2,100
	  itot=0
	do i=1,iout
	  occ(i,iadd)=occ(i,iadd-1)
	  if (occ(i,iadd).eq.1) itot=itot+1
	enddo

	do itrial=1,10000
	  r=rand()
	  i=int(r*rout)

	  if (occ(i,iadd).eq.9) goto 88
	  if (occ(i,iadd).eq.1) goto 88

	  occtot=0.
	  do j=1,6
	  if (occ(neighbor(i,j),iadd-1).eq.1) occtot=occtot+1
	  enddo
	  if (occtot.ge.1) then
	     occ(i,iadd)=1
	     itot=itot+1
 	     iaddout=iadd
             if (itot.gt.ngag/6) goto 55
	  endif

 88	continue
	enddo
	enddo

 55	continue

C write coordinates

	do i=1,iout
	   if (occ(i,iaddout).eq.1) then
 	   write(2,44) "ATOM",i," C  ","ICO","A",1,
     &             xout(i)*redge,yout(i)*redge,zout(i)*redge,
     &             rmin(i),2.-rmax(i)
	   endif
	enddo
 44	format(a4,2x,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3,2f6.2)

	end
