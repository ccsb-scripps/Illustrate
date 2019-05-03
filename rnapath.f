C rnapath.f
C Copyright 2019 David S. Goodsell, the Scripps Research Institute
C
C Licensed under the Apache License, Version 2.0 (the "License");
C you may not use this file except in compliance with the License.
C You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
C Unless required by applicable law or agreed to in writing, software
C distributed under the License is distributed on an "AS IS" BASIS,
C WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
C See the License for the specific language governing permissions and
C limitations under the License
C
C This work was supported in part by the US National Institutes of Health R01-GM120604
C and P50-GM103368 to the HIVE Center.
C
C Reads pruned quasisymmetrical lattice coordinates (from quasisymmetry.f)
C Generates two random walk chains extending from a random point on the lattice
C    ninterp points on each edge, ngenome total points in each chain
C
C reads from Unit 5:
C quasi_71.pdb      # input file name
C rnapath.pdb       # output file name
C 21,3058           # ninterp, ngenome
C 1234              # random seed
C
	character*80 line
	real*4 x(1000),y(1000),z(1000)
	integer*4 used(1000),neighbor(1000,6),n_neighbor(1000)
	integer*4 order(1000)
	real*4 xout(10000),yout(10000),zout(10000)

	read(5,*) line
	open(1,file=line)
	write(6,*) "input coords ",line
	read(5,*) line
	open(2,file=line)
	write(6,*) "output spherical coords ",line
	read(5,*) ninterp,ngenome
	write(6,*) "number of points to interpolate ",ninterp
	write(6,*) "total number of points ",ngenome
	read(5,*) rsphere
	write(6,*) "sphere radius: ",rsphere
	read(5,*) iseed
	do i=1,iseed
	r=rand()
	enddo

C ** read pruned quasisymmetrical lattice points
	i=1
 1	read(1,2,end=9) x(i),y(i),z(i)
 2	format(30x,3f8.3)
	i=i+1
	goto 1
 9	continue
	natom=i-1
	write(6,*) "natom ",natom
	ratom=float(natom)

C ** find neighboring points
	do ia=1,natom
	n_neighbor(ia)=0
	used(ia)=0
	do in=1,6
	neighbor(ia,in)=0
	enddo
	enddo

	do ia=1,natom
	do ib=1,natom
	rd=sqrt((x(ia)-x(ib))**2+ 
     &       (y(ia)-y(ib))**2+ 
     &       (z(ia)-z(ib))**2) 
	if ((ia.ne.ib).and.(rd.lt.8.0)) then
	 n_neighbor(ia)=n_neighbor(ia)+1
	 neighbor(ia,n_neighbor(ia))=ib
	endif
	enddo
	enddo
	write(6,*) "through initialization"
	write(6,*)


C ** pick random point for chain initiation
	order(1)=int(rand()*ratom)+1
	used(order(1))=1

C ** random walk in two directions as far as you can go
	do istep=2,natom
	do itry=1,20
	if (itry.eq.20) then
	  itot=istep-1
	  goto 100
	endif
	  iprev=order(istep-1)
          irandtry=int(rand()*float(n_neighbor(iprev)))+1
	inext=neighbor(iprev,irandtry)
	if (used(inext).eq.0) then
	used(inext)=1
	order(istep)=inext
	goto 7
	endif
	enddo
 7	continue
	enddo

 100	continue
	write(6,*) "initial path ",itot

C ** lengthen chain by lateral moves to neighboring vertices at random edges
	do itry=1,10000
	istep=int(rand()*float(itot-1))+1
	ia=order(istep)
	ib=order(istep+1)
	do in1=1,n_neighbor(ia)
	do in2=1,n_neighbor(ib)
	 if ((neighbor(ia,in1).eq.neighbor(ib,in2)).and.
     &       (used(neighbor(ia,in1)).eq.0)) then
	  itot=itot+1
	  do imove=itot,istep+2,-1
	    order(imove)=order(imove-1)
	  enddo
	  order(istep+1)=neighbor(ia,in1)
	  used(order(istep+1))=1
	  goto 200
	 endif
	enddo
	enddo
 200	continue
	enddo	
	write(6,*) "final path ",itot

C ** interpolate points along edges
	icount=0
	rscale=2.0
	do ia=1,itot-1
	rxa=x(order(ia))
	rya=y(order(ia))
	rza=z(order(ia))
	rxb=x(order(ia+1))
	ryb=y(order(ia+1))
	rzb=z(order(ia+1))
	ires=ia

	do interp=0,ninterp-1
	rinterp=float(interp)/float(ninterp)
	icount=icount+1

	if (icount.gt.ngenome*2+ninterp) goto 999

	xout(icount)=((rxb-rxa)*rinterp+rxa)*rscale
	yout(icount)=((ryb-rya)*rinterp+rya)*rscale
	zout(icount)=((rzb-rza)*rinterp+rza)*rscale
	enddo

	enddo

 999	continue
	write(6,*) "icount ",icount

C ** write coordinates
	icount=0
	rad=sqrt(xout(1)**2+yout(1)**2+zout(1)**2)
	rsrad=rsphere/rad
	do iout=ngenome,1,-1
	icount=icount+1
	write(2,88) "ATOM",icount," P  ","RNA","A",icount,
     &      xout(iout)*rsrad,yout(iout)*rsrad,zout(iout)*rsrad,1.,10.
	enddo
	icount=0
	do iout=ngenome+1,ngenome*2
	icount=icount+1
	write(2,88) "ATOM",iout," P  ","RNA","B",icount,
     &      xout(iout)*rsrad,yout(iout)*rsrad,zout(iout)*rsrad,1.,10.
	enddo

 88	format(a4,2x,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3,2f6.2)
	end
