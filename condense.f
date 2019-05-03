C inrelax.f
C Copyright 2019 David S. Goodsell
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
C Reads a collection of points and performs a simulated-annealing-like relaxation
C including user-defined constraints
C compile: gfortran inrelax.f -o inrelax
C run: inrelax
C reads from Unit 5:
C rnapath.pdb               # input coordinate file (here, from rnapath.f)
C 5UTR.const                # input constraint definition file
C NC.position               # input file for location of nucleocapsid
C rnapath_relax.pdb         # output coordinate file
C 1 3058                    # coordinate numbers in chain 1
C 3059 6116                 # coordinate numbers in chain 2
C -0.0001                   # step size for attraction/repulsion (-0.0001 for repulsion, 0.0010 for attraction)
C 1,38.                     # icellflag, if 1, constrain to cell radius
C 10000                     # number of steps
C 1234                      # random seed
C requires two input data files:
C 5UTR.const -- list of constraints
C NC.postion -- list of positions to place explicit nucleocapsid beads
C
	character*80 line,filename
	character*30 pdb(10000)
	character*4 type(5000)
	integer*4 iNC(5000),NCflag(10000)
	real*4 x(0:10000),y(0:10000),z(0:10000)
	integer*4 cspace(-60:60,-60:60,-60:60,200)
	integer*4 nspace(-60:60,-60:60,-60:60)
	integer*4 chain(10,2),start(10000),end(10000)
	real*4 rconstraint(1500),rconstraintsq(1500)
	real*4 radius(1500),radiussq(1500),radiuswRNAsq(1500)
	integer*4 nsphere(1500)
	integer*4 constraint(1500,2),imolconst(1500)
	integer*4 hist_constraint(-10:10),hist_sphsph(-10:10)
	integer*4 hist_clash(0:20),hist_bond(0:20)
	integer*4 hist_sphrna(-10:10),dhist(0:20)
	integer*4 irandatom(5000)
	real*4 rrx(5000),rry(5000),rrz(5000)
	real*4 rxrepulsion(10000),ryrepulsion(10000),rzrepulsion(10000)

	read(5,88) filename
	open(2,file=filename)
	write(6,*) "input lattice coordinate file ",filename
	read(5,88) filename
	open(3,file=filename)
	write(6,*) "constraint file ",filename
	read(5,88) filename
	open(4,file=filename)
	write(6,*) "NC position file ",filename
	read(5,88) filename
	write(6,*) "output coordinate file ",filename
	open(1,file=filename)
 88	format(a80)

C *** hardwired for 3 basepairs per bead
	rscale=0.34*3.
	write(6,*) "base pairs per bead ",3.

C *** manually define the start and stop numbers of two chains
	nchain=2
	do ichain=1,nchain
	read(5,*) chain(ichain,1),chain(ichain,2)
	write(6,*) "chain ends ",chain(ichain,1),chain(ichain,2)
	enddo
	do ia=1,10000
	start(ia)=0
	end(ia)=0
	enddo
	do ichain=1,nchain
	start(chain(ichain,1))=1
	end(chain(ichain,2))=1
	enddo

	read(5,*) rrepulsion_inp
	if (rrepulsion_inp.lt.0) then
	write(6,*) "use general repulsive term for RNA"
	else
	write(6,*) "use general attractive term for RNA and NC"
	endif

	read(5,*) icellflag,rcell
	if (icellflag.ne.0) then
	  write(6,*) "constrain to sphere of radius ",rcell
	  rcell=rcell/rscale
	endif
	rcellsq=rcell**2

	read(5,*) nstep
	write(6,*) "relaxation steps ",nstep

	read(5,*) irand
	write(6,*) "random seed ",irand
	do i=1,irand
	r=rand()
	enddo

C *** read bead positions
	i=1
 1	read(2,3,end=9) line
	if (line(1:4).eq."ATOM") then
         read(line,2) pdb(i),x(i),y(i),z(i)
	 i=i+1
	goto 1
	endif
 3	format(a80)
 2	format(a30,3f8.3,2f6.2)
	goto 1
 9	continue
	natom=i-1

	write(6,*) "number of input beads",natom

	write(6,*) "coordinates of start and end of chains"
	do ichain=1,nchain
	istart=chain(ichain,1)
	iend=chain(ichain,2)
	write(6,*) ichain,start(istart),x(istart),y(istart),z(istart)
	write(6,*) ichain,end(iend),x(iend),y(iend),z(iend)
	enddo

C *** read constraints
	ic=1
 77	read(3,76,end=79) imolconst(ic),constraint(ic,1),constraint(ic,2),
     &   rconstraint(ic),radius(ic),nsphere(ic)
 76	format(3i8,2f8.3,1x,i1,1x,a3)
	rconstraint(ic)=rconstraint(ic)/rscale
	rconstraintsq(ic)=rconstraint(ic)**2
	radius(ic)=radius(ic)/rscale
	radiussq(ic)=radius(ic)**2
C radius of RNA 0.5 + NC 1.0 = 1.5 nm
	rtmp=radius(ic)+1.5/rscale
	radiuswRNAsq(ic)=rtmp**2
	ic=ic+1
	goto 77
 79	nconstraint=ic-1
	write(6,*) "number of constraints ",nconstraint

C *** read NC positions
	nNC=1

	do in=1,10000
	NCflag(in)=0
	enddo

 345	read(4,346,end=349) iNC(nNC),type(nNC)
	do ifootprint=-1,1
	NCflag(iNC(nNC)+ifootprint)=1
	NCflag(iNC(nNC)+chain(2,1)+ifootprint)=1
	enddo
	nNC=nNC+1
	goto 345

 349	continue
	nNC=nNC-1

C *** build contact list
	do ix=-60,60
	do iy=-60,60
	do iz=-60,60
	nspace(ix,iy,iz)=0
	enddo
	enddo
	enddo

	do i=1,natom
	ix=int(x(i))
	iy=int(y(i))
	iz=int(z(i))
	nspace(ix,iy,iz)=nspace(ix,iy,iz)+1
	cspace(ix,iy,iz,nspace(ix,iy,iz))=i
	enddo

	write(6,*) "finished first contact list"

	iramp=0
	itoggle=0
	rnstep=float(nstep)
	irandcount=1

C *** this is the big step loop ***
	do istep=1,nstep

C *** iramp does things every 200 steps
	iramp=iramp+1
	if (iramp.eq.200) then
C *** histograms of how things are fitting to constraints
	write(6,*) "cycle ",istep
	write(6,555) "constr  ",(hist_constraint(ih),ih=-10,10)
	write(6,555) "sphsph  ",(hist_sphsph(ih),ih=-10,10)
	write(6,555) "sphRNA  ",(hist_sphrna(ih),ih=-10,10)
	write(6,556) "RNARNA ",(hist_clash(ih),ih=0,20)
	write(6,556) "bonds  ",(hist_bond(ih),ih=0,20)
C *** report overall size of nucleoid
	rxmid=0.
	rymid=0.
	rzmid=0.
	do ia=1,natom
	rxmid=rxmid+x(ia)
	rymid=rymid+y(ia)
	rzmid=rzmid+z(ia)
	enddo
	rxmid=rxmid/float(natom)
	rymid=rymid/float(natom)
	rzmid=rzmid/float(natom)
	rdtot=0.
	rdmax=-100.
	do ia=1,natom
	rd=sqrt((rxmid-x(ia))**2+
     &          (rymid-y(ia))**2+
     &          (rzmid-z(ia))**2)
	rdmax=max(rdmax,rd)
	rdtot=rdtot+rd
	enddo
	write(6,*) "average/max RNA radius: ",rdtot/float(natom),rdmax
	
 555	format(a8,i7,9i4,i7,3x,9i4,i7)
 556	format(a8,20i4,i7)
 78	format(a17,3f7.4,2x,2f9.6)
	iramp=0

C *** update contact list
	maxcontact=0
	do ix=-60,60
	do iy=-60,60
	do iz=-60,60
	nspace(ix,iy,iz)=0
	enddo
	enddo
	enddo

	do i=1,natom
	ix=x(i)
	iy=y(i)
	iz=z(i)
	nspace(ix,iy,iz)=nspace(ix,iy,iz)+1
	if (nspace(ix,iy,iz).gt.198) then
	   write(6,*) "there is a problem"
           write(6,*) "too many atoms in local neighborhood"
	   write(6,*) "nspace ",nspace(ix,iy,iz),ix,iy,iz
	   write(6,*) i,natom,x(i),y(i),z(i)
	endif
	cspace(ix,iy,iz,nspace(ix,iy,iz))=i
	enddo

	do id=1,10
	dhist(id)=0
	enddo
	do ix=-60,60
	do iy=-60,60
	do iz=-60,60
	id=min(nspace(ix,iy,iz),10)
	dhist(id)=dhist(id)+1
	enddo
	enddo
	enddo
	write(6,130) "dhist ",(dhist(id),id=1,10)
 130	format(a5,10i6)

C *** calculate simple mutual repulsion vector
C     only updated every 200 steps!
	do i=1,natom
	rxsum=0.
	rysum=0.
	rzsum=0.
	rsummax=-10.
	rxrepulsion(i)=0.
	ryrepulsion(i)=0.
	rzrepulsion(i)=0.

	do ia2=1,natom
	 if (abs(ia2-i).gt.4) then
	     rx=x(ia2)-x(i)
	     ry=y(ia2)-y(i)
	     rz=z(ia2)-z(i)
	     rd=sqrt(rx*rx+ry*ry+rz*rz)
	     rxsum=rxsum+rx/(rd*rd)
	     rysum=rysum+ry/(rd*rd)
	     rzsum=rzsum+rz/(rd*rd)
	 endif
	enddo

	  rxrepulsion(i)=-rxsum
	  ryrepulsion(i)=-rysum
	  rzrepulsion(i)=-rzsum

	rnorm=sqrt(rxsum**2+rysum**2+rzsum**2)
	if (rnorm.gt.rsummax) rsummax=rnorm
	enddo

C  normalize based on largest value for all atoms
	do i=1,natom
	rxrepulsion(i)=rxrepulsion(i)/rsummax
	ryrepulsion(i)=ryrepulsion(i)/rsummax
	rzrepulsion(i)=rzrepulsion(i)/rsummax
	enddo

C *** end of 1/200 iramp loop
	endif

C *** update displacements applied for each type of constraint
C bond constraints
	 rllow=(0.95+float(istep)*0.04/rnstep)**2
	 rlhigh=(1.05-float(istep)*0.04/rnstep)**2
	 rbond=0.030

C angle constraint
	 rdiag=(1.0+float(istep)*0.4/rnstep)**2
	 rdiagstep1=0.004
	 rdiagstep2=0.001

C constraint tolerances
	 rclow=-0.05+float(istep)*0.048/rnstep
	 rchigh=+0.05-float(istep)*0.048/rnstep
	 rcstep=0.030

C clash constraints
	rsphsph=0.004
	rsphrna=0.004

C radius of ssRNA, 0.5 nm, + NC, 1.5 nm
	rclash=(3.0+float(istep)*0.0/rnstep)**2
	rrnarna=0.004

C random step
	 rstep=0.04-float(istep)*0.038/rnstep

C general RNA-RNA attraction/repulsion
C positive for attraction, negative for repulsion
	 rrepulsion=rrepulsion_inp

	do ihist=-10,10
	hist_constraint(ihist)=0
	hist_sphsph(ihist)=0
	hist_sphrna(ihist)=0
	enddo
	do ihist=0,20
	hist_clash(ihist)=0
	hist_bond(ihist)=0
	enddo

C *** apply input constraints
	do ic=1,nconstraint

	ia1=constraint(ic,1)
	ia2=constraint(ic,2)

	if (ia1.ne.ia2) then
	rxv=x(ia1)-x(ia2)
	ryv=y(ia1)-y(ia2)
	rzv=z(ia1)-z(ia2)
	rdcsq=rxv**2+ryv**2+rzv**2
	rdc=sqrt(rdcsq)
	
	if (iramp.eq.199) then
	rdiff=rdc-rconstraint(ic)
	idiff=int(rdiff*10.)
	idiff=max(idiff,-10)
	idiff=min(idiff,10)
	hist_constraint(idiff)=hist_constraint(idiff)+1
	endif

	if (rdc.lt.rconstraint(ic)+rclow) then
	rxv=rxv/rdc*rcstep
	ryv=ryv/rdc*rcstep
	rzv=rzv/rdc*rcstep
	x(ia1)=x(ia1)+rxv
	y(ia1)=y(ia1)+ryv
	z(ia1)=z(ia1)+rzv
	x(ia2)=x(ia2)-rxv
	y(ia2)=y(ia2)-ryv
	z(ia2)=z(ia2)-rzv
	endif
	if (rdc.gt.rconstraint(ic)+rchigh) then
	rxv=rxv/rdc*rcstep
	ryv=ryv/rdc*rcstep
	rzv=rzv/rdc*rcstep
	x(ia1)=x(ia1)-rxv
	y(ia1)=y(ia1)-ryv
	z(ia1)=z(ia1)-rzv
	x(ia2)=x(ia2)+rxv
	y(ia2)=y(ia2)+ryv
	z(ia2)=z(ia2)+rzv
	endif

	endif

	enddo

C *** sphere-sphere overlap for constraints
	do ic=1,nconstraint
	if (radius(ic).ne.0) then

	ia1=constraint(ic,1)
	ia2=constraint(ic,2)
	rxc=(x(ia1)+x(ia2))/2.
	ryc=(y(ia1)+y(ia2))/2.
	rzc=(z(ia1)+z(ia2))/2.

	do ic2=1,nconstraint
	if ((ic.ne.ic2).and.(radius(ic2).ne.0.).and.
     &      (imolconst(ic).ne.imolconst(ic2))) then

	ib1=constraint(ic2,1)
	ib2=constraint(ic2,2)
	rxc2=(x(ib1)+x(ib2))/2.
	ryc2=(y(ib1)+y(ib2))/2.
	rzc2=(z(ib1)+z(ib2))/2.

	rxv=rxc2-rxc
	ryv=ryc2-ryc
	rzv=rzc2-rzc
	rdcsq=(rxv**2+ryv**2+rzv**2)
	rdc=sqrt(rdcsq)

	if (iramp.eq.199) then
	rdiff=rdc-radius(ic)-radius(ic2)
	idiff=int(rdiff*10.)
	idiff=max(idiff,-10)
	idiff=min(idiff,10)
	hist_sphsph(idiff)=hist_sphsph(idiff)+1
	endif

	if (rdc.lt.radius(ic)+radius(ic2)) then
	rxv=rxv/rdc*rsphsph
	ryv=ryv/rdc*rsphsph
	rzv=rzv/rdc*rsphsph
	x(ia1)=x(ia1)-rxv
	x(ia2)=x(ia2)-rxv
	y(ia1)=y(ia1)-ryv
	y(ia2)=y(ia2)-ryv
	z(ia1)=z(ia1)-rzv
	z(ia2)=z(ia2)-rzv
	x(ib1)=x(ib1)+rxv
	x(ib2)=x(ib2)+rxv
	y(ib1)=y(ib1)+ryv
	y(ib2)=y(ib2)+ryv
	z(ib1)=z(ib1)+rzv
	z(ib2)=z(ib2)+rzv
	endif

	endif	
	enddo
	endif	
	enddo

C *** small random moves
	ratom=float(natom)-6.

	do irand=1,2000

	ia=irandatom(irand)

	if (irandcount.eq.1) then
	irandatom(irand)=int(rand()*ratom)+3
	rrandx=rand()-0.5
	rrandy=rand()-0.5
	rrandz=rand()-0.5
	rd=sqrt(rrandx**2+rrandy**2+rrandz**2)
	rrandx=rrandx/rd
	rrandy=rrandy/rd
	rrandz=rrandz/rd

	ia=irandatom(irand)

	rxv=x(ia-1)-x(ia+1)
	ryv=y(ia-1)-y(ia+1)
	rzv=z(ia-1)-z(ia+1)
	rd=sqrt(rxv**2+ryv**2+rzv**2)
	rxv=rxv/rd
	ryv=ryv/rd
	rzv=rzv/rd
	rxcross=ryv*rrandz-rzv*rrandy
	rycross=rzv*rrandx-rxv*rrandz
	rzcross=rxv*rrandy-ryy*rrandx
	rd=sqrt(rxcross**2+rycross**2+rzcross**2)
	rrx(irand)=rxcross/rd*rstep
	rry(irand)=rycross/rd*rstep
	rrz(irand)=rzcross/rd*rstep

	endif

C *** apply random step to local area, not just one bead
	do ioff=-3,3
	roff=(4.-(float(abs(ioff))))/4.
	x(ia+ioff)=x(ia+ioff)+rrx(irand)*roff
	y(ia+ioff)=y(ia+ioff)+rry(irand)*roff
	z(ia+ioff)=z(ia+ioff)+rrz(irand)*roff
	enddo

	enddo
 876	format(3i6,9f9.5)
	irandcount=irandcount+1
	if (irandcount.gt.50) irandcount=1

C *** step through RNA ****
	do ia=1,natom

C *** input constraints as spherical clash with RNA
	rmin=999.
	do ic=1,nconstraint
	if (radius(ic).ne.0) then
	ia1=constraint(ic,1)
	ia2=constraint(ic,2)
	if ((ia.ne.ia1).and.(ia.ne.ia2)) then
	rxc=(x(ia1)+x(ia2))/2.
	ryc=(y(ia1)+y(ia2))/2.
	rzc=(z(ia1)+z(ia2))/2.
	rxv=x(ia)-rxc
	ryv=y(ia)-ryc
	rzv=z(ia)-rzc
	rdcsq=rxv**2+ryv**2+rzv**2

	if (iramp.eq.199) then
	rmin=min(rmin,sqrt(rdcsq)-sqrt(radiuswRNAsq(ic)))
	endif

	if (rdcsq.lt.radiuswRNAsq(ic)) then
	rdc=sqrt(rdcsq)
	x(ia)=x(ia)+rxv*rsphrna/rdc
	y(ia)=y(ia)+ryv*rsphrna/rdc
	z(ia)=z(ia)+rzv*rsphrna/rdc
	endif

	endif	
	endif	
	enddo

	if (iramp.eq.199) then
	idiff=int(rmin*2.)
	idiff=max(idiff,-10)
	idiff=min(idiff,10)
	hist_sphrna(idiff)=hist_sphrna(idiff)+1
	endif

C *** bond constraint
C     toggle alternates order of + and - moves
	if ((end(ia).eq.0).and.(itoggle.eq.1)) then
	itoggle=2
	rxv=x(ia)-x(ia+1)
	ryv=y(ia)-y(ia+1)
	rzv=z(ia)-z(ia+1)
	rd3sq=rxv**2+ryv**2+rzv**2
	if (rd3sq.lt.rllow) then
	x(ia)=x(ia)+rxv*rbond
	y(ia)=y(ia)+ryv*rbond
	z(ia)=z(ia)+rzv*rbond
	endif
	if (rd3sq.gt.rlhigh) then
	x(ia)=x(ia)-rxv*rbond
	y(ia)=y(ia)-ryv*rbond
	z(ia)=z(ia)-rzv*rbond
	endif
	endif

	if (start(ia).eq.0) then
	rxv=x(ia)-x(ia-1)
	ryv=y(ia)-y(ia-1)
	rzv=z(ia)-z(ia-1)
	rd3sq=rxv**2+ryv**2+rzv**2
	if (rd3sq.lt.rllow) then
	x(ia)=x(ia)+rxv*rbond
	y(ia)=y(ia)+ryv*rbond
	z(ia)=z(ia)+rzv*rbond
	endif
	if (rd3sq.gt.rlhigh) then
	x(ia)=x(ia)-rxv*rbond
	y(ia)=y(ia)-ryv*rbond
	z(ia)=z(ia)-rzv*rbond
	endif
	endif

	if ((end(ia).eq.0).and.(itoggle.eq.0)) then
	itoggle=1
	rxv=x(ia)-x(ia+1)
	ryv=y(ia)-y(ia+1)
	rzv=z(ia)-z(ia+1)
	rd3sq=rxv**2+ryv**2+rzv**2
	if (rd3sq.lt.rllow) then
	x(ia)=x(ia)+rxv*rbond
	y(ia)=y(ia)+ryv*rbond
	z(ia)=z(ia)+rzv*rbond
	endif
	if (rd3sq.gt.rlhigh) then
	x(ia)=x(ia)-rxv*rbond
	y(ia)=y(ia)-ryv*rbond
	z(ia)=z(ia)-rzv*rbond
	endif
	endif
	if (itoggle.eq.2) itoggle=0

	if (iramp.eq.199) then
	ibond=int(sqrt(rd3sq)*10.)
	ibond=min(ibond,20)
	hist_bond(ibond)=hist_bond(ibond)+1
	endif

C *** bond angles
	if ((start(ia).eq.0).and.(end(ia).eq.0)) then

	rxv=x(ia+1)-x(ia-1)
        ryv=y(ia+1)-y(ia-1)
        rzv=z(ia+1)-z(ia-1)
        rd3sq=rxv**2+ryv**2+rzv**2
        if (rd3sq.lt.rdiag) then
        rxa=(x(ia+1)+x(ia-1))/2.
        rya=(y(ia+1)+y(ia-1))/2.
        rza=(z(ia+1)+z(ia-1))/2.
        x(ia)=x(ia)+(rxa-x(ia))*rdiagstep1
        y(ia)=y(ia)+(rya-y(ia))*rdiagstep1
        z(ia)=z(ia)+(rza-z(ia))*rdiagstep1
        x(ia+1)=x(ia+1)-(rxa-x(ia))*rdiagstep2
        y(ia+1)=y(ia+1)-(rya-y(ia))*rdiagstep2
        z(ia+1)=z(ia+1)-(rza-z(ia))*rdiagstep2
        x(ia-1)=x(ia-1)-(rxa-x(ia))*rdiagstep2
        y(ia-1)=y(ia-1)-(rya-y(ia))*rdiagstep2
        z(ia-1)=z(ia-1)-(rza-z(ia))*rdiagstep2
        endif

	endif

c *** check for clash
	iclash=0
	rmin=999.
C--this search box needs to be large enough to have at least two strands
	do ixo=-2,2
	do iyo=-2,2
	do izo=-2,2
	ix=x(ia)+ixo
	iy=y(ia)+iyo
	iz=z(ia)+izo
	do icontact=1,nspace(ix,iy,iz)
	ic=cspace(ix,iy,iz,icontact)
	if ((ic.lt.ia-4).or.(ic.gt.ia+4)) then
	rcsq=(x(ia)-x(ic))**2+(y(ia)-y(ic))**2+(z(ia)-z(ic))**2
	if (rcsq.lt.rclash) iclash=1
	if (rcsq.lt.rmin) then
	 rvx=x(ia)-x(ic)
	 rvy=y(ia)-y(ic)
	 rvz=z(ia)-z(ic)
	 rmin=rcsq
	endif
	endif
	enddo
	enddo
	enddo
	enddo

	if (iramp.eq.199) then
	imin=int(sqrt(rmin)*10.)
	imin=min(imin,20)
	hist_clash(imin)=hist_clash(imin)+1
	endif

C *** move to resolve clash
	if (iclash.eq.1) then

	rnorm=sqrt(rvx**2+rvy**2+rvz**2)
	if (rnorm.eq.0.) rnorm=1.
	 x(ia)=x(ia)+rrnarna*rvx/rnorm
	 y(ia)=y(ia)+rrnarna*rvy/rnorm
	 z(ia)=z(ia)+rrnarna*rvz/rnorm
	
C this block adds a small motion to flanking regions with clash
	 if (start(ia).eq.0) then
	 x(ia-1)=x(ia-1)+rrnarna*rvx/rnorm/3.
	 y(ia-1)=y(ia-1)+rrnarna*rvy/rnorm/3.
	 z(ia-1)=z(ia-1)+rrnarna*rvz/rnorm/3.
	 endif
	 if (end(ia).eq.0) then
	 x(ia+1)=x(ia+1)+rrnarna*rvx/rnorm/3.
	 y(ia+1)=y(ia+1)+rrnarna*rvy/rnorm/3.
	 z(ia+1)=z(ia+1)+rrnarna*rvz/rnorm/3.
	 endif

	endif

C *** if no clash, apply simple mutual attraction/repulsion

	if (NCflag(ia).eq.1) then
C positions with NC
	x(ia)=x(ia)-rxrepulsion(ia)*rrepulsion
	y(ia)=y(ia)-ryrepulsion(ia)*rrepulsion
	z(ia)=z(ia)-rzrepulsion(ia)*rrepulsion
	else
C positions without NC
	x(ia)=x(ia)-rxrepulsion(ia)*rrepulsion/2.
	y(ia)=y(ia)-ryrepulsion(ia)*rrepulsion/2.
	z(ia)=z(ia)-rzrepulsion(ia)*rrepulsion/2.
	endif

C *** constrain to cell sphere
	if (icellflag.ne.0) then

	radsqcell=x(ia)**2+y(ia)**2+z(ia)**2
	rad=sqrt(radsqcell)
	if (radsqcell.gt.rcellsq) then
	  x(ia)=x(ia)-(x(ia)/rad)*0.004
	  y(ia)=y(ia)-(y(ia)/rad)*0.004
	  z(ia)=z(ia)-(z(ia)/rad)*0.004
	endif

	endif

	enddo

C *** end of the big step loop
	enddo

C *** write out the RNA!
	do i=1,natom
	write(1,2) pdb(i),x(i)*rscale,y(i)*rscale,z(i)*rscale,1.,0.5
	enddo

C *** write out crosslinking proteins 
	icount=natom+1
	do ic=1,nconstraint
	ia1=constraint(ic,1)
	ia2=constraint(ic,2)
	rxc=(x(ia1)+x(ia2))/2.
	ryc=(y(ia1)+y(ia2))/2.
	rzc=(z(ia1)+z(ia2))/2.

	if (radius(ic).ne.0) then
	if (nsphere(ic).eq.1) then
	write(1,700) "ATOM",icount,"C  ","IN ","C",ic,
     &               rxc*rscale,ryc*rscale,rzc*rscale,1.,radius(ic)
	icount=icount+1
	else
	rxc=x(ia1)+(x(ia2)-x(ia1))*0.25
	ryc=y(ia1)+(y(ia2)-y(ia1))*0.25
	rzc=z(ia1)+(z(ia2)-z(ia1))*0.25
	write(1,700) "ATOM",icount,"C  ","IN ","C",ic,
     &               rxc*rscale,ryc*rscale,rzc*rscale,1.,radius(ic)
	icount=icount+1
	rxc=x(ia1)+(x(ia2)-x(ia1))*0.75
	ryc=y(ia1)+(y(ia2)-y(ia1))*0.75
	rzc=z(ia1)+(z(ia2)-z(ia1))*0.75
	write(1,700) "ATOM",icount,"C  ","IN ","C",ic,
     &               rxc*rscale,ryc*rscale,rzc*rscale,1.,radius(ic)
	icount=icount+1
	endif
	endif

	enddo
 700	format(a4,i7,2x,a3,1x,a3,1x,a1,i4,4x,3f8.3,2f6.2)
	
 12	format(a6,i8)

C *** write NC coordinates

	do i=1,nNC
	if (type(i).eq."RAND") then
	type(i)="NCR "
	else
	type(i)="NCO "
	endif
	enddo
	
	do i=1,nNC
	write(1,700) "ATOM",icount,"N  ",type(i)(1:3),"D",iNC(i),
     &       x(iNC(i))*rscale,y(iNC(i))*rscale,z(iNC(i))*rscale,1.,1.
	icount=icount+1
	enddo
	do i=1,nNC
	iout=iNC(i)+natom/2
	write(1,700) "ATOM",icount,"N  ",type(i)(1:3),"E",iNC(i),
     &       x(iout)*rscale,y(iout)*rscale,z(iout)*rscale,1.,1.
	icount=icount+1
	enddo

 346	format(i8,2x,a4)

	end
