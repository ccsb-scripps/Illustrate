C makeconstraints.f
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
C Reads a RNA coordinate file, creates a constraint parameter file
C based on experimental localization of integrase and nucleocapsid
C
C reads from Unit 5:
C rnapath_relax.pdb         # input coordinate file
C inpeaks.dat               # data file with integrase sites
C ncpeaks.dat               # data file with nucleocapsid sites
C nucleoid.const            # output constraint parameter file
C 15.                       # max distance for crosslink
C 0,35,100                  # nINdimer, nINtetrame, nNC
C 1234                      # random seed
C
	character*80 filename,line
	real*4 x(6500),y(6500),z(6500)
	real*4 peakdist(200,200)
	integer*4 peaklist(200,200),npeaklist(200),peakcheck(200)
	integer*4 inpeaks(200),peaks(200),peakflag(200)
	integer*4 inncpeaks(200),ncpeaks(200)
	integer*4 const(100,4)
	integer*4 tmpflag(200)
	integer*4 seqflag(6500)
	integer*4 ip1max(100),ip1save(100)
	integer*4 ip2max(100),ip2save(100)
	integer*4 ip3max(100),ip3save(100)
	integer*4 ip4max(100),ip4save(100)

	do ip=1,200
	peakflag(ip)=0
	enddo

	read(5,*) filename
 	open(1,file=filename)
	write(6,*) "input coordinates ",filename
	read(5,*) filename
 	open(2,file=filename)
	write(6,*) "IN peak file ",filename
	read(5,*) filename
 	open(4,file=filename)
	write(6,*) "NC peak file ",filename
	read(5,*) filename
 	open(3,file=filename)
	write(6,*) "output constraint file ",filename
	read(5,*) filename
 	open(10,file=filename)
	write(6,*) "output NC position file ",filename
	read(5,*) rdmax
	write(6,*) "max distance for crosslink ",rdmax
	rdmaxsq=rdmax**2
	read(5,*) nDIM,nTET,nNC
	write(6,*) "number of integrase dimer: ",nDIM
	write(6,*) "number of integrase tetramer: ",nTET
	nIN=nDIM+nTET
	write(6,*) "number of nucleocapsid: ",nNC
	read(5,*) iseed
	write(6,*) "random seed: ",iseed

	do irand=1,iseed
	r=rand()
	enddo

C ** hardwired for NL4-3, 3 base per bead ***
	ngenome=3058

C ** read coordinates
	i=1
 1	read(1,3,end=9) line
 3	format(a80)
	if (line(18:20).eq."RNA") then
	read(line,2) x(i),y(i),z(i)
	i=i+1
	endif
	goto 1
 2	format(30x,3f8.3)
 9	continue
	write(6,*) "coordinates read",i-1

C *** read IN peaks
	ipeak=1
 10	read(2,*,end=11) inpeaks(ipeak)
	ipeak=ipeak+1
	goto 10
 11	npeak=ipeak-1

	write(6,*) "npeaks IN input ",npeak

	do ip=1,npeak
	peaks(ip)=inpeaks(ip)/3
	peaks(ip+npeak)=inpeaks(ip)/3 + ngenome
	enddo

	rpeak=float(npeak*2)
	rnpeak=float(npeak*2)
	npeak=npeak*2

C *** read NC peaks
	ipeak=1
 15	read(4,*,end=19) inncpeaks(ipeak)
	ipeak=ipeak+1
	goto 15
 19	ncpeak=ipeak-1

	write(6,*) "npeaks NC input ",ncpeak

	do ip=1,ncpeak
	ncpeaks(ip)=inncpeaks(ip)/3
cncpeaks(ip+ncpeak)=inncpeaks(ip)/3 + ngenome
	enddo

cncpeak=ncpeak*2

C---calculate distance between IN peaks
	do ip1=1,npeak
	do ip2=1,npeak
	rx1=x(peaks(ip1))
	ry1=y(peaks(ip1))
	rz1=z(peaks(ip1))
	rx2=x(peaks(ip2))
	ry2=y(peaks(ip2))
	rz2=z(peaks(ip2))

	peakdist(ip1,ip2)=(rx1-rx2)**2+(ry1-ry2)**2+(rz1-rz2)**2

	enddo
	enddo
	 
	write(6,*) "finished distance calculation"

C---create ordered IN peak list
	do ip1=1,npeak

	do iptest=1,npeak
	peakcheck(iptest)=0
	enddo

	do icount=1,npeak-1
	rmin=9999.

	do ip2=1,npeak
	if ((peakdist(ip1,ip2).le.rmin).and.
     &      (peakcheck(ip2).eq.0).and.(ip1.ne.ip2)) then
	   ipmin=ip2
	   rmin=peakdist(ip1,ip2)
	endif
	enddo

	peaklist(ip1,icount)=ipmin
	peakcheck(ipmin)=1
	enddo

	enddo
	
	write(6,*) "finished ordering of peak distance list"

C---place IN
	write(6,*) "picking IN less than ",rdmax
	rmaxall=9999.

	do itry=1,1000

	do ip=1,npeak
	tmpflag(ip)=0.
	enddo

	rmax=-9999.
	do ic=1,nIN

 22	r=rand()
	ip1=int(r*rnpeak)+1
	if (tmpflag(ip1).ne.0) goto 22
	tmpflag(ip1)=1

	do icount=1,npeak-1
	 ip2=peaklist(ip1,icount)
	 if (tmpflag(ip2).eq.0) then
	        tmpflag(ip2)=1
	        rmax=max(rmax,sqrt(peakdist(ip1,ip2)))
	        goto 23
	 endif
	enddo
 23	continue

	if (ic.le.nTET) then
	do icount=1,npeak-1
	 ip3=peaklist(ip1,icount)
	 if (tmpflag(ip3).eq.0) then
	        tmpflag(ip3)=1
	        rmax=max(rmax,sqrt(peakdist(ip1,ip3)))
	        goto 24
	 endif
	enddo
 24	continue
	do icount=1,npeak-1
	 ip4=peaklist(ip1,icount)
	 if (tmpflag(ip4).eq.0) then
	        tmpflag(ip4)=1
	        rmax=max(rmax,sqrt(peakdist(ip1,ip4)))
	        goto 25
	 endif
	enddo
 25	continue
	endif

	ip1save(ic)=ip1
	ip2save(ic)=ip2
	if (ic.le.nTET) then
	ip3save(ic)=ip3
	ip4save(ic)=ip4
	endif

C ** end of ic loop
	enddo

	if (rmax.lt.rmaxall) then
	do ic=1,nIN
	rmaxall=rmax
	ip1max(ic)=ip1save(ic)
	ip2max(ic)=ip2save(ic)
	if (ic.le.nTET) then
	ip3max(ic)=ip3save(ic)
	ip4max(ic)=ip4save(ic)
	endif
	enddo
	endif

C ** end of itry loop
	enddo

	do ic=1,nIN

	ip1=ip1max(ic)
	ip2=ip2max(ic)
	const(ic,1)=peaks(ip1)
	const(ic,2)=peaks(ip2)
	peakflag(ip1)=1
	peakflag(ip2)=1
         write(6,*) "IN dimer",ic,ip1,ip2,sqrt(peakdist(ip1,ip2))

	if (ic.le.nTET) then
	ip3=ip3max(ic)
	ip4=ip4max(ic)
	const(ic,3)=peaks(ip3)
	const(ic,4)=peaks(ip4)
	peakflag(ip3)=1
	peakflag(ip4)=1

         write(6,*) "IN tetra",ic,ip1,ip3,sqrt(peakdist(ip1,ip3))
         write(6,*) "IN tetra",ic,ip1,ip4,sqrt(peakdist(ip1,ip4))
	endif
	write(6,*)
	enddo
	write(6,*) "minimum rmax ",rmaxall
	   
C ** write NC constraints
	do iseq=1,3100
	seqflag(iseq)=0
	enddo

	iNCcount=0
	nNCtoplace=min(nNC,ncpeak)
	write(6,*) "number of specific NC to place ",nNCtoplace

	do inc=1,npeak
	seqflag(peaks(inc)-1)=1
	seqflag(peaks(inc))=1
	seqflag(peaks(inc)+1)=1
	enddo

	do inc=1,nNCtoplace
	if (seqflag(ncpeaks(inc)).eq.0) then
	  seqflag(ncpeaks(inc)-1)=1
	  seqflag(ncpeaks(inc))=2
	  seqflag(ncpeaks(inc)+1)=1
	  iNCcount=iNCcount+1
	endif
	enddo
	write(6,*) "specific NC added: ",iNCcount

C ** add randomly placed NC

	nrand=nNC-iNCcount
	icount=0
	rgenome=float(ngenome)
	do inc=1,nrand
	do itry=1,1000
C ** 120 is to keep random NC out of 5UTR
	rseq=rand()*(rgenome-124.)+120.
	iseq=int(rseq)
	if (seqflag(iseq).eq.0) then
	 icount=icount+1
	 seqflag(iseq-1)=1
	 seqflag(iseq)=3
	 seqflag(iseq+1)=1
	 goto 222
	endif
	enddo
 222	continue
	enddo
	write(6,*) "number of random NC placed ",icount
	write(6,*) ngenome
	do iseq=1,ngenome
	if (seqflag(iseq).eq.2) write(10,223) iseq,"PEAK"
	if (seqflag(iseq).eq.3) write(10,223) iseq,"RAND"
 223	format(i8,2x,a4)
	enddo

C ** write IN constraints
	do ic=1,nIN
	imol=ic+iNCcount+icount
	if (ic.gt.nTET) then
	write(3,77) imol,const(ic,1),const(ic,2),
     &              6.5,3.25,2,"IN dimer   "
	endif

	if (ic.le.nTET) then
	write(3,77) imol,const(ic,1),const(ic,3),
     &              9.1,4.55,2,"IN diagonal"
	write(3,77) imol,const(ic,2),const(ic,4),
     &              9.1,4.55,2,"IN diagonal"
	write(3,77) imol,const(ic,1),const(ic,2),
     &              6.5,0.00,0,"IN edge    "
	write(3,77) imol,const(ic,2),const(ic,3),
     &              6.5,0.00,0,"IN edge    "
	write(3,77) imol,const(ic,3),const(ic,4),
     &              6.5,0.00,0,"IN edge    "
	write(3,77) imol,const(ic,4),const(ic,1),
     &              6.5,0.00,0,"IN edge    "
	endif

 77	format(3i8,2f8.3,1x,i1,1x,a11)
	enddo

	end
