C--------------------------------------------------------------------
C                 ******************************
C                     I L L U S T R A T O R
C                   Biomolecular Illustration
C		  ******************************
C		        David S Goodsell
C--------------------------------------------------------------------
C
C		no anti-aliasing of edges
C	July 19 2007--conical shadows, smooth outlines with kernels,
C                     fixed  ppm output
C
C--------------------------------------------------------------------
		integer*4 scanline(9000)
c ***** FRAME BUFFER WITH COLORS *****
		real*4 pix(-10:3008,-10:3008,4)
c ***** FRAME BUFFER CONTAINING THE z HEIGHT *****
		real*4 zpix(-10:3008,-10:3008)
c ***** FRAME BUFFER CONTAINING THE ATOM NUMBER *****
		integer*4 atom(-10:3008,-10:3008)
		integer*4 n,ia,iatom
c ***** SHADING MAPS FOR THE EIGHT ATOM TYPES *****
		real*4 sphdat(0:32000,3)
		integer*4 numpix
c ***** PHONG SHADING PARAMETERS ****** 
		real*4 colortype(0:16,3)
		real*4 rfog(3)
		integer*4 iback(3),ifog(3)
c ***** ATOMIC INFORMATION *****
		real*4 coord(3,350000)
		integer*4 type(350000),res(0:350000),su(0:350000)
		real*4 radtype(16)
c ***** transformation matrices *****
		real*4 matrixin(4,4),rm(4,4)
c ***** STUFF FOR OUTLINES *****
		real*4 l_opacity,l_opacity_ave,g_opacity,opacity
		real*4 l(-1:1,-1:1),g
		real*4 l_low,l_high,g_low,g_high,l_low16,l_high16
		real*4 l_diff_min, l_diff_max
c ***** Conical shadows *****
		real*4 rtable(-51:51,-51:51),coneangle,pcone,rcone
c ***** ETC. *****
		real*4 x,y,z,rx,ry,rz,xp,yp,zp,d
		real*4 xn,yn,zn,ci,xs,xy,xz,cs
		character*3 comcode(14),command
		character*20 filename,inputfile
		integer*4 ixsize,iysize, idepth
		integer*4 ix,iy,iz
		integer illustrationflag
c ***** stuff for reading atoms
		integer*4 resrange(2,100),inptype(100),inpsu(100)
		character*6 atomdescriptor(100)
		character*10 descriptor(100)
		character*80 instring
c--------------------------------------------------------------------
		pi=3.141592
		rscale=1.

		illustrationflag=0

		data comcode/'rea','tra','xro','yro','zro','sca','cen',
     &		             'wor','map','cal','bac','ill','bon','mod'/
		psl=1.
		psh=1.
		lz=1.
		call clearmatrix(rm)
		l_low=1.
		l_high=10.
		g_low=350000.
		g_high=21000.
		xshmax=4000
		yshmax=4000
		su(0)=9999
		res(0)=9999
		xtran=0.
		ytran=0.
		ztran=0.
		l_diff_max=50.
		l_diff_min=1.
		ixsize=0
		iysize=0
		imodel=0
c ********************* READ CONTROL CARDS ***************************
 10		read (5,101,end=999) command
		icommand=10
		do icount=1,14
 100		if (command.eq.comcode(icount)) icommand=icount
		enddo
  		goto (1,2,3,4,5,6,7,8,10,111,11,12,13,14),icommand
		write (6,102) ' ***** invalid control card read: ',
     &  command, ' ***** '
		goto 10
 101		format (a3)
 102		format (a35,a3,a7)
c--------------------------------------------------------------------
C  *** read and classify atoms ***
c
 1	continue
		imodelcurrent=0
		read(5,113) inputfile
		open(1,file=inputfile,form='formatted',status='old')
c --- read atom descriptors ---
		ndes=0
C	default 50% gray
	do i=1,16
	do j=1,3
	colortype(i,j)=.5
	enddo
	enddo
 7020		read(5,7100) instring
 		if (instring(1:3).eq.'END') goto 7030
 		  ndes=ndes+1
 		  read(instring,21) atomdescriptor(ndes),descriptor(ndes)
 		  read(instring(18:80),*) (resrange(i,ndes),i=1,2),
     &                      inptype(ndes),inpsu(ndes),rr,rg,rb,rad
	           colortype(inptype(ndes),1)=rr
	           colortype(inptype(ndes),2)=rg
	           colortype(inptype(ndes),3)=rb
		   radtype(inptype(ndes))=rad
		goto 7020
 7030		write(6,*) ' atom descriptors: ',ndes
	do i=1,16
	write(6,*) "type, color ",i,(colortype(i,j),j=1,3)
	enddo
 21		format(a6,a10)
c --- read atoms and classify ---
		n=0
 7040		read(1,7100,end=7009) instring
		if (instring(1:5).eq.'MODEL') then
		  imodelcurrent=imodelcurrent+imodel
		  goto 7040
		endif
		if ((instring(1:4).ne.'ATOM').and.
     &      (instring(1:6).ne.'HETATM')) goto 7040
     	read(instring,200) ires
		do ides=1,ndes

 		 if (instring(1:6).ne.atomdescriptor(ides)(1:6)) goto 7050

 		 do ia=1,10
 		  if (descriptor(ides)(ia:ia).eq.'-') goto 7060
 		  if (instring(12+ia:12+ia).ne.descriptor(ides)(ia:ia)) goto 7050
 7060 	         continue
		enddo

 		 if ((ires.lt.resrange(1,ides)).or.
     &        (ires.gt.resrange(2,ides))) goto 7050

		 if (inptype(ides).eq.0) goto 7040

 		 n=n+1
		 read(instring,300) (coord(i,n),i=1,3)
		 type(n)=inptype(ides)
		 su(n)=inpsu(ides)+imodelcurrent
		 res(n)=ires
		 goto 7040
 7050 	continue
		enddo
		goto 7040
c --- done ---         
 7009		write(6,*)' atoms read: ', n, ' from: ',inputfile
		write(6,*)' '
c	 do j=1,n
c	 write(99,*) (coord(i,j),i=1,3),type(j),su(j),res(j)
c	 enddo
 7100	format(a80)
 200	format(22x,i4)
 300	format(30x,3f8.3)
		goto 10
c--------------------------------------------------------------------
 14		read(5,*) imodel
		write(6,*) "model increment", imodel
		goto 10
c--------------------------------------------------------------------
c       TRANSLATION 
 2		continue
		read(5,*) xtrani,ytrani,ztrani
		xtran=xtran+xtrani
		ytran=ytran+ytrani
		ztran=ztran+ztrani
 		write(6,*) 'translation vector : ',xtran,ytran,ztran
		write(6,*)
		goto 10
c--------------------------------------------------------------------
c		Z ROTATION
 5		call clearmatrix(matrixin)
		read(5,*) angle
		write(6,*) 'z rotation : ',angle
		angle=angle*3.141592/180.
		matrixin(1,1)=cos(angle)
		matrixin(1,2)=-sin(angle)
		matrixin(2,1)=sin(angle)
		matrixin(2,2)=cos(angle)
		call catenate(rm,matrixin)
		goto 10
c--------------------------------------------------------------------
c		Y ROTATION
 4		call clearmatrix(matrixin)
		read(5,*) angle
		write(6,*) 'y rotation : ',angle
		angle=angle*3.141592/180.
		matrixin(1,1)=cos(angle)
		matrixin(1,3)=sin(angle)
		matrixin(3,1)=-sin(angle)
		matrixin(3,3)=cos(angle)
		call catenate(rm,matrixin)
		goto 10
c--------------------------------------------------------------------
c		X ROTATION
 3		call clearmatrix(matrixin)
		read(5,*) angle
		write(6,*) 'x rotation : ',angle
		angle=angle*3.141592/180.
		matrixin(2,2)=cos(angle)
		matrixin(2,3)=-sin(angle)
		matrixin(3,2)=sin(angle)
		matrixin(3,3)=cos(angle)
		call catenate(rm,matrixin)
		goto 10
c--------------------------------------------------------------------
c       SCALE 
 6		continue
		read(5,*) rscalei
		rscale=rscale*rscalei
		write(6,*) 'scale factor : ',rscale
		write(6,*)
		goto 10
c--------------------------------------------------------------------
c		CENTERING of coordinates
c   three options: CEN(TER-OF-MASS) will place c-o-m on (0,0,0)
c                  EXT(ERNAL) will place input coordinates on (0,0,0)
c                  AUT(O-CENTER) will center the rotated coordinates
c                        in the frame, and place the uppermost atom
c                        just below the image plane
c    explicit TRANSLATIONs are applied after the centering
 7		continue
 		read(5,101) command
 		autocenter=0
 		if (command.eq.'aut') then
		   autocenter=1
		   goto 10
		endif

		if (command.eq.'ext') then
		write(6,*) 'reading external centering vector'
		read(5,*) xtranc,ytranc,ztranc
		write(6,*) 'external centering vector: ',xtranc,ytranc,ztranc
		endif

		if (command.eq.'cen') then
		xmin=10000.
		xmax=-10000.
 		ymin=10000.
		ymax=-10000.
 		zmin=10000.
		zmax=-10000.
		write(6,*) 'centering on center of mass'

		do ia=1,n
		xmin=MIN(xmin,coord(1,ia))
		xmax=MAX(xmax,coord(1,ia))
		ymin=MIN(ymin,coord(2,ia))
		ymax=MAX(ymax,coord(2,ia))
		zmin=MIN(zmin,coord(3,ia))
		zmax=MAX(zmax,coord(3,ia))
		enddo

		write(6,*) 'min coordinates : ',xmin,ymin,zmin
		write(6,*) 'max coordinates : ',xmax,ymax,zmax
		xtranc=-xmin-(xmax-xmin)/2.
		ytranc=-ymin-(ymax-ymin)/2.
		ztranc=-zmin-(zmax-zmin)/2.
		write(6,*) 'centering vector : ',xtranc,ytranc,ztranc
		endif

		do ia=1,n
		 coord(1,ia)=coord(1,ia)+xtranc
		 coord(2,ia)=coord(2,ia)+ytranc
 		 coord(3,ia)=coord(3,ia)+ztranc
		 if (su(ia).eq.0) coord(3,ia)=coord(3,ia)-(zmax-zmin)
		enddo

		goto 10
c--------------------------------------------------------------------
C read parameters that describe the environment
 8		continue
c
		read (5,*) (iback(i),i=1,3),(ifog(i),i=1,3),pfogh,pfogl
		if (pfogl.lt.0.) then
		   pfogl=-pfogl
		endif

		do i=1,3
		if (iback(i).gt.255) iback(i)=255
		if (ifog(i).gt.255) ifog(i)=255
 230		rfog(i)=float(ifog(i))/256.
	        colortype(0,i)=float(iback(i))/255.
		enddo

		if (pfogh.gt.1.) pfogh=1.
		if (pfogl.gt.1.) pfogl=1.

		write (6,204) ' background inten. :',(iback(i),i=1,3)
		write (6,204) ' fog intensity :    ',(ifog(i),i=1,3)
		write (6,202) ' upper fog percent :',pfogh*100.
		write (6,202) ' lower fog percent :',pfogl*100.
		pfogdiff=pfogh-pfogl
		read (5,*) icone,pcone,coneangle,rcone,pshadowmax
		if (icone.ne.0) write (6,202) ' draw conical shadows'
		read(5,*) ixsize,iysize
		if (ixsize.gt.3000) ixsize=3000
		if (iysize.gt.3000) iysize=3000
		write(6,*) 'input value for image size', ixsize,iysize
		goto 10
 201		format(4f12.5)
 202		format(1x,a20,1x,8f7.2)
 203    format(3(2f8.5,i8,f8.5),4f8.2)
 204    format(1x,a20,1x,8i7)
 205		format(6i8,2f10.4)
 206		format(4i8)
 207    format(/1x,a20/)
 208		format(8f10.5)
 209		format(f10.5,i8)
c--------------------------------------------------------------------
 11		continue
c		***** read backgrounds *****
C obsolete
		goto 10
c--------------------------------------------------------------------
 12		continue
c		***** read illustration parameters *****
		illustrationflag=1
		read(5,*) l_low,l_high,l_low16,l_high16,ikernel,
     &                    l_diff_min,l_diff_max
		read(5,*) r_low,r_high
		read(5,*) g_low,g_high,resdiff
		write(6,*) 'illustration parameters'
		write(6,*) 'l parameters: ',l_low,l_high
		write(6,*) 'g parameters: ',g_low,g_high
		goto 10
c--------------------------------------------------------------------
 13		continue
c		***** read bond parameters *****
C obsolete
		goto 10
c--------------------------------------------------------------------
 111    write (6,207) ' *begin calculation*'
c ***** Populate conical shadow table ****
		conemax=50.
		do i=-51,51
		do j=-51,51
		 rtable(i,j)=sqrt(float(i)**2+float(j)**2)
		 if (rtable(i,j).gt.conemax) rtable(i,j)=10000.
		enddo
		enddo
		 rtable(0,0)=10000.
c ***** SCALE RADII *****
		do i=1,16
		 radtype(i)=radtype(i)*rscale
		 if (radtype(i).gt.radius_max) radius_max=radtype(i)
 119		enddo
c ***** APPLY THE ROTATION MATRIX *****
 		do ia=1,n
		rx=coord(1,ia)*rm(1,1)+coord(2,ia)*rm(2,1)+
     &		   coord(3,ia)*rm(3,1)
		ry=coord(1,ia)*rm(1,2)+coord(2,ia)*rm(2,2)+
     &		   coord(3,ia)*rm(3,2)
		rz=coord(1,ia)*rm(1,3)+coord(2,ia)*rm(2,3)+
     &		   coord(3,ia)*rm(3,3)
		coord(1,ia)=rx
		coord(2,ia)=ry
		coord(3,ia)=rz
 150		continue
		enddo
c ***** APPLY AUTOCENTERING and AUTOCENTERING, if switched on *****
		if (autocenter.eq.1) then

		xmin=10000.
		xmax=-10000.
 		ymin=10000.
		ymax=-10000.
 		zmin=10000.
		zmax=-10000.
		write(6,*) 'automating centering'

		do ia=1,n
		 xmin=MIN(xmin,coord(1,ia))
		 xmax=MAX(xmax,coord(1,ia))
		 ymin=MIN(ymin,coord(2,ia))
		 ymax=MAX(ymax,coord(2,ia))
		 zmin=MIN(zmin,coord(3,ia))
		 zmax=MAX(zmax,coord(3,ia))
 117		continue
		enddo
		write(6,*) 'min coordinates : ',xmin,ymin,zmin
		write(6,*) 'max coordinates : ',xmax,ymax,zmax
		xtranc=-xmin-(xmax-xmin)/2.
		ytranc=-ymin-(ymax-ymin)/2.
		ztranc=-zmax-radius_max-1.
		write(6,*) 'centering vector : ',xtranc,ytranc,ztranc
c	if ((ixsize.le.0).or.(iysize.le.0)) then
		  write(6,*)
		  write(6,*) 'applying autosizing'
		  write(6,*) 'x and y frame width: ',-ixsize,-iysize
		  ixsize=-2.*ixsize+2.*radius_max+(xmax-xmin)*rscale
		  iysize=-2.*iysize+2.*radius_max+(ymax-ymin)*rscale
		  ixsize=min(ixsize,3000)
		  iysize=min(iysize,3000)
c	endif

		do ia=1,n
		 coord(1,ia)=coord(1,ia)+xtranc
		 coord(2,ia)=coord(2,ia)+ytranc
		 coord(3,ia)=coord(3,ia)+ztranc
		 if (su(ia).eq.0) coord(3,ia)=coord(3,ia)-(zmax-zmin)
 116		enddo

 		endif
c ***** OPEN OUTPUT FILES *****
		ixsize=int(ixsize/2)*2
		iysize=int(iysize/2)*2
		  write(6,*) 'xsize and ysize: ',ixsize,iysize
		  write(6,*)
		 read(5,113) filename
		 write(6,*) "output pnm filename: ",filename
		 open(8,file=filename,form='formatted')
		 write(8,1003) "P3"
 1003		 format(a2)
		 write(8,1004) iysize,ixsize
		 write(8,1004) 255
		 open(9,file="opacity.pnm",form='formatted')
		 write(9,1003) "P3"
		 write(9,1004) iysize,ixsize
		 write(9,1004) 255
 1004	format(2i5)
 113		format(a20)
c ***** APPLY TRANSLATION AND SCALING *****
		do ia=1,n
		 coord(1,ia)=(coord(1,ia)+xtran)*rscale
		 coord(2,ia)=(coord(2,ia)+ytran)*rscale
 		 coord(3,ia)=(coord(3,ia)+ztran)*rscale
 114		enddo
c ***** MAP SPHERICAL SURFACES OVER ATOMS *****
		do ix=1,ixsize
		do iy=1,iysize
 		  pix(ix,iy,1)=0.
 		  pix(ix,iy,2)=0.
 		  pix(ix,iy,3)=0.
 		  pix(ix,iy,4)=0.
		  atom(ix,iy)=0
 		  zpix(ix,iy)=-10000.
		enddo
		enddo

		if (n.gt.0) then

C ----- create the spherical shading map for atom types -----
		do irad=1,16
		ic=1
		irlim=int(radtype(irad))
		if (irlim.gt.100) then
		  write(6,*) 'atoms radius * scale > 100'
		  stop
		endif

		do ix=-irlim-1,irlim+1
		do iy=-irlim-1,irlim+1
		x=float(ix)
		y=float(iy)
		d=sqrt(x*x+y*y)
		if (d.gt.radtype(irad)) goto 350
		z=sqrt(radtype(irad)**2-d*d)
		sphdat(ic,1)=x
		sphdat(ic,2)=y
		sphdat(ic,3)=z
		ic=ic+1
 350		continue
		enddo
 352		enddo
		numpix=ic-1
c ----- then map spherical surface over atoms of the proper type ------
		icount=0
		do iatom=1,n
		if (type(iatom).ne.irad) goto 500
		icount=icount+1

		if (coord(3,iatom).lt.0.) then

		do ipix=1,numpix
		x=sphdat(ipix,1)+coord(1,iatom)+float(ixsize)/2.		
		y=sphdat(ipix,2)+coord(2,iatom)+float(iysize)/2.
		ix=int(x)
		iy=int(y)
		if ((x.gt.float(ixsize)).or.(x.lt.1.).or.
     &              (y.gt.float(iysize)).or.(y.lt.1.)) goto 510
		z=sphdat(ipix,3)+coord(3,iatom)		
		if (z.gt.zpix(ix,iy)) then
		 zpix(ix,iy)=z
		 atom(ix,iy)=iatom
		endif
 510		continue
		enddo

		endif

 500		continue
		enddo
		write(6,*) icount,' spheres added of type: ',irad

 340		enddo
		write(6,*) 'shading maps written into depth buffer'

 	endif

c***** CALCULATE SECOND DERIVATIVE OUTLINES ******
c ---- find maximum and minimum z levels ---
c		(note: I use a value of zpix=-10000. to distinguish background)
		zpix_max=-100000.
		zpix_min=100000.
		do ix=1,ixsize
		do iy=1,iysize
		if (zpix(ix,iy).gt.zpix_max) zpix_max=zpix(ix,iy)
		if ((zpix(ix,iy).ne.-10000.).and.
     &		    (zpix(ix,iy).lt.zpix_min)) zpix_min=zpix(ix,iy)
		enddo
		enddo
		zpix_max=min(zpix_max,0.)
		zpix_spread=zpix_max-zpix_min
		write(6,*) 'zpix_min,zpix_max ',zpix_min,zpix_max

c ***** PROCESSING OF THE IMAGE BEGINS HERE*****
		l_diff_min=l_diff_min*rscale
		l_diff_max=l_diff_max*rscale
		write(6,*) ' Pixel processing beginning '
		do ix=1,ixsize
		do iy=1,iysize
		zpix(ix,iy)=min(zpix(ix,iy),0.)
 2009		enddo
 1009		enddo

		do ix=1,ixsize
		do iy=1,iysize

c ***** SHADOW TESTING *****
		psh=1.
c		--- calculate conical shadows  ---
		pconetot=1.
		    if ((icone.ne.0).and.(atom(ix,iy).ne.0)) then	
		do i=-50,50,5
		do j=-50,50,5
		  if ((ix+i.gt.0).and.(ix+i.lt.ixsize).and.
     &                (iy+j.gt.0).and.(iy+j.lt.iysize)) then
		  rzdiff=zpix(ix+i,iy+j)-zpix(ix,iy)
		if (rzdiff.gt.rcone) then
		 if (rtable(i,j)*coneangle.lt.rzdiff+rcone) then
		pconetot=pconetot-pcone
		 endif
		 endif
		 endif
		enddo
		enddo
		pconetot=max(pconetot,pshadowmax)
		endif
c   		if (type(atom(ix,iy)).eq.16) psh=1.
c ***** CALCULATE THE FOG PERCENTAGE (DEPTH CUEING) *****
		pfh=pfogh-(zpix_max-zpix(ix,iy))/zpix_spread*pfogdiff
    		if (zpix(ix,iy).lt.zpix_min) pfh=1.
    		if (type(atom(ix,iy)).eq.16) pfh=1.
c ***** CALCULATE OUTLINES *****
		g_opacity=0.
		l_opacity=0.
		if (illustrationflag.ne.0) then
c ***** SUBUNIT OUTLINES *****
c ---- (this replaces original calculation of first derivatives)

		if ((ix.gt.1).and.(ix.lt.ixsize).and.
     &      (iy.gt.1).and.(iy.lt.iysize)) then
c		--- calculate subunit 'derivatives'  ---
		g=0.
		r=0.
		do i=-2,2
		do j=-2,2
		 if (abs(i*j).ne.4) then
		  if (su(atom(ix,iy)).ne.su(atom(ix+i,iy+j))) r=r+1.
                  if (abs(res(atom(ix,iy))-
     &                res(atom(ix+i,iy+j))).gt.resdiff) g=g+1.
		 endif
		enddo
		enddo
c		--- opacities are 1 for completely opaque ---
		g_opacity=min((g-g_low)/(g_high-g_low),1.)
		r_opacity=min((r-r_low)/(r_high-r_low),1.)
		g_opacity=max(g_opacity,r_opacity,0.)
		endif
c ***** SECOND DERIVATIVE OUTLINES *****
		if ((ix.gt.2).and.(ix.lt.ixsize-1).and.
     &      (iy.gt.2).and.(iy.lt.iysize-1)) then
		rl=0.
		l_opacity_ave=0.
		do ixl=-1,1
		do iyl=-1,1
		ixc=ix+ixl
		iyc=iy+iyl
	if (ikernel.eq.1) then
 	 l(ixl,iyl)=abs(1./3. * ( 
     &	 -0.8*zpix(ixc-1,iyc-1)-1.*zpix(ixc-1,iyc)-0.8*zpix(ixc-1,iyc+1)-
     &	 1.0* zpix(ixc,iyc-1)+7.2*zpix(ixc,iyc)-1.0*zpix(ixc,iyc+1)-
     &	 0.8* zpix(ixc+1,iyc-1)-1.*zpix(ixc+1,iyc)-0.8*zpix(ixc+1,iyc+1)
     &     ))
	endif
	if (ikernel.eq.2) then
		 l(ixl,iyl)=abs(1./3. * ( 
     &	-0.8*zpix(ixc-1,iyc-1)-1.0*zpix(ixc-1,iyc)-0.8*zpix(ixc-1,iyc+1)-
     &	 1.0*zpix(ixc,iyc-1)+8.8*zpix(ixc,iyc)-1.0*zpix(ixc,iyc+1)-
     &	 0.8*zpix(ixc+1,iyc-1)-1.0*zpix(ixc+1,iyc)-0.8*zpix(ixc+1,iyc+1)-
     &	 0.1*zpix(ixc+2,iyc-1)-0.2*zpix(ixc+2,iyc)-0.1*zpix(ixc+2,iyc+1)-
     &	 0.1*zpix(ixc-2,iyc-1)-0.2*zpix(ixc-2,iyc)-0.1*zpix(ixc-2,iyc+1)-
     &	 0.1*zpix(ixc-1,iyc+2)-0.2*zpix(ixc,iyc+2)-0.1*zpix(ixc+1,iyc+2)-
     &	 0.1*zpix(ixc-1,iyc-2)-0.2*zpix(ixc,iyc-2)-0.1*zpix(ixc+1,iyc-2)
     &     ))
	endif
	if (ikernel.eq.3) then
		do i=-1,1
		do j=-1,1
			rd=abs(zpix(ix,iy)-zpix(ix+i,iy+j))
			if (rd.gt.l_diff_min) then
			rd=(rd-l_diff_min)/(l_diff_max-l_diff_min)
			l(ixl,iyl)=l(ixl,iyl)+min(rd,1.)
			endif
		enddo
		enddo
	endif
	if (ikernel.eq.4) then
		do i=-2,2
		do j=-2,2
		 if (abs(i*j).ne.4) then
			rd=abs(zpix(ix,iy)-zpix(ix+i,iy+j))
			if (rd.gt.l_diff_min) then
			rd=abs(zpix(ix,iy)-zpix(ix+i,iy+j))
			rd=(rd-l_diff_min)/(l_diff_max-l_diff_min)
			l(ixl,iyl)=l(ixl,iyl)+min(rd,1.)
			endif
		 endif
		enddo
		enddo
	endif
		l(ixl,iyl)=min((l(ixl,iyl)-l_low)/(l_high-l_low),1.)
		l(ixl,iyl)=max(l(ixl,iyl),0.)
		if (l(ixl,iyl).gt.0.) rl=rl+1.
		l_opacity_ave=l_opacity_ave+l(ixl,iyl)
		enddo
		enddo
		if (rl.ge.6.) then
		l_opacity=l_opacity_ave/6.
		else
		l_opacity=l(0,0)
		endif
		l_opacity=min(l_opacity,1.)
		l_opacity=max(l_opacity,0.)
		endif
c ---- combine subunit outlines and derivative outlines ----
		l_opacity=max(l_opacity,g_opacity)
		endif
c ***** CALCULATE THE TOTAL PIXEL INTENSITY *****

	psh=min(psh,pconetot)
	ropacity=0.
	do icolor=1,3

	rcolor=
     &     pfh*(pconetot*(colortype(type(atom(ix,iy)),icolor)))+
     &     (1.-pfh)*rfog(icolor) 

	pix(ix,iy,icolor)=(1.-l_opacity)*rcolor
c ----calculate pixel opacity
	if (type(atom(ix,iy)).ne.0) ropacity=1.
	pix(ix,iy,4)=max(ropacity,l_opacity)

	enddo

 	enddo

c ***** output of a scan line *****
c ----- PPM format -----
	iscan=0
	do iout=1,iysize
	do ic=1,3
	iscan=iscan+1
	scanline(iscan)=int(pix(ix,iout,ic)*255.)
	scanline(iscan)=min(scanline(iscan),255)
	scanline(iscan)=max(scanline(iscan),0)
	enddo
	enddo
	if (iras.ne.0) write(8,1002) (scanline(if),if=1,iysize*3)
C -- write opacity
	iscan=0
	do iout=1,iysize
	do ic=1,3
	iscan=iscan+1
	scanline(iscan)=int(pix(ix,iout,4)*255.)
	scanline(iscan)=min(scanline(iscan),255)
	scanline(iscan)=max(scanline(iscan),0)
	enddo
	enddo
	if (iras.ne.0) write(9,1002) (scanline(if),if=1,iysize*3)
 1002	format(20i4)
c ----- diagnostic ------
		if (int(ix/20)*20.eq.int((float(ix)/20.)*20.)) then
		 write(6,6669) (int(pix(ix,iyo,1)*9.),iyo=1,iysize,20)
 6669	format(65i1)
		endif

 1000		enddo

 999		stop
		end
c--------------------------------------------------------------------
		subroutine catenate(m1,m2)
		real*4 m1(4,4),m2(4,4),m(4,4)
c--------------------------------------------------------------------
		m(1,1)=m1(1,1)*m2(1,1)+m1(2,1)*m2(1,2)+m1(3,1)*m2(1,3)+
     &		       m1(4,1)*m2(1,4)
		m(1,2)=m1(1,2)*m2(1,1)+m1(2,2)*m2(1,2)+m1(3,2)*m2(1,3)+
     &		       m1(4,2)*m2(1,4)
		m(1,3)=m1(1,3)*m2(1,1)+m1(2,3)*m2(1,2)+m1(3,3)*m2(1,3)+
     &		       m1(4,3)*m2(1,4)
		m(1,4)=m1(1,4)*m2(1,1)+m1(2,4)*m2(1,2)+m1(3,4)*m2(1,3)+
     &		       m1(4,4)*m2(1,4)
		m(2,1)=m1(1,1)*m2(2,1)+m1(2,1)*m2(2,2)+m1(3,1)*m2(2,3)+
     &		       m1(4,1)*m2(2,4)
		m(2,2)=m1(1,2)*m2(2,1)+m1(2,2)*m2(2,2)+m1(3,2)*m2(2,3)+
     &		       m1(4,2)*m2(2,4)
		m(2,3)=m1(1,3)*m2(2,1)+m1(2,3)*m2(2,2)+m1(3,3)*m2(2,3)+
     &		       m1(4,3)*m2(2,4)
		m(2,4)=m1(1,4)*m2(2,1)+m1(2,4)*m2(2,2)+m1(3,4)*m2(2,3)+
     &		       m1(4,4)*m2(2,4)
		m(3,1)=m1(1,1)*m2(3,1)+m1(2,1)*m2(3,2)+m1(3,1)*m2(3,3)+
     &		       m1(4,1)*m2(3,4)
		m(3,2)=m1(1,2)*m2(3,1)+m1(2,2)*m2(3,2)+m1(3,2)*m2(3,3)+
     &		       m1(4,2)*m2(3,4)
		m(3,3)=m1(1,3)*m2(3,1)+m1(2,3)*m2(3,2)+m1(3,3)*m2(3,3)+
     &		       m1(4,3)*m2(3,4)
		m(3,4)=m1(1,4)*m2(3,1)+m1(2,4)*m2(3,2)+m1(3,4)*m2(3,3)+
     &		       m1(4,4)*m2(3,4)
		m(4,1)=m1(1,1)*m2(4,1)+m1(2,1)*m2(4,2)+m1(3,1)*m2(4,3)+
     &		       m1(4,1)*m2(4,4)
		m(4,2)=m1(1,2)*m2(4,1)+m1(2,2)*m2(4,2)+m1(3,2)*m2(4,3)+
     &		       m1(4,2)*m2(4,4)
		m(4,3)=m1(1,3)*m2(4,1)+m1(2,3)*m2(4,2)+m1(3,3)*m2(4,3)+
     &		       m1(4,3)*m2(4,4)
		m(4,4)=m1(1,4)*m2(4,1)+m1(2,4)*m2(4,2)+m1(3,4)*m2(4,3)+
     &		       m1(4,4)*m2(4,4)
		do j=1,4
		do i=1,4
 100		m1(i,j)=m(i,j)
		enddo
		enddo
		return
		end
c--------------------------------------------------------------------
		subroutine clearmatrix(m)
		real*4 m(4,4)
c--------------------------------------------------------------------
		do i=1,4
		do j=1,4
		x=0.
		if (i.eq.j) x=1.
 100		m(i,j)=x
		enddo
		enddo
		return
		end
