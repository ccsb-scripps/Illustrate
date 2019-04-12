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
C                     fixed bond intersections, ppm output
C
C--------------------------------------------------------------------
		real*4 rframe(3,3000)
		integer*4 scanline(9000)
c ***** FRAME BUFFER WITH COLORS *****
		real*4 pix(-10:3008,-10:3008,3)
c ***** FRAME BUFFER CONTAINING THE z HEIGHT *****
		real*4 zpix(-10:3008,-10:3008)
c ***** FRAME BUFFER CONTAINING THE ATOM TYPE *****
		integer*4 atom(-10:3008,-10:3008)
		integer*4 n,ia,iatom,nclipped
c ***** FRAME BUFFERS FOR SHADOWS *****
		real*4 zshadow(4000,4000)
		integer*4 shadowatom(4000,4000)
c ***** SHADING MAPS FOR THE EIGHT ATOM TYPES *****
		real*4 sphdat(0:32000,4)
		real*4 colormap(-102:102,-102:102,3)
		integer*4 numpix,numshadowpix
c ***** PHONG SHADING PARAMETERS ****** 
		real*4 rp(16,3),wp(16,3),dif(16,3)
		integer*4 isp(16,3)
		real*4 rback(3),rbacki(3), rfog(3)
		integer*4 iback(3),ifog(3)
c ***** ATOMIC INFORMATION *****
		real*4 coord(3,350000)
		integer*4 type(350000),res(0:350000),su(0:350000)
		real*4 radtype(16)
		real*4 shclip(16,3)
c ***** transformation matrices *****
		real*4 matrixin(4,4),rm(4,4)
c ***** STUFF FOR OUTLINES *****
	real*4 l_opacity,l_opacity_ave,g_opacity,opacity
		real*4 l(-1:1,-1:1),g
	real*4 l_low,l_high,g_low,g_high,l_low16,l_high16
	real*4 l_diff_min, l_diff_max
c ***** BOND STUFF *****
		real*4 bondscale,bondradius,dxy,dxyz,dxn,dyn,dxyn
		real*4 xstart,xend,ystart,yend,dpixel,dxnp,dynp,dcircle
		real*4 x1,x2,y1,y2,z1,z2,xmid,ymid,zpixel
		integer*4 si,sj
		real*4 costh,bondlengthsq(16,16),halflength(16)
		integer*4 bondlist(10000,2)
c ***** Postscript STUFF *****
		character pshex(16)
		integer*4 psint(16)
c ***** ETC. *****
		real*4 lx,ly,lz,lr,lyrot,lzrot,pshadowmax,pshaa
		real*4 x,y,z,rx,ry,rz,xp,yp,zp,d
		real*4 xn,yn,zn,ci,xs,xy,xz,cs
		real*4 rtable(-51:51,-51:51),coneangle,pcone,rcone
		character*3 comcode(14),command
		character*20 filename,inputfile
		character*3 shortname
		integer*4 ixsize,iysize, idepth
		integer*4 ix,iy,iz
		real*4 dh(3)
		integer bondflag,illustrationflag,shadowflag,outsideflag
c--------------------------------------------------------------------
		pi=3.141592
		rscale=1.

		backflag=0
		bondflag=0
		illustrationflag=0
		shadowflag=0
 		outsideflag=0

		data comcode/'rea','tra','xro','yro','zro','sca','cen',
     &		             'wor','map','cal','bac','ill','bon','mod'/
		do i=1,16
		do j=1,3
		wp(i,j)=1.
		rp(i,j)=1.
 19		isp(i,j)=10
		enddo
		enddo
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
 1	        call readatoms(radtype,coord,type,su,res,n,inputfile,
     &                        imodel) 
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
 8		read (5,*) lx,ly,lz
		lr=sqrt(lx*lx+ly*ly+lz*lz)
		lx=lx/lr
		ly=ly/lr
		lz=lz/lr
        write (6,202) ' light vector :     ',lx,ly,lz
		do i=1,16
 		read (5,*) rp(i,1),wp(i,1),isp(i,1),dif(i,1),
     &               rp(i,2),wp(i,2),isp(i,2),dif(i,2),
     &               rp(i,3),wp(i,3),isp(i,3),dif(i,3)
		do j=1,3
		pertot=rp(i,j)+wp(i,j)+dif(i,j)
		 if (pertot.gt.1.) then
		   rp(i,j)=rp(i,j)/pertot
		   wp(i,j)=wp(i,j)/pertot
		   dif(i,j)=dif(i,j)/pertot
 		 endif
		enddo
		enddo
c
		do i=1,3
 220		read (5,*) (shclip(ia,i),ia=1,8)
		enddo
		do i=1,3
 221		read (5,*) (shclip(ia,i),ia=9,16)
		enddo
c
		read (5,*) (iback(i),i=1,3),(ifog(i),i=1,3),pfogh,pfogl
		ibackgradient=0
		if (pfogl.lt.0.) then
		   pfogl=-pfogl
		   ibackgradient=1
		   endif
		do i=1,3
		if (iback(i).gt.255) iback(i)=255
		if (ifog(i).gt.255) ifog(i)=255
		rback(i)=float(iback(i))/256.
		rbacki(i)=rback(i)
 230		rfog(i)=float(ifog(i))/256.
		enddo
		if (pfogh.gt.1.) pfogh=1.
		if (pfogl.gt.1.) pfogl=1.
		write (6,204) ' background inten. :',(iback(i),i=1,3)
		write (6,204) ' fog intensity :    ',(ifog(i),i=1,3)
		write (6,202) ' upper fog percent :',pfogh*100.
		write (6,202) ' lower fog percent :',pfogl*100.
		pfogdiff=pfogh-pfogl
		read (5,*) shadowflag,pshadowmax,pshaa
		if (shadowflag.ne.0) write (6,202) ' draw shadows        '
		read (5,*) icone,pcone,coneangle,rcone
		if (icone.ne.0) write (6,202) ' draw conical shadows'
		read(5,*) stereoangle,istereooffset
		write(6,202)  ' stereo rotation :  ',stereoangle
		write(6,204)  ' stereo offset :    ',istereooffset
		stereoangle=stereoangle/180.*pi
		rcosstereo=cos(stereoangle)
		rsinstereo=sin(stereoangle)
		rtanstereo=tan(stereoangle)
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
		backflag=1
		read(5,78) filename
 78 		format (a20)
c	call ppmopen(filename,iyback,ixback,idepth,'r')
		write(6,*) 'Background will be input'
		write(6,*) 'Size of background file : ',ixback,iyback
		do ix=1,ixback
c	call ppmread(scanline,iyback*3)
		do iy=1,iyback
		pix(ix,iy,1)=float(scanline(iy*3-0))/256.
		pix(ix,iy,2)=float(scanline(iy*3-1))/256.
		pix(ix,iy,3)=float(scanline(iy*3-2))/256.
		enddo
		enddo
c	call tifclose

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
		bondflag=1
		read(5,*) (halflength(i),i=1,8)
		read(5,*) (halflength(i),i=9,16)
		read(5,*) bondradius
	i=1
 820	read(7,*,end=821) bondlist(i,1),bondlist(i,2)
	bondlist(i+1,1)=bondlist(i,2)
	bondlist(i+1,2)=bondlist(i,1)
	i=i+2
	goto 820
 821	nbondlist=i-1
	write(6,*) "Number of Bonds: ",nbondlist/2
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
		bondradius=bondradius*rscale
c ****  CREATE AND SCALE BOND LENGTH ARRAY *****
		if (bondflag.ne.0) then
		write(6,*) "BondFlag: ",bondflag
		do i=1,16
		do j=1,16
		 bondlengthsq(i,j)=((halflength(i)+halflength(j))*rscale)**2
		enddo
		enddo
		write(6,*) ' bond length array'
		do i=1,16
		write(6,138) (sqrt(bondlengthsq(i,j))/rscale,j=1,16)
		enddo
 138		format(16f4.1)
		endif
c ***** ADD STEREO ROTATION TO THE ROTATION MATRIX *****
		call clearmatrix(matrixin)
		matrixin(2,2)=cos(stereoangle)
		matrixin(2,3)=-sin(stereoangle)
		matrixin(3,2)=sin(stereoangle)
		matrixin(3,3)=cos(stereoangle)
		call catenate(matrixin,rm)
		do i=1,3
		do j=1,3
 112		rm(i,j)=matrixin(i,j)
		enddo
		enddo
c ***** ROTATE LIGHT SOURCE *****
		lyrot=ly*cos(stereoangle)+lz*sin(stereoangle)
		lzrot=-ly*sin(stereoangle)+lz*cos(stereoangle)
		ly=lyrot
		lz=lzrot
c ***** ROTATE THE TRANSLATION VECTOR *****
		ytrrot=ytran*cos(stereoangle)+ztran*sin(stereoangle)
		ztrrot=-ytran*sin(stereoangle)+ztran*cos(stereoangle)
		ytran=ytrrot
		ztran=ztrrot
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
		read(5,*) iras,ips
		if (iras.ne.0) then
		 read(5,113) filename
		 write(6,*) "output pnm filename: ",filename
		 open(8,file=filename,form='formatted')
		 write(8,1003) "P3"
 1003		 format(a2)
		 write(8,1004) iysize,ixsize
		 write(8,1004) 255
 1004	format(2i5)
c         call tifopen(filename,iysize,ixsize,24,'w')
		endif
		if (ips.ne.0) then
		 		read(5,78) filename
		 		write(6,*) "output ps filename",filename
		 		open(7,file=filename,form='formatted',
     &		          status='new')
		endif
 113		format(a20)
c ***** APPLY TRANSLATION AND SCALING *****
		do ia=1,n
		 coord(1,ia)=(coord(1,ia)+xtran)*rscale
		 coord(2,ia)=(coord(2,ia)+ytran)*rscale
 		 coord(3,ia)=(coord(3,ia)+ztran)*rscale
 114		enddo
c ***** MAP SPHERICAL SURFACES OVER ATOMS *****
		if (backflag.eq.0) then
		  do ix=1,ixsize
		    if (ibackgradient.ne.0) then
		    rgrad=float(ix)/float(ixsize)
 		    rback(1)=(1.-rgrad)*rbacki(1)+rgrad*rfog(1)
 		    rback(2)=(1.-rgrad)*rbacki(2)+rgrad*rfog(2)
 		    rback(3)=(1.-rgrad)*rbacki(3)+rgrad*rfog(3)
		    endif
		  do iy=1,iysize
 		  pix(ix,iy,1)=rback(1)
 		  pix(ix,iy,2)=rback(2)
 1001		  pix(ix,iy,3)=rback(3)
		  enddo
		  enddo
		endif

C ***** MAP SPHERICAL SURFACES OVER ATOMS *****
		do ix=1,ixsize
		do iy=1,iysize
		  atom(ix,iy)=0
 400		zpix(ix,iy)=-10000.
		enddo
		enddo
		nclipped=0

		if (n.gt.0) then

		do irad=1,16
C ----- create the spherical shading map for one type -----
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
		sphdat(ic,4)=1.
		z=sqrt(radtype(irad)**2-d*d)
		sphdat(ic,1)=x
		sphdat(ic,2)=y
		sphdat(ic,3)=z
		xn=x/radtype(irad)
		yn=y/radtype(irad)
		zn=z/radtype(irad)
		ci=xn*lx+yn*ly+zn*lz
		if (ci.lt.0.) then
 		 colormap(ix,iy,1)=0.
 		 colormap(ix,iy,2)=0.
 		 colormap(ix,iy,3)=0.
		ic=ic+1
		 goto 350
		endif
		xs=2.*ci*xn-lx		
		ys=2.*ci*yn-ly		
		zs=2.*ci*zn-lz
		rs=sqrt(xs*xs+ys*ys+zs*zs)
		xs=xs/rs
		ys=ys/rs
		zs=zs/rs
		cs=xs*xn+ys*yn+zs*zn
		if (cs.eq.0) cs=.0001
		colormap(ix,iy,1)=rp(irad,1)*ci+wp(irad,1)*(cs**isp(irad,1))
		colormap(ix,iy,2)=rp(irad,2)*ci+wp(irad,2)*(cs**isp(irad,2))
		colormap(ix,iy,3)=rp(irad,3)*ci+wp(irad,3)*(cs**isp(irad,3))
		ic=ic+1
 350		continue
		enddo
 352		enddo
		numpix=ic-1
c ----- then map spherical surface over atoms of the proper type ------
		icount=0
		nclipped=0
		do iatom=1,n
		if (type(iatom).ne.irad) goto 500
		icount=icount+1
		dclip=coord(3,iatom)*rcosstereo+coord(2,iatom)*rsinstereo
		if (abs(dclip).lt.radtype(type(iatom))) then
		 nclipped=nclipped+1
		 coord(1,n+nclipped)=coord(1,iatom)
		 coord(2,n+nclipped)=coord(2,iatom)-dclip*rsinstereo
		 coord(3,n+nclipped)=600.
		radtype(type(n+nclipped))=
     &         sqrt(radtype(type(iatom))**2-dclip*dclip)+1.
		 type(n+nclipped)=type(iatom)

C INTEGER
		 rad=radtype(type(n+nclipped))
		 iradtype=int(rad)
		 do irx=-iradtype-1,iradtype+1
		 do iry=-iradtype-1,iradtype+1
		 ix=float(irx)+coord(1,n+nclipped)+float(ixsize)/2.
		 iy=float(iry)+coord(2,n+nclipped)+float(iysize)/2.
		 if ((ix.gt.ixsize).or.(ix.lt.1).or.(iy.gt.iysize).or.
     &   (iy.lt.1)) goto 410
		 dist=sqrt(x*x+y*y/rcosstereo/rcosstereo)
		 if (dist.gt.rad) goto 410
		 z=sqrt(rad**2-dist*dist)+601.
		 if (z.gt.zpix(ix,iy)) then
		   zpix(ix,iy)=z
		   atom(ix,iy)=n+nclipped
		   pix(ix,iy,1)=shclip(type(iatom),1)
		   pix(ix,iy,2)=shclip(type(iatom),2)
		   pix(ix,iy,3)=shclip(type(iatom),3)
		 endif
 410		continue
		enddo
		enddo

		endif

		if (coord(3,iatom).lt.-coord(2,iatom)*rtanstereo) then

		do ipix=1,numpix
		x=sphdat(ipix,1)+coord(1,iatom)+float(ixsize)/2.		
		y=sphdat(ipix,2)+coord(2,iatom)+float(iysize)/2.
C INTEGER
		ix=int(x)
		iy=int(y)
		if ((x.gt.float(ixsize)).or.(x.lt.1.).or.
     &          (y.gt.float(iysize)).or.
     &          (y.lt.1.)) goto 510
		z=sphdat(ipix,3)+coord(3,iatom)		
		if (z.ge.-(y-float(iysize)/2.)*rtanstereo) goto 510
		if (z.gt.zpix(ix,iy)) then
		 zpix(ix,iy)=z
		 atom(ix,iy)=iatom
C INTEGER
		ixsphdat=int(sphdat(ipix,1))
		iysphdat=int(sphdat(ipix,2))
		 pix(ix,iy,1)=colormap(ixsphdat,iysphdat,1)
		 pix(ix,iy,2)=colormap(ixsphdat,iysphdat,2)
		 pix(ix,iy,3)=colormap(ixsphdat,iysphdat,3)
		endif
 510		continue
		enddo
		endif
 500		continue
		enddo
		write(6,*) icount,' spheres added of type: ',irad

 340		enddo
		write(6,*) 'shading maps written into depth buffer'
		write(6,*) 'number of clipped atoms : ',nclipped

		endif
C ***** CALCULATE BOND CYLINDERS *****
		if (bondflag.ne.0) then
		bondscale=100./bondradius

		do ic=1,16
c ----- first calculate shading map for bond -----

C INTEGER
		do ibx=-100,100
		do iby=-100,100
		x=float(ibx)
		y=float(iby)

		d=sqrt(x*x+y*y)
		z=sqrt(max(0.,10000.-d*d))
		xn=x/100.
		yn=y/100.
		zn=z/100.
		ci=xn*lx+yn*ly+zn*lz
		if (ci.lt.0.) then
 		 colormap(ibx,iby,1)=0.
 		 colormap(ibx,iby,2)=0.
 		 colormap(ibx,iby,3)=0.
		 goto 392
		endif
		xs=2.*ci*xn-lx		
		ys=2.*ci*yn-ly		
		zs=2.*ci*zn-lz
		rs=sqrt(xs*xs+ys*ys+zs*zs)
		xs=xs/rs
		ys=ys/rs
		zs=zs/rs
		cs=xs*xn+ys*yn+zs*zn
		colormap(ibx,iby,1)=rp(ic,1)*ci+
     &   wp(ic,1)*(cs**isp(ic,1))
		colormap(ibx,iby,2)=rp(ic,2)*ci+
     &   wp(ic,2)*(cs**isp(ic,2))
  		colormap(ibx,iby,3)=rp(ic,3)*ci+
     &   wp(ic,3)*(cs**isp(ic,3))
 392		continue
		enddo
		enddo
		if (ic.eq.1) then
		do i=-100,100,10
		write(6,391) (int(colormap(i,j,1)*9.),j=-100,100,10)
		enddo
 391		format(1x,21i1)
		endif
c ----- then calculate half cylinders -----
		nbond=0
 	do ibond=1,nbondlist
	ia=bondlist(ibond,1)
	ib=bondlist(ibond,2)
		if (type(ia).ne.ic) goto 550
	write(6,*) ia,ib,ic,type(ia)

		rsq=(coord(1,ia)-coord(1,ib))**2+
     &      (coord(2,ia)-coord(2,ib))**2+
     &      (coord(3,ia)-coord(3,ib))**2
		nbond=nbond+1

		x1=coord(1,ia)+float(ixsize)/2.
		y1=coord(2,ia)+float(iysize)/2.
		z1=coord(3,ia)
		xb=coord(1,ib)+float(ixsize)/2.
		yb=coord(2,ib)+float(iysize)/2.
		zb=coord(3,ib)
		x2=x1*.49+xb*.51
		y2=y1*.49+yb*.51
		z2=z1*.49+zb*.51

		if (z1.gt.z2) then
		  xt=x1
		  yt=y1
		  zt=z1
		  x1=x2
		  y1=y2
		  z1=z2
		  x2=xt
		  y2=yt
		  z2=zt
		endif

		dx=x2-x1
		dy=y2-y1
		dz=z2-z1
		dxy=sqrt(dx*dx+dy*dy)
		dxyz=sqrt(dxy*dxy+dz*dz)
		dxn=dx/dxy
		dyn=dy/dxy
		dzn=dz/dxyz
		dxyn=dxy/dxyz

		xstart=min(x1,x2)-bondradius-1
		xend=max(x1,x2)+bondradius+1
		if ((xstart.gt.float(ixsize)).or.(xend.lt.0)) goto 550
		xstart=max(1.,xstart)
		xend=min(float(ixsize),xend)
		ystart=min(y1,y2)-bondradius-1
		yend=max(y1,y2)+bondradius+1
		if ((ystart.gt.float(iysize)).or.(yend.lt.0)) goto 550
		ystart=max(1.,ystart)
		yend=min(float(iysize),yend)

C INTEGER
		do ibx=int(xstart),int(xend)
		do iby=int(ystart),int(yend)
		x=float(ibx)
		y=float(iby)
		dpixel=sqrt((x-x1)**2+(y-y1)**2)
		if (dpixel.eq.0) then
		  dxnp=0.001
		  dynp=0.001
		else
		  dxnp=(x-x1)/dpixel
		  dynp=(y-y1)/dpixel
		endif
		costh=dxn*dxnp+dyn*dynp
		rpixel=dpixel*sqrt(1.-min(costh*costh,1.))

		if (rpixel.gt.bondradius) goto 560
		rcircle=sqrt(bondradius**2-rpixel**2)

		dcircle=rcircle*dz/dxyz
		zmid=z1+(costh*dpixel+dcircle)/dxy*dz
		if ((zmid.gt.z2).or.(zmid.lt.z1)) goto 560

		zpixel=zmid+abs(rcircle*dxyn)
		if (zpixel-1.0.lt.zpix(ibx,iby)) goto 560
		xmid=x1+(zmid-z1)/dz*dx
		ymid=y1+(zmid-z1)/dz*dy
		si=int((x-xmid)*bondscale)
		sj=int((y-ymid)*bondscale)

		atom(ibx,iby)=ia
		zpix(ibx,iby)=zpixel
		pix(ibx,iby,1)=colormap(si,sj,1)
		pix(ibx,iby,2)=colormap(si,sj,2)
		pix(ibx,iby,3)=colormap(si,sj,3)

 560		continue
 		enddo
 		enddo

 550		continue
		enddo
		write(6,*) nbond,' bonds added for type : ',ic

		enddo
		write(6,*)
		write(6,*) ' **** CYLINDRICAL BONDS ADDED TO DEPTH BUFFER *****'
		write(6,*)
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
c ***** MAP SHADOWS INTO SHADOW PLANE *****
		if (shadowflag.ne.0) then

		write(6,*) xshmax,yshmax
		do ix=1,int(xshmax)
		do iy=1,int(yshmax)
		shadowatom(ix,iy)=0
 600		zshadow(ix,iy)=-10000.
		enddo
		enddo
		write(6,*) xshmax,yshmax

		xsh_offset=xshmax-ixsize/2
		if (ly.lt.0) ysh_offset=yshmax-iysize/2
		if (ly.eq.0) ysh_offset=yshmax/2.
		if (ly.gt.0) ysh_offset=iysize/2
		write(6,*) 'shadow map size: ',xshmax,yshmax
		write(6,*) 'image is centered on: ',xsh_offset,ysh_offset

		do irad=1,16

c ----- first create shadow map for each atom type ----
		ic=1
		rshadow=radtype(irad)*0.9
C INTEGER
		do ishx=-int(rshadow*1.5),int(rshadow*1.5)
		do ishy=-int(rshadow*1.5),int(rshadow*1.5)
		x=int(ishx)
		y=int(ishy)

		xd=(ly*ly+lz*lz)*x-lx*ly*y
		yd=(lx*lx+lz*lz)*y-lx*ly*x
		zd=-lx*lz*x-ly*lz*y
		d=sqrt(xd*xd+yd*yd+zd*zd)

		if (d-rshadow.gt.0) goto 351

		rsha=(lx*lx+ly*ly)/(lz*lz)+1.
		rshb=2.*(lx*x+ly*y)/lz
		rshc=x*x+y*y-rshadow*rshadow
		z=(-rshb+sqrt(max(0.,rshb**2-4.*rsha*rshb)))/(2.*rsha)

		sphdat(ic,1)=x
		sphdat(ic,2)=y
		sphdat(ic,3)=z-     0.
c		** note kludge factor  ^^ ***
		ic=ic+1
 351		continue
		enddo
		enddo
		numshadowpix=ic-1
		write(6,*) 'numshadowpix : ',irad,numshadowpix

c ----- then map shadows into shadow plane -----

		do iatom=1,n
		if (type(iatom).ne.irad) goto 610

		if (coord(3,iatom).lt.-coord(2,iatom)*rtanstereo) then
     		rx=lx*coord(3,iatom)/lz
     		ry=ly*coord(3,iatom)/lz

		do ipix=1,numshadowpix
		z=coord(3,iatom)+sphdat(ipix,3)
		x=coord(1,iatom)+xsh_offset+sphdat(ipix,1)-rx
		y=coord(2,iatom)+ysh_offset+sphdat(ipix,2)-ry
		if ((x.gt.xshmax).or.(x.lt.1).or.(y.gt.yshmax).or.
     &		    (y.lt.1.)) goto 620
		if (z.ge.-(y-600.)*rtanstereo) goto 620
C INTEGER
		ishx=int(x)
		ishy=int(y)
		if (z.gt.zshadow(ishx,ishy)) then
		  zshadow(ishx,ishy)=z
		  shadowatom(ishx,ishy)=iatom
		endif
 620		continue
		enddo
		endif
 610		continue
		enddo

 341		continue
		enddo

cdo ixo=xsh_offset-ixsize/2,xsh_offset+ixsize/2,20
cwrite(6,6669) (int(shadowatom(ixo,iyo)),
c    &     iyo=ysh_offset-iysize/2,ysh_offset+iysize/2,20)
cenddo
 6669		format(65i1)
		write(6,*) 'Shadows mapped into shadow surface'
		endif
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

c **** DIFFUSE ILLUMINATION *****
		dh(1)=0.
		dh(2)=0.
		dh(3)=0.
		if (atom(ix,iy).ne.0) then
		  dh(1)=dif(type(atom(ix,iy)),1)
		  dh(2)=dif(type(atom(ix,iy)),2)
		  dh(3)=dif(type(atom(ix,iy)),3)
		endif
c ***** SHADOW TESTING *****
		psh=1.
		if ((shadowflag.ne.0).and.(atom(ix,iy).ne.0)) then
		  zp=zpix(ix,iy)
		  xp=float(ix)-float(ixsize)/2.+xsh_offset-lx*zp/lz
		  yp=float(iy)-float(iysize)/2.+ysh_offset-ly*zp/lz
C INTEGER
		  ixp=int(xp)
		  iyp=int(yp)
		  if ((xp.lt.0).or.(xp.gt.xshmax).or.(yp.lt.0).or.
     &        (yp.gt.yshmax)) then
	    if (outsideflag.eq.0) write(6,*) 'outside shadow map at: ',
     &          ix,iy,xp,yp
	   	outsideflag=1
		    psh=pshadowmax
          else
		  if ((atom(ix,iy).ne.shadowatom(ixp,iyp)).and.
     &        (zshadow(ixp,iyp).gt.zp+0.)) psh=pshadowmax
   		  endif
		endif
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
c	   pconexy=(10.-rtable(i,j))*coneangle
		pconetot=pconetot-pcone
c	psh=max(0,psh)
c	   psh=min(pcone,psh)
c	   goto 5000
		 endif
		 endif
		 endif
		enddo
		enddo
 5000		continue
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
c ----- diagnostic ------
c	if ((int(ix/20)*20.eq.int((float(ix)/20.)*20.)).and.
c    &			(iy.eq.iysize)) then
c	 write(6,6669) (int(pix(ix,iyo,1)*9.),iyo=1,iysize,20)
c	endif
c>>>NOTE>>>shadows applied to diffuse illumination!!
	psh=min(psh,pconetot)
		do icolor=1,3
         rframe(icolor,iy)=(1.-l_opacity)*
     &     (pfh*(psh*(pix(ix,iy,icolor)+dh(icolor)))+
     &     (1.-pfh)*rfog(icolor)) 
		   pix(ix,iy,icolor)=rframe(icolor,iy)
 5020		continue
		enddo

 2000		enddo

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
 1002	format(20i4)
c	if (iras.ne.0) call tifwf3(rframe)
c ----- diagnostic ------
		if (int(ix/20)*20.eq.int((float(ix)/20.)*20.)) then
		 write(6,6669) (int(pix(ix,iyo,1)*9.),iyo=1,iysize,20)
		endif
c ----- write z-level map ----
cif (zwrite.eq.1) write(77) (zpix(ix,iyo),iyo=1,iysize)

 1000		enddo
c	call tifclose

c ----- Postscript output -----
		if (ips.ne.0) then

		isc=4
		icx=100
		icy=150
		write(6,*) 'generating postscript'
		write(7,712) '%!'
		write(7,774) '/Times-Roman findfont 15 scalefont setfont'
		write(7,770) 100,100,'moveto   '
		write(7,779) '(',inputfile,') show'
		write(7,770) icx,icy,'translate'
		write(7,771) 3.84,3.84,'scale    '
		rsc=4./72.
		rxc=float(icx-ixsize/2)
		ryc=float(icy-iysize/2)
		iseed=int(100*rxc)

		lines=0
		guide_lines=1
		random_dot=0
		dither=0

		do ix=1,ixsize+8
		do iy=1,iysize+8
		  pixmin=min(pix(ix,iy,1),pix(ix,iy,2),pix(ix,iy,3))
		  pix(ix,iy,1)=1.

      if (random_dot.eq.1) then
C      rnum=ran(iseed)
		   if (pixmin.lt.rnum) pix(ix,iy,1)=0.
		  if (pixmin.gt.0.9) pix(ix,iy,1)=1.
		  if (pixmin.lt.0.1) pix(ix,iy,1)=0.
		  endif

		if (guide_lines.eq.1) then
		   if (pixmin.lt..75) then
		    rtestx=float(ix)/3.-float((ix)/3)
		    if (rtestx.eq.0) pix(ix,iy,1)=0.
		   endif
		   if (pixmin.lt..5) then
		    rtesty=float(iy)/3.-float((iy)/3)
		    if (rtesty.eq.0) pix(ix,iy,1)=0.
		   endif
		   if (pixmin.lt..25) then
		    rtesty=float(iy+1)/3.-float((iy+1)/3)
		    rtestx=float(ix+1)/3.-float((ix+1)/3)
		    if (rtesty.eq.0) pix(ix,iy,1)=0.
		    if (rtestx.eq.0) pix(ix,iy,1)=0.
		   endif
		   if (pixmin.lt..10) pix(ix,iy,1)=0.
		endif
		  if (lines.eq.1) then
		   if ((pixmin.gt..6).and.(pixmin.lt..7)) then
		    rtestx=float(ix)/3.-float((ix)/3)
		    rtesty=float(iy)/3.-float((iy)/3)
		    if ((rtestx.eq.0).and.(rtesty.eq.0)) pix(ix,iy,1)=0.
		   endif
		   if ((pixmin.gt..4).and.(pixmin.lt..6)) then
		    rtestx=float(ix)/3.-float((ix)/3)
		    rtesty=float(iy)/3.-float((iy)/3)
		    if ((rtestx.eq.0).and.(rtesty.eq.0)) pix(ix,iy,1)=0.
		    rtestx=float(ix+1)/3.-float((ix+1)/3)
		    rtesty=float(iy+1)/3.-float((iy+1)/3)
		    if ((rtestx.eq.0).and.(rtesty.eq.0)) pix(ix,iy,1)=0.
		    rtestx=float(ix+2)/3.-float((ix+2)/3)
		    rtesty=float(iy+2)/3.-float((iy+2)/3)
		    if ((rtestx.eq.0).and.(rtesty.eq.0)) pix(ix,iy,1)=0.
		   endif
		   if ((pixmin.gt..3).and.(pixmin.lt..4)) then
		    rtestx=float(ix)/3.-float((ix)/3)
		    rtesty=float(iy)/3.-float((iy)/3)
		    if ((rtestx.eq.0).or.(rtesty.eq.0)) pix(ix,iy,1)=0.
		   endif
		    if (pixmin.le.0.3) pix(ix,iy,1)=0.
		  endif

		  if (dither.eq.1) then
C      			rnum=ran(iseed)
C		pixmin=pixmin+(rnum)/20.
C		pixmin=min(1.0,pixmin)
		    rx=float(ix)/3.-float((ix)/3)
		    ry=float(iy)/3.-float((iy)/3)
		    if ((rx.eq.0).and.(ry.eq.0).and.
     &               (pixmin.le..9)) pix(ix,iy,1)=0.
		    rx=float(ix+1)/3.-float((ix+1)/3)
		    ry=float(iy)/3.-float((iy)/3)
		    if ((rx.eq.0).and.(ry.eq.0).and.
     &               (pixmin.le..8)) pix(ix,iy,1)=0.
		    rx=float(ix+1)/3.-float((ix+1)/3)
		    ry=float(iy+1)/3.-float((iy+1)/3)
		    if ((rx.eq.0).and.(ry.eq.0).and.
     &               (pixmin.le..7)) pix(ix,iy,1)=0.
		    rx=float(ix)/3.-float((ix)/3)
		    ry=float(iy+1)/3.-float((iy+1)/3)
		    if ((rx.eq.0).and.(ry.eq.0).and.
     &               (pixmin.le..6)) pix(ix,iy,1)=0.
		    rx=float(ix+2)/3.-float((ix+2)/3)
		    ry=float(iy)/3.-float((iy)/3)
		    if ((rx.eq.0).and.(ry.eq.0).and.
     &               (pixmin.le..5)) pix(ix,iy,1)=0.
		    rx=float(ix+2)/3.-float((ix+2)/3)
		    ry=float(iy+1)/3.-float((iy+1)/3)
		    if ((rx.eq.0).and.(ry.eq.0).and.
     &               (pixmin.le..4)) pix(ix,iy,1)=0.
		    rx=float(ix)/3.-float((ix)/3)
		    ry=float(iy+2)/3.-float((iy+2)/3)
		    if ((rx.eq.0).and.(ry.eq.0).and.
     &               (pixmin.le..3)) pix(ix,iy,1)=0.
		    rx=float(ix+1)/3.-float((ix+1)/3)
		    ry=float(iy+2)/3.-float((iy+2)/3)
		    if ((rx.eq.0).and.(ry.eq.0).and.
     &               (pixmin.le..2)) pix(ix,iy,1)=0.
		    rx=float(ix+2)/3.-float((ix+2)/3)
		    ry=float(iy+2)/3.-float((iy+2)/3)
		    if ((rx.eq.0).and.(ry.eq.0).and.
     &               (pixmin.le..1)) pix(ix,iy,1)=0.
		endif

		 		if ((ix.ge.ixsize).or.(iy.ge.iysize)) pix(ix,iy,1)=1.
		 		if (l_opacity.gt.0.1) pix(ix,iy,1)=0.
c if ((ix.le.1).and.(iy.le.iysize)) pix(ix,iy,1)=0.
c if ((iy.le.1).and.(ix.le.ixsize)) pix(ix,iy,1)=0.
c if ((ix.eq.ixsize).and.(iy.le.iysize)) pix(ix,iy,1)=0.
c if ((iy.eq.iysize).and.(ix.le.ixsize)) pix(ix,iy,1)=0.
 702		enddo
 701		enddo

		do ix=1,ixsize,8
		do iy=1,iysize,8
		 do i=0,7
		 psint(i*2+1)=pix(ix+i,iy+0,1)*8.+
     &		              pix(ix+i,iy+1,1)*4.+
     &		              pix(ix+i,iy+2,1)*2.+
     &		              pix(ix+i,iy+3,1)*1.
		 psint(i*2+2)=pix(ix+i,iy+4,1)*8.+
     &		              pix(ix+i,iy+5,1)*4.+
     &		              pix(ix+i,iy+6,1)*2.+
     &		              pix(ix+i,iy+7,1)*1.
 705		 enddo
		 iwhite=0
		do j=1,16
		 if (psint(j).eq.0) pshex(j)='0'
		 if (psint(j).eq.1) pshex(j)='1'
		 if (psint(j).eq.2) pshex(j)='2'
		 if (psint(j).eq.3) pshex(j)='3'
		 if (psint(j).eq.4) pshex(j)='4'
		 if (psint(j).eq.5) pshex(j)='5'
		 if (psint(j).eq.6) pshex(j)='6'
		 if (psint(j).eq.7) pshex(j)='7'
		 if (psint(j).eq.8) pshex(j)='8'
		 if (psint(j).eq.9) pshex(j)='9'
		 if (psint(j).eq.10) pshex(j)='A'
		 if (psint(j).eq.11) pshex(j)='B'
		 if (psint(j).eq.12) pshex(j)='C'
		 if (psint(j).eq.13) pshex(j)='D'
		 if (psint(j).eq.14) pshex(j)='E'
		 if (psint(j).eq.15) then
		    pshex(j)='F'
		    iwhite=iwhite+1
		 endif
 706		enddo

		 if (iwhite.lt.16) then
         write(7,771) float(iy)/8.-ryl,
     &               float(ixsize-ix)/8.-rxl,'translate'
		 rxl=float(ixsize-ix)/8.
		 ryl=float(iy)/8.
         write(7,773) '8 8 1 [8 0 0 -8 0 8] {<',(pshex(j),j=1,16),
     &		      '>} image'
		 endif

 704		 enddo
 703		 enddo
		 write(7,772) 'showpage'

 770		format(1x,2i5,1x,a9)
 771		format(1x,2f10.4,1x,a9)
 773		format(1x,a23,16a1,a8)
 772		format(1x,a8)
 712		format(a2)
 774		format(a42)
 779		format(1x,a1,a20,a6)

		endif


 999		stop
		end
c--------------------------------------------------------------------
	  subroutine readatoms(radtype,coord,type,su,res,n,inputfile,
     &                          imodel) 
c****** PDB FORMAT ! ********
		real*4 coord(3,350000)
		integer*4 type(350000),res(0:350000),su(0:350000)
		integer*4 resrange(2,100),inptype(100),inpsu(100)
		integer*4 n
		real*4 radtype(16)
		character*6 atomdescriptor(100)
		character*10 descriptor(100)
		character*20 inputfile,filetype
		character*80 instring
		integer*4 imodel,imodelcurrent
c--------------------------------------------------------------------
		imodelcurrent=0
		read(5,10) inputfile
 10		format(a20)
		open(1,file=inputfile,form='formatted',status='old')
c --- read atom radii ---
		read(5,*) (radtype(i),i=1,8)
		read(5,*) (radtype(i),i=9,16)
 111	format(8f10.5)
c --- read atom descriptors ---
		ndes=0
 20		read(5,100) instring
 		if (instring(1:3).eq.'END') goto 30
 		  ndes=ndes+1
 		  read(instring,21) atomdescriptor(ndes),descriptor(ndes)
 		  read(instring(18:80),*) (resrange(i,ndes),i=1,2),
     &                      inptype(ndes),inpsu(ndes)
		goto 20
 30		write(6,*) ' atom descriptors: ',ndes
C		do i=1,ndes
C		write(99,*) '.',atomdescriptor(i),'.',descriptor(i),
C    &		resrange(1,i),resrange(2,i),inptype(i),inpsu(i)
C		enddo
 21		format(a6,a10)
 22		format(4i8)
c --- read atoms and classify ---
		n=0
 40		read(1,100,end=9) instring
		if (instring(1:5).eq.'MODEL') then
		  imodelcurrent=imodelcurrent+imodel
		  goto 40
		endif
		if ((instring(1:4).ne.'ATOM').and.
     &      (instring(1:6).ne.'HETATM')) goto 40
     	read(instring,200) ires
C    	write(99,*) '.',instring(1:6),'.',
C    &   (instring(12+ia:12+ia),ia=1,10),ires
		do ides=1,ndes

 		 if (instring(1:6).ne.atomdescriptor(ides)(1:6)) goto 50

 		 do ia=1,10
 		  if (descriptor(ides)(ia:ia).eq.'-') goto 60
 		  if (instring(12+ia:12+ia).ne.descriptor(ides)(ia:ia)) goto 50
 60 	         continue
		enddo

 		 if ((ires.lt.resrange(1,ides)).or.
     &        (ires.gt.resrange(2,ides))) goto 50

		 if (inptype(ides).eq.0) goto 40

 		 n=n+1
		 read(instring,300) (coord(i,n),i=1,3)
		 type(n)=inptype(ides)
		 su(n)=inpsu(ides)+imodelcurrent
		 res(n)=ires
		 goto 40
 50 	continue
		enddo
		goto 40
c --- done ---         
 9		write(6,*)' atoms read: ', n, ' from: ',inputfile
		write(6,*)' '
c	 do j=1,n
c	 write(99,*) (coord(i,j),i=1,3),type(j),su(j),res(j)
c	 enddo
 100	format(a80)
 200	format(22x,i4)
 300	format(30x,3f8.3)
		return
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
