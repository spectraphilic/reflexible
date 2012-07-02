C     ################################################################
C     # FortFlex fortran routine for reading FLEXPART output
C     # Use with f2py, on linux:
C     #
C     # %f2py -c -m FortFlex FortFlex.f
C     # 
C     # This will create a FortFlex.so module that can be read by
C     # Python. It is used by PyFlex.readgridV8
C     #
C     # You will probably have to run this yourself as it is quite 
C     # dependent on libraries and machine variables. Also, will be
C     # different between 64bit and 32bit machines
C     #
C     # JFB, 02.07.2009
C     ################################################################
      
      
      subroutine readgrid(filegrid,
     +numxgrid,numygrid,numzgrid,numpointspec,nageclass,
     +grid,wetgrid,drygrid,itime,scaledepo,scaleconc,
     +decayconstant,npspec_int)
      
      parameter(smallnum=1.e-38)

      real decayconstant
      integer lage(0:nageclass),itime,numpointspec
Cf2py intent(out) itime
      real grid(numxgrid,numygrid,numzgrid,numpointspec,nageclass)
Cf2py intent(out) grid
      real wetgrid(numxgrid,numygrid,numpointspec,nageclass)
Cf2py intent(out) wetgrid
      real drygrid(numxgrid,numygrid,numpointspec,nageclass)
Cf2py intent(out) drygrid
c     creating a npspec option to read only one npspec
Cf2py intent(in) npspec_int  
      real npspec_int


***************** new sparse output
      integer fact,ii,ir,npspec,np
      real,allocatable,dimension (:) :: sparse_dump_r
      integer,allocatable,dimension (:) :: sparse_dump_i
      integer sp_count_i,sp_count_r,stat
***************** new sparse output

      character*500 filegrid
C      write(*,*) 'FortFlex Dynamic Memory Allocation'

c allocate memory for sparse matrix
      allocate(sparse_dump_r(numxgrid*numygrid*numzgrid),STAT=stat)
      if (stat.ne.0) print*,'***error allocating sparse_dump_r'
c      write(*,*) 'allocation 1'
      allocate(sparse_dump_i(numxgrid*numygrid*numzgrid),STAT=stat)
      if (stat.ne.0) print*,'***error allocating sparse_dump_i'
c      write(*,*) 'allocation 2'

c     print*,numxgrid,numygrid,numzgrid,nageclass,numpointspec

      open(10,file=filegrid,form='unformatted',status='old')
      read(10,end=99) itime 
c     write(*,*) itime

C Initialize age class fields
*****************************

      wetgrid=0.
      drygrid=0.
      grid=0.

C Loop about all age classes
****************************
      if (npspec_int.gt.0) then
          numpointspec=npspec_int
      endif

      do 38 np=1,numpointspec
      npspec=np
      if (npspec_int.gt.0) then
          npspec=1
      endif

      do 38 nage=1,nageclass

C      write(*,*) nage, np
c     print*,npspec

       do 63 ix=1,numxgrid
         do 63 jy=1,numygrid
           wetgrid(ix,jy,npspec,nage)=0.
           drygrid(ix,jy,npspec,nage)=0.
           do 63 kz=1,numzgrid
63           grid(ix,jy,kz,npspec,nage)=0.

C Read wet deposition
*********************

        fact=1
        read(10) sp_count_i
        read(10) (sparse_dump_i(ix),ix=1,sp_count_i)
        read(10) sp_count_r
        read(10) (sparse_dump_r(ix),ix=1,sp_count_r)
c       write(*,*) sp_count_i,sp_count_r
c       write(*,*) (sparse_dump_i(ix),ix=1,sp_count_i)
c       write(*,*) (sparse_dump_r(ix),ix=1,sp_count_r)

        ii=0
        do 32 ir=1,sp_count_r
          if ((sparse_dump_r(ir)*fact).gt.smallnum) then
             ii=ii+1


             n=sparse_dump_i(ii)
             fact=fact*(-1.)
          else
             pos=pos+1
          endif

          jy=n/numxgrid
          ix=n-numxgrid*jy
          wetgrid(ix+1,jy+1,npspec,nage)=sparse_dump_r(ir)*(-1)*fact

32      continue

c       print*,'wet grid read'

C Read dry deposition
*********************

        fact=1
        read(10) sp_count_i
        read(10) (sparse_dump_i(ix),ix=1,sp_count_i)
        read(10) sp_count_r
        read(10) (sparse_dump_r(ix),ix=1,sp_count_r)
c       write(*,*) sp_count_i,sp_count_r
c       write(*,*) (sparse_dump_i(ix),ix=1,sp_count_i)
c       write(*,*) (sparse_dump_r(ix),ix=1,sp_count_r)

        ii=0
        do 36 ir=1,sp_count_r
          if ((sparse_dump_r(ir)*fact).gt.smallnum) then
             ii=ii+1
             n=sparse_dump_i(ii)
             fact=fact*(-1.)
          else
             pos=pos+1
          endif

          jy=n/numxgrid
          ix=n-numxgrid*jy
          drygrid(ix+1,jy+1,npspec,nage)=sparse_dump_r(ir)*(-1)*fact

36      continue

c       print*,'dry grid read'

C Read concentrations
*********************

        fact=1
        read(10) sp_count_i
        read(10) (sparse_dump_i(ix),ix=1,sp_count_i)
        read(10) sp_count_r
        read(10) (sparse_dump_r(ix),ix=1,sp_count_r)
c       write(*,*) sp_count_i,sp_count_r
c       write(*,*) (sparse_dump_i(ix),ix=1,sp_count_i)
c       write(*,*) (sparse_dump_r(ix),ix=1,sp_count_r)

        ii=0
        do 47 ir=1,sp_count_r
          if ((sparse_dump_r(ir)*fact).gt.smallnum) then
             ii=ii+1
             n=sparse_dump_i(ii)
             fact=fact*(-1.)
          else
             n=n+1
          endif

          kz=n/(numxgrid*numygrid)
          jy=(n-kz*numxgrid*numygrid)/numxgrid
          ix=n-numxgrid*numygrid*kz-numxgrid*jy
          grid(ix+1,jy+1,kz,npspec,nage)=sparse_dump_r(ir)*(-1)*fact

47      continue

c       print*,'conc grid read'
c       print*,kz,n,numxgrid*numygrid

C Scale values as defined by user
*********************************

c        fac=exp(-1.*float(lage(nage)+lage(nage-1))/2./decayconstant)
c        if (lage(nage).gt.1e8) fac=1.0
c        print*,fac
        fac = 1.0
        do 41 ix=1,numxgrid
           do 41 jy=1,numygrid
             drygrid(ix,jy,npspec,nage)= 
     +             drygrid(ix,jy,npspec,nage)*scaledepo
             wetgrid(ix,jy,npspec,nage)= 
     +        wetgrid(ix,jy,npspec,nage)*scaledepo
             do 41 kz=1,numzgrid
              grid(ix,jy,kz,npspec,nage)=
     +             grid(ix,jy,kz,npspec,nage)*fac*scaleconc
41             continue


C End species loop, end age class loop
**************************************

138       continue
38      continue

      close(10)

99    continue

      if (allocated(sparse_dump_r)) deallocate(sparse_dump_r)
      if (allocated(sparse_dump_i)) deallocate(sparse_dump_i)

      return
      end

      
      
      
      subroutine sumgrid(zplot,grid,
     +numxgrid,numygrid,numzgrid,
     +numpoint,nageclass,
     +area,heightnn)
     
      real grid(numxgrid,numygrid,numzgrid,numpoint,nageclass)
Cf2py intent(in) grid      
      real zplot(numxgrid,numygrid,numzgrid,numpoint)
Cf2py intent(in,out) zplot
      real heightnn(numxgrid,numygrid,numzgrid)
Cf2py intent(in) heightnn    
      real area(numxgrid,numygrid)
Cf2py intent(in) areas

     
********************************************************
C Add contributions from this time step to gridded field
********************************************************
      do 203 k=1,numpoint
C      write(*,*) k
        do 201 ix=1,numxgrid
          do 201 jy=1,numygrid
              do 220 kz=1,numzgrid
              if (kz.eq.1) then
              contribution=grid(ix,jy,kz,k,nageclass)/area(ix,jy)
C     +       /(heightnn(ix,jy,kz))
              else
              contribution=grid(ix,jy,kz,k,nageclass)/area(ix,jy)
C     +        /(heightnn(ix,jy,kz)-heightnn(ix,jy,kz-1))
              endif
              zplot(ix,jy,kz,k)=zplot(ix,jy,kz,k)+contribution
220         continue
201         continue

203    continue 
       end


      subroutine readgrid_v6(filegrid,
     +numxgrid,numygrid,numzgrid,numpointspec,nageclass,
     +grid,wetgrid,drygrid,itime,scaledepo,scaleconc,
     +decayconstant)
      
C     +nzmax,numzgrid,maxspec,nspec,nspec2,maxageclass,nageclass,grid,
C     +lage,scaleconc,decayconstant)


      real decayconstant
      integer lage(0:nageclass),itime,numpointspec
Cf2py intent(out) itime
      real grid(numxgrid,numygrid,numzgrid,numpointspec,nageclass)
Cf2py intent(out) grid
      real wetgrid(numxgrid,numygrid,numpointspec,nageclass)
Cf2py intent(out) wetgrid
      real drygrid(numxgrid,numygrid,numpointspec,nageclass)
Cf2py intent(out) drygrid
***************** new sparse output
      real smallnum
      integer fact,ii,ir
      real sparse_dump_r(numxgrid*numygrid*numzgrid)
      integer sparse_dump_i(numxgrid*numygrid*numzgrid)
      integer sp_count_i,sp_count_r
***************** new sparse output


      character*500 filegrid

      open(10,file=filegrid,form='unformatted',status='old')
      read(10,end=99) itime 
      write(*,*) itime

C Loop about all species
************************

      do 38 k=1,numpointspec

C Loop about all age classes
****************************

        do 38 nage=1,nageclass

C Initialize age class fields
*****************************

        do 63 ix=1,numxgrid
          do 63 jy=1,numygrid
            do 63 kz=1,numzgrid
63            grid(ix,jy,kz,k,nage)=0.


C Read wet deposition
*********************

        fact=1
        read(10) sp_count_i
        read(10) (sparse_dump_i(ix),ix=1,sp_count_i)
        read(10) sp_count_r
        read(10) (sparse_dump_r(ix),ix=1,sp_count_r)


c       ii=0
c       do 32 ir=1,sp_count_r
c         if ((sparse_dump_r(ir)*fact).gt.smallnum) then
c            ii=ii+1
c            n=sparse_dump_i(ii)
c            fact=fact*(-1.)
c         else
c            pos=pos+1
c         endif

c         jy=n/numxgrid
c         ix=n-numxgrid*jy
c         wetgrid(ix+1,jy+1,k,nage)=sparse_dump_r(ir)

c2      continue

C Read dry deposition
*********************

        fact=1
        read(10) sp_count_i
        read(10) (sparse_dump_i(ix),ix=1,sp_count_i)
        read(10) sp_count_r
        read(10) (sparse_dump_r(ix),ix=1,sp_count_r)

c       ii=0
c       do 36 ir=1,sp_count_r
c         if ((sparse_dump_r(ir)*fact).gt.smallnum) then
c            ii=ii+1
c            n=sparse_dump_i(ii)
c            fact=fact*(-1.)
c         else
c            pos=pos+1
c         endif

c         jy=n/numxgrid
c         ix=n-numxgrid*jy
c         drygrid(ix+1,jy+1,k,nage)=sparse_dump_r(ir)

c6      continue


C Read concentrations
*********************

        fact=1
        read(10) sp_count_i
        read(10) (sparse_dump_i(ix),ix=1,sp_count_i)
        read(10) sp_count_r
        read(10) (sparse_dump_r(ix),ix=1,sp_count_r)
c       write(*,*) sp_count_i,sp_count_r
c       write(*,*) (sparse_dump_i(ix),ix=1,sp_count_i)
c       write(*,*) (sparse_dump_r(ix),ix=1,sp_count_r)

        ii=0
        do 47 ir=1,sp_count_r
          if ((sparse_dump_r(ir)*fact).gt.smallnum) then
             ii=ii+1
             n=sparse_dump_i(ii)
             fact=fact*(-1.)
          else
             n=n+1
          endif

          kz=n/(numxgrid*numygrid)
          jy=(n-kz*numxgrid*numygrid)/numxgrid
          ix=n-numxgrid*numygrid*kz-numxgrid*jy
          grid(ix+1,jy+1,kz,k,nage)=abs(sparse_dump_r(ir))

c         if (kz.eq.1) write(96,*) k,ix,jy,grid(ix+1,jy+1,kz,k,nage)

47      continue



C Sum up all age classes to total fields and scale values as defined by user
****************************************************************************

          do 41 ix=1,numxgrid
            do 41 jy=1,numygrid
              do 41 kz=1,numzgrid
c               grid(ix,jy,kz,k,nage)=grid(ix,jy,kz,k,nage)*
c    +          exp(-1.*float(lage(nage)+lage(nage-1))/2./
c    +          decayconstant)
                grid(ix,jy,kz,k,nage)=grid(ix,jy,kz,k,nage)*scaleconc
41              continue


C End species loop, end age class loop
**************************************

138       continue
        if (k.eq.lspec2) goto 99
38      continue

      close(10)

99    continue


      end

       
      subroutine readheader(filename,nxmax,numxgrid,nymax,numygrid,
     +nzmax,numzgrid,outlon0,outlat0,dxout,dyout,outheight,ibdate,
     +ibtime,loutstep,maxspec,nspec,maxageclass,nageclass,lage,
     +ireleasestart,ireleaseend,maxpoint,numpoint,xpoint,ypoint,
     +zpoint1,zpoint2,heightnn,area)

      parameter(pi=3.14159265,r_earth=6.371e6,pih=pi/180.)

      character filename*150,compoint(maxpoint)*45,species(maxspec)*7
Cf2py intent(out) compoint
Cf2py intent(out) species
Cf2py intent(in) filename 
      real outheight(nzmax),heightnn(nxmax,nymax,0:nzmax)
Cf2py intent(out) outheight
Cf2py intent(out) heightnn
      integer ireleasestart(maxpoint),ireleaseend(maxpoint)
Cf2py intent(out) ireleasestart
Cf2py intent(out) ireleaseend
      integer npart(maxpoint),kind(maxpoint),lage(0:maxageclass)
Cf2py intent(out) npart
Cf2py intent(out) kind
Cf2py intent(out) lage
      real xpoint(maxpoint),ypoint(maxpoint),zpoint1(maxpoint)
Cf2py intent(out) xpoint
Cf2py intent(out) ypoint
Cf2py intent(out) zpoint1
      real zpoint2(maxpoint),xmass(maxpoint,maxspec)
Cf2py intent(out) zpoint2
Cf2py intent(out) xmass
      real oro(nxmax,nymax),area(nxmax,nymax)
Cf2py intent(out) oro
Cf2py intent(out) area


      open(10,file=filename,form='unformatted',status='old')
      read(10) ibdate,ibtime
Cf2py intent(out) ibdate
Cf2py intent(out) ibtime
      write(*,*) ibdate,ibtime
      read(10) loutstep,loutaver,loutsample
Cf2py intent(out) loutstep
Cf2py intent(out) loutaver
Cf2py intent(out) loutsample
      read(10) outlon0,outlat0,numxgrid,numygrid,dxout,dyout
Cf2py intent(out) outlon0
Cf2py intent(out) outlat0
Cf2py intent(out) numxgrid
Cf2py intent(out) numygrid
Cf2py intent(out) dxout
Cf2py intent(out) dyout
      write(*,*) outlon0,outlat0,numxgrid,numygrid,dxout,dyout
      read(10) numzgrid,(outheight(i),i=1,numzgrid)
Cf2py intent(out) numzgrid
      read(10) jjjjmmdd,ihmmss
Cf2py intent(out) jjjjmmdd
Cf2py intent(out) hhmmss
      write(*,*) jjjjmmdd,ihmss
      read(10) nspec
Cf2py intent(out) nspec
      nspec=nspec/3
      do 8 n=1,nspec
        read(10) numzgrid,species(n)
        read(10) numzgrid,species(n)
8       read(10) numzgrid,species(n)


      read(10) numpoint
Cf2py intent(out) numpoint
      do 13 i=1,numpoint
        read(10) ireleasestart(i),ireleaseend(i)

        read(10) xpoint(i),ypoint(i),xp2,yp2,zpoint1(i),zpoint2(i)
        read(10) npart(i),kind(i)
        read(10) compoint(i)
        do 13 j=1,nspec
          read(10)
          read(10)
13        read(10) xmass(i,j)
      read(10) method
Cf2py intent(out) method
      read(10) nageclass,(lage(i),i=1,nageclass)
Cf2py intent(out) nageclass
      lage(0)=0
      write(*,*) (lage(i),i=0,nageclass)


      do 130 ix=1,numxgrid
130     read(10) (oro(ix,jy),jy=1,numygrid)

      close(10)
      write(*,*) (outheight(i),i=1,numzgrid)

      if (loutstep.lt.0) nspec=numpoint


c Calculate height, which is outheight plus topography
******************************************************

      do 150 ix=1,numxgrid
        do 150 jy=1,numygrid
          if (ltopo.eq.1) then
            heightnn (ix,jy,0) = oro(ix,jy)
          else
            heightnn (ix,jy,0) = 0.
          endif
          do 150 i=1,numzgrid
            if (ltopo.eq.1) then
              heightnn (ix,jy,i) = outheight(i) + oro(ix,jy)
            else
              heightnn (ix,jy,i) = outheight(i)
            endif
150         continue


C Determine area of each output grid box
****************************************

      do 140 jy=1,numygrid
        ylata=outlat0+(float(jy-1)+0.5)*dyout
        ylatp=ylata+0.5*dyout
        ylatm=ylata-0.5*dyout
        if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
          hzone=dyout*r_earth*pih
        else
          cosfact=cos(ylata*pih)*r_earth
          cosfactp=cos(ylatp*pih)*r_earth
          cosfactm=cos(ylatm*pih)*r_earth
          if (cosfactp.lt.cosfactm) then
            hzone=sqrt(r_earth**2-cosfactp**2)-
     +      sqrt(r_earth**2-cosfactm**2)
          else
            hzone=sqrt(r_earth**2-cosfactm**2)-
     +      sqrt(r_earth**2-cosfactp**2)
          endif
        endif
        gridarea=2.*pi*r_earth*hzone*dxout/360.
        do 140 ix=1,numxgrid
140       area(ix,jy)=gridarea


      end

      
