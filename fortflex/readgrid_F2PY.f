      subroutine readgrid(filegrid,nxmax,numxgrid,nymax,numygrid,
     +nzmax,numzgrid,nspec,nageclass,scaleconc,decayconstant)

      real grid(nxmax,nymax,nzmax,maxspec,maxageclass)
Cf2py intent(out) grid

      integer lage(0:maxageclass)

***************** new sparse output
      real smallnum
      integer fact,ii,ir
      real sparse_dump_r(nxmax*nymax*nzmax)
      integer sparse_dump_i(nxmax*nymax*nzmax)
      integer sp_count_i,sp_count_r
***************** new sparse output


      character*250 filegrid
Cf2py intent(in) filegrid
Cf2py intent(in) nxmax
Cf2py intent(in) numxgrid
Cf2py intent(in) nymax
Cf2py intent(in) numygrid
Cf2py intent(in)  nzmax
Cf2py intent(in) numzgrid
Cf2py intent(in) scaleconc

      open(10,file=filegrid,form='unformatted',status='old')
      read(10,end=99) itime 
      write(*,*) itime

C Loop about all species
************************

      do 38 k=1,nspec

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
