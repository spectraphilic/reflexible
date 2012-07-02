      subroutine reademissions(nemission,outlon0,outlat0,numxgrid,
     +numygrid,dxout,dyout,nxmax,nymax,nzmax,
     +emissions,index_cont)
C   
C     This file is modified by jfb to create a module for reading in the
C     flexpart emissions files.
C     
C     It was used in combination with gen_flexpart_netcdf_emission.py
C     to create netcdf emissions for the plotting routines
C
C
C     This file can be compiled into a python module:
C     f2py -c -m readflexemissions reademissions_f2py.f
C
C      real area(nxmax,nymax)
C      real heightnn(nxmax,nymax,0:nzmax)
Cf2py intent(in) nemission
Cf2py intent(in) nxmax
Cf2py intent(in) nymax
Cf2py intent(in) outlon0
Cf2py intent(in) outlat0
Cf2py intent(in) numxgrid
Cf2py intent(in) numygrid
Cf2py intent(in) dxout
Cf2py intent(in) dyout
Cf2py intent(in) nxmax
Cf2py intent(in) nymax
Cf2py intent(in) nzmax
C     intent(in) area
C     intent(in) heightnn
Cf2py intent(out) emissions
Cf2py intent(out) index_cont
      real emissions(nxmax,nymax)
      integer index_cont(nxmax,nymax)
      parameter(pi=3.14159265,r_earth=6.371e6,pih=pi/180.)
      parameter(maxemissions=4,nxem=360,nyem=180)
      real edgar_emissions(nxem,nyem,maxemissions)
      real emission_tami(360,180)
      integer indcont(nxem,nyem)

      parameter(maxpoint=1000000)
      real xlonl(maxpoint),xlonr(maxpoint),ylatl(maxpoint)
      real ylatr(maxpoint),zl(maxpoint),zr(maxpoint)
      real emission_highres(maxpoint,maxemissions)

      double precision juldate,julyear
      character*40 edgarname

      real weightmolar(0:maxemissions)
      data weightmolar/28.97,28.,46.,64.,12.01/

C Determine fraction of emissions released per second
*****************************************************

      julyear=juldate(19960101,0)-juldate(19950101,0)
      fract=1./julyear/86400.


*******************************
C Read EDGAR emission inventory
*******************************

C Give indices to the different continents and set emissions to zero
********************************************************************

      do 710 i=1,nxem
        do 710 j=1,nyem
          indcont(i,j)=0
          if ((j.gt.102).and.(i.gt.10).and.(i.le.130)) then
            indcont(i,j)=1       ! N-America
            if ((i.gt.40).and.(i.le.130).and.(j.gt.114).and.(j.le.142))
     +      indcont(i,j)=7       ! Use American regional inventory
            if ((i.ge.81).and.(i.le.82).and.(j.ge.109).and.(j.le.110))
     +      indcont(i,j)=8       ! Use Mexico City inventory
          else if ((j.gt.126).and.(i.gt.155).and.(i.le.240)) then
            indcont(i,j)=2       ! Europe
            if (((i.le.207).or.(j.gt.132)).and.
     +      ((i.le.232).or.(j.gt.142))) then
              indcont(i,j)=2       ! Europe
            else
              indcont(i,j)=3       ! Asia
            endif
          else if (((j.gt.80).and.(i.gt.240).and.(i.le.350)).or.
     +    ((j.gt.105).and.(j.le.126).and.
     +    (i.gt.214).and.(i.le.240))) then
            indcont(i,j)=3       ! Asia
          else
            indcont(i,j)=0       ! No continent
          endif

          if ((j.lt.102).and.(i.gt.95).and.(i.le.150)) then
            indcont(i,j)=4       ! S-America
          else if ((j.le.126).and.(i.gt.160).and.(i.le.232)) then
            if ((i.gt.192).and.(j.gt.112)) then
c             indcont(i,j)=0
              continue
            else if ((i.gt.214).and.(j.gt.110)) then
c             indcont(i,j)=0
              continue
            else if ((i.gt.222).and.(j.gt.104)) then
c             indcont(i,j)=0
              continue
            else
              indcont(i,j)=5       ! Africa
            endif
          else if ((i.gt.290).and.(j.gt.50).and.(j.lt.80)) then
            indcont(i,j)=6       ! Australia
          else
c           indcont(i,j)=0
            continue
          endif
          do 710 n=1,maxemissions
710         edgar_emissions(i,j,n)=0.


C Read EDGAR inventory
**********************

        n=nemission
        if (n.eq.1) then
          open(10,file='CATEGORIES.co')
        else if (n.eq.2) then
          open(10,file='CATEGORIES.nox')
        else if (n.eq.3) then
          open(10,file='CATEGORIES.so2')
        else if (n.eq.4) then
          open(10,file='CATEGORIES.bc')
        endif

200       read(10,'(a)',end=199) edgarname       !which categories are to be used

          open(20,file='edgarv32ft_2000/'//edgarname)

          if (n.le.3) then    ! take EDGAR inventory

            do 11 i=1,12
11            read(20,*)

100           read(20,*,end=99) igrid,jgrid,emcell
              i=igrid+181
              j=jgrid+91
              edgar_emissions(i,j,n)=edgar_emissions(i,j,n)+emcell
              goto 100
99          continue
            close(20)

          else if (n.eq.4) then ! take Tami Bond's BC inventory

            do 12 i=1,9
12            read(20,*)
            read(20,*) ((emission_tami(i,j),i=1,360),j=1,180)
            close(20)
C Convert emissions from Gg to kg
            do 13 i=1,360
              do 13 j=1,180
13              edgar_emissions(i,j,n)=edgar_emissions(i,j,n)+
     +          emission_tami(i,j)*1.e6

          endif

          goto 200
199     continue
        close(10)

C Convert kg/year -> kg/s
*************************

      do 25 j=1,nyem
        do 25 i=1,nxem
            edgar_emissions(i,j,nemission)=
     +      edgar_emissions(i,j,nemission)*fract
25          continue


      write(*,*) 'End reading EDGAR emission data'

C Read regional North American inventory
****************************************

        open(20,file='american_emission_inventory.dat')
        l=1
300      read(20,*,end=299) xlonl(l),xlonr(l),ylatl(l),ylatr(l),
     +    zl(l),zr(l),(emission_highres(l,n),n=1,maxemissions)
C Convert kg/year -> kg/s
          do 60 n=1,maxemissions
60          emission_highres(l,n)=emission_highres(l,n)*fract
          l=l+1
          goto 300
299     nhighres=l-1
        close(20)

C Read Mexico City inventory
****************************

        open(20,file='mexico_city_inventory.dat')
        l=nhighres+1
310      read(20,*,end=309) xlonl(l),xlonr(l),ylatl(l),ylatr(l),
     +    zl(l),zr(l),(emission_highres(l,n),n=1,maxemissions)
C Convert kg/year -> kg/s
          do 66 n=1,maxemissions
66          emission_highres(l,n)=emission_highres(l,n)*fract
          l=l+1
          goto 310
309     nhighres=l-1
        close(20)

C Read European inventory
*************************

        open(20,file='european_emission_inventory.dat')
        l=nhighres+1
320      read(20,*,end=319) xlonl(l),xlonr(l),ylatl(l),ylatr(l),
     +    zl(l),zr(l),(emission_highres(l,n),n=1,3)
C Convert kg/year -> kg/s
          do 67 n=1,3
67          emission_highres(l,n)=emission_highres(l,n)*fract
          l=l+1
          goto 320
319     nhighres=l-1
        close(20)


C Now attribute emissions to model grid cells
*********************************************

      do 10 ix=1,numxgrid
        xl=outlon0+float(ix-1)*dxout
        xr=outlon0+float(ix)*dxout
        xm=(xl+xr)/2.
        do 10 jy=1,numygrid
          emissions(ix,jy)=0.

          yl=outlat0+float(jy-1)*dyout
          yr=outlat0+float(jy)*dyout
          ym=(yl+yr)/2.

          indx=int(xm+180.)+1
          indx=mod(indx-1,360)+1
          indy=int(ym+90.)+1
          index_cont(ix,jy)=indcont(indx,indy)
          neur=0
          if ((indx.ge.190).and.(indx.le.229).and.(indy.ge.128).and.
     +    (indy.le.161)) neur=1     ! use EMEP European inventory
          if ((nemission.gt.3).or.((index_cont(ix,jy).le.6)
     +    .and.(neur.eq.0))) then
            emissions(ix,jy)=edgar_emissions(indx,indy,nemission)*
     +      dxout*dyout
          else        ! use N. American, European, and Mexico City inventory
            do 20 l=1,nhighres
             if ((xlonl(l).eq.xlonr(l)).and.(ylatl(l).eq.ylatr(l))) then   ! point sources
                if ((xlonl(l).ge.xl).and.(xlonl(l).lt.xr).and.
     +              (ylatl(l).ge.yl).and.(ylatl(l).lt.yr))
     +          emissions(ix,jy)=emissions(ix,jy)+
     +          emission_highres(l,nemission)
              else              ! area sources
                xll=max(xl,xlonl(l))
                yll=max(yl,ylatl(l))
                xrr=min(xr,xlonr(l))
                yrr=min(yr,ylatr(l))
                dxl=xrr-xll
                dyl=yrr-yll
                if ((dxl.gt.0.).and.(dyl.gt.0.)) then
                  areafract=dxl*dyl/(xlonr(l)-xlonl(l))/
     +            (ylatr(l)-ylatl(l))
                  emissions(ix,jy)=emissions(ix,jy)+
     +            emission_highres(l,nemission)*areafract
                endif
              endif
20            continue
          endif
10        continue

C Convert kg/s -> kg/m3/s
C Also account for molar weight to convert to volume mixing ratio
C      do 61 ix=1,numxgrid
C        do 61 jy=1,numygrid
C61        emissions(ix,jy)=emissions(ix,jy)/
C     +    (heightnn(ix,jy,1)-heightnn(ix,jy,0))/area(ix,jy)*
C     +    weightmolar(0)/weightmolar(nemission)


      return
      end

      FUNCTION JULDATE(YYYYMMDD,HHMISS)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                             *
*     Calculates the Julian date                                              *
*                                                                             *
*     AUTHOR: Andreas Stohl (15 October 1993)                                 *
*                                                                             *
*     Variables:                                                              *
*     DD             Day                                                      *
*     HH             Hour                                                     *
*     HHMISS         Hour, minute + second                                    *
*     JA,JM,JY       help variables                                           *
*     JULDATE        Julian Date                                              *
*     JULDAY         help variable                                            *
*     MI             Minute                                                   *
*     MM             Month                                                    *
*     SS             Second                                                   *
*     YYYY           Year                                                     *
*     YYYYMMDDHH     Date and Time                                            *
*                                                                             *
*     Constants:                                                              *
*     IGREG          help constant                                            *
*                                                                             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      IMPLICIT NONE

      INTEGER YYYYMMDD,YYYY,MM,DD,HH,MI,SS,HHMISS
      INTEGER JULDAY,JY,JM,JA,IGREG
      DOUBLE PRECISION JULDATE
      PARAMETER (IGREG=15+31*(10+12*1582))

      YYYY=YYYYMMDD/10000
      MM=(YYYYMMDD-10000*YYYY)/100
      DD=YYYYMMDD-10000*YYYY-100*MM
      HH=HHMISS/10000
      MI=(HHMISS-10000*HH)/100
      SS=HHMISS-10000*HH-100*MI

      IF (YYYY.EQ.0) PAUSE 'There is no Year Zero.'
      IF (YYYY.LT.0) YYYY=YYYY+1
      IF (MM.GT.2) THEN
        JY=YYYY
        JM=MM+1
      ELSE
        JY=YYYY-1
        JM=MM+13
      ENDIF
      JULDAY=INT(365.25*JY)+INT(30.6001*JM)+DD+1720995
      IF (DD+31*(MM+12*YYYY).GE.IGREG) THEN
        JA=INT(0.01*JY)
        JULDAY=JULDAY+2-JA+INT(0.25*JA)
      ENDIF

      JULDATE=DBLE(FLOAT(JULDAY))+DBLE(FLOAT(HH)/24.)+
     +DBLE(FLOAT(MI)/1440.)+DBLE(FLOAT(SS)/86400.)

      END
