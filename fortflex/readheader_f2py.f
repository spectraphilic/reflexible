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
