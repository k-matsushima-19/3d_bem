c
c     merge.f
c
c     2002.06.29.v0
c     2002.06.29.v1 for fast execution introduce mmm,nnn
c     2003.07.17.v2 for EFIE(v0.1-).
c     2003.09.17.v3 for EFIE(v0.2-).
c     2004.08.23.v4 for EFIE(v0.2.1-): writing format was removed
c     2005.02.05.v4 for EFIE(v0.2.1-): writing format was added
c     2005.06.06.v4 for EFIE(v0.2.2-): ibc was removed
*
      program main
      implicit none
      integer max,iunit,maxmlt
      parameter(max=1 000 000,iunit=6,maxmlt=20)
*
      integer i,j,k,l,nb,nnd,nod(3,max),icheck(max),list(max),nlist,
     &     nndnew,nodnew(3,max),mmm(max),nnn(maxmlt,max),itmp,
*050606     &     nd,ibc(max),iblng(2,max)
     &     nd,iblng(2,max)
      real*8 x(3,max),dist,distmin,eps,tmp,xnew(3,max)
c      parameter(eps=1.0e-10)
c      parameter(eps=1.0e-7)
      parameter(eps=1.0e-18)
c
CCC      call loadoff(nb,nnd,nod,x)
CCC      do i=1,nb
CCC         do j=1,3
CCC            nod(j,i)=nod(j,i)+1
CCC         enddo
CCC      enddo
*050606      call loaddat(nb,nnd,nd,nod,ibc,iblng,x)
      call loaddat(nb,nnd,nd,nod,iblng,x)

      if((nnd.gt.max).or.(nb.gt.max))then ! check
         write(0,*)'max is too small; stop'
         stop
      endif

      distmin=1.0d+100
      do i=1,nb
         tmp=dist(x,nod(1,i),nod(2,i))
         if(distmin.gt.tmp)distmin=tmp
         tmp=dist(x,nod(2,i),nod(3,i))
         if(distmin.gt.tmp)distmin=tmp
         tmp=dist(x,nod(3,i),nod(1,i))
         if(distmin.gt.tmp)distmin=tmp
      enddo
      write(0,*)'distmin=',distmin
      if(eps.gt.distmin)then
         write(0,*)' eps is too big; stop'
         stop
      endif

      do i=1,nnd
         mmm(i)=0
      enddo
      do i=1,nnd
         if(mod(i,1000).eq.0)then
            write(0,2)i
 2          format(i7,$)
         endif
         do j=1,nb
            do k=1,3
               if(i.eq.nod(k,j))then
                  mmm(i)=mmm(i)+1
                  if(mmm(i).gt.maxmlt)stop
                  nnn(mmm(i),i)=j
               endif
            enddo
         enddo
      enddo
      write(0,*)


      do i=1,nnd
         icheck(i)=0
      enddo

      nndnew=0
      do i=1,nnd
         if(mod(i,1000).eq.0)then
            write(0,1)i
 1          format('[',i8,']',$)
         endif
         if(icheck(i).eq.0)then !consider a node 'i' which have not checked yet
            nndnew=nndnew+1     !the node is defined as 'nndnew'th one
            do j=1,3
               xnew(j,nndnew)=x(j,i)
            enddo

            do j=1,nnd          !obtain all nodes existing at the same point of node 'i'
               list(j)=0
            enddo
            nlist=0
            do j=1,nnd
               if(icheck(j).eq.0)then !consider only unchecked nodes
                  if(dist(x,i,j).lt.eps)then !node 'j' and 'i' shares the same point
                     nlist=nlist+1 !add 'j' to the list
                     list(nlist)=j
                     icheck(j)=1
                  endif
               endif
            enddo
*v0            do j=1,nlist        !change conectivities which contain nodes listed in 'list'
*v0               do k=1,nb
*v0                  do l=1,3
*v0                     if(list(j).eq.nod(l,k))then
*v0                        nodnew(l,k)=nndnew
*v0                     endif
*v0                  enddo
*v0               enddo
*v0            enddo
            do j=1,nlist
               do k=1,mmm(list(j))
                  itmp=nnn(k,list(j))
                  do l=1,3
                     if(list(j).eq.nod(l,itmp))then
                        nodnew(l,itmp)=nndnew
                     endif
                  enddo
               enddo
            enddo
         endif
      enddo
      write(0,*)
      write(0,*)'nndnew=',nndnew

CCC      do i=1,nb
CCC         do j=1,3
CCC            nodnew(j,i)=nodnew(j,i)-1
CCC         enddo
CCC      enddo
CCC      call wrtpch(nb,nndnew,nodnew,xnew)
*050606      call wrtdat(nb,nndnew,nd,nodnew,ibc,iblng,xnew)
      call wrtdat(nb,nndnew,nd,nodnew,iblng,xnew)
c
      stop
      end

c***********************************************************************
      subroutine loadoff(nb,nnd,nod,x)
c***********************************************************************
      implicit none
      integer nb,nnd,nod(3,*)
      real*8 x(3,*)
*
      character coff*3,infile*64 
      integer i,j,idummy,itmp
      real*8 dummy
c      
      write(0,*)'# OFF file '
      read(*,'(a64)')infile
      open(2,file=infile,status='old')
      read(2,*)coff
      if(coff.ne.'OFF')then
         write(0,*)'bad format; stop'
         stop
      endif
      read(2,*)nnd,nb,idummy
      do i=1,nnd
         read(2,*)(x(j,i),j=1,3)
      enddo
      do i=1,nb
         read(2,*)itmp,(nod(j,i),j=1,3),(dummy,j=1,3)
      enddo
      close(2)
c
      return
      end

c***********************************************************************
      real*8 function dist(x,i,j)
c***********************************************************************
      implicit none
      integer i,j
      real*8 x(3,*)
c
      dist=(x(1,i)-x(1,j))**2+(x(2,i)-x(2,j))**2+(x(3,i)-x(3,j))**2
c
      return
      end

c***********************************************************************
      subroutine wrtpch(nb,nnd,nod,x)
c***********************************************************************
      implicit none
      integer nb,nnd,nod(3,*)
      real*8 x(3,*)
*
      integer i,j,iunit
      parameter(iunit=6)
c      
      write(iunit,10)nnd
      do i=1,nnd
*NOFORMAT         write(iunit,20)(x(j,i),j=1,3)
*050205         write(iunit,*)(x(j,i),j=1,3)
         write(iunit,20)(x(j,i),j=1,3)
      enddo
      write(iunit,10)nb
      do i=1,nb
         write(iunit,30)(nod(j,i),j=1,3)
      enddo
 10   format(i8)
*050205 20   format(3e15.7)
 20   format(3e24.15e3)
 30   format(3i8)
c
      return
      end

c***********************************************************************
*050606      subroutine loaddat(nb,nnd,nd,nod,ibc,iblng,x)
      subroutine loaddat(nb,nnd,nd,nod,iblng,x)
c***********************************************************************
      implicit none
*050606      integer nb,nnd,nd,nod(3,*),ibc(*),iblng(2,*)
      integer nb,nnd,nd,nod(3,*),iblng(2,*)
      real*8 x(3,*)
*
      character infile*64 
      integer i,j,itmp
c      
      write(0,*)'# dat file '
      read(*,'(a64)')infile
c$$$      write(*,'(a64)')infile
c$$$      write(*,*) "aaa"
c$$$      stop
      open(2,file=infile,status='old')
      read(2,*)nb,nnd,nd
      do i=1,nb
         read(2,*)itmp,(nod(j,itmp),j=1,3),
*050606     &        ibc(itmp),(iblng(j,itmp),j=1,2)
     &        (iblng(j,itmp),j=1,2)
      enddo
      do i=1,nnd
         read(2,*)itmp,(x(j,itmp),j=1,3)
      enddo
      close(2)
c
      return
      end

c***********************************************************************
*050606      subroutine wrtdat(nb,nnd,nd,nod,ibc,iblng,x)
      subroutine wrtdat(nb,nnd,nd,nod,iblng,x)
c***********************************************************************
      implicit none
*050606      integer nb,nnd,nd,nod(3,*),ibc(*),iblng(2,*)
      integer nb,nnd,nd,nod(3,*),iblng(2,*)
      real*8 x(3,*)
*
      integer i,j,iunit
      parameter(iunit=6)
c      
      write(iunit,10)nb,nnd,nd
      do i=1,nb
         write(iunit,20)i,(nod(j,i),j=1,3),
*050606     &        ibc(i),(iblng(j,i),j=1,2)
     &        (iblng(j,i),j=1,2)
      enddo
      do i=1,nnd
         write(iunit,30)i,(x(j,i),j=1,3)
      enddo
 10   format(3i8)
 20   format(7i8)
 30   format(i8,3e24.15e3)
c
      return
      end
