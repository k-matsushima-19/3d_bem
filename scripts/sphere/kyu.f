c
c     kyu.f
c     
c     1999.10. 3 for difu3.f(v0.0-)
c     1999.10. 3 for pot3.f(v0.4-)
c     2000. 1.22 for fmel.f(v0.5-)
c     2000.11.13 for dyna3d.f(v0.0-)
c     2002. 1.15 for fmel.f
c     2002. 3. 4 added iu1,iu2,iu3,it1,it2,it3
c     2003. 5.21 for efie3d.f(v0.0-)
c     2003. 9.17 for efie3d.f(v0.2-)
c     --- outer domain is '1(+)', while inner domain is '2(-)'. Then iflag must be '0'.
c     2004. 5.17 for hel3d.f(v0.9.1-), 'session' was added.
c     2004. 6.19 based on kyu2.f. 'mknode2' and 'mkedge' were introduced.
c     2005. 2. 5 for efie3d.f(v0.2.1-)
c     2005. 6. 6 for efie3d.f(v0.2.2-): the column for ibc was removed

      program kyu
      implicit none
      integer id,i,ie,in,ih,nod(3,20),iflag
      real*8 rad,x(3,12),center(3)

*     session

      write(0,*)'iflag(0:inward,1:outward), id, rad, center'
      read(*,*)iflag,id,rad,center(1),center(2),center(3)
      write(0,*)'iflag=',iflag,' id=',id,' rad=',rad,
     &     ' center=',(center(i),i=1,3)

*     header

      write(*,10)20*id**2,20*(id+2)*(id+1)/2,2
 10   format(3i8)
      write(0,*)' nb=',20*id**2
      write(0,*)'nnd=',20*(id+2)*(id+1)/2

      if(iflag.eq.1)then
         write(0,*)'OUTWARD NORMAL TO SPHERE'
      else
         write(0,*)'INWARD NORMAL TO SPHERE'
      endif

*     body

      ie=1
      ih=1
      do i=1,20
         call mkelem(rad,id,ie,ih,iflag)
      enddo
      call lddelta20(nod,x)
      in=1
      do i=1,20
*040618         call mknode(rad,id,in,nod,x,i)
         call mknode2(rad,id,in,nod,x,i,center)
      enddo
      call wrtpara(id,rad,iflag)
c
      stop
      end

c***********************************************************************
      subroutine mkelem(rad,id,ie,ih,iflag)
c***********************************************************************
      implicit none
      integer id,ie,ih,iflag
      real*8 rad
*
      integer i,j,k,l
      real*8 tmp(3)
c
      do i=1,id
         do j=1,i
            if(iflag.eq.1)then  ! INTERIOR
               write(*,30)ie,
     &              (ih  )+(j-1), ! INNER
     &              (ih+i)+(j-1),
     &              (ih+i)+(j  ),
*050606     &              1,1,2
     &              1,2
            else                ! EXTERIOR
               write(*,30)ie,
     &              (ih  )+(j-1),
     &              (ih+i)+(j  ),
     &              (ih+i)+(j-1),
*050606     &              1,1,2
     &              1,2
            endif
            ie=ie+1
            if(j.ne.i)then
               if(iflag.eq.1)then ! OUTWARD
                  write(*,30)ie,
     &                 (ih  )+(j-1),
     &                 (ih+i)+(j  ),
     &                 (ih  )+(j  ),
*050606     &                 1,1,2
     &                 1,2
               else             ! INWARD
                  write(*,30)ie,
     &                 (ih  )+(j-1),
     &                 (ih  )+(j  ),
     &                 (ih+i)+(j  ),
*050606     &                 1,1,2
     &                 1,2
               endif
               ie=ie+1
            endif
         enddo
         ih=ih+i
      enddo
      ih=ih+id+1
 30   format(10i10)
c
      return
      end

c***********************************************************************
      subroutine mknode(rad,id,in,nod,x,icase)
c***********************************************************************
      implicit none
      integer id,in,icase,nod(3,20)
      real*8 rad,x(3,12)
*
      integer i,j,k
      real*8 vec12(3),vec23(3),tmp1,tmp2,hen,p(3)
c
      do i=1,3
         vec12(i)=x(i,nod(2,icase))-x(i,nod(1,icase))
         vec23(i)=x(i,nod(3,icase))-x(i,nod(2,icase))
      enddo
      tmp1=dsqrt(dabs(vec12(1)**2+vec12(2)**2+vec12(3)**2))
      tmp2=dsqrt(dabs(vec23(1)**2+vec23(2)**2+vec23(3)**2))
      do i=1,3
         vec12(i)=vec12(i)/tmp1
         vec23(i)=vec23(i)/tmp2
      enddo

      hen=tmp1/dble(id)
      do i=0,id
         do j=0,i
            do k=1,3
               p(k)=x(k,nod(1,icase))+dble(i)*hen*vec12(k)
     &           +dble(j)*hen*vec23(k)
            enddo
            call henkan(p(1),p(2),p(3),rad)
            write(*,50)in,(p(k),k=1,3)
            in=in+1
         enddo
      enddo
*050205 50   format(i7,3e15.7)
 50   format(i8,3e24.15e3)
c
      return
      end

c***********************************************************************
      subroutine henkan(x,y,z,rad)
c***********************************************************************
      implicit none
      real*8 x,y,z,rad
*
      real*8 tmp
c
      tmp=dsqrt(dabs(x**2+y**2+z**2))
      x=x/tmp*rad
      y=y/tmp*rad
      z=z/tmp*rad
c     
      return
      end
      
c***********************************************************************
      subroutine lddelta20(nod,x)
c***********************************************************************
* Refer to /usr/local/Geomview/data/geom/polyhedra/icosahedron 
*
      implicit none
      integer nod(3,20)
      real*8 x(3,12)
*     
      integer i,j
c
      x(1, 1)=  0.0000000E+00
      x(2, 1)=  0.0000000E+00
      x(3, 1)=  0.2000000E+01
      x(1, 2)=  0.1788854E+01
      x(2, 2)=  0.0000000E+00
      x(3, 2)=  0.8944270E+00
      x(1, 3)=  0.5527860E+00
      x(2, 3)=  0.1701302E+01
      x(3, 3)=  0.8944270E+00
      x(1, 4)= -0.1447214E+01
      x(2, 4)=  0.1051462E+01
      x(3, 4)=  0.8944270E+00
      x(1, 5)= -0.1447214E+01
      x(2, 5)= -0.1051462E+01
      x(3, 5)=  0.8944270E+00
      x(1, 6)=  0.5527860E+00
      x(2, 6)= -0.1701302E+01
      x(3, 6)=  0.8944270E+00
      x(1, 7)=  0.1447214E+01
      x(2, 7)=  0.1051462E+01
      x(3, 7)= -0.8944270E+00
      x(1, 8)= -0.5527860E+00
      x(2, 8)=  0.1701302E+01
      x(3, 8)= -0.8944270E+00
      x(1, 9)= -0.1788854E+01
      x(2, 9)=  0.0000000E+00
      x(3, 9)= -0.8944270E+00
      x(1,10)= -0.5527860E+00
      x(2,10)= -0.1701302E+01
      x(3,10)= -0.8944270E+00
      x(1,11)=  0.1447214E+01
      x(2,11)= -0.1051462E+01
      x(3,11)= -0.8944270E+00
      x(1,12)=  0.0000000E+00
      x(2,12)=  0.0000000E+00
      x(3,12)= -0.2000000E+01
      do i=1,12
         do j=1,3
            x(j,i)=x(j,i)*0.50d0
         enddo
      enddo
      nod(1, 1)= 3
      nod(2, 1)= 1
      nod(3, 1)= 2
      nod(1, 2)= 4
      nod(2, 2)= 1
      nod(3, 2)= 3
      nod(1, 3)= 5
      nod(2, 3)= 1
      nod(3, 3)= 4
      nod(1, 4)= 6
      nod(2, 4)= 1
      nod(3, 4)= 5
      nod(1, 5)= 2
      nod(2, 5)= 1
      nod(3, 5)= 6
      nod(1, 6)= 3
      nod(2, 6)= 2
      nod(3, 6)= 7
      nod(1, 7)= 8
      nod(2, 7)= 3
      nod(3, 7)= 7
      nod(1, 8)= 4
      nod(2, 8)= 3
      nod(3, 8)= 8
      nod(1, 9)= 9
      nod(2, 9)= 4
      nod(3, 9)= 8
      nod(1,10)= 5
      nod(2,10)= 4
      nod(3,10)= 9
      nod(1,11)=10
      nod(2,11)= 5
      nod(3,11)= 9
      nod(1,12)= 6
      nod(2,12)= 5
      nod(3,12)=10
      nod(1,13)=11
      nod(2,13)= 6
      nod(3,13)=10
      nod(1,14)= 7
      nod(2,14)= 2
      nod(3,14)=11
      nod(1,15)= 2
      nod(2,15)= 6
      nod(3,15)=11
      nod(1,16)= 7
      nod(2,16)=12
      nod(3,16)= 8
      nod(1,17)= 8
      nod(2,17)=12
      nod(3,17)= 9
      nod(1,18)= 9
      nod(2,18)=12
      nod(3,18)=10
      nod(1,19)=10
      nod(2,19)=12
      nod(3,19)=11
      nod(1,20)=11
      nod(2,20)=12
      nod(3,20)= 7
c
      return
      end
         
c***********************************************************************
      subroutine wrtpara(id,rad,iflag)
c***********************************************************************
      implicit none
      integer id,iflag
      real*8 rad
c
c$$$      write(*,10)
c$$$      write(*,20)'id',id
c$$$      write(*,30)'rad',rad
c$$$      write(*,20)'iflag',iflag
*050205      write(*,20)'iu',iu
*050205      write(*,20)'iq',iq
*040618 10   format('# created by kyu2.f')
 10   format('# created by kyu.f')
 20   format('# ',a8,'=',i8)
 30   format('# ',a8,'=',e15.7)
c
      return
      end

c***********************************************************************
      subroutine mknode2(rad,id,in,nod,x,icase,center)
c***********************************************************************
      implicit none
      integer id,in,icase,nod(3,20)
      real*8 rad,x(3,12),center(3)
*
      integer i,j,k
      real*8 v1(3),v2(3),v3(3),u1(3),u2(3),p(3)
c
      do i=1,3
         v1(i)=x(i,nod(1,icase))
         v2(i)=x(i,nod(2,icase))
         v3(i)=x(i,nod(3,icase))
      enddo

      do i=0,id
         call getedge(id,rad,i,v1,v2,u1)
         call getedge(id,rad,i,v1,v3,u2)
         do j=0,i
            call getedge(i,rad,j,u1,u2,p)
            write(*,50)in,(p(k)+center(k),k=1,3)
            in=in+1
         enddo
      enddo
*050205 50   format(i7,3e15.7)
 50   format(i8,3e24.15e3)
c
      return
      end

c***********************************************************************
      subroutine getedge(nmax,rad,n,s,e,v)
c***********************************************************************
      implicit none
      integer nmax,n
      real*8 rad,s(3),e(3),v(3)
*
      integer i
      real*8 ds,de,dn,vs(3),ve(3),vn(3),vp(3),angle,rot,cr,sr,eps
      parameter(eps=1.0d-12)
c
      if(nmax.eq.0)then         ! s=e=v
         do i=1,3
            v(i)=s(i)
         enddo
      else
         ds=dsqrt(s(1)**2+s(2)**2+s(3)**2)
         de=dsqrt(e(1)**2+e(2)**2+e(3)**2)
         do i=1,3
            vs(i)=s(i)/ds
            ve(i)=e(i)/de
         enddo
         vn(1)=vs(2)*ve(3)-vs(3)*ve(2)
         vn(2)=vs(3)*ve(1)-vs(1)*ve(3)
         vn(3)=vs(1)*ve(2)-vs(2)*ve(1)
         dn=dsqrt(vn(1)**2+vn(2)**2+vn(3)**2)
         do i=1,3
            vn(i)=vn(i)/dn
         enddo
         vp(1)=vn(2)*vs(3)-vn(3)*vs(2)
         vp(2)=vn(3)*vs(1)-vn(1)*vs(3)
         vp(3)=vn(1)*vs(2)-vn(2)*vs(1)
         angle=dacos(vs(1)*ve(1)+vs(2)*ve(2)+vs(3)*ve(3))
         rot=angle*dble(n)/dble(nmax)
         cr=dcos(rot)
         sr=dsin(rot)
         do i=1,3
            v(i)=rad*(cr*vs(i)+sr*vp(i))
         enddo
      endif
c
      return
      end
