module basic
implicit none
private
public :: datack, bisec, vec12, orthg

contains
!************* data check of pairs of sites ************
subroutine datack(ipair,ibond,n)
implicit none
integer, intent(in) :: ipair(ibond*2), ibond, n
integer :: k, isite1, isite2
do k=1,ibond
  isite1=ipair(k*2-1)
  isite2=ipair(k*2  )
  if((isite1.le.0).or.(isite2.le.0).or.(isite1.gt.n).or.(isite2.gt.n))then
            print *,' #(E03)# Incorrect data in ipair'
            print *,'         Location :  ',k*2-1,k*2
            stop
  end if
enddo
end subroutine datack
!!
!********* eigenvalues by the bisection method **********
!
subroutine bisec(alpha,beta,ndim,E,ne,eps)
!
!    alpha  @ diagonal element
!    beta   @ subdiagonal element
!    ndim   @ matrix dimension
!    E      # eigenvalues
!    ne     @ number of eigenvalues to calculate
!    eps    @ limit of error
!
implicit none
integer, intent(in) :: ne, ndim
double precision, intent(in) :: eps
double precision, intent(inout) :: E(:)
double precision, intent(inout) :: alpha(:),beta(:)
double precision :: b2(8000)
!
integer :: k, i, j, numneg, ipass
double precision :: a, c, g, range, epsabs, b
!
if(ndim > 8000)then
          print *,' #(E04)# ndim given to bisec exceeds 8000'
          stop
end if
if( (ne>ndim) .or. (ne<0) )then
          print *,' #(E05)# ne given to bisec out of range',ne
          stop
end if
!!
!*** initial bound
range=dabs(alpha(1))+dabs(beta(1))
do k=2,ndim-1
  range=dmax1(range,dabs(beta(k-1))+dabs(alpha(k))+dabs(beta(k)))
enddo
range=dmax1(range,dabs(beta(ndim-1))+dabs(alpha(ndim)))
range=-range
!!
b2(1)=0.d0
do i=2,ndim
  b2(i)=beta(i-1)**2
enddo
!!
epsabs=dabs(range)*eps
do i=1,ne
  E(i)=-range
enddo
b=range
!!
!*** bisection method
do k=1,ne
  a=E(k)
  do j=1,100
    c=(a+b)/2.d0
    if(dabs(a-b)<epsabs) cycle
    numneg=0
    g=1.d0
    ipass=0
    do i=1,ndim
      if(ipass==0)then
        g=c-alpha(i)-b2(i)/g
      else if(ipass.eq.1)then
        ipass=2
      else
        g=c-alpha(i)
        ipass=0
      end if
!!
      if(ipass==0)then
        if(g<0.d0) numneg=numneg+1
        if(dabs(g)<dabs(b2(i)*epsabs*eps)) ipass=1
      end if
    enddo
    numneg=ndim-numneg
    if(numneg<k)then
      b=c
    else
      a=c
      do i=k,min(numneg,ne)
        E(i)=c
      enddo
    end if
  enddo
enddo
end subroutine bisec
!!
!*** eigenvector of a tridiagonal matrix by inverse iteration ***
!                 for the large/medium routines
!
subroutine vec12(alpha,beta,v,E,ndim,nvec,di,bl,bu,bv,cm,lex)
!
!    E(4)       @  4 lowest eigenvalues
!    ndim       @  matrix dimension
!    nvec       @  number of vectors to calculate
!    di - lex      working areas
!
implicit none
integer, intent(in) :: ndim, nvec
double precision, intent(inout) :: alpha(:),beta(:),v(:,:)
double precision, intent(inout) :: E(:)
double precision, intent(inout) :: di(:),bl(:),bu(:),bv(:),cm(:)
integer, intent(inout) :: lex(:)
integer :: k, j, i, l, m, km
double precision :: s, prd, dnorm
!
do k=1,nvec
  do j=1,ndim
    di(j)=E(k)-alpha(j)
    bl(j)=-beta(j)
    bu(j)=-beta(j)
  enddo 
!
!*** LU decomposition
  do j=1,ndim-1
    if(dabs(di(j))>dabs(bl(j)))then
!--- non pivoting
      lex(j)=0
      if(abs(di(j)).lt.1.d-13) di(j)=1.d-13
      cm(j+1)=bl(j)/di(j)
      di(j+1)=di(j+1)-cm(j+1)*bu(j)
      bv(j)=0.d0
    else
!--- pivoting
      lex(j)=1
      cm(j+1)=di(j)/bl(j)
      di(j)=bl(j)
      s=bu(j)
      bu(j)=di(j+1)
      bv(j)=bu(j+1)
      di(j+1)=s-cm(j+1)*bu(j)
      bu(j+1)= -cm(j+1)*bv(j)
    end if
  enddo
  if(dabs(di(ndim))<1.d-13)  di(ndim)=1.d-13
!
!--- initial vector
  do j=1,ndim
    v(j,k)=1.d0/(dble(j)*5.d0)
  enddo
!
!*** degeneracy check up
  if(k==1)then
    km=k
  else if(dabs(E(k)-E(km))>1.d-13)then
    km=k
  else
    do i=km,k-1
      prd=0.d0
      do j=1,ndim
        prd=prd+v(j,i)*v(j,k)
      enddo
      do j=1,ndim
        v(j,k)=v(j,k)-prd*v(j,i)
      enddo
    enddo
  end if
!
!*** inverse iteration
  do l=1,k-km+3
    if((l.ne.1).or.(k.ne.km))then
!--- forward substitution
      do j=1,ndim-1
        if(lex(j)==0)then
          v(j+1,k)=v(j+1,k)-cm(j+1)*v(j,k)
        else
          s=v(j,k)
          v(j,k)=v(j+1,k)
          v(j+1,k)=s-cm(j+1)*v(j,k)
        end if
      enddo
    end if
!--- backward substitution
    do j=ndim,1,-1
      s=v(j,k)
      if(j.le.ndim-1)s=s-bu(j)*v(j+1,k)
      if(j.le.ndim-2)s=s-bv(j)*v(j+2,k)
      v(j,k)=s/di(j)
    enddo
!
!*** normalization
    dnorm=0.d0
    do j=1,ndim
      dnorm=dnorm+v(j,k)**2
    enddo
    if(dnorm>1.d-13)dnorm=1.d0/dsqrt(dnorm)
    do j=1,ndim
      v(j,k)=v(j,k)*dnorm
    enddo
  enddo
enddo
end subroutine vec12
!
!
!!
!******* Orthogonalization of the eigenvectors ************
!
subroutine orthg(idim,ideclr,ev,norm,idgn,numvec)
!
!   idim    @  matrix dimension
!   ideclr  @  declared array size in the main program
!   ev      @# vectors to be orthogonalized / orthogonalized vectors
!   norm(j) #  norm of the j-th vector returned
!   idgn    #  degree of degenearcy
!   numvec  @  number of vectors to be checked
!
implicit none
integer, intent(in) :: idim, ideclr
integer, intent(inout) :: idgn, numvec
double precision, intent(inout) :: ev(ideclr,numvec)
integer, intent(inout) :: norm(numvec)

integer :: i, j, k, l, m
double precision :: dnorm, prjct, vnorm, prd
	  
if(numvec.le.1)then
         print *,' #(W03)# Number of vectors is less than 2 in orthg'
         return
end if
do i=1,numvec
  dnorm=0.0d0
  do j=1,idim
    dnorm=dnorm+ev(j,i)**2
  enddo
  if(dnorm<1.0d-20)then
           print *,' #(W04)# Null vector given to orthg. Location is',i
           return
  end if
  dnorm=1.0d0/dsqrt(dnorm)
  do j=1,idim
    ev(j,i)=ev(j,i)*dnorm
  enddo
enddo

idgn=numvec
norm(1)=1
!*** orthogonalization
do i=2,numvec
  norm(i)=1
  do j=1,i-1
    prjct=0.0d0
    do l=1,idim
      prjct=prjct+ev(l,i)*ev(l,j)
    enddo
    do l=1,idim
      ev(l,i)=ev(l,i)-prjct*ev(l,j)
    enddo
  enddo
  vnorm=0.0d0
  do l=1,idim
    vnorm=vnorm+ev(l,i)**2
  enddo
  if(vnorm>1.0d-15)then
    vnorm=1.0d0/dsqrt(vnorm)
    do l=1,idim
      ev(l,i)=ev(l,i)*vnorm
    enddo
  else
    do l=1,idim
      ev(l,i)=0.0d0
    enddo
    idgn=idgn-1
    norm(i)=0
  end if
enddo
!*** check orthogonality
do i=2,numvec
  do j=1,i-1
    prd=0.0d0
    do l=1,idim
      prd=prd+ev(l,i)*ev(l,j)
    enddo
    if(dabs(prd)<1.0d-13) exit
        print 200,i,j
200     format(' #(W05)# Non-orthogonal vectors at',2i4)
        print 210,prd
210     format('         Overlap : ',d14.7)
        print *,'         Unsuccessful orthogonalization'
        return
  enddo
enddo
end subroutine orthg
end module basic
