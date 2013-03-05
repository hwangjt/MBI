subroutine evaluatePt(d1, d2, na, nx, nf, nC, nP, ks, ms, t, C, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) d1, d2, na, nx, nf, nC, nP, ks, ms, t, C
  !f2py intent(out) P
  !f2py depend(nx) ks, ms
  !f2py depend(nP,nx) t
  !f2py depend(nC, nf) C
  !f2py depend(nP,nf) P

  !Input
  integer, intent(in) ::  d1, d2, na, nx, nf, nC, nP, ks(nx), ms(nx)
  double precision, intent(in) ::  t(nP,nx), C(nC,nf)

  !Output
  double precision, intent(out) ::  P(nP,nf)

  !Working
  integer k, m
  integer iC
  integer i, ix, ia, ik, i0, rem, iP
  double precision Ba(nP,na)
  integer Bj(nP,na)
  double precision, allocatable, dimension(:) ::  d, B

  Ba(:,:) = 1.0
  Bj(:,:) = 1
  do ix=1,nx
     k = ks(ix)
     m = ms(ix)
     allocate(d(k+m))
     allocate(B(k))
     do iP=1,nP
        call knotopen(k, k+m, d)
        if ((d1 .ne. ix) .and. (d2 .ne. ix)) then
           call basis(k, k+m, t(iP,ix), d, B, i0)
        else if ((d1 .eq. ix) .and. (d2 .ne. ix)) then
           call basis1(k, k+m, t(iP,ix), d, B, i0)
        else if ((d1 .ne. ix) .and. (d2 .eq. ix)) then
           call basis1(k, k+m, t(iP,ix), d, B, i0)
        else if ((d1 .eq. ix) .and. (d2 .eq. ix)) then
           call basis2(k, k+m, t(iP,ix), d, B, i0)
        end if
        do ia=1,na
           rem = ia - 1
           do i=nx,ix,-1
              rem = mod(rem,product(ks(:i)))
           end do
           ik = rem/product(ks(:ix-1)) + 1
           iC = (ik + i0 - 1)*product(ms(:ix-1))
           Ba(iP,ia) = Ba(iP,ia) * B(ik)
           Bj(iP,ia) = Bj(iP,ia) + iC
        end do
     end do
     deallocate(d)
     deallocate(B)
  end do

  P(:,:) = 0.0
  do ia=1,na
     do iP=1,nP
        P(iP,:) = P(iP,:) + Ba(iP,ia)*C(Bj(iP,ia),:)
     end do
  end do

end subroutine evaluatePt




subroutine inverseMap(k, m, nP, P0, C, t)

  implicit none

  !Fortran-python interface directive
  !f2py intent(in) k, m, nP, P0, C
  !f2py intent(out) t
  !f2py depend(nP) P0
  !f2py depend(m) C
  !f2py depend(nP) t

  !Input
  integer, intent(in) ::  k, m, nP
  double precision, intent(in) ::  P0(nP), C(m)

  !Output
  double precision, intent(out) ::  t(nP)

  !Working
  integer i, j, i0, i1, iP
  double precision f, dfdx, B(k), B1(k), d(k+m), x

  call knotopen(k, k+m, d)

  do iP=1,nP
     x = (P0(iP)-C(1))/(C(m)-C(1))

     call basis(k, k+m, x, d, B, i0)
     call basis1(k, k+m, x, d, B1, i1)
     f = -P0(iP)
     dfdx = 0
     do i=1,k
        f = f + B(i)*C(i0+i)
        dfdx = dfdx + B1(i)*C(i1+i)
     end do
     
     do j=1,100    
        !print *, j, f
        if (abs(f) .lt. 1e-15) then
           exit
        end if
        x = x - f/dfdx
        
        call basis(k, k+m, x, d, B, i0)
        call basis1(k, k+m, x, d, B1, i1)
        f = -P0(iP)
        dfdx = 0
        do i=1,k
           f = f + B(i)*C(i0+i)
           dfdx = dfdx + B1(i)*C(i1+i)
        end do
     end do
     t(iP) = x
  end do

end subroutine inverseMap
