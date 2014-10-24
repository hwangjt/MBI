subroutine evaluate(d1, d2, nx, nf, nC, nCx1, nCx2, nP, ks, ms, t, C, Cx1, Cx2, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) d1, d2, nx, nf, nC, nCx1, nCx2, nP, ks, ms, t, C, Cx1, Cx2
  !f2py intent(out) P
  !f2py depend(nx) ks, ms
  !f2py depend(nP,nx) t
  !f2py depend(nC,nf) C
  !f2py depend(nCx1) Cx1
  !f2py depend(nCx2) Cx2
  !f2py depend(nP,nf) P

  !Input
  integer, intent(in) ::  d1, d2, nx, nf, nC, nCx1, nCx2, nP, ks(nx), ms(nx)
  double precision, intent(in) ::  t(nP,nx), C(nC,nf), Cx1(nCx1,1), Cx2(nCx2,1)

  !Output
  double precision, intent(out) ::  P(nP,nf)

  !Working
  integer iP, na, na1, na2, k1(1), k2(1), m1(1), m2(1)
  double precision t1(nP,1), t2(nP,1)
  double precision f(nP,nf), dfdu(nP,nf), d2fdu2(nP,nf), d2fdudv(nP,nf)
  double precision dxdu(nP,1), dxdv(nP,1), d2xdu2(nP,1)

  na = product(ks)
  if (d1 .gt. 0) then
     na1 = ks(d1)
     k1(1) = ks(d1)
     m1(1) = ms(d1)
     t1(:,1) = t(:,d1)
  end if
  if (d2 .gt. 0) then
     na2 = ks(d2)
     k2(1) = ks(d2)
     m2(1) = ms(d2)
     t2(:,1) = t(:,d2)
  end if

  if ((d1 .eq. 0) .and. (d2 .eq. 0)) then
     call evaluatePts(0, 0, na, nx, nf, nC, nP, ks, ms, t, C, f)
     P(:,:) = f(:,:)
  else if ((d1 .ne. 0) .and. (d2 .eq. 0)) then
     call evaluatePts(d1, 0,  na, nx, nf,   nC, nP, ks, ms,  t,   C, dfdu)
     call evaluatePts( 1, 0, na1,  1,  1, nCx1, nP, k1, m1, t1, Cx1, dxdu)
     do iP=1,nP
        P(iP,:) = dfdu(iP,:)/dxdu(iP,1)
     end do
  else if ((d1 .ne. 0) .and. (d1 .eq. d2)) then
     call evaluatePts(d1,  0,  na, nx, nf,   nC, nP, ks, ms,  t,   C, dfdu)
     call evaluatePts(d1, d1,  na, nx, nf,   nC, nP, ks, ms,  t,   C, d2fdu2)
     call evaluatePts( 1,  0, na1,  1,  1, nCx1, nP, k1, m1, t1, Cx1, dxdu)
     call evaluatePts( 1,  1, na1,  1,  1, nCx1, nP, k1, m1, t1, Cx1, d2xdu2)
     do iP=1,nP
        P(iP,:) = (d2fdu2(iP,:)*dxdu(iP,1) - dfdu(iP,:)*d2xdu2(iP,1))/dxdu(iP,1)**3
     end do
  else if ((d1 .ne. 0) .and. (d2 .ne. 0)) then
     call evaluatePts(d1, d2,  na, nx, nf,   nC, nP, ks, ms,  t,   C, d2fdudv)
     call evaluatePts( 1,  0, na1,  1,  1, nCx1, nP, k1, m1, t1, Cx1, dxdu)
     call evaluatePts( 1,  0, na2,  1,  1, nCx2, nP, k2, m2, t2, Cx2, dxdv)
     do iP=1,nP
        P(iP,:) = d2fdudv(iP,:)/dxdu(iP,1)/dxdv(iP,1)
     end do
  else
     P(:,:) = 0.0
  end if
  
end subroutine evaluate



subroutine evaluatePts(d1, d2, na, nx, nf, nC, nP, ks, ms, t, C, P)

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

end subroutine evaluatePts




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
     dfdx = 0.0
     do i=1,k
        f = f + B(i)*C(i0+i)
        dfdx = dfdx + B1(i)*C(i1+1)
     end do

     do j=1,100
        !print *, j, x, f
        if (abs(f) .lt. 1e-15) then
           exit
        end if

        x = x - f/dfdx

        if (x .lt. 0) then
           x = 0.0
        else if (x .gt. 1) then
           x = 1.0
        end if

        call basis(k, k+m, x, d, B, i0)
        call basis1(k, k+m, x, d, B1, i1)
        f = -P0(iP)
        dfdx = 0.0
        do i=1,k
           f = f + B(i)*C(i0+i)
           dfdx = dfdx + B1(i) * C(i1+i)
        end do

     end do
     t(iP) = x
  end do

end subroutine inverseMap
