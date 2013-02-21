subroutine evaluatePt(d1, d2, nx, nf, nC, ks, ms, t, C, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) d1, d2, nx, nf, nC, ks, ms, t, C
  !f2py intent(out) P
  !f2py depend(nx) ks, ms
  !f2py depend(nx) t
  !f2py depend(nC, nf) C
  !f2py depend(nf) P

  !Input
  integer, intent(in) ::  d1, d2, nx, nf, nC, ks(nx), ms(nx)
  double precision, intent(in) ::  t(nx), C(nC,nf)

  !Output
  double precision, intent(out) ::  P(nf)

  !Working
  integer k, m, na
  integer iC
  integer i, ix, ia, ik, i0, rem
  double precision, allocatable, dimension(:) ::  d, B
  double precision, allocatable, dimension(:) ::  Ba
  integer, allocatable, dimension(:) ::  Bj

  na = product(ks)

  allocate(Ba(na))
  allocate(Bj(na))

  Ba(:) = 1.0
  Bj(:) = 1
  do ix=1,nx
     k = ks(ix)
     m = ms(ix)
     allocate(d(k+m))
     allocate(B(k))
     call knotopen(k, k+m, d)
     if ((d1 .ne. ix) .and. (d2 .ne. ix)) then
        call basis(k, k+m, t(ix), d, B, i0)
     else if ((d1 .eq. ix) .and. (d2 .ne. ix)) then
        call basis1(k, k+m, t(ix), d, B, i0)
     else if ((d1 .ne. ix) .and. (d2 .eq. ix)) then
        call basis1(k, k+m, t(ix), d, B, i0)
     else if ((d1 .eq. ix) .and. (d2 .eq. ix)) then
        call basis2(k, k+m, t(ix), d, B, i0)
     end if
     do ia=1,na
        rem = ia - 1
        do i=nx,ix,-1
           rem = mod(rem,product(ks(:i)))
        end do
        ik = rem/product(ks(:ix-1)) + 1
        iC = (ik + i0 - 1)*product(ms(:ix-1))
        Ba(ia) = Ba(ia) * B(ik)
        Bj(ia) = Bj(ia) + iC
     end do
     deallocate(d)
     deallocate(B)
  end do

  P(:) = 0.0
  do ia=1,na
     P = P + Ba(ia)*C(Bj(ia),:)
  end do

  deallocate(Ba)
  deallocate(Bj)

end subroutine evaluatePt




subroutine inverseMap(k, m, P0, C, x)

  implicit none

  !Fortran-python interface directive
  !f2py intent(in) k, m, P0, C
  !f2py intent(out) x
  !f2py depend(m) C

  !Input
  integer, intent(in) ::  k, m
  double precision, intent(in) ::  P0, C(m)

  !Output
  double precision, intent(out) ::  x

  !Working
  integer i, j, i0, i1
  double precision f, dfdx, B(k), B1(k), d(k+m)

  call knotopen(k, k+m, d)

  x = (P0-C(1))/(C(m)-C(1))

  call basis(k, k+m, x, d, B, i0)
  call basis1(k, k+m, x, d, B1, i1)
  f = -P0
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
     f = -P0
     dfdx = 0
     do i=1,k
        f = f + B(i)*C(i0+i)
        dfdx = dfdx + B1(i)*C(i1+i)
     end do
  end do

end subroutine inverseMap
