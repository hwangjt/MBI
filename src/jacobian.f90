subroutine computeJacobian(d1, d2, nx, nP, nB, ks, ms, t, Ba, Bi, Bj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) d1, d2, nx, nP, nB, ks, ms, t
  !f2py intent(out) Ba, Bi, Bj
  !f2py depend(nx) ks, ms
  !f2py depend(nP,nx) t
  !f2py depend(nB) Ba, Bi, Bj

  !Input
  integer, intent(in) ::  d1, d2, nx, nP, nB, ks(nx), ms(nx)
  double precision, intent(in) ::  t(nP,nx)

  !Output
  double precision, intent(out) ::  Ba(nB)
  integer, intent(out) ::  Bi(nB), Bj(nB)

  !Working
  integer k, m, na
  integer iP, iB, iC
  integer i, ix, ia, ik, i0, rem
  double precision, allocatable, dimension(:) ::  d, B

  na = product(ks)

  do iP=1,nP
     Bi((iP-1)*na+1:iP*na) = iP - 1
  end do

  Ba(:) = 1.0
  Bj(:) = 0
  do ix=1,nx
     k = ks(ix)
     m = ms(ix)
     allocate(d(k+m))
     allocate(B(k))
     call knotopen(k, k+m, d)
     do iP=1,nP
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
           iB = (iP-1)*na + ia
           iC = (ik + i0 - 1)*product(ms(:ix-1))
           Ba(iB) = Ba(iB) * B(ik)
           Bj(iB) = Bj(iB) + iC
        end do
     end do
     deallocate(d)
     deallocate(B)
  end do

end subroutine computeJacobian
