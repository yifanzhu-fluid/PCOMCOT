subroutine solveDynamicPressure(GP, LP, iLayer, LocalData, LocalDataLength)
! Subroutine solveDynamicPressure solves the Poisson-type equation of dynamic pressure.
! The coefficients of 5-point equations are stored in PL, PR, PB, PT and PC.
! Solution of dynamic pressure is stored in Q.
! This solver uses the iterative algorithm Bi-CGSTAB to solve the sparse equation.
! Incomplete LU preconditioning (ILU) is used to accelerate the solver
use mpi
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer, LocalDatalength
real*8      ::  LocalData(LocalDatalength)
integer*4   ::  irank, nsize, master
integer*4   ::  nstartx, nendx, nstarty, nendy
real*8      ::  D_Inverse(LP(iLayer)%NX, LP(iLayer)%NY)
real*8      ::  r(LP(iLayer)%NX, LP(iLayer)%NY), rc0(LP(iLayer)%NX, LP(iLayer)%NY)
real*8      ::  v(LP(iLayer)%NX, LP(iLayer)%NY),   p(LP(iLayer)%NX, LP(iLayer)%NY)
real*8      ::  s(LP(iLayer)%NX, LP(iLayer)%NY),   t(LP(iLayer)%NX, LP(iLayer)%NY)
real*8      ::  sc(LP(iLayer)%NX, LP(iLayer)%NY), tc(LP(iLayer)%NX, LP(iLayer)%NY)
real*8      ::  y(LP(iLayer)%NX, LP(iLayer)%NY),   z(LP(iLayer)%NX, LP(iLayer)%NY)
real*8      ::  rho0, rho, alpha, omega, beta
real*8      ::  rc0_v, tc_sc, tc_tc
real*8      ::  Min_InnerProduct, tol, r_max, b_max, r_norm, b_norm
integer*4   ::  i, j, errorcode, ierror, iter, itmax
real*8      ::  CPUTime1, CPUTime2

irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master 
tol = 1.0e-7; Min_InnerProduct = 1.0e-30; itmax = 100
call CPU_TIME(CPUTime1)

if(irank.lt.LP(iLayer)%nsize) then
    nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
    nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
    nstarty = LP(iLayer)%PartitionInfo(irank+1,3) 
    nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
else
    nstartx = 0; nendx = -1; nstarty = 0; nendy = -1
endif
!LP(iLayer)%Q = 0.0d0

call ILUpreconditioner(GP, LP, iLayer, D_Inverse)
call BiCGSTAB_INIT(GP, LP, iLayer, r, rc0, rho0, alpha, omega, v, p)

!/// check if the initial guess of dynami pressure is already accurate enough ///!
r_max = 0.0d0; b_max = 0.0d0
do j = nstarty,nendy
do i = nstartx,nendx
    b_max = MAX(b_max, ABS(LP(iLayer)%PQ(i,j)))
    r_max = MAX(r_max, ABS(r(i,j)))
enddo
enddo
call MPI_ALLREDUCE(b_max, b_norm, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
call MPI_ALLREDUCE(r_max, r_norm, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)

if(b_norm.le.tol) then
    LP(iLayer)%Q = 0.0d0; return
elseif(r_norm/b_norm.le.tol) then
    return
endif

iter = 0
do
    iter = iter + 1
    if(iter.gt.itmax) then
        if(irank.eq.master) then
            write(*,*)
            write(*,'(a,f4.1)') 'ERROR: Bi-CGSTAB algorithm does not converge when solving dispersion term!'
            write(*,*)
            call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
        endif
    endif

    call innerProduct(rc0, r, rho, LP(iLayer)%NX, LP(iLayer)%NY, nstartx, nendx, nstarty, nendy)
    beta = (rho/rho0)*(alpha/omega)
    do j = nstarty,nendy
    do i = nstartx,nendx
        p(i,j) = r(i,j) + beta*(p(i,j)-omega*v(i,j))
    enddo
    enddo
    call solveILUprecondition(GP, LP, iLayer, D_Inverse, p, y)
    call productAx(GP, LP, iLayer, y, v, LocalData, LocalDataLength)
    call innerProduct(rc0, v, rc0_v, LP(iLayer)%NX, LP(iLayer)%NY, nstartx, nendx, nstarty, nendy)

    !/// check if the value of rc0 needs modifying to ensure convergence ///!
    if(ABS(rho).le.Min_InnerProduct.or.ABS(rc0_v).le.Min_InnerProduct) then
        do j = nstarty,nendy
        do i = nstartx,nendx
            call random_number(rc0(i,j))
            rc0(i,j) = 1.0/rc0_v*rc0(i,j)
        enddo
        enddo
        v = 0.0d0; p = 0.0d0
        rho0 = 1.0d0; alpha = 1.0d0; omega = 1.0d0
        cycle
    endif

    alpha = rho/rc0_v
    do j = nstarty,nendy
    do i = nstartx,nendx
        s(i,j) = r(i,j) - alpha*v(i,j)
    enddo
    enddo
    call solveILUprecondition(GP, LP, iLayer, D_Inverse, s, z)
    call productAx(GP, LP, iLayer, z, t, LocalData, LocalDataLength)
    call solveLowerTriangular(GP, LP, iLayer, D_Inverse, t, tc)
    call solveLowerTriangular(GP, LP, iLayer, D_Inverse, s, sc)
    call innerProduct(tc, sc, tc_sc, LP(iLayer)%NX, LP(iLayer)%NY, nstartx, nendx, nstarty, nendy)
    call innerProduct(tc, tc, tc_tc, LP(iLayer)%NX, LP(iLayer)%NY, nstartx, nendx, nstarty, nendy)

    omega = tc_sc/tc_tc; r_max = 0.0d0
    do j = nstarty,nendy
    do i = nstartx,nendx
        LP(iLayer)%Q(i,j) = LP(iLayer)%Q(i,j) + alpha*y(i,j) + omega*z(i,j)
        r(i,j) = s(i,j) - omega*t(i,j)
        r_max = MAX(r_max, ABS(r(i,j)))
    enddo
    enddo
    rho0 = rho
    call MPI_ALLREDUCE(r_max, r_norm, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
    if(r_norm/b_norm.le.tol) exit
enddo
call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1 

end subroutine solveDynamicPressure



subroutine BiCGSTAB_INIT(GP, LP, iLayer, r, rc0, rho0, alpha, omega, v, p)

use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer
real*8      ::  r(LP(iLayer)%NX,LP(iLayer)%NY), rc0(LP(iLayer)%NX,LP(iLayer)%NY)
real*8      ::  v(LP(iLayer)%NX,LP(iLayer)%NY),   p(LP(iLayer)%NX,LP(iLayer)%NY)
real*8      ::  rho0, alpha, omega
integer*4   ::  irank, nstartx, nendx, nstarty, nendy, i, j

irank = GP%irank

if(irank.lt.LP(iLayer)%nsize) then

nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)

do j = nstarty,nendy
do i = nstartx,nendx
    r(i,j) = LP(iLayer)%PQ(i,j) - LP(iLayer)%PC(i,j)*LP(iLayer)%Q(i,j)
    if(i.ge.2)                  r(i,j) = r(i,j) - LP(iLayer)%PL(i,j)*LP(iLayer)%Q(i-1,j)
    if(i.le.LP(iLayer)%NX-1)    r(i,j) = r(i,j) - LP(iLayer)%PR(i,j)*LP(iLayer)%Q(i+1,j)
    if(j.ge.2)                  r(i,j) = r(i,j) - LP(iLayer)%PB(i,j)*LP(iLayer)%Q(i,j-1)
    if(j.le.LP(iLayer)%NY-1)    r(i,j) = r(i,j) - LP(iLayer)%PT(i,j)*LP(iLayer)%Q(i,j+1)
    rc0(i,j) = r(i,j)
enddo
enddo

endif !if:this node is used by this layer

v = 0.0d0; p = 0.0d0
rho0 = 1.0d0; alpha = 1.0d0; omega = 1.0d0


end subroutine BiCGSTAB_INIT



subroutine ILUpreconditioner(GP, LP, iLayer, D_Inverse)
! Subroutine ILUpreconditioner applies partitioned ILU preconditioning to the original linear equations of dynamic pressure.
! The coefficient matrix is factorized into a lower and an upper triangular matrixes.
! The lower triangular matrix is L = D+L, and the upper matrix is U = I+U/D.
! The inverse of D is stored distributedly in D_Inverse to save time.
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer
real*8      ::  D_Inverse(LP(iLayer)%NX, LP(iLayer)%NY)
integer*4   ::  irank, nstartx, nendx, nstarty, nendy, i, j
real*8      ::  d0, d1, d2, dmin

irank = GP%irank; dmin = 1.0e-6; 

if(irank.lt.LP(iLayer)%nsize) then

nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)

do j = nstarty,nendy
do i = nstartx,nendx
    if(i.ge.nstartx+1) then
        d1 = LP(iLayer)%PL(i,j)*D_Inverse(i-1,j)*LP(iLayer)%PR(i-1,j)
    else
        d1 = 0.0d0
    endif
    if(j.ge.nstarty+1) then
        d2 = LP(iLayer)%PB(i,j)*D_Inverse(i,j-1)*LP(iLayer)%PT(i,j-1)
    else
        d2 = 0.0d0
    endif
    d0 = LP(iLayer)%PC(i,j) - d1 - d2
    if(ABS(d0).gt.dmin) then
        D_Inverse(i,j) = 1.0/d0
    elseif(d0.ge.0.0) then
        D_Inverse(i,j) = 1.0/dmin
    else
        D_Inverse(i,j) = -1.0/dmin
    endif
enddo
enddo

endif !if:this node is used by this layer

end subroutine ILUpreconditioner



subroutine solveILUprecondition(GP, LP, iLayer, D_Inverse, y, x)
! Subroutine solveILUprecondition solves the preconditioning eqution LUx=y.
! First, the lower triangular equation Lz = y is solved.
! Then, the upper triangular equation Ux = z is solved.
! For parallelization, the solving process is implemented independently at each node.
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer
real*8, intent(in)  ::  D_Inverse(LP(iLayer)%NX, LP(iLayer)%NY), y(LP(iLayer)%NX, LP(iLayer)%NY)
real*8, intent(out) ::  x(LP(iLayer)%NX, LP(iLayer)%NY)
integer*4   ::  irank, nstartx, nendx, nstarty, nendy, i, j
real*8      ::  lx1, lx2, ux1, ux2

irank = GP%irank

if(irank.lt.LP(iLayer)%nsize) then

nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)

!/// solve lower triangular equation Lx=y ///!
do j = nstarty,nendy
do i = nstartx,nendx
    if(i.ge.nstartx+1) then
        lx1 = LP(iLayer)%PL(i,j)*x(i-1,j)
    else
        lx1 = 0.0d0
    endif
    if(j.ge.nstarty+1) then
        lx2 = LP(iLayer)%PB(i,j)*x(i,j-1)
    else
        lx2 = 0.0d0
    endif
    x(i,j) = D_Inverse(i,j)*(y(i,j) - lx1 - lx2)
enddo
enddo

!/// solve upper triangular equation Ux=x ///!
do j = nendy,nstarty,-1
do i = nendx,nstartx,-1
    if(i.le.nendx-1) then
        ux1 = LP(iLayer)%PR(i,j)*x(i+1,j)
    else
        ux1 = 0.0d0
    endif
    if(j.le.nendy-1) then
        ux2 = LP(iLayer)%PT(i,j)*x(i,j+1)
    else
        ux2 = 0.0d0
    endif
    x(i,j) = x(i,j) - D_Inverse(i,j)*(ux1 + ux2)
enddo
enddo

endif !if:this node is used by this layer


end subroutine solveILUprecondition



subroutine solveLowerTriangular(GP, LP, iLayer, D_Inverse, y, x)
! Subroutine solveLowerTriangular solves equation (L+D)x=y independently on each node.
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer
real*8, intent(in)  ::  D_Inverse(LP(iLayer)%NX, LP(iLayer)%NY), y(LP(iLayer)%NX, LP(iLayer)%NY)
real*8, intent(out) ::  x(LP(iLayer)%NX, LP(iLayer)%NY)
integer*4   ::  irank, nstartx, nendx, nstarty, nendy, i, j
real*8      ::  lx1, lx2

irank = GP%irank

if(irank.lt.LP(iLayer)%nsize) then

nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)

do j = nstarty,nendy
do i = nstartx,nendx
    if(i.ge.nstartx+1) then
        lx1 = LP(iLayer)%PL(i,j)*x(i-1,j)
    else
        lx1 = 0.0d0
    endif
    if(j.ge.nstarty+1) then
        lx2 = LP(iLayer)%PB(i,j)*x(i,j-1)
    else
        lx2 = 0.0d0
    endif
    x(i,j) = D_Inverse(i,j)*(y(i,j) - lx1 - lx2)
enddo
enddo

endif !if: this node is used by this layer

end subroutine solveLowerTriangular



subroutine solveUpperTriangular(GP, LP, iLayer, D_Inverse, y, x)
! Subroutine solveUpperTriangular solve equation (I+U/D)x=y independently on each node.
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer
real*8, intent(in)  ::  D_Inverse(LP(iLayer)%NX, LP(iLayer)%NY), y(LP(iLayer)%NX, LP(iLayer)%NY)
real*8, intent(out) ::  x(LP(iLayer)%NX, LP(iLayer)%NY)
integer*4   ::  irank, nstartx, nendx, nstarty, nendy, i, j
real*8      ::  ux1, ux2

irank = GP%irank

if(irank.lt.LP(iLayer)%nsize) then

nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)

do j = nendy,nstarty,-1
do i = nendx,nstartx,-1
    if(i.le.nendx-1) then
        ux1 = LP(iLayer)%PR(i,j)*x(i+1,j)
    else
        ux1 = 0.0d0
    endif
    if(j.le.nendy-1) then
        ux2 = LP(iLayer)%PT(i,j)*x(i,j+1)
    else
        ux2 = 0.0d0
    endif
    x(i,j) = y(i,j) - D_Inverse(i,j)*(ux1 + ux2)
enddo
enddo

endif !if: this node is used by this layer

end subroutine solveUpperTriangular



subroutine innerProduct(a,b,p,m,n,nstartx,nendx,nstarty,nendy)
! Subroutine innerProduct calculates inner product of two matrix-form vectors a and b. 
! a(m,n) and b(m,n) are stored distributedly in different processes.
! nstartx, nendx, nstarty and nendy are the range of index of part of a,b in this process.
! The result of inner product is obtained by summing up parts in all processes by global reduce and stored in p.
use mpi
implicit NONE
integer*4, intent(in)   ::  m, n, nstartx, nendx, nstarty,nendy
real*8, intent(in)      ::  a(m,n), b(m,n)
real*8, intent(out)     ::  p
real*8                  ::  p_part
integer*4               ::  i, j, ierror

p_part = 0.0d0
do j = nstarty,nendy
do i = nstartx,nendx
    p_part = p_part + a(i,j)*b(i,j)
enddo
enddo
call MPI_ALLREDUCE(p_part, p, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)

end subroutine innerProduct



subroutine productAx(GP, LP, iLayer, x, y, LocalData, LocalDataLength)
! Subroutine productAx calculates product of coefficient matrix A and a vector x. 
! A is composed of PL, PR, PB, PT and PC. x is a vector in the form of matrix, and needs synchronizing before calculation.
! The result of product is stored distributedly in y.
use mpi
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer, LocalDataLength
real*8      ::  LocalData(LocalDataLength)
real*8      ::  x(LP(iLayer)%NX, LP(iLayer)%NY), y(LP(iLayer)%NX, LP(iLayer)%NY)
integer*4   ::  irank, nsize, master
integer*4   ::  ierror, istatus(MPI_STATUS_SIZE)
integer*4   ::  iBC, i, j
integer*4   ::  nstartx, nendx, nstarty, nendy

irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master

!/// synchronize the values of matrix-form vector x on the boundaries between subdomains ///!
do iBC = 1,LP(iLayer)%BoundarySendRecvCount
    if(irank.eq.LP(iLayer)%BoundarySendRecv(iBC, 1)) then
        nstartx = LP(iLayer)%BoundarySendRecv(iBC,3)
        nendx   = LP(iLayer)%BoundarySendRecv(iBC,4)
        nstarty = LP(iLayer)%BoundarySendRecv(iBC,5)
        nendy   = LP(iLayer)%BoundarySendRecv(iBC,6)
        if(nstartx.le.nendx.and.nstarty.le.nendy) then
            do j = nstarty,nendy
            do i = nstartx,nendx
                LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = x(i,j)
            enddo
            enddo
            call MPI_SEND(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1),MPI_DOUBLE_PRECISION, &
                LP(iLayer)%BoundarySendRecv(iBC,2),2022,MPI_COMM_WORLD,ierror)
        endif
    elseif(irank.eq.LP(iLayer)%BoundarySendRecv(iBC,2)) then
        nstartx = LP(iLayer)%BoundarySendRecv(iBC,3)
        nendx   = LP(iLayer)%BoundarySendRecv(iBC,4)
        nstarty = LP(iLayer)%BoundarySendRecv(iBC,5)
        nendy   = LP(iLayer)%BoundarySendRecv(iBC,6)
        if(nstartx.le.nendx.and.nstarty.le.nendy) then
            call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1),MPI_DOUBLE_PRECISION, &
                LP(iLayer)%BoundarySendRecv(iBC,1),2022,MPI_COMM_WORLD,istatus,ierror)
            do j = nstarty, nendy
            do i = nstartx, nendx
                x(i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
            enddo
            enddo
        endif
    endif
enddo

!/// compute part of A*x on each node ///!
if(irank.lt.LP(iLayer)%nsize) then

nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)

do j = nstarty,nendy
do i = nstartx,nendx
    y(i,j) = LP(iLayer)%PC(i,j)*x(i,j)
    if(i.ge.2)                  y(i,j) = y(i,j) + LP(iLayer)%PL(i,j)*x(i-1,j)
    if(i.le.LP(iLayer)%NX-1)    y(i,j) = y(i,j) + LP(iLayer)%PR(i,j)*x(i+1,j)
    if(j.ge.2)                  y(i,j) = y(i,j) + LP(iLayer)%PB(i,j)*x(i,j-1)
    if(j.le.LP(iLayer)%NY-1)    y(i,j) = y(i,j) + LP(iLayer)%PT(i,j)*x(i,j+1)
enddo
enddo

endif !if:this node is used by this layer

end subroutine productAx

                

