subroutine addDispersion(GP, LP, iLayer, iStep, iTimeStep, LocalData, LocalDataLength)

use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(LayerParameters)    ::  LP(100)
integer*4   ::  iLayer, iStep, iTimeStep, LocalDataLength
real*8      ::  LocalData(LocalDataLength)      

call bottomVerticalVelocity(GP, LP, iLayer)

if(iTimeStep.gt.1) then ! We do not estimate dynamic pressure at the first major step due to lack of inital condition
    call checkZeroDynamicPressure(GP, LP, iLayer)
    call constructDynamicPressureEquations(GP, LP, iLayer)
    call solveDynamicPressure(GP, LP, iLayer, LocalData, LocalDataLength)
    call bcastComputeDomainBoundaryDynamicPressure(GP, LP, iLayer, LocalData, LocalDataLength)
    if(GP%CoordinatesType.eq.1) then
        call dispersionCartesian(GP, LP, iLayer)
    elseif(GP%CoordinatesType.eq.0) then
        call dispersionSpherical(GP, LP, iLayer)
    endif
else
    LP(iLayer)%Q = 0.0d0
endif

call surfaceVerticalVelocity(GP, LP, iLayer, iTimeStep)

end subroutine addDispersion



subroutine dispersionCartesian(GP, LP, iLayer)

use VariableDefination
implicit NONE
type(GlobalParameters)    ::  GP
type(LayerParameters)     ::  LP(100)
integer*4   ::  iLayer, irank
integer*4   ::  nstartx, nendx, nstarty, nendy
integer*4   ::  istart, iend, jstart, jend
integer*4   ::  i, j, k
real*8      ::  alpha, beta
real*8      ::  q1, q2, h1, h2, z1, z2
real*8      ::  CPUTime1, CPUTime2

irank = GP%irank
call CPU_TIME(CPUTime1)

!/// use different governing equations based on rate of change of water depth ///!
if(LP(iLayer)%DepthVariability.eq.0) then
    beta = 2.0/3.0; alpha = 0.5
else
    beta = 0.5; alpha = 1.0
endif

if(irank.lt.LP(iLayer)%nsize) then
   
nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)

!///calculate dispersive flux component M///!
if(LP(iLayer)%Level.eq.1) then !calculate all grids in top layer
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX-1,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY,nendy)
elseif(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids in child layer
    istart = MAX(GP%nGHOST+1,nstartx); iend = MIN(LP(iLayer)%NX-GP%nGHOST-1,nendx)
    jstart = MAX(GP%nGHOST+1,nstarty); jend = MIN(LP(iLayer)%NY-GP%nGHOST,nendy)
endif
do j = jstart,jend
do i = istart,iend
    !/// exclude the cell where water cannot flow freely, because dynamic pressure is zero on both sides ///!
    if(LP(iLayer)%Sign_M(i,j).ne.0) cycle
    !/// modify the value of M estimated from shallow water equations ///!
    q1 = LP(iLayer)%Q(i,j);   q2 = LP(iLayer)%Q(i+1,j)
    h1 = LP(iLayer)%H(2,i,j); h2 = LP(iLayer)%H(2,i+1,j)
    z1 = LP(iLayer)%Z(i,j);   z2 = LP(iLayer)%Z(i+1,j)
    LP(iLayer)%M(2,i,j) = LP(iLayer)%M(2,i,j) - beta*LP(iLayer)%CC1*LP(iLayer)%D_M(i,j)*(q2-q1) &
        -0.5*beta*LP(iLayer)%CC1*(q1+q2)*((h2-alpha*z2)-(h1-alpha*z1))
enddo
enddo

!///calculate dispersive flux component N///!
if(LP(iLayer)%Level.eq.1) then !calculate all grids in top layer 
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY-1,nendy)
elseif(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids on child layers
    istart = MAX(GP%nGHOST+1,nstartx); iend = MIN(LP(iLayer)%NX-GP%nGHOST,nendx)
    jstart = MAX(GP%nGHOST+1,nstarty); jend = MIN(LP(iLayer)%NY-GP%nGHOST-1,nendy)
endif
do j = jstart,jend
do i = istart,iend
    !/// exclude the cell where water cannot flow freely, because dynamic pressure is zero on both sides ///!
    if(LP(iLayer)%Sign_N(i,j).ne.0) cycle
    !///modify the value of N estimated from shallow water equations///!
    q1 = LP(iLayer)%Q(i,j);   q2 = LP(iLayer)%Q(i,j+1)
    h1 = LP(iLayer)%H(2,i,j); h2 = LP(iLayer)%H(2,i,j+1)
    z1 = LP(iLayer)%Z(i,j);   z2 = LP(iLayer)%Z(i,j+1)
    LP(iLayer)%N(2,i,j) = LP(iLayer)%N(2,i,j) - beta*LP(iLayer)%CC2*LP(iLayer)%D_N(i,j)*(q2-q1) &
        -0.5*beta*LP(iLayer)%CC2*(q1+q2)*((h2-alpha*z2)-(h1-alpha*z1))
enddo
enddo

!//// set zero water flux around dry cells ////!
do k = 1,LP(iLayer)%PermanentDryCellsCount
    i = LP(iLayer)%PermanentDryCells(k,1)
    j = LP(iLayer)%PermanentDryCells(k,2)
    if(i.le.LP(iLayer)%NX-1) LP(iLayer)%M(2,i,j) = 0.0d0
    if(i.ge.2) LP(iLayer)%M(2,i-1,j) = 0.0d0
    if(j.le.LP(iLayer)%NY-1) LP(iLayer)%N(2,i,j) = 0.0d0
    if(j.ge.2) LP(iLayer)%N(2,i,j-1) = 0.0d0
enddo

endif !if: this node is used by this layer

call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1

end subroutine dispersionCartesian



subroutine dispersionSpherical(GP, LP, iLayer)

use VariableDefination
implicit NONE
type(GlobalParameters)    ::  GP
type(LayerParameters)     ::  LP(100)
integer*4   ::  iLayer, irank
integer*4   ::  nstartx, nendx, nstarty, nendy
integer*4   ::  istart, iend, jstart, jend
integer*4   ::  i, j, k
real*8      ::  alpha, beta
real*8      ::  q1, q2, h1, h2, z1, z2
real*8      ::  CPUTime1, CPUTime2
    
irank = GP%irank
call CPU_TIME(CPUTime1)

!/// use different governing equations based on rate of change of water depth ///!
if(LP(iLayer)%DepthVariability.eq.0) then
    beta = 2.0/3.0; alpha = 0.5
else
    beta = 0.5; alpha = 1.0
endif
    
if(irank.lt.LP(iLayer)%nsize) then
       
nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
    
!///calculate dispersive flux component M///!
if(LP(iLayer)%Level.eq.1) then !calculate all grids in top layer
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX-1,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY,nendy)
elseif(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids in child layer
    istart = MAX(GP%nGHOST+1,nstartx); iend = MIN(LP(iLayer)%NX-GP%nGHOST-1,nendx)
    jstart = MAX(GP%nGHOST+1,nstarty); jend = MIN(LP(iLayer)%NY-GP%nGHOST,nendy)
endif
do j = jstart,jend
do i = istart,iend
    !/// exclude the cell where water cannot flow freely, because dynamic pressure is zero on both sides ///!
    if(LP(iLayer)%Sign_M(i,j).ne.0) cycle
    !/// modify the value of M estimated from shallow water equations ///!
    q1 = LP(iLayer)%Q(i,j);   q2 = LP(iLayer)%Q(i+1,j)
    h1 = LP(iLayer)%H(2,i,j); h2 = LP(iLayer)%H(2,i+1,j)
    z1 = LP(iLayer)%Z(i,j);   z2 = LP(iLayer)%Z(i+1,j)
    LP(iLayer)%M(2,i,j) = LP(iLayer)%M(2,i,j) - beta*LP(iLayer)%CS1(j)*LP(iLayer)%D_M(i,j)*(q2-q1) &
            -0.5*beta*LP(iLayer)%CS1(j)*(q1+q2)*((h2-alpha*z2)-(h1-alpha*z1))
enddo
enddo
    
!///calculate dispersive flux component N///!
if(LP(iLayer)%Level.eq.1) then !calculate all grids in top layer 
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY-1,nendy)
elseif(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids on child layers
    istart = MAX(GP%nGHOST+1,nstartx); iend = MIN(LP(iLayer)%NX-GP%nGHOST,nendx)
    jstart = MAX(GP%nGHOST+1,nstarty); jend = MIN(LP(iLayer)%NY-GP%nGHOST-1,nendy)
endif
do j = jstart,jend
do i = istart,iend
    !/// exclude the cell where water cannot flow freely, because dynamic pressure is zero on both sides ///!
    if(LP(iLayer)%Sign_N(i,j).ne.0) cycle
    !///modify the value of N estimated from shallow water equations///!
    q1 = LP(iLayer)%Q(i,j);   q2 = LP(iLayer)%Q(i,j+1)
    h1 = LP(iLayer)%H(2,i,j); h2 = LP(iLayer)%H(2,i,j+1)
    z1 = LP(iLayer)%Z(i,j);   z2 = LP(iLayer)%Z(i,j+1)
    LP(iLayer)%N(2,i,j) = LP(iLayer)%N(2,i,j) - beta*LP(iLayer)%CS5*LP(iLayer)%D_N(i,j)*(q2-q1) &
            -0.5*beta*LP(iLayer)%CS5*(q1+q2)*((h2-alpha*z2)-(h1-alpha*z1))
enddo
enddo
    
!//// set zero water flux around dry cells ////!
do k = 1,LP(iLayer)%PermanentDryCellsCount
    i = LP(iLayer)%PermanentDryCells(k,1)
    j = LP(iLayer)%PermanentDryCells(k,2)
    if(i.le.LP(iLayer)%NX-1) LP(iLayer)%M(2,i,j) = 0.0d0
    if(i.ge.2) LP(iLayer)%M(2,i-1,j) = 0.0d0
    if(j.le.LP(iLayer)%NY-1) LP(iLayer)%N(2,i,j) = 0.0d0
    if(j.ge.2) LP(iLayer)%N(2,i,j-1) = 0.0d0
enddo
    
endif !if: this node is used by this layer
    
call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1
    
end subroutine dispersionSpherical



subroutine checkZeroDynamicPressure(GP, LP, iLayer)

use VariableDefination
implicit NONE
type(GlobalParameters)    ::  GP
type(LayerParameters)     ::  LP(100)
integer*4   ::  iLayer, irank
integer*4   ::  nstartx, nendx, nstarty, nendy
integer*4   ::  i, j
real*8      ::  CPUTime1, CPUTime2

irank = GP%irank
call CPU_TIME(CPUTime1) 

if(irank.lt.LP(iLayer)%nsize) then

nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)

!/// check if dynamic pressure at each cell is zero ///!
do j = nstarty,nendy
do i = nstartx,nendx
    !/// Dynamic pressure in outermost cells of the first layer is set zero as boundary condition ///! 
    if(LP(iLayer)%Level.eq.1.and.(i.eq.1.or.i.eq.LP(iLayer)%NX.or.j.eq.1.or.j.eq.LP(iLayer)%NY)) then
        LP(iLayer)%Has_Q(i,j) = 0
    !/// Dynamic pressure is zero in too shallow area or on the wet-dry boundary ///!        
    elseif(LP(iLayer)%H(2,i,j)+LP(iLayer)%Z(i,j).le.GP%MinDispersionDepth) then
        LP(iLayer)%Has_Q(i,j) = 0
    elseif(i.ge.2.and.(LP(iLayer)%D_M(i-1,j).le.GP%MinDispersionDepth.or.LP(iLayer)%Sign_M(i-1,j).ne.0)) then
        LP(iLayer)%Has_Q(i,j) = 0
    elseif(i.le.LP(iLayer)%NX-1.and.(LP(iLayer)%D_M(i,j).le.GP%MinDispersionDepth.or.LP(iLayer)%Sign_M(i,j).ne.0)) then
        LP(iLayer)%Has_Q(i,j) = 0
    elseif(j.ge.2.and.(LP(iLayer)%D_N(i,j-1).le.GP%MinDispersionDepth.or.LP(iLayer)%Sign_N(i,j-1).ne.0)) then
        LP(iLayer)%Has_Q(i,j) = 0
    elseif(j.le.LP(iLayer)%NY-1.and.(LP(iLayer)%D_N(i,j).le.GP%MinDispersionDepth.or.LP(iLayer)%Sign_N(i,j).ne.0)) then
        LP(iLayer)%Has_Q(i,j) = 0
    else
        LP(iLayer)%Has_Q(i,j) = 1
    endif
enddo
enddo

endif ! if this node is used by this layer

call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1 

end subroutine checkZeroDynamicPressure



subroutine constructDynamicPressureEquations(GP, LP, iLayer)

use VariableDefination
implicit NONE
type(GlobalParameters)    ::  GP
type(LayerParameters)     ::  LP(100)
integer*4   ::  iLayer, irank
integer*4   ::  nstartx, nendx, nstarty, nendy, i, j
real*8      ::  A(LP(iLayer)%NX-1, LP(iLayer)%NY), B(LP(iLayer)%NX, LP(iLayer)%NY-1)
real*8      ::  alpha, beta
real*8      ::  h1, h2, z1, z2, D, u1, u2, v1, v2, ws1, wb1, wb2
real*8      ::  csy1, csy2
real*8      ::  CPUTime1, CPUTime2

irank = GP%irank
call CPU_TIME(CPUTime1)

!/// use different governing equations based on rate of change of water depth ///!
if(LP(iLayer)%DepthVariability.eq.0) then
    beta = 2.0/3.0; alpha = 0.5
else
    beta = 0.5; alpha = 1.0
endif

if(irank.lt.LP(iLayer)%nsize) then

nstartx = MAX(LP(iLayer)%PartitionInfo(irank+1,1)-GP%nRowBoundary,1)
nendx   = MIN(LP(iLayer)%PartitionInfo(irank+1,2)+GP%nRowBoundary,LP(iLayer)%NX)
nstarty = MAX(LP(iLayer)%PartitionInfo(irank+1,3)-GP%nRowBoundary,1)
nendy   = MIN(LP(iLayer)%PartitionInfo(irank+1,4)+GP%nRowBoundary,LP(iLayer)%NY)

!/// calculate coefficients A and B before constructing the matrix ///! 
do j = nstarty, nendy
do i = nstartx, nendx
    if(i.le.nendx-1) then !this cell has A(i,j)
    if(LP(iLayer)%Sign_M(i,j).ne.0.or.LP(iLayer)%D_M(i,j).le.GP%MinDispersionDepth) then
        A(i,j) = 0.0d0
    else
        h1 = LP(iLayer)%H(2,i,j); h2 = LP(iLayer)%H(2,i+1,j)
        z1 = LP(iLayer)%Z(i,j);   z2 = LP(iLayer)%Z(i+1,j); D = LP(iLayer)%D_M(i,j)
        A(i,j) = ((h2-alpha*z2) - (h1-alpha*z1))/D
    endif
    endif

    if(j.le.nendy-1) then !this cell has B(i,j)
    if(LP(iLayer)%Sign_N(i,j).ne.0.or.LP(iLayer)%D_N(i,j).le.GP%MinDispersionDepth) then
        B(i,j) = 0.0d0
    else
        h1 = LP(iLayer)%H(2,i,j); h2 = LP(iLayer)%H(2,i,j+1)
        z1 = LP(iLayer)%Z(i,j);   z2 = LP(iLayer)%Z(i,j+1); D = LP(iLayer)%D_N(i,j)
        B(i,j) = ((h2-alpha*z2) - (h1-alpha*z1))/D
    endif
    endif
enddo
enddo

!/// calculate the value of coefficients in the linear equations of dynamic pressure
nstartx = LP(iLayer)%PartitionInfo(irank+1,1); nendx = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3); nendy = LP(iLayer)%PartitionInfo(irank+1,4)
do j = nstarty, nendy
do i = nstartx, nendx
    if(LP(iLayer)%Has_Q(i,j).eq.0) then ! dynamic pressure is forced to be zero (boundary of top layer, some inner cells)
        LP(iLayer)%PL(i,j) = 0.0d0; LP(iLayer)%PR(i,j) = 0.0d0
        LP(iLayer)%PB(i,j) = 0.0d0; LP(iLayer)%PT(i,j) = 0.0d0
        LP(iLayer)%PC(i,j) = 1.0d0; LP(iLayer)%PQ(i,j) = 0.0d0
    elseif(LP(iLayer)%Level.gt.1.and.(i.le.GP%nGHOST.or.i.ge.LP(iLayer)%NX-GP%nGHOST+1 &
            .or.j.le.GP%nGHOST.or.j.ge.LP(iLayer)%NY-GP%nGHOST+1)) then ! dynamic pressure is forced to be the value interpolated from parent layer
        LP(iLayer)%PL(i,j) = 0.0d0; LP(iLayer)%PR(i,j) = 0.0d0                
        LP(iLayer)%PB(i,j) = 0.0d0; LP(iLayer)%PT(i,j) = 0.0d0
        LP(iLayer)%PC(i,j) = 1.0d0; LP(iLayer)%PQ(i,j) = LP(iLayer)%Q(i,j)
    else
        D = LP(iLayer)%H(2,i,j) + LP(iLayer)%Z(i,j)
        u1 = LP(iLayer)%M(2,i-1,j)/LP(iLayer)%D_M(i-1,j); u2 = LP(iLayer)%M(2,i,j)/LP(iLayer)%D_M(i,j)
        v1 = LP(iLayer)%N(2,i,j-1)/LP(iLayer)%D_N(i,j-1); v2 = LP(iLayer)%N(2,i,j)/LP(iLayer)%D_N(i,j)
        ws1 = LP(iLayer)%WS(1,i,j); wb1 = LP(iLayer)%WB(1,i,j); wb2 = LP(iLayer)%WB(2,i,j)
        if(GP%CoordinatesType.eq.1) then
            LP(iLayer)%PL(i,j) = beta*LP(iLayer)%CCD1*(-1+0.5*A(i-1,j))
            LP(iLayer)%PR(i,j) = beta*LP(iLayer)%CCD1*(-1-0.5*A(i,j))
            LP(iLayer)%PB(i,j) = beta*LP(iLayer)%CCD2*(-1+0.5*B(i,j-1))
            LP(iLayer)%PT(i,j) = beta*LP(iLayer)%CCD2*(-1-0.5*B(i,j))
            LP(iLayer)%PC(i,j) = beta*LP(iLayer)%CCD1*(1+0.5*A(i-1,j)+1-0.5*A(i,j)) + &
                beta*LP(iLayer)%CCD2*(1+0.5*B(i,j-1)+1-0.5*B(i,j)) + 2.0*LP(iLayer)%dt/D**2
            LP(iLayer)%PQ(i,j) = -LP(iLayer)%CCD3*(u2-u1) - LP(iLayer)%CCD4*(v2-v1) - (ws1+wb1-2.0*wb2)/D
        elseif(GP%CoordinatesType.eq.0) then
            csy1 = LP(iLayer)%CSY(j-1); csy2 = LP(iLayer)%CSY(j)
            LP(iLayer)%PL(i,j) = beta*LP(iLayer)%CSD1(j)*(-1+0.5*A(i-1,j))
            LP(iLayer)%PR(i,j) = beta*LP(iLayer)%CSD1(j)*(-1-0.5*A(i,j))
            LP(iLayer)%PB(i,j) = beta*LP(iLayer)%CSD2(j)*(-1+0.5*B(i,j-1))*csy1
            LP(iLayer)%PT(i,j) = beta*LP(iLayer)%CSD2(j)*(-1-0.5*B(i,j))*csy2
            LP(iLayer)%PC(i,j) = beta*LP(iLayer)%CSD1(j)*(1+0.5*A(i-1,j)+1-0.5*A(i,j)) + &
                beta*LP(iLayer)%CSD2(j)*((1+0.5*B(i,j-1))*csy1+(1-0.5*B(i,j))*csy2) + 2.0*LP(iLayer)%dt/D**2
            LP(iLayer)%PQ(i,j) = -LP(iLayer)%CSD3(j)*(u2-u1) - LP(iLayer)%CSD4(j)*(v2*csy2-v1*csy1) - (ws1+wb1-2.0*wb2)/D
        endif
    endif
enddo
enddo

endif !if this rank is used

call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1

end subroutine constructDynamicPressureEquations



subroutine bottomVerticalVelocity(GP, LP, iLayer)

use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer, irank
integer*4   ::  nstartx, nendx, nstarty, nendy
integer*4   ::  i, j
real*8      ::  m1, m2, n1, n2, z0, z1, z2, z3, z4, zl, zr, zb, zt
real*8      ::  D0, U0, V0, Dlimit  
real*8      ::  CPUTime1, CPUTime2

irank = GP%irank
call CPU_TIME(CPUTime1); Dlimit = 1e-3

if(irank.lt.LP(iLayer)%nsize) then

nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
if(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids on child layers
    nstartx = MAX(nstartx,GP%nGHOST+1)
    nendx   = MIN(nendx,LP(iLayer)%NX-GP%nGHOST)
    nstarty = MAX(nstarty,GP%nGHOST+1)
    nendy   = MIN(nendy,LP(iLayer)%NY-GP%nGHOST)
endif

!/// calculate bottom vertical velocity at wet cells ///!
do j = nstarty, nendy
do i = nstartx, nendx
    z0 = LP(iLayer)%Z(i,j); D0 = LP(iLayer)%Z(i,j) + LP(iLayer)%H(2,i,j)
    if(D0.le.0.0) then !exclude dry cells
        LP(iLayer)%WB(2,i,j) = 0.0d0; cycle
    endif
    U0 = 0.0d0; V0 = 0.0d0
    if(i.ge.2) then
        z1 = LP(iLayer)%Z(i-1,j); m1 = LP(iLayer)%M(1,i-1,j)   
    else
        z1 = z0; m1 = 0.0d0
    endif
    if(i.le.LP(iLayer)%NX-1) then
        z2 = LP(iLayer)%Z(i+1,j); m2 = LP(iLayer)%M(1,i,j)   
    else
        z2 = z0; m2 = 0.0d0
    endif
    if(j.ge.2) then
        z3 = LP(iLayer)%Z(i,j-1); n1 = LP(iLayer)%N(1,i,j-1)
    else
        z3 = z0; n1 = 0.0d0
    endif
    if(j.le.LP(iLayer)%NY-1) then
        z4 = LP(iLayer)%Z(i,j+1); n2 = LP(iLayer)%N(1,i,j)
    else
        z4 = z0; n2 = 0.0d0
    endif
    if(D0.gt.Dlimit) then
        U0 = 0.5*(m1+m2)/D0; V0 = 0.5*(n1+n2)/D0
    endif

    !/// estimate bottom slope with upwind scheme ///!
    if(U0.ge.0.0) then
        zl = z1; zr = z0
    else
        zl = z0; zr = z2
    endif
    if(V0.ge.0.0) then
        zb = z3; zt = z0
    else
        zb = z0; zt = z4
    endif

    if(GP%CoordinatesType.eq.1) then
        LP(iLayer)%WB(2,i,j) = -U0*LP(iLayer)%CCD3*(zr-zl)-V0*LP(iLayer)%CCD4*(zt-zb) 
    elseif(GP%CoordinatesType.eq.0) then
        LP(iLayer)%WB(2,i,j) = -U0*LP(iLayer)%CSD3(j)*(zr-zl)-V0*LP(iLayer)%CSD5*(zt-zb) 
    endif
enddo
enddo

endif !if: irank is used by this layer

call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1

end subroutine bottomVerticalVelocity



subroutine surfaceVerticalVelocity(GP, LP, iLayer, iTimeStep)

use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer, iTimeStep, irank
integer*4   ::  nstartx, nendx, nstarty, nendy
integer*4   ::  i, j
real*8      ::  D0, DM1, DM2, DN1, DN2, Dlimit
real*8      ::  m1, m2, n1, n2, u1, u2, v1, v2
real*8      ::  csy1, csy2 
real*8      ::  CPUTime1, CPUTime2

irank = GP%irank
call CPU_TIME(CPUTime1); Dlimit = 1e-3

if(irank.lt.LP(iLayer)%nsize) then

nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
if(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids on child layers
    nstartx = MAX(nstartx, GP%nGHOST+1)
    nendx   = MIN(nendx,LP(iLayer)%NX-GP%nGHOST)
    nstarty = MAX(nstarty, GP%nGHOST+1)
    nendy   = MIN(nendy, LP(iLayer)%NY-GP%nGHOST)
endif

!/// Dynamic pressure cannot be calculated at first major step due to unknown initial value of surface vertical velocity ///!
!/// ws is estimated from 3D continuity equation at first major step, and obtained directly from dynamic pressure at subsequent steps ///!
if(iTimeStep.gt.1) then 
do j = nstarty, nendy
do i = nstartx, nendx
    D0 = LP(iLayer)%Z(i,j) + LP(iLayer)%H(2,i,j)
    if(D0.le.0.0) then
        LP(iLayer)%WS(2,i,j) = 0.0d0; cycle
    endif
    LP(iLayer)%WS(2,i,j) = LP(iLayer)%WS(1,i,j) + LP(iLayer)%WB(1,i,j) - LP(iLayer)%WB(2,i,j)
    if(D0.gt.Dlimit) LP(iLayer)%WS(2,i,j) = LP(iLayer)%WS(2,i,j) + 2.0*LP(iLayer)%dt*LP(iLayer)%Q(i,j)/D0
enddo
enddo

elseif(iTimeStep.eq.1) then 
do j = nstarty, nendy
do i = nstartx, nendx
    D0 = LP(iLayer)%Z(i,j) + LP(iLayer)%H(2,i,j)
    if(D0.le.0.0) then
        LP(iLayer)%WS(2,i,j) = 0.0d0; cycle
    endif
    u1 = 0.0d0; u2 = 0.0d0; v1 = 0.0d0; v2 = 0.0d0
    if(i.ge.2) then
        m1 = LP(iLayer)%M(2,i-1,j); DM1 = LP(iLayer)%D_M(i-1,j)
    else
        m1 = 0.0d0; DM1 = D0
    endif
    if(i.le.LP(iLayer)%NX-1) then
        m2 = LP(iLayer)%M(2,i,j); DM2 = LP(iLayer)%D_M(i,j)
    else
        m2 = 0.0d0; DM2 = D0
    endif
    if(j.ge.2) then
        n1 = LP(iLayer)%N(2,i,j-1); DN1 = LP(iLayer)%D_N(i,j-1)
    else
        n1 = 0.0d0; DN1 = D0
    endif
    if(j.le.LP(iLayer)%NY-1) then
        n2 = LP(iLayer)%N(2,i,j); DN2 = LP(iLayer)%D_N(i,j)
    else
        n2 = 0.0d0; DN2 = D0
    endif
    if(GP%CoordinatesType.eq.0) then
        if(j.ge.2) then
            csy1 = LP(iLayer)%CSY(j-1)
        else
            csy1 = COS((LP(iLayer)%ymin-0.5*LP(iLayer)%dy)*GP%PI/180.0)
        endif
        if(j.le.LP(iLayer)%NY-1) then
            csy2 = LP(iLayer)%CSY(j)
        else
            csy2 = COS((LP(iLayer)%ymax+0.5*LP(iLayer)%dy)*GP%PI/180.0)
        endif
    endif

    if(DM1.gt.Dlimit) u1 = m1/DM1; if(DM2.gt.Dlimit) u2 = m2/DM2
    if(DN1.gt.Dlimit) v1 = n1/DN1; if(DN2.gt.Dlimit) v2 = n2/DN2
    if(GP%CoordinatesType.eq.1) then
        LP(iLayer)%WS(2,i,j) = LP(iLayer)%WB(2,i,j) - D0*(LP(iLayer)%CCD3*(u2-u1)+LP(iLayer)%CCD4*(v2-v1))
    elseif(GP%CoordinatesType.eq.0) then
        LP(iLayer)%WS(2,i,j) = LP(iLayer)%WB(2,i,j) - D0*(LP(iLayer)%CSD3(j)*(u2-u1)+LP(iLayer)%CSD4(j)*(v2*csy2-v1*csy1))
    endif
enddo
enddo
endif !if calculation dynamic pressure has started

endif !if: irank is used by this layer

call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1

end subroutine surfaceVerticalVelocity




subroutine getLayerBoundaryDynamicPressureFromParent(GP, LP, iLayer, Localdata, LocalDataLength)

use mpi
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer
integer*4   ::  LocalDataLength
real*8      ::  LocalData(LocalDataLength)
integer*4   ::  irank, nsize, master
integer*4   ::  ierror, istatus(MPI_STATUS_SIZE)
integer*4   ::  pLayer
integer*4   ::  iPC, iNodeFrom, iNodeTo, iBoundary, iHMN
integer*4   ::  istart, iend, jstart, jend, i, j
real*8      ::  x, y, val
real*8      ::  CPUTime1, CPUTime2     
    
irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
    
!/// For each major step iTimeStep, store initial value of dynamic pressure in Q0, and final value in QF
pLayer = LP(iLayer)%Parent; LP(iLayer)%Q0 = LP(iLayer)%Q; LP(iLayer)%QF = 0.0d0
if(LP(pLayer)%Dispersion.eq.0) return
    
call CPU_TIME(CPUTime1)
    
do iPC = 1,LP(iLayer)%ParentToChildSendRecvCount
    iNodeFrom = LP(iLayer)%ParentToChildSendRecv(iPC,1)
    iNodeTo   = LP(iLayer)%ParentToChildSendRecv(iPC,2)
    iBoundary = LP(iLayer)%ParentToChildSendRecv(iPC,3)
    iHMN      = LP(iLayer)%ParentToChildSendRecv(iPC,4)
    istart    = LP(iLayer)%ParentToChildSendRecv(iPC,5)
    iend      = LP(iLayer)%ParentToChildSendRecv(iPC,6)
    jstart    = LP(iLayer)%ParentToChildSendRecv(iPC,7)
    jend      = LP(iLayer)%ParentToChildSendRecv(iPC,8)
    
    !/// use the ParentToChildSendRecv table of H for Q ///!
    if(iHMN.eq.1) then
        !/// interpolate final value of Q on iNodeFrom ///!
        if(irank.eq.iNodeFrom) then
            do j = jstart,jend
            do i = istart,iend
                x = LP(iLayer)%X(i); y = LP(iLayer)%Y(j)
                call interpData(GP,LP,pLayer,4,x,y,val)
                LP(iLayer)%QF(i,j) = val
                if(iNodeFrom.eq.iNodeTo) then
                    if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit.or. &
                        LP(iLayer)%Z(i,j)+LP(iLayer)%HF(i,j).le.0.0) then
                            LP(iLayer)%QF(i,j) = 0.0d0
                    endif
                endif
            enddo
            enddo
        endif
    
        !/// send interpolated dynamic pressure value from iNodeFrom to iNodeTo ///!
        if(irank.eq.iNodeFrom.and.iNodeTo.ne.iNodeFrom) then
            do j = jstart,jend
            do i = istart,iend
                LocalData((i-istart)*(jend-jstart+1)+j-jstart+1) = LP(iLayer)%QF(i,j)
            enddo
            enddo
            call MPI_SEND(LocalData,(iend-istart+1)*(jend-jstart+1),MPI_DOUBLE_PRECISION,&
                iNodeTo,2022,MPI_COMM_WORLD,ierror)
        endif
    
        if(irank.eq.iNodeTo.and.iNodeTo.ne.iNodeFrom) then
            call MPI_RECV(LocalData,(iend-istart+1)*(jend-jstart+1),MPI_DOUBLE_PRECISION,&
                iNodeFrom,2022,MPI_COMM_WORLD,istatus,ierror)
            do j = jstart,jend
            do i = istart,iend
                LP(iLayer)%QF(i,j) = LocalData((i-istart)*(jend-jstart+1)+j-jstart+1)
                if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit.or. &
                    LP(iLayer)%Z(i,j)+LP(iLayer)%HF(i,j).le.0.0) then
                        LP(iLayer)%QF(i,j) = 0.0d0
                endif
            enddo
            enddo
        endif
    endif !if this SendRecv table is for H
enddo !loop for all ParenToChildSendRecv
    
call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,4) = GP%CPUTime(irank+1,4)+CPUTime2-CPUTime1
    
end subroutine getLayerBoundaryDynamicPressureFromParent
    
    

    
subroutine getLayerBoundaryDynamicPressureAtFineTimeStep(GP, LP, iLayer, iStep)
    
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer, iStep, irank
integer*4   ::  iPC, iNodeFrom, iNodeTo, iBoundary, iHMN
integer*4   ::  istart, iend, jstart, jend, i, j
real*8      ::  t1, t2, t, deltat
    
irank = GP%irank
if(irank.lt.LP(iLayer)%nsize) then
    
t1 = 0.0d0; t2 = GP%dt
deltat = 1.0d0/(t2-t1); t = iStep*LP(iLayer)%dt
    
do iPC = 1,LP(iLayer)%ParentToChildSendRecvCount
    iNodeFrom = LP(iLayer)%ParentToChildSendRecv(iPC,1)
    iNodeTo   = LP(iLayer)%ParentToChildSendRecv(iPC,2)
    iBoundary = LP(iLayer)%ParentToChildSendRecv(iPC,3)
    iHMN      = LP(iLayer)%ParentToChildSendRecv(iPC,4)
    istart    = LP(iLayer)%ParentToChildSendRecv(iPC,5)
    iend      = LP(iLayer)%ParentToChildSendRecv(iPC,6)
    jstart    = LP(iLayer)%ParentToChildSendRecv(iPC,7)
    jend      = LP(iLayer)%ParentToChildSendRecv(iPC,8)
    
    if(iHMN.eq.1.and.irank.eq.iNodeTo) then
    do j = jstart,jend
    do i = istart,iend
        LP(iLayer)%Q(i,j) = LP(iLayer)%Q0(i,j) + &
            (LP(iLayer)%QF(i,j)-LP(iLayer)%Q0(i,j))*deltat*(t-t1)
        if(LP(iLayer)%Z(i,j)+LP(iLayer)%H(2,i,j).le.GP%MinDispersionDepth) LP(iLayer)%Q(i,j) = 0.0d0
    enddo
    enddo
    endif !if this node has the boundary of dynamic pressure
    
enddo ! loop for all ParentToChildSendRecv
endif ! if this node is used by this layer
    
end subroutine getLayerBoundaryDynamicPressureAtFineTimeStep



subroutine bcastComputeDomainBoundaryDynamicPressure(GP, LP, iLayer, LocalData, LocalDataLength)

use mpi
use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(LayerParameters)    ::  LP(100)
integer*4  ::  iLayer, LocalDataLength
real*8     ::  LocalData(LocalDataLength)
integer*4  ::  irank, nsize, master
integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
integer*4  ::  iBC, i, j
integer*4  ::  nstartx, nendx, nstarty, nendy
real*8     ::  CPUTime1, CPUTime2

irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
call CPU_TIME(CPUTime1)

do iBC = 1,LP(iLayer)%BoundarySendRecvCount
    if(irank.eq.LP(iLayer)%BoundarySendRecv(iBC, 1)) then
        nstartx = LP(iLayer)%BoundarySendRecv(iBC,3)
        nendx   = LP(iLayer)%BoundarySendRecv(iBC,4)
        nstarty = LP(iLayer)%BoundarySendRecv(iBC,5)
        nendy   = LP(iLayer)%BoundarySendRecv(iBC,6)
        if(nstartx.le.nendx.and.nstarty.le.nendy) then
        do j = nstarty, nendy
        do i = nstartx, nendx
            LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%Q(i,j)
        enddo
        enddo
        call MPI_SEND(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1),MPI_DOUBLE_PRECISION, &
            LP(iLayer)%BoundarySendRecv(iBC,2),iLayer,MPI_COMM_WORLD,ierror)
        endif
    elseif(irank.eq.LP(iLayer)%BoundarySendRecv(iBC,2)) then
        nstartx = LP(iLayer)%BoundarySendRecv(iBC,3)
        nendx   = LP(iLayer)%BoundarySendRecv(iBC,4)
        nstarty = LP(iLayer)%BoundarySendRecv(iBC,5)
        nendy   = LP(iLayer)%BoundarySendRecv(iBC,6)
        if(nstartx.le.nendx.and.nstarty.le.nendy) then
        call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1),MPI_DOUBLE_PRECISION, &
            LP(iLayer)%BoundarySendRecv(iBC,1),iLayer,MPI_COMM_WORLD,istatus,ierror)
        do j = nstarty, nendy
        do i = nstartx, nendx
            LP(iLayer)%Q(i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
        enddo
        enddo
        endif
    endif
enddo
call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,4) = GP%CPUTime(irank+1,4)+CPUTime2-CPUTime1

end subroutine bcastComputeDomainBoundaryDynamicPressure

