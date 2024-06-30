subroutine momentumNonlinearCartesian(GP, LP, iLayer)

use VariableDefination
implicit NONE
type(GlobalParameters)    ::  GP
type(LayerParameters)     ::  LP(100)
integer*4   ::  iLayer, irank
integer*4   ::  nstartx, nendx, nstarty, nendy
integer*4   ::  istart, iend, jstart, jend
integer*4   ::  i, j, k, FluxSign
real*8      ::  z1, z2, h1, h2, MM, NN, Cf
real*8      ::  D0, D1, D2, D3, D5, Dlimit
real*8      ::  m0, m1, m2, m3, m4, m5, m6, m7, m8
real*8      ::  n0, n1, n2, n3, n4, n5, n6, n7, n8
real*8      ::  MU1, MU2, NU1, NU2, NV1, NV2, MV1, MV2, Fx, Fy, phi
real*8      ::  CPUTime1, CPUTime2

irank = GP%irank; call CPU_TIME(CPUTime1)

Dlimit = 1e-3 
if(LP(iLayer)%FluxCenter.eq.0) then
    phi = GP%centerWeighting0
else
    phi = GP%centerWeighting1
endif

if(irank.lt.LP(iLayer)%nsize) then
    
nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)

!****************** calculate flux component M ******************!
if(LP(iLayer)%Level.eq.1) then !calculate all grids in top layer
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX-1,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY,nendy)
elseif(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids in child layer
    istart = MAX(GP%nGHOST+1,nstartx); iend = MIN(LP(iLayer)%NX-GP%nGHOST-1,nendx)
    jstart = MAX(GP%nGHOST+1,nstarty); jend = MIN(LP(iLayer)%NY-GP%nGHOST,nendy)
endif
do j = jstart,jend
do i = istart,iend
    z1 = LP(iLayer)%Z(i,j);   h1 = LP(iLayer)%H(2,i,j)
    z2 = LP(iLayer)%Z(i+1,j); h2 = LP(iLayer)%H(2,i+1,j)

    !/// exclude the cell where water cannot flow ///!
    call checkFluxDirection(GP%PermanentDryLimit, z1, z2, h1, h2, FluxSign)
    if(FluxSign.eq.999) then !water cannot flow between there two cells
        LP(iLayer)%M(2,i,j) = 0.0d0; cycle
    endif

    !/// calculate the value of M by solving the nonlinear equation ///!
    MU1 = 0.0d0; MU2 = 0.0d0; NU1 = 0.0d0; NU2 = 0.0d0
    m0 = LP(iLayer)%M(1,i,j); D0 = LP(iLayer)%D_M(i,j)
    if(i.eq.1) then
        m1 = 0.0d0; D1 = D0
    else
        m1 = LP(iLayer)%M(1,i-1,j); D1 = LP(iLayer)%D_M(i-1,j)
    endif
    if(i.eq.LP(iLayer)%NX-1) then
        m2 = 0.0d0; D2 =  D0
    else
        m2 = LP(iLayer)%M(1,i+1,j); D2 = LP(iLayer)%D_M(i+1,j)
    endif
    if(j.eq.1) then
        m3 = 0.0d0; D3 = D0; n3 = 0.0d0; n4 = 0.0d0
    else
        m3 = LP(iLayer)%M(1,i,j-1); D3 = LP(iLayer)%D_M(i,j-1)
        n3 = LP(iLayer)%N(1,i,j-1); n4 = LP(iLayer)%N(1,i+1,j-1)
    endif
    if(j.eq.LP(iLayer)%NY) then
        m5 = 0.0d0; D5 = D0; n0 = 0.0d0; n2 = 0.0d0
    else
        m5 = LP(iLayer)%M(1,i,j+1); D5 = LP(iLayer)%D_M(i,j+1)
        n0 = LP(iLayer)%N(1,i,j); n2 = LP(iLayer)%N(1,i+1,j)
    endif
    !/// extra terms n5, n6, n7, n8 ///!
    if(j.le.2) then
        n7 = 0.0d0; n8 = 0.0d0
    else
        n7 = LP(iLayer)%N(1,i,j-2); n8 = LP(iLayer)%N(1,i+1,j-2)
    endif
    if(j.ge.LP(iLayer)%NY-1) then
        n5 = 0.0d0; n6 = 0.0d0
    else
        n5 = LP(iLayer)%N(1,i,j+1); n6 = LP(iLayer)%N(1,i+1,j+1)
    endif

    !/// upwind scheme ///!
    if(m0.ge.0.0) then
        if(D1.gt.Dlimit) MU1 = m1**2/D1
        if(D0.gt.Dlimit) MU2 = m0**2/D0
    else
        if(D0.gt.Dlimit) MU1 = m0**2/D0
        if(D2.gt.Dlimit) MU2 = m2**2/D2
    endif
    if(n0+n2+n3+n4.ge.0.0) then
        if(D3.gt.Dlimit) NU1 = 0.25*m3*(n3+n4+n7+n8)/D3
        if(D0.gt.Dlimit) NU2 = 0.25*m0*(n0+n2+n3+n4)/D0
    else
        if(D0.gt.Dlimit) NU1 = 0.25*m0*(n0+n2+n3+n4)/D0
        if(D5.gt.Dlimit) NU2 = 0.25*m5*(n0+n2+n5+n6)/D5
    endif
    LP(iLayer)%M(2,i,j) = phi*m0 + 0.5*(1.0-phi)*(m1+m2) - LP(iLayer)%CC3*D0*(h2-h1) &
                        - LP(iLayer)%CC1*(MU2-MU1) - LP(iLayer)%CC2*(NU2-NU1)

    !/// bottom friction ///!
    if(D0.gt.GP%FrictionDepthLimit) then
        MM = m0**2; NN = (0.25*(n0+n2+n3+n4))**2
        if(GP%BoundaryConditionType.eq.2.and.iLayer.eq.GP%TopLayer) then
            Cf = GP%MANNING + 0.5*(LP(iLayer)%SpongeMANNING(i,j)+LP(iLayer)%SpongeMANNING(i+1,j))
        else
            Cf = GP%MANNING
        endif
        Fx = GP%GRAV*Cf**2/(D0**2.33)*SQRT(MM+NN)*m0
        if(ABS(LP(iLayer)%dt*Fx).gt.ABS(m0)) Fx = m0/LP(iLayer)%dt
        LP(iLayer)%M(2,i,j) = LP(iLayer)%M(2,i,j) - LP(iLayer)%dt*Fx
    endif

    !/// check if the direction of M is right ///!
    if(FluxSign.ne.0) then
        if(LP(iLayer)%M(2,i,j)*FluxSign.lt.0.0) LP(iLayer)%M(2,i,j) = 0.0d0
    endif

enddo
enddo

!********************* calculate flux component N *********************!    
if(LP(iLayer)%Level.eq.1) then !calculate all grids in top layer 
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY-1,nendy)
elseif(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids on child layers
    istart = MAX(GP%nGHOST+1,nstartx); iend = MIN(LP(iLayer)%NX-GP%nGHOST,nendx)
    jstart = MAX(GP%nGHOST+1,nstarty); jend = MIN(LP(iLayer)%NY-GP%nGHOST-1,nendy)
endif
do j = jstart,jend
do i = istart,iend
    z1 = LP(iLayer)%Z(i,j);   h1 = LP(iLayer)%H(2,i,j)
    z2 = LP(iLayer)%Z(i,j+1); h2 = LP(iLayer)%H(2,i,j+1)

    !/// exclude the cell where water cannot flow ///!
    call checkFluxDirection(GP%PermanentDryLimit, z1, z2, h1, h2, FluxSign)
    if(FluxSign.eq.999) then !water cannot flow between these two cells
        LP(iLayer)%N(2,i,j) = 0.0d0; cycle
    endif

    !/// calculate the value of N by solving the nonlinear equation ///!
    NV1 = 0.0d0; NV2 = 0.0d0; MV1 = 0.0d0; MV2 = 0.0d0
    n0 = LP(iLayer)%N(1,i,j); D0 = LP(iLayer)%D_N(i,j)
    if(j.eq.1) then
        n3 = 0.0d0; D3 = D0
    else
        n3 = LP(iLayer)%N(1,i,j-1); D3 = LP(iLayer)%D_N(i,j-1)
    endif
    if(j.eq.LP(iLayer)%NY-1) then
        n5 = 0.0d0; D5 = D0
    else
        n5 = LP(iLayer)%N(1,i,j+1); D5 = LP(iLayer)%D_N(i,j+1)
    endif
    if(i.eq.1) then
        n1 = 0.0d0; D1 = D0; m1 = 0.0d0; m4 = 0.0d0
    else
        n1 = LP(iLayer)%N(1,i-1,j); D1 = LP(iLayer)%D_N(i-1,j)
        m1 = LP(iLayer)%M(1,i-1,j); m4 = LP(iLayer)%M(1,i-1,j+1)
    endif
    if(i.eq.LP(iLayer)%NX) then
        n2 = 0.0d0; D2 = D0; m0 = 0.0d0; m5 = 0.0d0
    else
        n2 = LP(iLayer)%N(1,i+1,j); D2 = LP(ilayer)%D_N(i+1,j)
        m0 = LP(iLayer)%M(1,i,j); m5 = LP(iLayer)%M(1,i,j+1)
    endif
    !/// extra terms m2, m6, m7, m8///!
    if(i.le.2) then
        m7 = 0.0d0; m8 = 0.0d0
    else
        m7 = LP(iLayer)%M(1,i-2,j); m8 = LP(iLayer)%M(1,i-2,j+1)
    endif
    if(i.ge.LP(iLayer)%NX-1) then
        m2 = 0.0d0; m6 = 0.0d0
    else
        m2 = LP(iLayer)%M(1,i+1,j); m6 = LP(iLayer)%M(1,i+1,j+1)
    endif

    !/// upwind scheme ///!
    if(n0.ge.0.0) then
        if(D3.gt.Dlimit) NV1 = n3**2/D3
        if(D0.gt.Dlimit) NV2 = n0**2/D0
    else
        if(D0.gt.Dlimit) NV1 = n0**2/D0
        if(D5.gt.Dlimit) NV2 = n5**2/D5
    endif
    if(m0+m1+m4+m5.ge.0.0) then
        if(D1.gt.Dlimit) MV1 = 0.25*n1*(m1+m4+m7+m8)/D1
        if(D0.gt.Dlimit) MV2 = 0.25*n0*(m0+m1+m4+m5)/D0
    else
        if(D0.gt.Dlimit) MV1 = 0.25*n0*(m0+m1+m4+m5)/D0
        if(D2.gt.Dlimit) MV2 = 0.25*n2*(m0+m2+m5+m6)/D2
    endif                        
    LP(iLayer)%N(2,i,j) = phi*n0 + 0.5*(1-phi)*(n3+n5) - LP(iLayer)%CC4*D0*(h2-h1) &
                        - LP(iLayer)%CC2*(NV2-NV1) - LP(iLayer)%CC1*(MV2-MV1)

    !/// bottom friction ///!
    if(D0.gt.GP%FrictionDepthLimit) then 
        MM = (0.25*(m0+m1+m4+m5))**2; NN = n0**2
        if(GP%BoundaryConditionType.eq.2.and.iLayer.eq.GP%TopLayer) then
            Cf = GP%MANNING + 0.5*(LP(iLayer)%SpongeMANNING(i,j)+LP(iLayer)%SpongeMANNING(i,j+1))
        else
            Cf = GP%MANNING
        endif
        Fy = GP%GRAV*Cf**2/(D0**2.33)*SQRT(MM+NN)*n0
        if(ABS(LP(iLayer)%dt*Fy).gt.ABS(n0)) Fy = n0/LP(iLayer)%dt
        LP(iLayer)%N(2,i,j) = LP(iLayer)%N(2,i,j) - LP(iLayer)%dt*Fy
    endif
            
    !/// check if the direction of N is right ///!
    if(FluxSign.ne.0) then
        if(LP(iLayer)%N(2,i,j)*FluxSign.lt.0.0) LP(iLayer)%N(2,i,j) = 0.0d0
    endif

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

endif !if irank is used 
call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1


end subroutine momentumNonlinearCartesian




subroutine momentumNonlinearSpherical(GP, LP, iLayer)

use VariableDefination
implicit NONE
type(GlobalParameters)    ::  GP
type(LayerParameters)     ::  LP(100)
integer*4   ::  iLayer, irank
integer*4   ::  nstartx, nendx, nstarty, nendy
integer*4   ::  istart, iend, jstart, jend
integer*4   ::  i, j, k, FluxSign
real*8      ::  z1, z2, h1, h2, MM, NN, Cf
real*8      ::  D0, D1, D2, D3, D5, Dlimit
real*8      ::  m0, m1, m2, m3, m4, m5, m6, m7, m8
real*8      ::  n0, n1, n2, n3, n4, n5, n6, n7, n8
real*8      ::  MU1, MU2, NU1, NU2, NV1, NV2, MV1, MV2 
real*8      ::  NU0, NV0, MU0, Fx, Fy, flux, phi
real*8      ::  CPUTime1, CPUTime2
    
irank = GP%irank; call CPU_TIME(CPUTime1)

Dlimit = 1e-3
if(LP(iLayer)%FluxCenter.eq.0) then
    phi = GP%centerWeighting0
else
    phi = GP%centerWeighting1
endif
    
if(irank.lt.LP(iLayer)%nsize) then
        
nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
    
!****************** calculate flux component M ******************!
if(LP(iLayer)%Level.eq.1) then !calculate all grids in top layer
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX-1,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY,nendy)
elseif(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids in child layer
    istart = MAX(GP%nGHOST+1,nstartx); iend = MIN(LP(iLayer)%NX-GP%nGHOST-1,nendx)
    jstart = MAX(GP%nGHOST+1,nstarty); jend = MIN(LP(iLayer)%NY-GP%nGHOST,nendy)
endif
do j = jstart,jend
do i = istart,iend
    z1 = LP(iLayer)%Z(i,j);   h1 = LP(iLayer)%H(2,i,j)
    z2 = LP(iLayer)%Z(i+1,j); h2 = LP(iLayer)%H(2,i+1,j)

    !/// exclude the cell where water cannot flow ///!
    call checkFluxDirection(GP%PermanentDryLimit, z1, z2, h1, h2, FluxSign)
    if(FluxSign.eq.999) then !water cannot flow between these two cells
        LP(iLayer)%M(2,i,j) = 0.0d0; cycle
    endif
    
    !/// calculate the value of M by solving the nonlinear equation ///!
    MU1 = 0.0d0; MU2 = 0.0d0; NU1 = 0.0d0; NU2 = 0.0d0; NU0 = 0.0d0
    m0 = LP(iLayer)%M(1,i,j); D0 = LP(iLayer)%D_M(i,j)
    if(i.eq.1) then
        m1 = 0.0d0; D1 = D0
    else
        m1 = LP(iLayer)%M(1,i-1,j); D1 = LP(iLayer)%D_M(i-1,j)
    endif
    if(i.eq.LP(iLayer)%NX-1) then
        m2 = 0.0d0; D2 =  D0
    else
        m2 = LP(iLayer)%M(1,i+1,j); D2 = LP(iLayer)%D_M(i+1,j)
    endif
    if(j.eq.1) then
        m3 = 0.0d0; D3 = D0; n3 = 0.0d0; n4 = 0.0d0
    else
        m3 = LP(iLayer)%M(1,i,j-1); D3 = LP(iLayer)%D_M(i,j-1)
        n3 = LP(iLayer)%N(1,i,j-1); n4 = LP(iLayer)%N(1,i+1,j-1)
    endif
    if(j.eq.LP(iLayer)%NY) then
        m5 = 0.0d0; D5 = D0; n0 = 0.0d0; n2 = 0.0d0
    else
        m5 = LP(iLayer)%M(1,i,j+1); D5 = LP(iLayer)%D_M(i,j+1)
        n0 = LP(iLayer)%N(1,i,j); n2 = LP(iLayer)%N(1,i+1,j)
    endif
    !/// extra terms n5, n6, n7, n8 ///!
    if(j.le.2) then
        n7 = 0.0d0; n8 = 0.0d0
    else
        n7 = LP(iLayer)%N(1,i,j-2); n8 = LP(iLayer)%N(1,i+1,j-2)
    endif
    if(j.ge.LP(iLayer)%NY-1) then
        n5 = 0.0d0; n6 = 0.0d0
    else
        n5 = LP(iLayer)%N(1,i,j+1); n6 = LP(iLayer)%N(1,i+1,j+1)
    endif

    !/// upwind scheme ///!
    if(m0.ge.0.0) then
        if(D1.gt.Dlimit) MU1 = m1**2/D1
        if(D0.gt.Dlimit) MU2 = m0**2/D0
    else
        if(D0.gt.Dlimit) MU1 = m0**2/D0
        if(D2.gt.Dlimit) MU2 = m2**2/D2
    endif
    if(n0+n2+n3+n4.ge.0.0) then
        if(D3.gt.Dlimit) NU1 = 0.25*m3*(n3+n4+n7+n8)/D3
        if(D0.gt.Dlimit) NU2 = 0.25*m0*(n0+n2+n3+n4)/D0
    else
        if(D0.gt.Dlimit) NU1 = 0.25*m0*(n0+n2+n3+n4)/D0
        if(D5.gt.Dlimit) NU2 = 0.25*m5*(n0+n2+n5+n6)/D5
    endif
    if(D0.gt.Dlimit) NU0 = 0.25*m0*(n0+n2+n3+n4)/D0
    LP(iLayer)%M(2,i,j) = phi*m0 + 0.5*(1.0-phi)*(m1+m2) - LP(iLayer)%CS3(j)*D0*(h2-h1) &
        -LP(iLayer)%CS1(j)*(MU2-MU1) - LP(iLayer)%CS5*(NU2-NU1) + 2.0*LP(iLayer)%CS6(j)*NU0

    !/// Coriolis force ///!
    flux = 0.25*(n0+n2+n3+n4)
    LP(iLayer)%M(2,i,j) = LP(iLayer)%M(2,i,j) + LP(iLayer)%CPX(j)*flux

    !/// bottom friction ///!
    if(D0.gt.GP%FrictionDepthLimit) then
        MM = m0**2; NN = (0.25*(n0+n2+n3+n4))**2
        if(GP%BoundaryConditionType.eq.2.and.iLayer.eq.GP%TopLayer) then
            Cf = GP%MANNING + 0.5*(LP(iLayer)%SpongeMANNING(i,j)+LP(iLayer)%SpongeMANNING(i+1,j))
        else
            Cf = GP%MANNING
        endif
        Fx = GP%GRAV*Cf**2/(D0**2.33)*SQRT(MM+NN)*m0
        if(ABS(LP(iLayer)%dt*Fx).gt.ABS(m0)) Fx = m0/LP(iLayer)%dt
        LP(iLayer)%M(2,i,j) = LP(iLayer)%M(2,i,j) - LP(iLayer)%dt*Fx
    endif           
    
    !/// check if the direction of M is right ///!
    if(FluxSign.ne.0) then
        if(LP(iLayer)%M(2,i,j)*FluxSign.lt.0.0) LP(iLayer)%M(2,i,j) = 0.0d0
    endif
    
enddo
enddo

!********************* calculate flux component N *********************!
if(LP(iLayer)%Level.eq.1) then !calculate all grids in top layer 
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY-1,nendy)
elseif(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids on child layers
    istart = MAX(GP%nGHOST+1,nstartx); iend = MIN(LP(iLayer)%NX-GP%nGHOST,nendx)
    jstart = MAX(GP%nGHOST+1,nstarty); jend = MIN(LP(iLayer)%NY-GP%nGHOST-1,nendy)
endif
do j = jstart,jend
do i = istart,iend
    z1 = LP(iLayer)%Z(i,j);   h1 = LP(iLayer)%H(2,i,j)
    z2 = LP(iLayer)%Z(i,j+1); h2 = LP(iLayer)%H(2,i,j+1)

    !/// exclude the cell where water cannot flow ///!
    call checkFluxDirection(GP%PermanentDryLimit, z1, z2, h1, h2, FluxSign)
    if(FluxSign.eq.999) then !water cannot flow between these two cells
        LP(iLayer)%N(2,i,j) = 0.0d0; cycle
    endif
        
    !/// calculate the value of N by solving the nonlinear equation ///!
    NV1 = 0.0d0; NV2 = 0.0d0; MV1 = 0.0d0; MV2 = 0.0d0; NV0 = 0.0d0; MU0 = 0.0d0
    n0 = LP(iLayer)%N(1,i,j); D0 = LP(iLayer)%D_N(i,j)
    if(j.eq.1) then
        n3 = 0.0d0; D3 = D0
    else
        n3 = LP(iLayer)%N(1,i,j-1); D3 = LP(iLayer)%D_N(i,j-1)
    endif
    if(j.eq.LP(iLayer)%NY-1) then
        n5 = 0.0d0; D5 = D0
    else
        n5 = LP(iLayer)%N(1,i,j+1); D5 = LP(iLayer)%D_N(i,j+1)
    endif
    if(i.eq.1) then
        n1 = 0.0d0; D1 = D0; m1 = 0.0d0; m4 = 0.0d0
    else
        n1 = LP(iLayer)%N(1,i-1,j); D1 = LP(iLayer)%D_N(i-1,j)
        m1 = LP(iLayer)%M(1,i-1,j); m4 = LP(iLayer)%M(1,i-1,j+1)
    endif
    if(i.eq.LP(iLayer)%NX) then
        n2 = 0.0d0; D2 = D0; m0 = 0.0d0; m5 = 0.0d0
    else
        n2 = LP(iLayer)%N(1,i+1,j); D2 = LP(ilayer)%D_N(i+1,j)
        m0 = LP(iLayer)%M(1,i,j); m5 = LP(iLayer)%M(1,i,j+1)
    endif
    !/// extra terms m2, m6,  m7, m8 ///!
    if(i.le.2) then
        m7 = 0.0d0; m8 = 0.0d0
    else
        m7 = LP(iLayer)%M(1,i-2,j); m8 = LP(iLayer)%M(1,i-2,j+1)
    endif
    if(i.ge.LP(iLayer)%NX-1) then
        m2 = 0.0d0; m6 = 0.0d0
    else
        m2 = LP(iLayer)%M(1,i+1,j); m6 = LP(iLayer)%M(1,i+1,j+1)
    endif

    !/// upwind scheme ///!
    if(n0.ge.0.0) then
        if(D3.gt.Dlimit) NV1 = n3**2/D3
        if(D0.gt.Dlimit) NV2 = n0**2/D0
    else
        if(D0.gt.Dlimit) NV1 = n0**2/D0
        if(D5.gt.Dlimit) NV2 = n5**2/D5
    endif
    if(m0+m1+m4+m5.ge.0.0) then
        if(D1.gt.Dlimit) MV1 = 0.25*n1*(m1+m4+m7+m8)/D1
        if(D0.gt.Dlimit) MV2 = 0.25*n0*(m0+m1+m4+m5)/D0
    else
        if(D0.gt.Dlimit) MV1 = 0.25*n0*(m0+m1+m4+m5)/D0
        if(D2.gt.Dlimit) MV2 = 0.25*n2*(m0+m2+m5+m6)/D2
    endif
    if(D0.gt.Dlimit) then
        MU0 = (0.25*(m0+m1+m4+m5))**2/D0; NV0 = n0**2/D0
    endif

    LP(iLayer)%N(2,i,j) = phi*n0 + 0.5*(1.0-phi)*(n3+n5) - LP(iLayer)%CS4*D0*(h2-h1) &
        -LP(iLayer)%CS5*(NV2-NV1) - LP(iLayer)%CS7(j)*(MV2-MV1) - LP(iLayer)%CS8(j)*(MU0-NV0)

    !/// Coriolis force ///!
    flux = 0.25*(m0+m1+m4+m5)
    LP(iLayer)%N(2,i,j) = LP(iLayer)%N(2,i,j) - LP(iLayer)%CPY(j)*flux

    !/// bottom friction ///!
    if(D0.gt.GP%FrictionDepthLimit) then 
        MM = (0.25*(m0+m1+m4+m5))**2; NN = n0**2
        if(GP%BoundaryConditionType.eq.2.and.iLayer.eq.GP%TopLayer) then
            Cf = GP%MANNING + 0.5*(LP(iLayer)%SpongeMANNING(i,j)+LP(iLayer)%SpongeMANNING(i,j+1))
        else
            Cf = GP%MANNING
        endif
        Fy = GP%GRAV*Cf**2/(D0**2.33)*SQRT(MM+NN)*n0
        if(ABS(LP(iLayer)%dt*Fy).gt.ABS(n0)) Fy = n0/LP(iLayer)%dt
        LP(iLayer)%N(2,i,j) = LP(iLayer)%N(2,i,j) - LP(iLayer)%dt*Fy
    endif
            
    !/// check if the direction of N is right ///!
    if(FluxSign.ne.0) then
        if(LP(iLayer)%N(2,i,j)*FluxSign.lt.0.0) LP(iLayer)%N(2,i,j) = 0.0d0
    endif

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
    
endif !if irank is used 
call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1

    
end subroutine momentumNonlinearSpherical