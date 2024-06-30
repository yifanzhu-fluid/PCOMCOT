
    subroutine mass(GP, LP, iLayer, iStep, iTimeStep, LocalData, LocalDataLength)

    use VariableDefination
    implicit NONE
    type(GlobalParameters)    ::  GP
    type(LayerParameters)     ::  LP(100)
    integer*4   ::  iLayer, iStep, iTimeStep, LocalDataLength
    real*8      ::  LocalData(LocalDataLength)

    if(GP%CoordinatesType.eq.1) then
        call massCartesian(GP, LP, iLayer)
    elseif(GP%CoordinatesType.eq.0) then
        call massSpherical(GP, LP, iLayer)
    endif
    !///filter the wave height if needed///!
    if(MOD(iTimeStep,LP(iLayer)%NDTFilterData).eq.0.and.iStep.eq.LP(iLayer)%nStepsPerTimeStep) then
        call bcastComputeDomainBoundaryValue(GP, LP, iLayer, LocalData, LocalDataLength)
        call smoothWaveHeight(GP, LP, iLayer)
    endif

    end subroutine mass



    subroutine momentum(GP, LP, iLayer, iStep, iTimeStep, LocalData, LocalDataLength)

    use VariableDefination
    implicit NONE
    type(GlobalParameters)    ::  GP
    type(LayerParameters)     ::  LP(100)
    integer*4   ::  iLayer, iStep, iTimeStep, LocalDataLength
    real*8      ::  LocalData(LocalDataLength)

    call reconstructFlowDepth(GP, LP, iLayer)
    if(GP%CoordinatesType.eq.1) then
        if(LP(iLayer)%Nonlinearity.eq.0) then
            call momentumLinearCartesian(GP, LP, iLayer)
        elseif(LP(iLayer)%Nonlinearity.eq.1) then
            call momentumNonlinearCartesian(GP, LP, iLayer)
        endif
    elseif(GP%CoordinatesType.eq.0) then
        if(LP(iLayer)%Nonlinearity.eq.0) then
            call momentumLinearSpherical(GP, LP, iLayer)
        elseif(LP(iLayer)%Nonlinearity.eq.1) then
            call momentumNonlinearSpherical(GP, LP, iLayer)
        endif
    endif
    !///add dispersive terms to shallow water solution///!
    if(LP(iLayer)%Dispersion.eq.1) then
        call bcastComputeDomainBoundaryFlux(GP, LP, iLayer, LocalData, LocalDataLength)
        call addDispersion(GP, LP, iLayer, iStep, iTimeStep, LocalData, LocalDataLength)
    endif
    !///add energy dissipation for wave breaking///!
    if(LP(iLayer)%Breaking.eq.1) then
        call addBreaking(GP, LP, iLayer, iStep, iTimeStep, LocalData, LocalDataLength)
    endif
    !///supress unreasonable volume flux with Froude number cap///!
    call bcastComputeDomainBoundaryFlux(GP, LP, iLayer, LocalData, LocalDataLength)
    call froudeCap(GP, LP, iLayer)
    !///filter volume flux if needed///!
    if(MOD(iTimeStep,LP(iLayer)%NDTFilterData).eq.0.and.iStep.eq.LP(iLayer)%nStepsPerTimeStep) then
        call bcastComputeDomainBoundaryFlux(GP, LP, iLayer, LocalData, LocalDataLength)
        call smoothFlux(GP, LP, iLayer)
    endif 

    end subroutine momentum



    subroutine massCartesian(GP, LP, iLayer)

    use VariableDefination
    implicit NONE
    type(GlobalParameters)    ::  GP
    type(LayerParameters)     ::  LP(100)
    integer*4  ::  iLayer, irank
    integer*4  ::  nstartx, nendx, nstarty, nendy
    integer*4  ::  i, j, k
    real*8     ::  m1, m2, n1, n2
    real*8     ::  CPUTime1, CPUTime2

    irank = GP%irank
    call CPU_TIME(CPUTime1)

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

    !//// solve the mass equation ////!
    do j = nstarty, nendy
    do i = nstartx, nendx
        !/// exclude permanent dry cells from equations ///! 
        if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit) then
            LP(iLayer)%H(2,i,j) = 0.0d0; cycle
        endif
        if(i.eq.1) then
            m1 = 0.0d0
        else
            m1 = LP(iLayer)%M(1,i-1,j)
        endif
        if(i.eq.LP(iLayer)%NX) then
            m2 = 0.0d0
        else
            m2 = LP(iLayer)%M(1,i,j)
        endif
        if(j.eq.1) then
            n1 = 0.0d0
        else
            n1 = LP(iLayer)%N(1,i,j-1)
        endif
        if(j.eq.LP(iLayer)%NY) then
            n2 = 0.0d0
        else
            n2 = LP(iLayer)%N(1,i,j)
        endif
        LP(iLayer)%H(2,i,j) = LP(iLayer)%H(1,i,j) - &
            LP(iLayer)%CC1*(m2-m1) - LP(iLayer)%CC2*(n2-n1)
    enddo
    enddo
 
    !//// calculate inundation in flooding zone ///!
    do k = 1,LP(iLayer)%FloodingCellsCount
        i = LP(iLayer)%FloodingCells(k,1)
        j = LP(iLayer)%FloodingCells(k,2)

        !/// wet -> dry ///!
        if(LP(iLayer)%Z(i,j)+LP(iLayer)%H(1,i,j).gt.0.0.and.&
          LP(iLayer)%H(2,i,j)-LP(iLayer)%H(1,i,j).lt.0.0.and.&
            LP(iLayer)%Z(i,j)+LP(iLayer)%H(2,i,j).le.GP%MinWaterDepth) then
            if(LP(iLayer)%Z(i,j).gt.0.0) then
                LP(iLayer)%H(2,i,j) = -LP(iLayer)%Z(i,j)
            else
                LP(iLayer)%H(2,i,j) = 0.0d0
            endif
        endif

        !/// dry -> wet/dry ///
        if(LP(iLayer)%Z(i,j)+LP(iLayer)%H(1,i,j).le.0.0) then
            if(LP(iLayer)%H(2,i,j)-LP(iLayer)%H(1,i,j).gt.0.0) then ! dry -> wet
                LP(iLayer)%H(2,i,j) = LP(iLayer)%H(2,i,j)-LP(iLayer)%H(1,i,j)-LP(iLayer)%Z(i,j)
            else ! dry ->dry
                if(LP(iLayer)%Z(i,j).gt.0.0) then
                    LP(iLayer)%H(2,i,j) = -LP(iLayer)%Z(i,j)
                else
                    LP(iLayer)%H(2,i,j) = 0.0d0
                endif
            endif
        endif
    enddo

    !/// store maximum and minimum water height in this layer ///!
    nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
    nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
    nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
    nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
    do j = nstarty, nendy
    do i = nstartx, nendx
        LP(iLayer)%Hmax(i,j) = MAX(LP(iLayer)%Hmax(i,j), LP(iLayer)%H(2,i,j))
        LP(iLayer)%Hmin(i,j) = MIN(LP(iLayer)%Hmin(i,j), LP(iLayer)%H(2,i,j))
    enddo
    enddo

    endif !if: irank is used

    call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1

    end subroutine massCartesian




    subroutine massSpherical(GP, LP, iLayer)

    use VariableDefination
    implicit NONE
    type(GlobalParameters)    ::  GP
    type(LayerParameters)     ::  LP(100)
    integer*4  ::  iLayer, irank
    integer*4  ::  nstartx, nendx, nstarty, nendy
    integer*4  ::  i, j, k
    real*8     ::  m1, m2, n1, n2, tmp1, tmp2
    real*8     ::  CPUTime1, CPUTime2

    irank = GP%irank
    call CPU_TIME(CPUTime1)

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

    !//// solve the mass equation ////!
    do j = nstarty,nendy
    do i = nstartx,nendx
        !/// exclude permanent dry cells from equations ///! 
        if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit) then
            LP(iLayer)%H(2,i,j) = 0.0; cycle
        endif
        if(i.eq.1) then
            m1 = 0.0d0
        else
            m1 = LP(iLayer)%M(1,i-1,j)
        endif
        if(i.eq.LP(iLayer)%NX) then
            m2 = 0.0d0
        else
            m2 = LP(iLayer)%M(1,i,j)
        endif
        if(j.eq.1) then
            n1 = 0.0d0; tmp1 = COS((LP(iLayer)%ymin-0.5*LP(iLayer)%dy)*GP%PI/180.0)
        else
            n1 = LP(iLayer)%N(1,i,j-1); tmp1 = LP(iLayer)%CSY(j-1)
        endif
        if(j.eq.LP(iLayer)%NY) then
            n2 = 0.0d0; tmp2 = COS((LP(iLayer)%ymax+0.5*LP(iLayer)%dy)*GP%PI/180.0)
        else
            n2 = LP(iLayer)%N(1,i,j); tmp2 = LP(iLayer)%CSY(j)
        endif
        LP(iLayer)%H(2,i,j) = LP(iLayer)%H(1,i,j) - LP(iLayer)%CS1(j)*(m2-m1) - &
            LP(iLayer)%CS2(j)*(n2*tmp2-n1*tmp1)
    enddo
    enddo   
     
    !//// calculate inundation in flooding zone ///!
    do k = 1,LP(iLayer)%FloodingCellsCount
        i = LP(iLayer)%FloodingCells(k,1)
        j = LP(iLayer)%FloodingCells(k,2)

        !/// wet -> dry ///!
        if(LP(iLayer)%Z(i,j)+LP(iLayer)%H(1,i,j).gt.0.0.and.&
          LP(iLayer)%H(2,i,j)-LP(iLayer)%H(1,i,j).lt.0.0.and.&
            LP(iLayer)%Z(i,j)+LP(iLayer)%H(2,i,j).le.GP%MinWaterDepth) then
            if(LP(iLayer)%Z(i,j).gt.0.0) then
                LP(iLayer)%H(2,i,j) = -LP(iLayer)%Z(i,j)
            else
                LP(iLayer)%H(2,i,j) = 0.0d0
            endif
        endif

        !/// dry -> wet/dry ///
        if(LP(iLayer)%Z(i,j)+LP(iLayer)%H(1,i,j).le.0.0) then
            if(LP(iLayer)%H(2,i,j)-LP(iLayer)%H(1,i,j).gt.0.0) then ! dry -> wet
                LP(iLayer)%H(2,i,j) = LP(iLayer)%H(2,i,j)-LP(iLayer)%H(1,i,j)-LP(iLayer)%Z(i,j)
            else ! dry ->dry
                if(LP(iLayer)%Z(i,j).gt.0.0) then
                    LP(iLayer)%H(2,i,j) = -LP(iLayer)%Z(i,j)
                else
                    LP(iLayer)%H(2,i,j) = 0.0d0
                endif
            endif
        endif
    enddo

    !/// store maximum water height in this layer ///!
    nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
    nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
    nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
    nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
    do j = nstarty, nendy
    do i = nstartx, nendx
        LP(iLayer)%Hmax(i,j) = MAX(LP(iLayer)%Hmax(i,j), LP(iLayer)%H(2,i,j))
        LP(iLayer)%Hmin(i,j) = MIN(LP(iLayer)%Hmin(i,j), LP(iLayer)%H(2,i,j))
    enddo
    enddo

    endif !if: irank is used

    call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1

    end subroutine massSpherical



    subroutine momentumLinearCartesian(GP, LP, iLayer)

    use VariableDefination
    implicit NONE
    type(GlobalParameters)    ::  GP
    type(LayerParameters)     ::  LP(100)
    integer*4   ::  iLayer, irank
    integer*4   ::  nstartx, nendx, nstarty, nendy
    integer*4   ::  istart, iend, jstart, jend
    integer*4   ::  i, j, k, FluxSign
    real*8      ::  z1, z2, h1, h2
    real*8      ::  m0, m1, m2, m3, m4
    real*8      ::  n0, n1, n2, n3, n4
    real*8      ::  D0, MM, NN, Cf, Fx, Fy, phi
    real*8      ::  CPUTime1, CPUTime2

    irank = GP%irank; call CPU_TIME(CPUTime1)

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
    
    !///calculate flux component M///!
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
        if(FluxSign.eq.999) then
            LP(iLayer)%M(2,i,j) = 0.0d0; cycle
        endif 
        !/// calculate the value of M ///!
        m0 = LP(iLayer)%M(1,i,j)
        if(i.eq.1) then
            m1 = 0.0d0
        else
            m1 = LP(iLayer)%M(1,i-1,j)
        endif
        if(i.eq.LP(iLayer)%NX-1) then
            m2 = 0.0d0
        else
            m2 = LP(iLayer)%M(1,i+1,j)
        endif
        LP(iLayer)%M(2,i,j) = phi*m0 + 0.5*(1.0-phi)*(m1+m2) - LP(iLayer)%CC3*LP(iLayer)%D_M(i,j)*(h2-h1)

        !/// bottom friction ///!
        D0 = LP(iLayer)%D_M(i,j)
        if(j.ge.2) then
            n1 = LP(iLayer)%N(1,i,j-1); n2 = LP(iLayer)%N(1,i+1,j-1)
        else
            n1 = 0.0d0; n2 = 0.0d0
        endif
        if(j.le.LP(iLayer)%NY-1) then
            n3 = LP(iLayer)%N(1,i,j); n4 = LP(iLayer)%N(1,i+1,j)
        else
            n3 = 0.0d0; n4 = 0.0d0
        endif
        if(D0.gt.GP%FrictionDepthLimit) then
            MM = m0**2; NN = (0.25*(n1+n2+n3+n4))**2
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

    !///calculate flux component N///!
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
        if(FluxSign.eq.999) then
            LP(iLayer)%N(2,i,j) = 0.0d0; cycle
        endif
        !/// calculate the value of N ///!
        n0 = LP(iLayer)%N(1,i,j)
        if(j.eq.1) then
            n1 = 0.0d0
        else
            n1 = LP(iLayer)%N(1,i,j-1)
        endif
        if(j.eq.LP(iLayer)%NY-1) then
            n2 = 0.0d0
        else
            n2 = LP(iLayer)%N(1,i,j+1)
        endif
        LP(iLayer)%N(2,i,j) = phi*n0 + 0.5*(1.0-phi)*(n1+n2) - LP(iLayer)%CC4*LP(iLayer)%D_N(i,j)*(h2-h1)

        !/// bottom friction ///!
        D0 = LP(iLayer)%D_N(i,j)
        if(i.ge.2) then
            m1 = LP(iLayer)%M(1,i-1,j); m2 = LP(iLayer)%M(1,i-1,j+1)
        else
            m1 = 0.0d0; m2 = 0.0d0
        endif
        if(i.le.LP(iLayer)%NX-1) then
            m3 = LP(iLayer)%M(1,i,j); m4 = LP(iLayer)%M(1,i,j+1)
        else
            m3 = 0.0d0; m4 = 0.0d0
        endif
        if(D0.gt.GP%FrictionDepthLimit) then
            MM = (0.25*(m1+m2+m3+m4))**2; NN = n0**2
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

    end subroutine momentumLinearCartesian



    subroutine momentumLinearSpherical(GP, LP, iLayer)

    use VariableDefination
    implicit NONE
    type(GlobalParameters)    ::  GP
    type(LayerParameters)     ::  LP(100)
    integer*4   ::  iLayer, irank
    integer*4   ::  nstartx, nendx, nstarty, nendy
    integer*4   ::  istart, iend, jstart, jend
    integer*4   ::  i, j, k, FluxSign
    real*8      ::  z1, z2, h1, h2, flux
    real*8      ::  m0, m1, m2, m3, m4
    real*8      ::  n0, n1, n2, n3, n4
    real*8      ::  D0, MM, NN, Cf, Fx, Fy, phi
    real*8      ::  CPUTime1, CPUTime2

    irank = GP%irank; call CPU_TIME(CPUTime1)

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

    !/// calculate flux component M ///!
    if(LP(iLayer)%Level.eq.1) then !calculate all grids on top layer
        istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX-1,nendx)
        jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY,nendy)
    elseif(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids on child layers
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
        !/// calculate the value of M ///!
        m0 = LP(iLayer)%M(1,i,j)
        if(i.eq.1) then
            m1 = 0.0d0
        else
            m1 = LP(iLayer)%M(1,i-1,j)
        endif
        if(i.eq.LP(iLayer)%NX-1) then
            m2 = 0.0d0
        else
            m2 = LP(iLayer)%M(1,i+1,j)
        endif
        LP(iLayer)%M(2,i,j) = phi*m0 + 0.5*(1.0-phi)*(m1+m2) - LP(iLayer)%CS3(j)*LP(iLayer)%D_M(i,j)*(h2-h1)

        !/// Coriolis force ///!
        if(j.eq.1) then
            n1 = 0.0d0; n2 = 0.0d0
        else
            n1 = LP(iLayer)%N(1,i,j-1); n2 = LP(iLayer)%N(1,i+1,j-1)
        endif
        if(j.eq.LP(iLayer)%NY) then
            n3 = 0.0d0; n4 = 0.0d0
        else
            n3 = LP(iLayer)%N(1,i,j); n4 = LP(iLayer)%N(1,i+1,j)
        endif
        flux = 0.25*(n1+n2+n3+n4)
        LP(iLayer)%M(2,i,j) = LP(iLayer)%M(2,i,j) + LP(iLayer)%CPX(j)*flux

        !/// bottom friction ///!
        D0 = LP(iLayer)%D_M(i,j)
        if(D0.gt.GP%FrictionDepthLimit) then
            MM = m0**2; NN = (0.25*(n1+n2+n3+n4))**2
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

    !/// calculate flux component N ///!
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
        !/// calculate the value of N ///!
        n0 = LP(iLayer)%N(1,i,j)
        if(j.eq.1) then
            n1 = 0.0d0
        else
            n1 = LP(iLayer)%N(1,i,j-1)
        endif
        if(j.eq.LP(iLayer)%NY-1) then
            n2 = 0.0d0
        else
            n2 = LP(iLayer)%N(1,i,j+1)
        endif
        LP(iLayer)%N(2,i,j) = phi*n0 + 0.5*(1.0-phi)*(n1+n2) - LP(iLayer)%CS4*LP(iLayer)%D_N(i,j)*(h2-h1)

        !/// Coriolis force ///!
        if(i.eq.1) then
            m1 = 0.0d0; m2 = 0.0d0
        else
            m1 = LP(iLayer)%M(1,i-1,j); m2 = LP(iLayer)%M(1,i-1,j+1)
        endif
        if(i.eq.LP(iLayer)%NX) then
            m3 = 0.0d0; m4 = 0.0d0
        else
            m3 = LP(iLayer)%M(1,i,j); m4 = LP(iLayer)%M(1,i,j+1)
        endif
        flux = 0.25*(m1+m2+m3+m4)
        LP(iLayer)%N(2,i,j) = LP(iLayer)%N(2,i,j) - LP(iLayer)%CPY(j)*flux

        !/// bottom friction ///!
        D0 = LP(iLayer)%D_N(i,j)
        if(D0.gt.GP%FrictionDepthLimit) then
            MM = (0.25*(m1+m2+m3+m4))**2; NN = n0**2
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

    !//// set zero water flux at boundaries of dry cells ////!
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

    end subroutine momentumLinearSpherical



    subroutine reconstructFlowDepth(GP, LP, iLayer)

    use VariableDefination
    implicit NONE
    type(GlobalParameters)    ::  GP
    type(LayerParameters)     ::  LP(100)
    integer*4   ::  iLayer, irank
    integer*4   ::  nstartx, nendx, nstarty, nendy
    integer*4   ::  i, j, s
    real*8      ::  z1, z2, h1, h2, m0, n0
    real*8      ::  CPUTime1, CPUTime2

    irank = GP%irank
    call CPU_TIME(CPUTime1)

    if(irank.lt.LP(iLayer)%nsize) then

    nstartx = MAX(LP(iLayer)%PartitionInfo(irank+1,1)-GP%nRowBoundary,1)
    nendx   = MIN(LP(iLayer)%PartitionInfo(irank+1,2)+GP%nRowBoundary,LP(iLayer)%NX)
    nstarty = MAX(LP(iLayer)%PartitionInfo(irank+1,3)-GP%nRowBoundary,1)
    nendy   = MIN(LP(iLayer)%PartitionInfo(irank+1,4)+GP%nRowBoundary,LP(iLayer)%NY)

    !/// estimate flow depth for flux component M ///!
    do j = nstarty,nendy
    do i = nstartx,nendx-1
        z1 = LP(iLayer)%Z(i,j);   h1 = LP(iLayer)%H(2,i,j)
        z2 = LP(iLayer)%Z(i+1,j); h2 = LP(iLayer)%H(2,i+1,j)
        m0 = LP(iLayer)%M(1,i,j)
        call checkFluxDirection(GP%PermanentDryLimit,z1,z2,h1,h2,s)
        if(s.eq.999) then
            LP(iLayer)%D_M(i,j) = 0.0d0
        elseif(s.eq.0) then
            if(m0.ge.0.0) then
                LP(iLayer)%D_M(i,j) = (z1+z2)/2.0 + h1
            else
                LP(iLayer)%D_M(i,j) = (z1+z2)/2.0 + h2
            endif
        else
            LP(iLayer)%D_M(i,j) = (z1+z2)/2.0 + MAX(h1,h2)
        endif
        LP(iLayer)%Sign_M(i,j) = s
    enddo
    enddo

    !/// estimate flow depth for flux component N ///!
    do j = nstarty, nendy-1
    do i = nstartx, nendx
        z1 = LP(iLayer)%Z(i,j);   h1 = LP(iLayer)%H(2,i,j)
        z2 = LP(iLayer)%Z(i,j+1); h2 = LP(iLayer)%H(2,i,j+1)
        n0 = LP(iLayer)%N(1,i,j)
        call checkFluxDirection(GP%PermanentDryLimit,z1,z2,h1,h2,s)
        if(s.eq.999) then
            LP(iLayer)%D_N(i,j) = 0.0d0
        elseif(s.eq.0) then
            if(n0.ge.0.0) then
                LP(iLayer)%D_N(i,j) = (z1+z2)/2.0 + h1
            else
                LP(iLayer)%D_N(i,j) = (z1+z2)/2.0 + h2
            endif
        else
            LP(iLayer)%D_N(i,j) = (z1+z2)/2.0 + MAX(h1,h2)
        endif
        LP(iLayer)%Sign_N(i,j) = s
    enddo
    enddo

    endif !if irank is used
    call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1 

    end subroutine reconstructFlowDepth



    subroutine checkFluxDirection(dryLimit,z1,z2,h1,h2,FluxSign)
    ! This subroutine checks the reasonble situations of flux M/N
    ! FluxSign can be {-2,-1,0,1,2,999}, with the same definition as Sign_M/N in LayerParameters
    ! dryLimit      : permanent dry limit defined in GlobalParametes
    ! z1,z2,h1,h2   : water depth and elevation on both sides of flux M/N
    ! 0 :water can flow in both directions; 999:water cannot flow
    ! 1 :water can only flow in positive direction, and from higher to lower
    ! 2 :water can only flow in positive direction, and from lower to higher
    !-1 :water can only flow in negative direction, and from higher to lower
    !-2 :water can only flow in negative direction, and from lower to higher
    ! If water can flow in only one direction, values of h1 and h2 are modified for momentum equation
            
    implicit NONE
    real*8, intent(in)      ::  dryLimit, z1, z2
    real*8, intent(inout)   ::  h1, h2          
    integer*4, intent(out)  ::  FluxSign
        
    if(z1.le.-dryLimit.or.z2.le.-dryLimit) then
        FluxSign = 999; return
    endif
        
    if(z1+h1.le.0.0.or.z2+h2.le.0.0) then
        if(z1+h1.le.0.0.and.z2+h2.le.0.0) then
            FluxSign = 999; return
        elseif(z1.gt.z2.and.z1+h1.gt.0.0.and.z2+h1.le.0.0) then
            FluxSign = 999; return
        elseif(z2.gt.z1.and.z2+h2.gt.0.0.and.z1+h2.le.0.0) then
            FluxSign = 999; return
        endif
    endif
        
    if(z1.le.z2.and.z1+h1.gt.0.0.and.z1+h2.le.0.0) then
        FluxSign = 1; h2 = -z1
    elseif(z1.ge.z2.and.z2+h2.le.0.0.and.z2+h1.gt.0.0) then
        FluxSign = 2; h2 = -z2
    elseif(z2.le.z1.and.z2+h2.gt.0.0.and.z2+h1.le.0.0) then
        FluxSign = -1; h1 = -z2
    elseif(z2.ge.z1.and.z1+h1.le.0.0.and.z1+h2.gt.0.0) then
        FluxSign = -2; h1 = -z1
    else
        FluxSign = 0
    endif
        
    end subroutine checkFluxDirection



    subroutine froudeCap(GP, LP, iLayer)

    use VariableDefination
    implicit NONE
    type(GlobalParameters)    ::  GP
    type(LayerParameters)     ::  LP(100)
    integer*4   ::  iLayer, irank
    integer*4   ::  nstartx, nendx, nstarty, nendy
    integer*4   ::  istart, iend, jstart, jend, i, j
    real*8      ::  m0, n0, u0, v0, Utotal, Ulimit, u_cap, v_cap, theta
    real*8      ::  Ma(LP(iLayer)%NX-1,LP(iLayer)%NY), Na(LP(iLayer)%NX,LP(iLayer)%NY-1)
    real*8      ::  CPUTime1, CPUTime2
        
    irank = GP%irank
    call CPU_TIME(CPUTime1)
                
    if(irank.lt.LP(iLayer)%nsize) then
        
    Ma = LP(iLayer)%M(2,:,:); Na = LP(iLayer)%N(2,:,:)              
    nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
    nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
    nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
    nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
        
    !/// apply Froude cap to flux component M///!
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX-1,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY,nendy)
    do j = jstart,jend
    do i = istart,iend
        !///leave the value on the outermost cells unchanged///!
        if(i.eq.1.or.i.eq.LP(iLayer)%NX-1.or.j.eq.1.or.j.eq.LP(iLayer)%NY) cycle
        !///leave the value on the wet-dry boundary unchanged///!
        if(LP(iLayer)%Sign_M(i,j).ne.0) cycle
        m0 = Ma(i,j)
        n0 = 0.25*(Na(i,j-1)+Na(i+1,j-1)+Na(i,j)+Na(i+1,j)) 
        u0 = m0/MAX(LP(iLayer)%D_M(i,j),GP%MinWaterDepth)
        v0 = n0/MAX(LP(iLayer)%D_M(i,j),GP%MinWaterDepth)
        Utotal = SQRT(u0*u0+v0*v0)
        Ulimit = GP%FrCap*SQRT(GP%GRAV*MAX(LP(iLayer)%D_M(i,j),GP%MinWaterDepth))
        if(Utotal.gt.Ulimit) then
            theta = ATAN2(v0,u0)
            u_cap = Ulimit*COS(theta)
            LP(iLayer)%M(2,i,j) = u_cap*MAX(LP(iLayer)%D_M(i,j),GP%MinWaterDepth)
        endif
    enddo
    enddo
        
    !///apply Froude cap to flux component N///!
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY-1,nendy)
    do j = jstart,jend
    do i = istart,iend
        !///leave the value on the outermost cells unchanged///!
        if(i.eq.1.or.i.eq.LP(iLayer)%NX.or.j.eq.1.or.j.eq.LP(iLayer)%NY-1) cycle
        !///leave the value on the wet-dry boundary unchanged///!
        if(LP(iLayer)%Sign_N(i,j).ne.0) cycle
        m0 = 0.25*(Ma(i-1,j)+Ma(i,j)+Ma(i-1,j+1)+Ma(i,j+1))
        n0 = Na(i,j)
        u0 = m0/MAX(LP(iLayer)%D_N(i,j),GP%MinWaterDepth)
        v0 = n0/MAX(LP(iLayer)%D_N(i,j),GP%MinWaterDepth)
        Utotal = SQRT(u0*u0+v0*v0)
        Ulimit = GP%FrCap*SQRT(GP%GRAV*MAX(LP(iLayer)%D_N(i,j),GP%MinWaterDepth))
        if(Utotal.gt.Ulimit) then
            theta = ATAN2(v0,u0)
            v_cap = Ulimit*SIN(theta)
            LP(iLayer)%N(2,i,j) = v_cap*MAX(LP(iLayer)%D_N(i,j),GP%MinWaterDepth)
        endif
    enddo
    enddo
        
    endif !if: this node is used by this layer
        
    call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1
        
    end subroutine froudeCap



    subroutine smoothWaveHeight(GP, LP, iLayer)

    use VariableDefination
    implicit NONE
    type(GlobalParameters)    ::  GP
    type(LayerParameters)     ::  LP(100)
    integer*4  ::  iLayer, irank
    integer*4  ::  nstartx, nendx, nstarty, nendy
    integer*4  ::  i, j
    real*8     ::  h0, h1, h2, h3, h4, h5, h6, h7, h8
    real*8     ::  z0, z1, z2, z3, z4, z5, z6, z7, z8
    real*8     ::  d0, d1, d2, d3, d4, d5, d6, d7, d8
    real*8     ::  HS(LP(iLayer)%NX,LP(iLayer)%NY)
    real*8     ::  CPUTime1, CPUTime2
        
    irank = GP%irank
    call CPU_TIME(CPUTime1)
        
    if(irank.lt.LP(iLayer)%nsize) then
        
    nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
    nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
    nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
    nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
        
    !/// apply 9-point filter to wave height ///!
    do j = nstarty,nendy
    do i = nstartx,nendx
        !///leave the value in the outermost cells unchanged///!
        if(i.eq.1.or.i.eq.LP(iLayer)%NX.or.j.eq.1.or.j.eq.LP(iLayer)%NY) then
            HS(i,j) = LP(iLayer)%H(2,i,j); cycle
        endif
        h0 = LP(iLayer)%H(2,i,j);     z0 = LP(iLayer)%Z(i,j)   
        h1 = LP(iLayer)%H(2,i-1,j);   z1 = LP(iLayer)%Z(i-1,j)
        h2 = LP(iLayer)%H(2,i+1,j);   z2 = LP(iLayer)%Z(i+1,j) 
        h3 = LP(iLayer)%H(2,i,j-1);   z3 = LP(iLayer)%Z(i,j-1)
        h4 = LP(iLayer)%H(2,i,j+1);   z4 = LP(iLayer)%Z(i,j+1)
        h5 = LP(iLayer)%H(2,i-1,j-1); z5 = LP(iLayer)%Z(i-1,j-1)
        h6 = LP(iLayer)%H(2,i+1,j-1); z6 = LP(iLayer)%Z(i+1,j-1)
        h7 = LP(iLayer)%H(2,i-1,j+1); z7 = LP(iLayer)%Z(i-1,j+1)
        h8 = LP(iLayer)%H(2,i+1,j+1); z8 = LP(iLayer)%Z(i+1,j+1)
        d0 = z0+h0; d1 = z1+h1; d2 = z2+h2; d3 = z3+h3; d4 = z4+h4
        d5 = z5+h5; d6 = z6+h6; d7 = z7+h7; d8 = z8+h8
        if(d0.le.0.0.or.d1.le.0.0.or.d2.le.0.0.or.d3.le.0.0.or.d4.le.0.0.or. &
            d5.le.0.0.or.d6.le.0.0.or.d7.le.0.0.or.d8.le.0.0) then
            HS(i,j) = h0
        else
            HS(i,j) = h0 + (h1+h2+h3+h4-4.0*h0)/8.0 + (h5+h6+h7+h8-4.0*h0)/16.0
        endif
    enddo
    enddo
        
    !///give the filtered value back to H///!
    do j = nstarty,nendy
    do i = nstartx,nendx
        LP(iLayer)%H(2,i,j) = HS(i,j)
        if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit) then
            LP(iLayer)%H(2,i,j) = 0.0d0
        elseif(LP(iLayer)%H(2,i,j)+LP(iLayer)%Z(i,j).le.0.0) then
            if(LP(iLayer)%Z(i,j).gt.0.0) then
                LP(iLayer)%H(2,i,j) = -LP(iLayer)%Z(i,j)
            else
                LP(iLayer)%H(2,i,j) = 0.0d0
            endif
        endif
    enddo
    enddo
        
    endif !if: this node is used by this layer
        
    call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1
        
    end subroutine smoothWaveHeight
        
        
        
    subroutine smoothFlux(GP, LP, iLayer)
        
    use VariableDefination
    implicit NONE
    type(GlobalParameters)    ::  GP
    type(LayerParameters)     ::  LP(100)
    integer*4   ::  iLayer, irank
    integer*4   ::  nstartx, nendx, nstarty, nendy
    integer*4   ::  istart, iend, jstart, jend, i, j
    real*8      ::  m0, m1, m2, m3, m4, m5, m6, m7, m8
    real*8      ::  n0, n1, n2, n3, n4, n5, n6, n7, n8
    real*8      ::  MS(LP(iLayer)%NX-1,LP(iLayer)%NY), NS(LP(iLayer)%NX,LP(iLayer)%NY-1)
    real*8      ::  CPUTime1, CPUTime2
        
    irank = GP%irank
    call CPU_TIME(CPUTime1)
        
    if(irank.lt.LP(iLayer)%nsize) then
           
    nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
    nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
    nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
    nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
        
    !///apply 9-point filter to flux component M///!
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX-1,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY,nendy)
    do j = jstart,jend
    do i = istart,iend
        !///leave the value in the outermost cells unchanged///!
        if(i.eq.1.or.i.eq.LP(iLayer)%NX-1.or.j.eq.1.or.j.eq.LP(iLayer)%NY) then
            MS(i,j) = LP(iLayer)%M(2,i,j); cycle
        endif
        !///leave the value on the wet-dry boundary unchanged///!
        if(LP(iLayer)%Sign_M(i,j).ne.0) then
            MS(i,j) = LP(iLayer)%M(2,i,j); cycle
        endif
        m0 = LP(iLayer)%M(2,i,j);     m1 = LP(iLayer)%M(2,i-1,j)
        m2 = LP(iLayer)%M(2,i+1,j);   m3 = LP(iLayer)%M(2,i,j-1) 
        m4 = LP(iLayer)%M(2,i,j+1);   m5 = LP(iLayer)%M(2,i-1,j-1)  
        m6 = LP(iLayer)%M(2,i+1,j-1); m7 = LP(iLayer)%M(2,i-1,j+1) 
        m8 = LP(iLayer)%M(2,i+1,j+1)
        MS(i,j) = m0 + (m1+m2+m3+m4-4.0*m0)/8.0 + (m5+m6+m7+m8-4.0*m0)/16.0
    enddo
    enddo
    !///give the filtered value back to M///!
    do j = jstart,jend
    do i = istart,iend
        if(LP(iLayer)%Sign_M(i,j).eq.999) then
            LP(iLayer)%M(2,i,j) = 0.0d0; cycle
        endif
        LP(iLayer)%M(2,i,j) = MS(i,j)
        if(LP(iLayer)%Sign_M(i,j).ne.0) then
            if(LP(iLayer)%M(2,i,j)*LP(iLayer)%Sign_M(i,j).le.0.0) LP(iLayer)%M(2,i,j) = 0.0d0
        endif
    enddo
    enddo
        
    !///apply 9-point filter to flux component N///!
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY-1,nendy)
    do j = jstart,jend
    do i = istart,iend
        !///leave the value in the outermost cells unchanged///!
        if(i.eq.1.or.i.eq.LP(iLayer)%NX.or.j.eq.1.or.j.eq.LP(iLayer)%NY-1) then
            NS(i,j) = LP(iLayer)%N(2,i,j); cycle
        endif
        !///leave the value on the wet-dry boundary unchanged///!
        if(LP(iLayer)%Sign_N(i,j).ne.0) then
            NS(i,j) = LP(iLayer)%N(2,i,j); cycle
        endif
        n0 = LP(iLayer)%N(2,i,j);     n1 = LP(iLayer)%N(2,i-1,j)
        n2 = LP(iLayer)%N(2,i+1,j);   n3 = LP(iLayer)%N(2,i,j-1) 
        n4 = LP(iLayer)%N(2,i,j+1);   n5 = LP(iLayer)%N(2,i-1,j-1)  
        n6 = LP(iLayer)%N(2,i+1,j-1); n7 = LP(iLayer)%N(2,i-1,j+1) 
        n8 = LP(iLayer)%N(2,i+1,j+1)
        NS(i,j) = n0 + (n1+n2+n3+n4-4.0*n0)/8.0 + (n5+n6+n7+n8-4.0*n0)/16.0
    enddo
    enddo
    !///give the filtered value back to N///!
    do j = jstart,jend
    do i = istart,iend
        if(LP(iLayer)%Sign_N(i,j).eq.999) then
            LP(iLayer)%N(2,i,j) = 0.0d0; cycle
        endif
        LP(iLayer)%N(2,i,j) = NS(i,j)
        if(LP(iLayer)%Sign_N(i,j).ne.0) then
            if(LP(iLayer)%N(2,i,j)*LP(iLayer)%Sign_N(i,j).le.0.0) LP(iLayer)%N(2,i,j) = 0.0d0
        endif
    enddo
    enddo
        
    endif !if: this node is used by this layer
        
    call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1
        
    end subroutine smoothFlux

