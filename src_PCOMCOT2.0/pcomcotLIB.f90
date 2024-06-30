subroutine readConfig(GP)

use mpi 
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
integer*4               ::  ios, i, indx1, indx2, errorcode, ierror
real*8                  ::  version
character(999)          ::  s

open(666,file=TRIM(ADJUSTL(GP%COMCOTParametersFileName)),status='old',form='formatted',action='read',iostat=ios)
if(ios.ne.0) then
    write(*,*)
    write(*,*) 'ERROR: cann''t find file ',TRIM(ADJUSTL(GP%COMCOTParametersFileName))
    write(*,*)
    call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
endif

do i = 1,3
    read(666, '(a)') s
enddo
indx1 = INDEX(s, '(v') + 2
indx2 = INDEX(s, ')') - 1
read(s(indx1:indx2),*) version
if(ABS(GP%Version-version).gt.0.05) then
    write(*,*)
    write(*,'(a,f4.1)') 'ERROR: pcomcot.ctl version should be ', GP%Version
    write(*,*)
    call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
endif

write(*,*)
write(*,*) '##################################################'
write(*,*) '#                                                #'
write(*,'(a,f4.1,a)') ' #         PCOMCOT Simulation Version',GP%Version,'         #'
write(*,*) '#                                                #'
write(*,*) '##################################################'
read(666, '(4/)')
read(666, '(49x,i30)')      GP%PurposeCalculation
read(666, '(49x,i30)')      GP%InitialConditionType
read(666, '(49x,i30)')      GP%CoordinatesType
read(666, '(49x,f30.6)')    GP%TotalTime
read(666, '(49x,f30.6)')    GP%dt
read(666, '(49x,f30.6)')    GP%DTSaveData
read(666, '(49x,i30)')      GP%SaveFlux
read(666, '(49x,i30)')      GP%SaveDynamicPressure
read(666, '(49x,i30)')      GP%MinGridsPerNode
read(666, '(49x,i30)')      GP%FeedbackToParentLayer
read(666, '(3/)')
read(666, '(49x,i30)')      GP%HorizontalMotion
read(666, '(49x,i30)')      GP%KajiuraFilter
read(666, '(49x,i30)')      GP%UseAverageDepth
read(666, '(49x,f30.6)')    GP%MinKajiuraDepth
read(666, '(3/)')
read(666, '(49x,i30)')      GP%Nonlinearity
read(666, '(49x,i30)')      GP%Dispersion
read(666, '(49x,i30)')      GP%DepthVariability
read(666, '(49x,f30.6)')    GP%MinDispersionDepth
read(666, '(49x,i30)')      GP%Breaking
read(666, '(49x,i30)')      GP%FluxCenter
read(666, '(49x,f30.6)')    GP%FrCap
read(666, '(49x,f30.6)')    GP%DTFilterData
read(666, '(3/)')
read(666, '(49x,i30)')      GP%BoundaryConditionType
read(666, '(49x,f30.6)')    GP%SpongeWidthX
read(666, '(49x,f30.6)')    GP%SpongeWidthY
read(666, '(49x,f30.6)')    GP%MaxSpongeMANNING
read(666, '(49x,f30.6)')    GP%SpongeDampingA
read(666, '(49x,f30.6)')    GP%SpongeDampingR
read(666, '(3/)')
read(666, '(49x,f30.6)')    GP%PermanentDryLimit
read(666, '(49x,f30.6)')    GP%MinWaterDepth
read(666, '(49x,f30.6)')    GP%FrictionDepthLimit
read(666, '(49x,f30.6)')    GP%MANNING 
read(666, '(3/)')
read(666, '(49x,f30.6)')    GP%SourceStartX
read(666, '(49x,f30.6)')    GP%SourceEndX
read(666, '(49x,f30.6)')    GP%SourceDX
read(666, '(49x,f30.6)')    GP%SourceStartY
read(666, '(49x,f30.6)')    GP%SourceEndY
read(666, '(49x,f30.6)')    GP%SourceDY
read(666, '(49x,i30)')      GP%SourceBasisFunctionType
read(666, '(49x,f30.6)')    GP%GaussianRatio
read(666, '(49x,f30.6)')    GP%SigmoidCoefficient
close(666)

GP%MinKajiuraDepth = MAX(GP%MinKajiuraDepth,5e-2)
GP%MinDispersionDepth = MAX(GP%MinDispersionDepth, 5e-2)
GP%FrCap = MAX(GP%FrCap,0.0d0)
GP%MaxSpongeMANNING = MAX(GP%MaxSpongeMANNING,0.0d0)
GP%SpongeDampingA = MAX(GP%SpongeDampingA, 1.0d0)
GP%PermanentDryLimit = MAX(GP%PermanentDryLimit, 0.0d0)
GP%MinWaterDepth = MAX(GP%MinWaterDepth, 5e-3)
GP%FrictionDepthLimit = MAX(GP%FrictionDepthLimit, 5e-3)
GP%MANNING = MAX(GP%MANNING, 0.0d0)
GP%NDTFilterData = MAX(1, NINT(GP%DTFilterData/GP%dt))

if(GP%BoundaryConditionType.ne.1.and.GP%BoundaryConditionType.ne.2) then
    write(*,*) "ERROR: only wall and sponge boundary are supported."
    call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
endif

end subroutine ReadConfig



subroutine checkFiles(GP)

use mpi      
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
integer*4               ::  iLayer, ios, errorcode, ierror
character(999)          ::  fn1, fn2, s
logical                 ::  fe, fe2

write(*,*)
write(*,*) 'checking files that are needed ...'
write(*,*)

write(*,'(a)',advance='no') '         Bathymetry Data File:  '
GP%NumLayers = 0
do iLayer = 1,99
    write(fn1,'(a,i2.2,a)') TRIM(ADJUSTL(GP%BathymetryFilePrefix)),iLayer,'.nf'
    inquire(file=fn1, exist=fe)
    write(fn2,'(a,i2.2,a)') TRIM(ADJUSTL(GP%BathymetryFilePrefix)),iLayer,'.xyz'
    inquire(file=fn2, exist=fe2)
    if(fe.or.fe2) then
    GP%NumLayers = GP%NumLayers+1
    if(GP%NumLayers.ne.1.and.MOD(GP%NumLayers-1,5).eq.0) then
        write(*,*)
        write(*,'(a)',advance='no') '                                '
    elseif(GP%NumLayers.ne.1) then
        write(*,'(a)',advance='no') ', '
    endif
    if(fe) then
        write(*,'(a11)',advance='no') ADJUSTL(fn1)
    else
        write(*,'(a11)',advance='no') ADJUSTL(fn2)
    endif
    endif
enddo
write(*,*)
if(GP%NumLayers.eq.0) then
    write(*,*) 'ERROR: cann''t find any bathymetry data files:  layerXX.nf or layerXX.xyz'
    call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
endif

if(GP%PurposeCalculation.eq.1.and.GP%InitialConditionType.eq.0) then
    fn1 = TRIM(ADJUSTL(GP%InitialElevationFilePrefix))//'.nf'
    inquire(file=fn1, exist=fe)
    fn2 = TRIM(ADJUSTL(GP%InitialElevationFilePrefix))//'.xyz'
    inquire(file=fn2, exist=fe2)
    if((.not.fe).and.(.not.fe2)) then
        write(*,*) 'ERROR: cann''t find initial water elevation file:  ',TRIM(ADJUSTL(fn1)),' or ',TRIM(ADJUSTL(fn2))
        call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
    endif
    if(fe) then
        GP%InitialElevationFileName = fn1
        GP%InitialElevationFileFormat = 1
    elseif(fe2) then
        GP%InitialElevationFileName = fn2
        GP%InitialElevationFileFormat = 2
    endif
    write(*,*) '        Initial Water Elevation Data File:  ',TRIM(ADJUSTL(GP%InitialElevationFileName))
    fn1 = TRIM(ADJUSTL(GP%InitialMFileName)) 
    inquire(file=fn1,exist=fe)
    fn2 = TRIM(ADJUSTL(GP%InitialNFileName))
    inquire(file=fn2,exist=fe2)
    if(fe.or.fe2) then
        write(*,'(a)',advance='no') '         Initial Flux Data File:  '
        if(fe) then
            write(*,'(a,a)',advance='no') TRIM(ADJUSTL(GP%InitialMFileName)), '  '
        endif
        if(fe2) then
            write(*,'(a)',advance='no') TRIM(ADJUSTL(GP%InitialNFileName))
        endif
        write(*,*)
    endif
endif
      
if((GP%PurposeCalculation.eq.1.and.GP%InitialConditionType.eq.1).or.(GP%PurposeCalculation.eq.2)) then
    inquire(file=TRIM(ADJUSTL(GP%FaultParametersFileName)),exist=fe)
    if(.not.fe) then
        write(*,*) 'ERROR: cann''t find muli-fault information file:  ', TRIM(ADJUSTL(GP%FaultParametersFileName))
        call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
    endif
    write(*,*) '        Muli-Fault Information File:  ',TRIM(ADJUSTL(GP%FaultParametersFileName))
endif

inquire(file=TRIM(ADJUSTL(GP%StationFileName)),exist=fe)
if(.not.fe) then
    write(*,*) 'ERROR: cann''t find stations information file:  ',TRIM(ADJUSTL(GP%StationFileName))
    call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
endif
write(*,*) '        Stations coordinates file:  ',TRIM(ADJUSTL(GP%StationFileName))

inquire(file=TRIM(ADJUSTL(GP%ResultPath)),exist=fe)
if(.not.fe) then
    ierror = system('mkdir '//TRIM(ADJUSTL(GP%ResultPath)))
endif


end subroutine checkFiles



subroutine getBathymetryDataSize(GP, LP)

use mpi      
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4               ::  iLayer, ierror, iLayerCheck
character(999)          ::  fn1, fn2
logical                 ::  fe, fe2
real*8     ::  x1,x2,x,y1,y
integer*4  ::  i, ios, endoffile


write(*,*)
write(*,*) 'getting bathymetry data size ...'
write(*,*)

iLayer = 0
do iLayerCheck = 1,99
    write(fn1,'(a,i2.2,a)') TRIM(ADJUSTL(GP%BathymetryFilePrefix)),iLayerCheck,'.nf'
    inquire(file=fn1, exist=fe)
    write(fn2,'(a,i2.2,a)') TRIM(ADJUSTL(GP%BathymetryFilePrefix)),iLayerCheck,'.xyz'
    inquire(file=fn2, exist=fe2)
    if(fe.or.fe2) then
    iLayer = iLayer+1
    if(fe) then
        LP(iLayer)%BathymetryFileName = fn1
        LP(iLayer)%BathymetryFileFormat = 1
    elseif(fe2) then
        LP(iLayer)%BathymetryFileName = fn2
        LP(iLayer)%BathymetryFileFormat = 2
    endif
    if(LP(iLayer)%BathymetryFileFormat.eq.1) then
        call getBathymetryDataSizeNetCDF(GP, LP, iLayer)
    else
        open(23,file=TRIM(ADJUSTL(LP(iLayer)%BathymetryFileName)),status='old',form='formatted',action='read',iostat=ios)
        read(23,*) x1,y1
        if(iLayer.eq.1) then
            GP%StartEastWest = 0
            if((GP%CoordinatesType.eq.0).and.(x1.ge.0).and.(x1.le.180))then
                GP%StartEastWest = 1
            elseif((GP%CoordinatesType.eq.0).and.(x1.gt.180.or.x1.lt.0.0))then
                GP%StartEastWest = 2
            endif
        endif
        if((GP%StartEastWest.eq.1).and.(x1.lt.0)) x1 = x1+360
        if((GP%StartEastWest.eq.2).and.(x1.gt.180)) x1 = x1-360
        LP(iLayer)%xmin = x1; LP(iLayer)%ymin = y1
        read(23,*) x2
        if((GP%StartEastWest.eq.1).and.(x2.lt.0)) x2 = x2+360
        if((GP%StartEastWest.eq.2).and.(x2.gt.180)) x2 = x2-360
        LP(iLayer)%dx = x2-x1
        LP(iLayer)%NX = 2; LP(iLayer)%NY = 1
        do
            read(23,*,iostat=ios) x,y
            if((GP%StartEastWest.eq.1).and.(x.lt.0)) x = x+360
            if((GP%StartEastWest.eq.2).and.(x.gt.180)) x = x-360
            if(ABS(x-x1).gt.0.1*LP(iLayer)%dx) then
                LP(iLayer)%NX = LP(iLayer)%NX+1;  LP(iLayer)%xmax = x
            else
                LP(iLayer)%NY = LP(iLayer)%NY+1;  LP(iLayer)%dy = y-y1
                exit
            endif
        enddo
        do i = 1,LP(iLayer)%NX-1
            read(23,*) x, LP(iLayer)%ymax
        enddo
        do
            endoffile = 0
            do i = 1,LP(iLayer)%NX
                read(23,*,iostat=ios) x, LP(iLayer)%ymax
                if(ios.ne.0) then
                    endoffile = 1
                    exit
                endif
            enddo
            if(endoffile.eq.1) exit
            LP(iLayer)%NY = LP(iLayer)%NY + 1
        enddo
        close(23)
    endif

    write(*,'(a,a11)',advance='no') '        ',ADJUSTL(LP(iLayer)%BathymetryFileName)
    write(*,'(a,f12.3,a,f12.3,a,f15.6,a,i5)')                        &
        ':  xmin: ',LP(iLayer)%xmin,'    xmax: ',LP(iLayer)%xmax, &
        '    dx:',LP(iLayer)%dx,'    nx: ',LP(iLayer)%NX
    write(*,'(a,a,f12.3,a,f12.3,a,f15.6,a,i5)')                      &
        '                   ',                                    &
        '   ymin: ',LP(iLayer)%ymin,'    ymax: ',LP(iLayer)%ymax, &
        '    dy:',LP(iLayer)%dy,'    ny: ',LP(iLayer)%NY
    if(LP(iLayer)%BathymetryFileFormat.eq.2)                      &
    write(*,'(a)') '        NOTE: xyz format, bathymetry data must vary with x first.'
    write(*,*)
    endif
enddo

end subroutine getBathymetryDataSize



subroutine determineLayerDependency(GP, LP)

use mpi
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4               ::  errorcode, ierror
integer*4               ::  i, j, k, NumMotherLayers, NumTopLayers, IsOverlap
real*8                  ::  p1(4,2), p2(4,2)

GP%NumLayerLevels = 0
do i = 1,GP%NumLayers
    NumMotherLayers = 0
    do j = 1,GP%NumLayers
        if(j.ne.i) then
            if((LP(i)%xmin.ge.LP(j)%xmin).and.(LP(i)%xmax.le.LP(j)%xmax).and. &
                (LP(i)%ymin.ge.LP(j)%ymin).and.(LP(i)%ymax.le.LP(j)%ymax)) &
                NumMotherLayers = NumMotherLayers+1
        endif
    enddo
    LP(i)%Level = NumMotherLayers+1
    GP%NumLayerLevels = MAX(GP%NumLayerLevels, LP(i)%Level)
enddo
      
GP%TopLayer = 0; NumTopLayers = 0
do i = 1,GP%NumLayers
    if(LP(i)%Level.eq.1) then
        GP%TopLayer = i
        NumTopLayers = NumTopLayers+1
    endif
enddo
if(NumTopLayers.gt.1) then
    write(*,'(a)',advance='no') 'ERROR: more than 1 top layer found: '
    do i = 1,GP%NumLayers
        if(LP(i)%Level.eq.1) write(*,'(a,a)',advance='no') TRIM(ADJUSTL(LP(i)%BathymetryFileName)),'  '
    enddo
    write(*,*)
    call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
endif

do i = 1,GP%NumLayers
    do j = 1,GP%NumLayers
        if(i.ne.j.and.LP(i)%Level.eq.LP(j)%Level) then
            IsOverlap = 0
            p1(1,1)=LP(i)%xmin; p1(1,2)=LP(i)%ymin; p1(2,1)=LP(i)%xmax; p1(2,2)=LP(i)%ymin
            p1(3,1)=LP(i)%xmax; p1(3,2)=LP(i)%ymax; p1(4,1)=LP(i)%xmin; p1(4,2)=LP(i)%ymax
            p2(1,1)=LP(j)%xmin; p2(1,2)=LP(j)%ymin; p2(2,1)=LP(j)%xmax; p2(2,2)=LP(j)%ymin
            p2(3,1)=LP(j)%xmax; p2(3,2)=LP(j)%ymax; p2(4,1)=LP(j)%xmin; p2(4,2)=LP(j)%ymax
            do k = 1,4
                if(p1(k,1).ge.p2(1,1).and.p1(k,1).le.p2(3,1).and. &
                    p1(k,2).ge.p2(1,2).and.p1(k,2).le.p2(3,2)) then
                    IsOverlap = 1
                    exit
                endif
            enddo
            do k = 1,4
                if(p2(k,1).ge.p1(1,1).and.p2(k,1).le.p1(3,1).and. &
                    p2(k,2).ge.p1(1,2).and.p2(k,2).le.p1(3,2)) then
                    IsOverlap = 1
                    exit
                endif
            enddo
            if(IsOverlap.eq.1) then
                write(*,*) "Warning: layers overlaped: ",               &
                    TRIM(ADJUSTL(LP(i)%BathymetryFileName)), " and ", &
                    TRIM(ADJUSTL(LP(j)%BathymetryFileName))
            endif
        endif
    enddo
enddo

do i = 1,GP%NumLayers
    LP(i)%Parent = 0
    do j = 1,GP%NumLayers
        if(LP(i)%Level-LP(j)%Level.eq.1) then
            if((LP(i)%xmin.ge.LP(j)%xmin).and.(LP(i)%xmax.le.LP(j)%xmax).and. &
                (LP(i)%ymin.ge.LP(j)%ymin).and.(LP(i)%ymax.le.LP(j)%ymax)) then
                LP(i)%Parent = j; exit
            endif
        endif
    enddo
enddo

do i = 1,GP%NumLayers
    write(*,'(a,a11)',advance='no') '        ',ADJUSTL(LP(i)%BathymetryFileName)
    write(*,'(a,i3,a)',advance='no') ':  Level: ',LP(i)%Level,'    Parent: '
    if(LP(i)%Parent.ne.0) then
        write(*,'(a11)') ADJUSTL(LP(LP(i)%Parent)%BathymetryFileName)
    else
        write(*,*)
    endif
    enddo
write(*,*)
do i = 1,GP%NumLayers
if(LP(i)%Level.gt.1) then
    if(LP(i)%dx.gt.LP(LP(i)%Parent)%dx) &
        write(*,'(a,a,a,f7.4,a,a,a,f7.4)')                           &
        'WARNING: ',TRIM(ADJUSTL(LP(i)%BathymetryFileName)),         &
        ' dx: ',LP(i)%dx,' larger than Parent ',                     &
        TRIM(ADJUSTL(LP(LP(i)%Parent)%BathymetryFileName)),' dx: ', LP(LP(i)%Parent)%dx
    if(LP(i)%dy.gt.LP(LP(i)%Parent)%dy) &
        write(*,'(a,a,a,f7.4,a,a,a,f7.4)')                           &
        'WARNING: ',TRIM(ADJUSTL(LP(i)%BathymetryFileName)),         &
        ' dy: ',LP(i)%dy,' larger than Parent ',                     &
        TRIM(ADJUSTL(LP(LP(i)%Parent)%BathymetryFileName)),' dy: ', LP(LP(i)%Parent)%dy
endif
enddo

end subroutine determineLayerDependency



subroutine readLayerConfig(GP, LP)

use mpi
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer, pLayer, i, ios, errorcode, ierror, hasWrong
integer*4   ::  layerID, indx1, indx2
integer*4   ::  nonlinear, dispersive, variableDepth, breakingWave, fluxCentered, ndt
real*8      ::  version, T0, dt
character(999)  ::  ctlfile, s, fn1, fn2

write(*,*)
write(*,*) 'setting simulation parameters of each layer...'
write(*,*)
!///Parameters of each layer are the same as global parameters by default ///!
do iLayer = 1,GP%NumLayers
    LP(iLayer)%Nonlinearity = GP%Nonlinearity
    LP(iLayer)%Dispersion = GP%Dispersion
    LP(iLayer)%DepthVariability = GP%DepthVariability
    LP(iLayer)%Breaking = GP%Breaking
    LP(iLayer)%FluxCenter = GP%FluxCenter
    LP(iLayer)%ComputingStartTime = 0.0d0
    LP(iLayer)%DTFilterData = GP%DTFilterData
    LP(iLayer)%NDTFilterData = GP%NDTFilterData
enddo

ctlfile = GP%COMCOTLayerCtlFileName
open(999,file=TRIM(ADJUSTL(ctlfile)),status='old',action='read',form='formatted',iostat=ios)
if(ios.eq.0) then
do i = 1,3
    read(999,'(a)') s
enddo
indx1 = INDEX(s,'(v') + 2
indx2 = INDEX(s,')') -1
read(s(indx1:indx2),*) version
if(ABS(GP%Version-version).gt.0.05) then
    write(*,*)
    write(*,'(a,f4.1)') 'ERROR: layers.ctl version should be ', GP%Version
    write(*,*)
    close(999)
    call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
endif
read(999,*)
do 
    read(999,'(a)',iostat=ios) s
    if(ios.ne.0) exit
    read(999,'(a)',iostat=ios) s
    if(ios.ne.0) exit
    read(999,'(a)',iostat=ios) s
    if(ios.ne.0.or.INDEX(s,'Layer Name').eq.0) exit
    read(s,'(49x,i30)') layerID
    read(999,'(a)')
    read(999,'(49x,i30)')   nonlinear
    read(999,'(49x,i30)')   dispersive
    read(999,'(49x,i30)')   variableDepth
    read(999,'(49x,i30)')   breakingWave
    read(999,'(49x,i30)')   fluxCentered
    read(999,'(49x,f30.6)') T0
    read(999,'(49x,f30.6)') dt
    write(fn1,'(a,i2.2,a)') TRIM(ADJUSTL(GP%BathymetryFilePrefix)),layerID,'.nf'
    write(fn2,'(a,i2.2,a)') TRIM(ADJUSTL(GP%BathymetryFilePrefix)),layerID,'.xyz'
    do iLayer = 1,GP%NumLayers
    if(TRIM(ADJUSTL(LP(iLayer)%BathymetryFileName)).eq.TRIM(ADJUSTL(fn1)).or.&
        TRIM(ADJUSTL(LP(iLayer)%BathymetryFileName)).eq.TRIM(ADJUSTL(fn2))) then
        LP(iLayer)%Nonlinearity = nonlinear
        LP(iLayer)%Dispersion = dispersive
        LP(iLayer)%DepthVariability = variableDepth
        LP(iLayer)%Breaking = breakingWave
        LP(iLayer)%FluxCenter = fluxCentered
        LP(iLayer)%ComputingStartTime = MAX(0.0d0, T0)
        LP(iLayer)%DTFilterData = dt
        LP(iLayer)%NDTFilterData = MAX(1, NINT(LP(iLayer)%DTFilterData/GP%dt))
        exit
    endif
    enddo
enddo
close(999)
endif !if:layers.ctl exists 
    
hasWrong = 0
do iLayer = 1,GP%NumLayers
    pLayer = LP(iLayer)%Parent
    write(*,'(a,a)',advance='no') '        ', TRIM(ADJUSTL(LP(iLayer)%BathymetryFileName))
    write(*,'(a,i2,a)',advance='no') ':   Level: ', LP(iLayer)%Level, '    Parent: '
    if(pLayer.ge.1) then
        write(*,'(a11)') TRIM(ADJUSTL(LP(pLayer)%BathymetryFileName))
    else
        write(*,*)
    endif
    write(*,'(a)',advance='no') '        Governing Equation: '
    if(LP(iLayer)%Nonlinearity.eq.0) then
        write(*,'(a)',advance='no') 'linear, '
    elseif(LP(iLayer)%Nonlinearity.eq.1) then
        write(*,'(a)',advance='no') 'nonlinear, '
    endif
    if(LP(iLayer)%Dispersion.eq.0) then
        write(*,'(a)',advance='no') 'non-dispersive, '
    elseif(LP(iLayer)%Dispersion.eq.1) then
        write(*,'(a)',advance='no') 'dispersive, '
    endif
    if(LP(iLayer)%Breaking.eq.0) then
        write(*,'(a)') 'non-breaking waves'
    elseif(LP(iLayer)%Breaking.eq.1) then
        write(*,'(a)') 'breaking waves'
    endif
    write(*,'(a)',advance='no') '        Scheme for LSWEs: '
    if(LP(iLayer)%FluxCenter.eq.0) then
        write(*,'(a)',advance='no') 'forward time-centered space (no diffusion)'
    elseif(LP(iLayer)%FluxCenter.eq.1) then
        write(*,'(a)',advance='no') 'flux-centered (extra diffusion)'
    endif

    write(s,'(f9.1)') LP(iLayer)%ComputingStartTime
    write(*,'(a,a)') '    Starting Time(seconds): ', TRIM(ADJUSTL(s))
    if(pLayer.ge.1.and.LP(iLayer)%ComputingStartTime.lt.LP(pLayer)%ComputingStartTime) then
        write(*,'(a)') 'ERROR: '//TRIM(ADJUSTL(LP(iLayer)%BathymetryFileName))// 'starts earliear than its parent layer!'
        hasWrong = 1
    endif
    write(*,*)
enddo

if(hasWrong.eq.1) call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)


end subroutine readLayerConfig



subroutine readBathymetry(GP, LP)

use mpi    
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4  ::  iLayer, i, j, ios, iLayerLevel, iFlag
real*8     ::  x, y, z, xdistance, ydistance
character(999)  ::  s

write(*,*)
write(*,*) 'reading bathymetry data ...'
write(*,*)
write(*,'(a)',advance='no') '        '

do iLayer = 1,GP%NumLayers
    if(iLayer.ne.1.and.MOD(iLayer-1,5).eq.0) then
        write(*,*)
        write(*,'(a)',advance='no') '        '
    elseif(iLayer.ne.1) then
        write(*,'(a)',advance='no') ', '
    endif
    write(*,'(a11)',advance='no') ADJUSTL(LP(iLayer)%BathymetryFileName)
    LP(iLayer)%zmax = -2.0e10; LP(iLayer)%zmin = 2.0e10
    if(LP(iLayer)%BathymetryFileFormat.eq.1) then
        call readBathymetryNetCDF(GP, LP, iLayer)
    else
        open(23,file=TRIM(ADJUSTL(LP(iLayer)%BathymetryFileName)), &
            status='old',form='formatted',action='read',iostat=ios)
        do j = 1,LP(iLayer)%NY
            do i = 1,LP(iLayer)%NX
                read(23,*) x,y,LP(iLayer)%Z(i,j)
                LP(iLayer)%zmax = MAX(LP(iLayer)%zmax, LP(iLayer)%Z(i,j))
                LP(iLayer)%zmin = MIN(LP(iLayer)%zmin, LP(iLayer)%Z(i,j))
                if((GP%StartEastWest.eq.1).and.(x.lt.0)) x = x+360
                if((GP%StartEastWest.eq.2).and.(x.gt.180)) x = x-360
                if(j.eq.1) LP(iLayer)%X(i) = x
                if(i.eq.1) LP(iLayer)%Y(j) = y
            enddo
        enddo
        close(23)
    endif
enddo

do iLayerLevel = 1,GP%NumLayerLevels
do iLayer = 1,GP%NumLayers
if(LP(iLayer)%Level.eq.iLayerLevel) then

    if(iLayerLevel.ge.2) then
        iFlag = 0
        do j = 1,LP(iLayer)%NY
        do i = 1,LP(iLayer)%NX
            if(LP(iLayer)%Z(i,j).ne.LP(iLayer)%Z(i,j)) then
                call interpData(GP, LP, LP(iLayer)%Parent, iFlag, LP(iLayer)%X(i), LP(iLayer)%Y(j), z)
                LP(iLayer)%Z(i,j) = z
            endif
        enddo
        enddo
    endif

    do j = 1,LP(iLayer)%NY
        do i = 1,LP(iLayer)%NX
            if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit) then
                call removeNarrowWater(GP,LP,iLayer,i,j)
            endif
        enddo
    enddo
            
    write(s,'(a,i2.2,a)') TRIM(ADJUSTL(GP%ResultPath))//'_xcoordinate',iLayer,'.dat'
    open(23,file=TRIM(ADJUSTL(s)),form='formatted',status='replace')
    do i = 1,LP(iLayer)%NX
        write(23,'(f15.5)',advance='no') LP(iLayer)%X(i)
    enddo
    close(23)
    write(s,'(a,i2.2,a)') TRIM(ADJUSTL(GP%ResultPath))//'_ycoordinate',iLayer,'.dat'
    open(23,file=TRIM(ADJUSTL(s)),form='formatted',status='replace')
    do i = 1,LP(iLayer)%NY
        write(23,'(f15.5)',advance='no') LP(iLayer)%Y(i)
    enddo
    close(23)
    write(s,'(a,i2.2,a)') TRIM(ADJUSTL(GP%ResultPath))//'_bathymetry',iLayer,'.dat'
    open(23,file=TRIM(ADJUSTL(s)),form='unformatted',status='replace')
    write(23) LP(iLayer)%NX, LP(iLayer)%NY
    do j = 1,LP(iLayer)%NY
        write(23) (LP(iLayer)%Z(i,j), i=1,LP(iLayer)%NX)
    enddo
    close(23)
      
endif
enddo
enddo
write(*,*)

end subroutine readBathymetry



recursive subroutine removeNarrowWater(GP, LP, iLayer, m, n)
      
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4  ::  iLayer, m, n, i, j

if(m.le.GP%nRowBathymetryBoundary) then
    do i=m-1, 1, -1
        if(LP(iLayer)%Z(i,n).gt.-GP%PermanentDryLimit) then
            LP(iLayer)%Z(i,n) = -GP%PermanentDryLimit
            call removeNarrowWater(GP,LP,iLayer,i,n)
        endif
    enddo
endif
if(m.ge.LP(iLayer)%NX-GP%nRowBathymetryBoundary+1) then
    do i=m+1, LP(iLayer)%NX
        if(LP(iLayer)%Z(i,n).gt.-GP%PermanentDryLimit) then
            LP(iLayer)%Z(i,n) = -GP%PermanentDryLimit
            call removeNarrowWater(GP,LP,iLayer,i,n)
        endif
    enddo
endif
if(n.le.GP%nRowBathymetryBoundary) then
    do i=n-1, 1, -1
        if(LP(iLayer)%Z(m,i).gt.-GP%PermanentDryLimit) then
            LP(iLayer)%Z(m,i) = -GP%PermanentDryLimit
            call removeNarrowWater(GP,LP,iLayer,m,i)
        endif
    enddo
endif
if(n.ge.LP(iLayer)%NY-GP%nRowBathymetryBoundary+1) then
    do i=n+1, LP(iLayer)%NY
        if(LP(iLayer)%Z(m,i).gt.-GP%PermanentDryLimit) then
            LP(iLayer)%Z(m,i) = -GP%PermanentDryLimit
            call removeNarrowWater(GP,LP,iLayer,m,i)
        endif
    enddo
endif
if(m.ge.GP%nRowBathymetryBoundary+1) then
    do i=m-1, m-GP%nRowBathymetry, -1
        if(LP(iLayer)%Z(i,n).le.-GP%PermanentDryLimit) then
            do j=m-1, i+1, -1
                if(LP(iLayer)%Z(j,n).gt.-GP%PermanentDryLimit) then
                    LP(iLayer)%Z(j,n) = -GP%PermanentDryLimit
                    call removeNarrowWater(GP,LP,iLayer,j,n)
                endif
            enddo
        endif
    enddo
endif
if(m.le.LP(iLayer)%NX-GP%nRowBathymetryBoundary) then
    do i=m+1, m+GP%nRowBathymetry
        if(LP(iLayer)%Z(i,n).le.-GP%PermanentDryLimit) then
            do j=m+1, i-1
                if(LP(iLayer)%Z(j,n).gt.-GP%PermanentDryLimit) then
                    LP(iLayer)%Z(j,n) = -GP%PermanentDryLimit
                    call removeNarrowWater(GP,LP,iLayer,j,n)
                endif
            enddo
         endif
    enddo
endif
if(n.ge.GP%nRowBathymetryBoundary+1) then
    do i=n-1, n-GP%nRowBathymetry, -1
        if(LP(iLayer)%Z(m,i).le.-GP%PermanentDryLimit) then
            do j=n-1, i+1, -1
                if(LP(iLayer)%Z(m,j).gt.-GP%PermanentDryLimit) then
                    LP(iLayer)%Z(m,j) = -GP%PermanentDryLimit
                    call removeNarrowWater(GP,LP,iLayer,m,j)
                endif
            enddo
        endif
    enddo
endif
if(n.le.LP(iLayer)%NY-GP%nRowBathymetryBoundary) then
    do i=n+1, n+GP%nRowBathymetry
        if(LP(iLayer)%Z(m,i).le.-GP%PermanentDryLimit) then
            do j=n+1, i-1
                if(LP(iLayer)%Z(m,j).gt.-GP%PermanentDryLimit) then
                    LP(iLayer)%Z(m,j) = -GP%PermanentDryLimit
                    call removeNarrowWater(GP,LP,iLayer,m,j)
                endif
            enddo
        endif
    enddo
endif

end subroutine removeNarrowWater



subroutine cflCheck(GP, LP)
      
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4               ::  iLayer, i, j
real*8                  ::  zmax, dmin, cr

write(*,*)
write(*,*) 'check CFL and determine dt for each layer ...'
zmax = LP(GP%TopLayer)%zmax
if(GP%CoordinatesType.eq.0) then !spehrical coordinate
    dmin = MIN(COS(LP(GP%TopLayer)%ymin*GP%PI/180.0),COS(LP(GP%TopLayer)%ymax*GP%PI/180.0))*GP%R_EARTH
    dmin = MIN(GP%R_EARTH*LP(GP%TopLayer)%dy, dmin*LP(GP%TopLayer)%dx)*GP%PI/180.0
else
    dmin = MIN(LP(GP%TopLayer)%DX, LP(GP%TopLayer)%DY)
endif
cr = GP%dt*SQRT(GP%GRAV*zmax)/dmin
write(*,*)
write(*,'(a,f6.3)',advance='no') '        CFL number:',cr
write(*,*) '  for linear shallow water equations ...'
if(cr.gt.0.5) then
    write(*,*) 'WARNING: CFL number is greater than 0.5 ...'
endif
write(*,*) '       determine dt for each layer according to cfl number:'
write(*,*)
do iLayer = 1,GP%NumLayers
    if(iLayer.eq.GP%TopLayer) then
        LP(iLayer)%dt = GP%dt
        LP(iLayer)%nStepsPerTimeStep = 1
    else
        zmax = LP(iLayer)%zmax
        if(GP%CoordinatesType.eq.0) then !spehrical coordinate
            dmin = MIN(COS(LP(iLayer)%ymin*GP%PI/180.0),COS(LP(iLayer)%ymax*GP%PI/180.0))*GP%R_EARTH
            dmin = MIN(GP%R_EARTH*LP(iLayer)%dy, dmin*LP(iLayer)%dx)*GP%PI/180.0
        else
            dmin = MIN(LP(iLayer)%DX, LP(iLayer)%DY)
        endif
        LP(iLayer)%dt = cr*dmin/SQRT(GP%GRAV*zmax)
        i = 1
        do
            if(LP(iLayer)%dt*i.ge.GP%dt.or.ABS(LP(iLayer)%dt*i-GP%dt).lt.GP%dt*1e-4) then
                LP(iLayer)%nStepsPerTimeStep = i
                LP(iLayer)%dt = GP%dt/i
                exit
            endif
            i = i + 1
        enddo
    endif
    write(*,'(a,a11,a,f7.4)') '        ',ADJUSTL(LP(iLayer)%BathymetryFileName),':  dt = ',LP(iLayer)%dt
enddo
write(*,*)

end subroutine cflCheck



subroutine readStations(GP, LP, SP)
      
use VariableDefination
implicit NONE
type(GlobalParameters)    ::  GP
type(LayerParameters)     ::  LP(100)
type(StationParameters)   ::  SP(999)
character(999)            ::  s
integer*4  ::  i, i0, j, ios
real*8     ::  x, y

write(*,*) 'reading location of stations from file ', &
    TRIM(ADJUSTL(GP%StationFileName)),' ...'
open(23,file=TRIM(ADJUSTL(GP%StationFileName)), &
    status='old',form='formatted',action='read',iostat=ios)
i = 0; i0 = 0
do
    read(23,'(a)',iostat=ios) s
    if(ios.ne.0) exit
    if(LEN(TRIM(ADJUSTL(s))).eq.0) exit
    i = i+1
    read(s,*) x,y
    if((GP%StartEastWest.eq.1).and.(x.lt.0.0)) x=x+360.0
    if((GP%StartEastWest.eq.2).and.(x.gt.180.0)) x=x-360.0
    if((x.ge.LP(GP%TopLayer)%xmin+LP(GP%TopLayer)%dx).and. &
        (x.le.LP(GP%TopLayer)%xmax-LP(GP%TopLayer)%dx).and. &
        (y.ge.LP(GP%TopLayer)%ymin+LP(GP%TopLayer)%dy).and. &
        (y.le.LP(GP%TopLayer)%ymax-LP(GP%TopLayer)%dy)) then
        i0 = i0 + 1
        SP(i0)%X = x; SP(i0)%Y = y
    else
        if(i0.eq.i-1)  write(*,'(a)',advance='no') &
            '        Station(s) not in computing domain: #'
        write(*,'(i4)',advance='no') i
    endif
enddo
close(23)
if(i0.ne.GP%NumStations) then
    write(*,*)
    GP%NumStations = i0
endif
do i = 1,GP%NumStations
    SP(i)%nLayer = 0
    do j = 1,GP%NumLayers
        if(SP(i)%X.ge.LP(j)%xmin.and.SP(i)%X.le.LP(j)%xmax.and.  &
            SP(i)%Y.ge.LP(j)%ymin.and.SP(i)%Y.le.LP(j)%ymax) then
            if(SP(i)%nLayer.eq.0) then
                SP(i)%nLayer = j
            elseif(LP(j)%Level.gt.LP(SP(i)%nLayer)%Level) then
                SP(i)%nLayer = j
            endif
        endif
    enddo
    write(*,'(a,i3,a,f12.4,a,f12.4,a,i3)')  &
        '        station ',i,':    x:',SP(i)%x,',    y:',SP(i)%y,',    layer:',SP(i)%nLayer
enddo
write(*,*)

end subroutine readStations



subroutine readFaultParameters(GP, FP)
      
use VariableDefination
implicit NONE
type(GlobalParameters)    ::  GP
type(FaultParameters)     ::  FP(4000)
type(FaultParameters)     ::  FPtmp
integer*4                 ::  i, j, k, ios
character(999)            ::  s

write(*,*) 'reading fault parameters from file: ', &
    TRIM(ADJUSTL(GP%FaultParametersFileName))
open(23,file=TRIM(ADJUSTL(GP%FaultParametersFileName)), &
    status='old',form='formatted',action='read',iostat=ios)
read(23, '(3/)')
i = 1
do
    read(23, '(a)', iostat=ios) s
    if(ios.ne.0) exit
    read(23, '(a)', iostat=ios) s
    if(ios.ne.0) exit
    read(23, '(a)', iostat=ios) s
    if(ios.ne.0) exit
    if(INDEX(s,'Parameters for Fault Segment').eq.0) exit
    read(23, '(a)', iostat=ios) s
    read(23, '(49x,f30.6)') FP(i)%T0
    if((FP(i)%T0).lt.0.0.or.FP(i)%T0.gt.GP%TotalTime) then
        FP(i)%NT = -1
    else
        FP(i)%NT = NINT(FP(i)%T0/GP%dt)
    endif
    read(23, '(49x,f30.6)') FP(i)%Depth
    read(23, '(49x,f30.6)') FP(i)%Length
    read(23, '(49x,f30.6)') FP(i)%Width
    read(23, '(49x,f30.6)') FP(i)%Slip
    read(23, '(49x,f30.6)') FP(i)%Rake
    read(23, '(49x,f30.6)') FP(i)%Strike
    read(23, '(49x,f30.6)') FP(i)%Dip
    read(23, '(49x,f30.6)') FP(i)%Y0
    read(23, '(49x,f30.6)') FP(i)%X0
    if(GP%StartEastWest.eq.1.and.FP(i)%X0.lt.0.0) FP(i)%X0 = FP(i)%X0+360.0
    if(GP%StartEastWest.eq.2.and.FP(i)%X0.gt.180.0) FP(i)%X0 = FP(i)%X0-360.0
    i=i+1
enddo
GP%NumFaults = i-1
close(23)
write(*,*)
write(*,'(a,i4,a)') '        There are ',GP%NumFaults,'  fault segments.'
write(*,*)

! if(GP%PurposeCalculation.eq.1) then
! do i = 1,GP%NumFaults-1
! do j = i+1,GP%NumFaults
!     if(FP(i)%NT.gt.FP(j)%NT) then
!         FPtmp = FP(i); FP(i) = FP(j); FP(j) = FPtmp
!     endif
! enddo
! enddo
! endif

end subroutine readFaultParameters



subroutine computeGlobalParameters(GP, LP, SP)
      
use VariableDefination
implicit NONE
type(GlobalParameters)    ::  GP
type(LayerParameters)     ::  LP(100)
type(StationParameters)   ::  SP(999)

write(*,*) 'computing global paramters ...'
write(*,*)

GP%TotalTimeSteps = MAX(1,NINT(GP%TotalTime/GP%dt))
GP%NDTSaveData = MAX(1,NINT(GP%DTSaveData/GP%dt))
      
GP%nCalculations = 1
if(GP%PurposeCalculation.eq.2) then
    GP%nCalculations = GP%NumFaults
elseif(GP%PurposeCalculation.eq.3) then
    if(GP%CoordinatesType.eq.0) then
        if((GP%StartEastWest.eq.1).and.(GP%SourceStartX.lt.0)) &
            GP%SourceStartX = GP%SourceStartX+360
        if((GP%StartEastWest.eq.2).and.(GP%SourceStartX.gt.180)) &
            GP%SourceStartX = GP%SourceStartX-360
        if((GP%StartEastWest.eq.1).and.(GP%SourceEndX.lt.0)) &
            GP%SourceEndX = GP%SourceEndX+360
        if((GP%StartEastWest.eq.2).and.(GP%SourceEndX.gt.180)) &
            GP%SourceEndX = GP%SourceEndX-360
        GP%SourceDX = GP%SourceDX/(GP%R_EARTH*COS((GP%SourceStartY+GP%SourceEndY)/2.0*GP%PI/180.0))*180.0/GP%PI
        GP%SourceDY = GP%SourceDY/GP%R_EARTH*180.0/GP%PI
    endif
    GP%SourceNX = NINT((GP%SourceEndX-GP%SourceStartX)/GP%SourceDX)
    GP%SourceNY = NINT((GP%SourceEndY-GP%SourceStartY)/GP%SourceDY)
    GP%nCalculations = GP%SourceNX*GP%SourceNY
endif
      
end subroutine computeGlobalParameters



subroutine allocateLayerVariables(GP, LP, SP)

use mpi    
use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(LayerParameters)    ::  LP(100)
type(StationParameters)  ::  SP(999)
integer*4                ::  iLayer, iSta
integer*4                ::  irank, master

irank = GP%irank; master = GP%master
do iLayer = 1,GP%NumLayers
    if(irank.ne.master) then ! bathymetry variables on master node has been allocated before
        ALLOCATE(LP(iLayer)%X(LP(iLayer)%NX))
        ALLOCATE(LP(iLayer)%Y(LP(iLayer)%NY))
        ALLOCATE(LP(iLayer)%Z(LP(iLayer)%NX,LP(iLayer)%NY))
    endif
    ALLOCATE(LP(iLayer)%H(2, LP(iLayer)%NX, LP(iLayer)%NY))
    ALLOCATE(LP(iLayer)%M(2, LP(iLayer)%NX-1, LP(iLayer)%NY))
    ALLOCATE(LP(iLayer)%N(2, LP(iLayer)%NX, LP(iLayer)%NY-1))
    ALLOCATE(LP(iLayer)%D_M(LP(iLayer)%NX-1, LP(iLayer)%NY))
    ALLOCATE(LP(iLayer)%D_N(LP(iLayer)%NX, LP(iLayer)%NY-1))
    ALLOCATE(LP(iLayer)%Hmax(LP(iLayer)%NX, LP(iLayer)%NY))
    ALLOCATE(LP(iLayer)%Hmin(LP(iLayer)%NX, LP(iLayer)%NY))
    if(LP(iLayer)%Level.gt.1) then !these variables are only used by child layers
        ALLOCATE(LP(iLayer)%H0(LP(iLayer)%NX, LP(iLayer)%NY))
        ALLOCATE(LP(iLayer)%M0(LP(iLayer)%NX-1, LP(iLayer)%NY))
        ALLOCATE(LP(iLayer)%N0(LP(iLayer)%NX, LP(iLayer)%NY-1))
        ALLOCATE(LP(iLayer)%HF(LP(iLayer)%NX, LP(iLayer)%NY))
        ALLOCATE(LP(iLayer)%MF(LP(iLayer)%NX-1, LP(iLayer)%NY))
        ALLOCATE(LP(iLayer)%NF(LP(iLayer)%NX, LP(iLayer)%NY-1))
    endif
    ALLOCATE(LP(iLayer)%Sign_M(LP(iLayer)%NX-1, LP(iLayer)%NY))
    ALLOCATE(LP(iLayer)%Sign_N(LP(iLayer)%NX, LP(iLayer)%NY-1))

    ALLOCATE(LP(iLayer)%PermanentDryCells((LP(iLayer)%MaxNX+10)*(LP(iLayer)%MaxNY+10),2))
    ALLOCATE(LP(iLayer)%FloodingCells((LP(iLayer)%MaxNX+10)*(LP(iLayer)%MaxNY+10),2))
    ALLOCATE(LP(iLayer)%FloodingCellsWaterHeight((LP(iLayer)%MaxNX+10)*(LP(iLayer)%MaxNY+10)))
    
    if(GP%KajiuraFilter.eq.1) then
        ALLOCATE(LP(iLayer)%HK(LP(iLayer)%NX, LP(iLayer)%NY))
    endif

    if(GP%BoundaryConditionType.eq.2.and.iLayer.eq.GP%TopLayer) then
        ALLOCATE(LP(iLayer)%SpongeMANNING(LP(iLayer)%NX, LP(iLayer)%NY))
    endif

    if(GP%CoordinatesType.eq.0) then
        ALLOCATE(LP(iLayer)%CPX(LP(iLayer)%NY))
        ALLOCATE(LP(iLayer)%CPY(LP(iLayer)%NY-1))
        ALLOCATE(LP(iLayer)%CSY(LP(iLayer)%NY-1))
        ALLOCATE(LP(iLayer)%CS1(LP(iLayer)%NY))
        ALLOCATE(LP(iLayer)%CS2(LP(iLayer)%NY))
        ALLOCATE(LP(iLayer)%CS3(LP(iLayer)%NY))
        if(LP(iLayer)%Nonlinearity.eq.1) then
            ALLOCATE(LP(iLayer)%CS6(LP(iLayer)%NY))
            ALLOCATE(LP(iLayer)%CS7(LP(iLayer)%NY-1))
            ALLOCATE(LP(iLayer)%CS8(LP(iLayer)%NY-1))
        endif
    endif

    if(LP(iLayer)%Dispersion.eq.1) then
        ALLOCATE(LP(iLayer)%Q(LP(iLayer)%NX,LP(iLayer)%NY))
        if(LP(iLayer)%Level.gt.1) then
            ALLOCATE(LP(iLayer)%Q0(LP(iLayer)%NX,LP(iLayer)%NY))
            ALLOCATE(LP(iLayer)%QF(LP(iLayer)%NX,LP(iLayer)%NY))
        endif
        ALLOCATE(LP(iLayer)%Has_Q(LP(iLayer)%NX, LP(iLayer)%NY))
        ALLOCATE(LP(iLayer)%WB(2,LP(iLayer)%NX,LP(iLayer)%NY), LP(iLayer)%WS(2,LP(iLayer)%NX,LP(iLayer)%NY))
        ALLOCATE(LP(iLayer)%PL(LP(iLayer)%NX,LP(iLayer)%NY), LP(iLayer)%PR(LP(iLayer)%NX,LP(iLayer)%NY))
        ALLOCATE(LP(iLayer)%PB(LP(iLayer)%NX,LP(iLayer)%NY), LP(iLayer)%PT(LP(iLayer)%NX,LP(iLayer)%NY))
        ALLOCATE(LP(iLayer)%PC(LP(iLayer)%NX,LP(iLayer)%NY), LP(iLayer)%PQ(LP(iLayer)%NX,LP(iLayer)%NY))
        if(GP%CoordinatesType.eq.0) then
            ALLOCATE(LP(iLayer)%CSD1(LP(iLayer)%NY), LP(iLayer)%CSD2(LP(iLayer)%NY))
            ALLOCATE(LP(iLayer)%CSD3(LP(iLayer)%NY), LP(iLayer)%CSD4(LP(iLayer)%NY))
        endif
    endif

    if(LP(iLayer)%Breaking.eq.1) then
        ALLOCATE(LP(iLayer)%BrkAge(LP(iLayer)%NX, LP(iLayer)%NY))
        ALLOCATE(LP(iLayer)%Brkv(LP(iLayer)%NX, LP(iLayer)%NY))
        if(LP(iLayer)%Level.gt.1) then
            ALLOCATE(LP(iLayer)%BrkAge0(LP(iLayer)%NX, LP(iLayer)%NY))
            ALLOCATE(LP(iLayer)%BrkAgeF(LP(iLayer)%NX, LP(iLayer)%NY))
            ALLOCATE(LP(iLayer)%Brkv0(LP(iLayer)%NX, LP(iLayer)%NY))
            ALLOCATE(LP(iLayer)%BrkvF(LP(iLayer)%NX, LP(iLayer)%NY))
        endif
        if(GP%CoordinatesType.eq.0) then
            ALLOCATE(LP(iLayer)%CSB1(LP(iLayer)%NY), LP(iLayer)%CSB2(LP(iLayer)%NY))
            ALLOCATE(LP(iLayer)%CSB3(LP(iLayer)%NY-1), LP(iLayer)%CSB4(LP(iLayer)%NY-1))
        endif
    endif

enddo

do iSta = 1,GP%NumStations
    ALLOCATE(SP(iSta)%H(GP%TotalTimeSteps+1))
    if(GP%SaveFlux.eq.1) then
        ALLOCATE(SP(iSta)%M(GP%TotalTimeSteps+1))
        ALLOCATE(SP(iSta)%N(GP%TotalTimeSteps+1))
    endif
    if(GP%SaveDynamicPressure.eq.1.and.LP(SP(iSta)%nLayer)%Dispersion.eq.1) then
        ALLOCATE(SP(iSta)%Q(GP%TotalTimeSteps+1))
    endif
enddo

end subroutine allocateLayerVariables


subroutine computeParameters(GP, LP, SP, FP)

use mpi    
use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(LayerParameters)    ::  LP(100)
type(StationParameters)  ::  SP(999)
type(FaultParameters)    ::  FP(4000)
integer*4  ::  irank, nsize, master
integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
integer*4  ::  jstart, jend, nstarty, nendy, iLayer, j, i, k, istart, iend
integer*4  ::  Iwidth, Jwidth, iin, jin
real*8     ::  manningx, manningy

irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
if(irank.eq.master) then
    write(*,*) 'computing parameters on computing nodes (e.g., Coriolis) ...'
    write(*,*)
endif

do iLayer = 1,GP%NumLayers
if(irank.lt.LP(iLayer)%nsize) then

    istart = MAX(LP(iLayer)%PartitionInfo(irank+1,1)-GP%nRowBoundary,1)
    iend   = MIN(LP(iLayer)%PartitionInfo(irank+1,2)+GP%nRowBoundary,LP(iLayer)%NX)
    jstart = MAX(LP(iLayer)%PartitionInfo(irank+1,3)-GP%nRowBoundary,1)
    jend   = MIN(LP(iLayer)%PartitionInfo(irank+1,4)+GP%nRowBoundary,LP(iLayer)%NY)
    !/// Sponge Damping coefficients ///!
    if(GP%BoundaryConditionType.eq.2.and.iLayer.eq.GP%TopLayer) then
        if(GP%CoordinatesType.eq.1) then
            Iwidth = MAX(NINT(GP%SpongeWidthX/LP(iLayer)%dx), 0)
            Jwidth = MAX(NINT(GP%SpongeWidthY/LP(iLayer)%dy), 0)
        elseif(GP%CoordinatesType.eq.0) then
            Iwidth = NINT(GP%SpongeWidthX/GP%R_EARTH/LP(iLayer)%dx/GP%PI*180.0 &
                /COS(0.5*(LP(iLayer)%ymin+LP(iLayer)%ymax)*GP%PI/180.0))
            Jwidth = NINT(GP%SpongeWidthY/GP%R_EARTH/LP(iLayer)%dy/GP%PI*180.0)
            Iwidth = MAX(Iwidth, 0); Jwidth = MAX(Jwidth, 0)
        endif
        do j = jstart,jend
        do i = istart,iend
            iin = MIN(i-1, LP(iLayer)%NX-i)
            jin = MIN(j-1, LP(iLayer)%NY-j)
            if(Iwidth.ge.2.and.iin.le.Iwidth-1) then
                manningx = GP%MaxSpongeMANNING*(1-TANH(10.0*iin/(Iwidth-1)))
            else
                manningx = 0.0d0
            endif
            if(Jwidth.ge.2.and.jin.le.Jwidth-1) then
                manningy = GP%MaxSpongeMANNING*(1-TANH(10.0*jin/(Jwidth-1)))
            else
                manningy = 0.0d0
            endif
            LP(iLayer)%SpongeMANNING(i,j) = MAX(manningx, manningy)
        enddo
        enddo
    endif

    !////// Permanent Dry Cells //////!
    LP(iLayer)%PermanentDryCellsCount = 0
    do j = jstart,jend
    do i = istart,iend
        if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit) then
            LP(iLayer)%PermanentDryCellsCount = LP(iLayer)%PermanentDryCellsCount+1
            LP(iLayer)%PermanentDryCells(LP(iLayer)%PermanentDryCellsCount,1) = i
            LP(iLayer)%PermanentDryCells(LP(iLayer)%PermanentDryCellsCount,2) = j
        endif
    enddo
    enddo

    !///// Assign Flux Signs to Permanent Wet and Dry Cells
    LP(iLayer)%Sign_M = 0; LP(iLayer)%Sign_N = 0
    do k = 1, LP(iLayer)%PermanentDryCellsCount
        i = LP(iLayer)%PermanentDryCells(k,1)
        j = LP(iLayer)%PermanentDryCells(k,2)
        if(i.le.LP(iLayer)%NX-1) LP(iLayer)%Sign_M(i,j) = 999
        if(i.ge.2) LP(iLayer)%Sign_M(i-1,j) = 999
        if(j.le.LP(iLayer)%NY-1) LP(iLayer)%Sign_N(i,j) = 999
        if(j.ge.2) LP(iLayer)%Sign_N(i,j-1) = 999
    enddo

    !////// Flooding Cells //////!
    LP(iLayer)%FloodingCellsCount = 0
    do j = jstart,jend
    do i = istart,iend
        if(LP(iLayer)%Z(i,j).lt.GP%PermanentDryLimit.and.LP(iLayer)%Z(i,j).gt.-GP%PermanentDryLimit) then
            LP(iLayer)%FloodingCellsCount = LP(iLayer)%FloodingCellsCount+1
            LP(iLayer)%FloodingCells(LP(iLayer)%FloodingCellsCount,1) = i
            LP(iLayer)%FloodingCells(LP(iLayer)%FloodingCellsCount,2) = j
        endif
    enddo
    enddo
                      
    !////// Coriolis parameters //////!
    if(GP%CoordinatesType.eq.0) then
        nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
        nendy = LP(iLayer)%PartitionInfo(irank+1,4)
        jstart = MAX(1,nstarty)
        jend = MIN(LP(iLayer)%NY,nendy)
        do j = jstart,jend
            LP(iLayer)%CPX(j)=2.0*GP%Omega*SIN(LP(iLayer)%Y(j)*GP%PI/180.0)*LP(iLayer)%dt
        enddo
        jend = MIN(LP(iLayer)%NY-1,nendy);
        do j = jstart,jend
            LP(iLayer)%CPY(j)=2.0*GP%Omega*SIN((LP(iLayer)%Y(j)+0.5*LP(iLayer)%dy)*GP%PI/180.0)*LP(iLayer)%dt
        enddo
    endif

    !//// mass/momentum equations parameters //////!
    if(GP%CoordinatesType.eq.1) then
        LP(iLayer)%CC1 = LP(iLayer)%dt/LP(iLayer)%dx
        LP(iLayer)%CC2 = LP(iLayer)%dt/LP(iLayer)%dy
        LP(iLayer)%CC3 = LP(iLayer)%CC1*GP%GRAV
        LP(iLayer)%CC4 = LP(iLayer)%CC2*GP%GRAV
    elseif(GP%CoordinatesType.eq.0) then
        nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
        nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
        jstart  = MAX(nstarty-GP%nRowBoundary,1)
        jend    = MIN(nendy+GP%nRowBoundary,LP(iLayer)%NY)
        do j = jstart,jend
            LP(iLayer)%CS1(j) = LP(iLayer)%dt/GP%R_EARTH/LP(iLayer)%dx/GP%PI*180.0/COS(LP(iLayer)%Y(j)*GP%PI/180.0)
            LP(iLayer)%CS2(j) = LP(iLayer)%CS1(j)*LP(iLayer)%dx/LP(iLayer)%dy
            LP(iLayer)%CS3(j) = LP(iLayer)%CS1(j)*GP%GRAV 
        enddo
            LP(iLayer)%CS4    = LP(iLayer)%dt*GP%GRAV/GP%R_EARTH/LP(iLayer)%dy/GP%PI*180.0
        do j = jstart, MIN(jend,LP(iLayer)%NY-1)
            LP(iLayer)%CSY(j) = COS((LP(iLayer)%Y(j)+0.5*LP(iLayer)%dy)*GP%PI/180.0)
        enddo
        if(LP(iLayer)%Nonlinearity.eq.1) then
            LP(iLayer)%CS5    = LP(iLayer)%dt/GP%R_EARTH/LP(iLayer)%dy/GP%PI*180.0
            do j = jstart,jend
                LP(iLayer)%CS6(j) = LP(iLayer)%dt*TAN(LP(iLayer)%Y(j)*GP%PI/180.0)/GP%R_EARTH
            enddo
            do j = jstart,MIN(jend,LP(iLayer)%NY-1)
                LP(iLayer)%CS7(j) = LP(iLayer)%dt/GP%R_EARTH/LP(iLayer)%dx/GP%PI*180.0/LP(iLayer)%CSY(j)
                LP(iLayer)%CS8(j) = LP(iLayer)%dt*TAN((LP(iLayer)%Y(j)+0.5*LP(iLayer)%dy)*GP%PI/180.0)/GP%R_EARTH
            enddo
        endif   
    endif

    !/// dispersion parameters ///!
    if(LP(iLayer)%Dispersion.eq.1) then
        if(GP%CoordinatesType.eq.1) then
            LP(iLayer)%CCD1 = LP(iLayer)%dt/LP(iLayer)%dx**2
            LP(iLayer)%CCD2 = LP(iLayer)%dt/LP(iLayer)%dy**2
            LP(iLayer)%CCD3 = 1.0/LP(iLayer)%dx
            LP(iLayer)%CCD4 = 1.0/LP(iLayer)%dy
        elseif(GP%CoordinatesType.eq.0) then
            nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
            nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
            jstart  = MAX(nstarty-GP%nRowBoundary,1)
            jend    = MIN(nendy+GP%nRowBoundary,LP(iLayer)%NY)
            do j = jstart,jend
                LP(iLayer)%CSD1(j) = LP(iLayer)%dt/ &
                    (GP%R_EARTH*LP(iLayer)%dx*GP%PI/180.0*COS(LP(iLayer)%Y(j)*GP%PI/180.0))**2
                LP(iLayer)%CSD2(j) = LP(iLayer)%dt/ &
                    (GP%R_EARTH*LP(iLayer)%dy*GP%PI/180.0)**2/COS(LP(iLayer)%Y(j)*GP%PI/180.0)
                LP(iLayer)%CSD3(j) = 1.0/(GP%R_EARTH*LP(iLayer)%dx*GP%PI/180.0*COS(LP(iLayer)%Y(j)*GP%PI/180.0))
                LP(iLayer)%CSD4(j) = 1.0/(GP%R_EARTH*LP(iLayer)%dy*GP%PI/180.0*COS(LP(iLayer)%Y(j)*GP%PI/180.0))
            enddo
            LP(iLayer)%CSD5 = 1.0/GP%R_EARTH/LP(iLayer)%dy/GP%PI*180.0
            LP(iLayer)%CS5  = LP(iLayer)%dt/GP%R_EARTH/LP(iLayer)%dy/GP%PI*180.0
        endif
    endif

    !/// breaking parameters ///!
    if(LP(iLayer)%Breaking.eq.1) then
        if(GP%CoordinatesType.eq.1) then
            LP(iLayer)%CCB1 = LP(iLayer)%dt/LP(iLayer)%dx**2
            LP(iLayer)%CCB2 = LP(iLayer)%dt/LP(iLayer)%dy**2
        elseif(GP%CoordinatesType.eq.0) then
            nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
            nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
            jstart  = MAX(nstarty-GP%nRowBoundary,1)
            jend    = MIN(nendy+GP%nRowBoundary,LP(iLayer)%NY)
            do j = jstart,jend
                LP(iLayer)%CSB1(j) = LP(iLayer)%dt/ &
                    (GP%R_EARTH*LP(iLayer)%dx*GP%PI/180.0*COS(LP(iLayer)%Y(j)*GP%PI/180.0))**2
                LP(iLayer)%CSB2(j) = LP(iLayer)%dt/ &
                    (GP%R_EARTH*LP(iLayer)%dy*GP%PI/180.0)**2/COS(LP(iLayer)%Y(j)*GP%PI/180.0)
            enddo
            do j = jstart,MIN(jend,LP(iLayer)%NY-1)
                LP(iLayer)%CSB3(j) = LP(iLayer)%dt/ &
                    (GP%R_EARTH*LP(iLayer)%dx*GP%PI/180.0*LP(iLayer)%CSY(j))**2
                LP(iLayer)%CSB4(j) = LP(iLayer)%dt/ &
                    (GP%R_EARTH*LP(iLayer)%dy*GP%PI/180.0)**2/LP(iLayer)%CSY(j)
            enddo
        endif
    endif

endif !if:this node is used by this layer
enddo !loop for all layers

end subroutine computeParameters



subroutine removeStationFiles(GP, LP, SP)

use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(LayerParameters)    ::  LP(100)
type(StationParameters)  ::  SP(999)
integer*4                ::  i
character(999)           ::  s

do i = 1,GP%NumStations
    write(s,'(a,i4.4,a4)') TRIM(ADJUSTL(GP%ResultPath))//'Station',i,'.dat'
    open(1000+i,file=s,status='replace',form='unformatted')
    write(1000+i) GP%PurposeCalculation, GP%nCalculations, GP%TotalTimeSteps+1
    close(1000+i)
enddo
if(GP%SaveFlux.eq.1) then
    do i = 1,GP%NumStations
        write(s,'(a,i4.4,a6)') TRIM(ADJUSTL(GP%ResultPath))//'Station',i,'_M.dat'
        open(1000+i,file=s,status='replace',form='unformatted')
        write(1000+i) GP%PurposeCalculation, GP%nCalculations, GP%TotalTimeSteps+1
        close(1000+i)
    enddo
    do i = 1,GP%NumStations
        write(s,'(a,i4.4,a6)') TRIM(ADJUSTL(GP%ResultPath))//'Station',i,'_N.dat'
        open(1000+i,file=s,status='replace',form='unformatted')
        write(1000+i) GP%PurposeCalculation, GP%nCalculations, GP%TotalTimeSteps+1
        close(1000+i)
    enddo
endif
if(GP%SaveDynamicPressure.eq.1) then
    do i = 1,GP%NumStations
        if(LP(SP(i)%nLayer)%Dispersion.eq.1) then
            write(s,'(a,i4.4,a6)') TRIM(ADJUSTL(GP%ResultPath))//'Station',i,'_Q.dat'
            open(1000+i,file=s,status='replace',form='unformatted')
            write(1000+i) GP%PurposeCalculation, GP%nCalculations, GP%TotalTimeSteps+1
            close(1000+i)
        endif
    enddo
endif

end subroutine removeStationFiles



subroutine modifyFaultParameters(GP, FP, iCal)

use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(FaultParameters)    ::  FP(4000)
integer*4                ::  iCal, iFault
integer*4                ::  irank, nsize, master

irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
if(irank.eq.master) then
    write(*,*)
    write(*,*) '***************************************************'
    write(*,'(a,i4,a,i4)') ' Green function for subfault ',iCal,' of ',GP%NumFaults
    write(*,*) '***************************************************'
    write(*,*)
endif
do iFault = 1,GP%Numfaults
    if(iFault.eq.iCal) then
        FP(iFault)%T0 = 0
        FP(iFault)%NT = 0
        FP(iFault)%Slip = 1.0
    else
        FP(iFault)%T0 = GP%TotalTime*2
        FP(iFault)%NT = GP%TotalTimeSteps*2
    endif
enddo

end subroutine modifyFaultParameters



subroutine getInitialCondition(GP, LP, SP, FP, iCal, iTimeStep, LocalData, LocalDataLength)

use mpi
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
type(StationParameters) ::  SP(999)
type(FaultParameters)   ::  FP(4000)
integer*4   ::  iCal, iTimeStep
integer*4   ::  LocalDataLength
real*8      ::  LocalData(LocalDataLength)
integer*4   ::  irank, nsize, master
integer*4   ::  errorcode, ierror, ios, istatus(MPI_STATUS_SIZE)
integer*4   ::  iLayer, iLayerLevel, iNode, iFault, iFlag, parent
integer*4   ::  iHMN, i, j, ix, iy, idx
integer*4   ::  istartx, iendx, istarty, iendy
integer*4   ::  nRowBoundary, nRB, s, nFaultsNow
real*8      ::  x, y, val, z1, z2, h1, h2
real*8      ::  h, h_horizontal, u1, u2 
real*8      ::  dx, dy, zx, zy, z_slope
real*8      ::  depth,xlength,xwidth
real*8      ::  slip,strike,dip,rake,x0,y0,lx2,ly2,x_left,x_right,y_left,y_right,lx,ly
real*8      ::  kajiuraRmax=6.0d0
character(999)  ::  fn
logical         ::  fe

irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
nRowBoundary = GP%nRowBoundary; nRB= MAX(GP%nRowBoundary, GP%nRowBoundaryFlux)

if(GP%PurposeCalculation.eq.1.and.GP%InitialConditionType.eq.0.and.iTimeStep.eq.0) then

do iHMN = 1,3 !read initial value of water elevation, flux M and flux N
    if(iHMN.eq.1) fn = TRIM(ADJUSTL(GP%InitialElevationFileName))
    if(iHMN.eq.2) fn = TRIM(ADJUSTL(GP%InitialMFileName))
    if(iHMN.eq.3) fn = TRIM(ADJUSTL(GP%InitialNFileName))
    if(iHMN.ne.1) then !only initial elevation is required, initial flux is optional
        inquire(file=fn, exist=fe)
        if(.not.fe) cycle
    endif
    !///read initial H/M/N of top layer on the master node///!
    if(irank.eq.master) then
        if(iHMN.eq.1) then
            write(*,*) 'reading initial water elevation from file ',TRIM(ADJUSTL(GP%InitialElevationFileName)),' ...'
        elseif(iHMN.eq.2) then
            write(*,*) 'reading initial flux component M from file ',TRIM(ADJUSTL(GP%InitialMFileName)),' ...'
            write(*,*) 'ATTENTION: M is defined at (x+0.5dx,y), and its dimension is (NX-1,NY)!'
        elseif(iHMN.eq.3) then
            write(*,*) 'reading initial flux component N from file ',TRIM(ADJUSTL(GP%InitialNFileName)),' ...'
            write(*,*) 'ATTENTION: N is defined at (x,y+0.5dy), and its dimension is (NX,NY-1)!'
        endif
        write(*,*)

        if(iHMN.eq.1.and.GP%InitialElevationFileFormat.eq.1) then
            call readInitialElevationNetCDF(GP, LP)
        else
            open(23,file=TRIM(ADJUSTL(fn)),status='old',action='read',form='formatted')
            if(iHMN.eq.1) then
                iendx = LP(GP%TopLayer)%NX; iendy = LP(GP%TopLayer)%NY
            elseif(iHMN.eq.2) then
                iendx = LP(GP%TopLayer)%NX-1; iendy = LP(GP%TopLayer)%NY
            elseif(iHMN.eq.3) then
                iendx = LP(GP%TopLayer)%NX; iendy = LP(GP%TopLayer)%NY-1
            endif
            do j = 1,iendy
            do i = 1,iendx
                read(23,*,iostat=ios) x, y, val
                if(ios.ne.0) then
                    write(*,*) 'ERROR: something wrong when reading the initial file: ', TRIM(ADJUSTL(fn))
                    call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
                endif
                if(iHMN.eq.2.and.ABS(x-LP(GP%TopLayer)%X(i)-0.5*LP(GP%TopLayer)%dx).gt.0.1*LP(GP%TopLayer)%dx) then
                    write(*,*) 'ERROR: the file for initial flux is not given properly: ', TRIM(ADJUSTL(fn))
                    call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
                endif
                if(iHMN.eq.3.and.ABS(y-LP(GP%TopLayer)%Y(j)-0.5*LP(GP%TopLayer)%dy).gt.0.1*LP(GP%TopLayer)%dy) then
                    write(*,*) 'ERROR: the file for initial flux is not given properly: ', TRIM(ADJUSTL(fn))
                    call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
                endif
                if(iHMN.eq.1) LP(GP%TopLayer)%H(1,i,j) = val
                if(iHMN.eq.2) LP(GP%TopLayer)%M(1,i,j) = val
                if(iHMN.eq.3) LP(GP%TopLayer)%N(1,i,j) = val
            enddo
            enddo
            close(23)
        endif
    endif
    !///apply Kajiura filter to the surface elevation///!
    if(iHMN.eq.1.and.GP%KajiuraFilter.eq.1) then
        if(irank.eq.master) then
            LP(GP%TopLayer)%HK = LP(GP%TopLayer)%H(1,:,:)
            do iNode = 0,LP(GP%TopLayer)%nsize-1
                if(iNode.ne.master) then
                    istartx = LP(GP%TopLayer)%PartitionInfo(iNode+1,1)
                    iendx   = LP(GP%TopLayer)%PartitionInfo(iNode+1,2)
                    istarty = LP(GP%TopLayer)%PartitionInfo(iNode+1,3)
                    iendy   = LP(GP%TopLayer)%PartitionInfo(iNode+1,4)
                    do i = istartx,iendx
                    do j = istarty,iendy
                        LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1)) = LP(GP%TopLayer)%HK(i,j)
                    enddo
                    enddo
                    call MPI_SEND(LocalData,(iendx-istartx+1)*(iendy-istarty+1), &
                        MPI_DOUBLE_PRECISION,iNode,2023,MPI_COMM_WORLD,ierror)
                endif
            enddo
        elseif(irank.lt.LP(GP%TopLayer)%nsize) then
            istartx = LP(GP%TopLayer)%PartitionInfo(irank+1,1)
            iendx   = LP(GP%TopLayer)%PartitionInfo(irank+1,2)
            istarty = LP(GP%TopLayer)%PartitionInfo(irank+1,3)
            iendy   = LP(GP%TopLayer)%PartitionInfo(irank+1,4)
            call MPI_RECV(LocalData,(iendx-istartx+1)*(iendy-istarty+1),  &
                MPI_DOUBLE_PRECISION,master,2023,MPI_COMM_WORLD,istatus,ierror)
            do i = istartx,iendx
            do j = istarty,iendy
                LP(GP%TopLayer)%HK(i,j) = LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1))
            enddo
            enddo
        endif
        call smoothElevation_Kajiura(GP, LP, GP%TopLayer, kajiuraRmax, LocalData, LocalDataLength)
        if(irank.eq.master) LP(GP%TopLayer)%H(1,:,:) = LP(GP%TopLayer)%HK  
    endif
    !///modify the surface elevation of the top layer, and then interpolate it into all child layers///!
    if(irank.eq.master) then
        if(iHMN.eq.1) then
            iendx = LP(GP%TopLayer)%NX; iendy = LP(GP%TopLayer)%NY
        elseif(iHMN.eq.2) then
            iendx = LP(GP%TopLayer)%NX-1; iendy = LP(GP%TopLayer)%NY
        elseif(iHMN.eq.3) then
            iendx = LP(GP%TopLayer)%NX; iendy = LP(GP%TopLayer)%NY-1
        endif
        do j = 1,iendy
        do i = 1,iendx
            if(iHMN.eq.1) then
                if(LP(GP%TopLayer)%Z(i,j).le.-GP%PermanentDryLimit) then 
                    LP(GP%TopLayer)%H(1,i,j) = 0.0d0
                elseif(LP(GP%TopLayer)%Z(i,j)+LP(GP%TopLayer)%H(1,i,j).le.0.0) then
                    if(LP(GP%TopLayer)%Z(i,j).gt.0.0) then
                        LP(GP%TopLayer)%H(1,i,j) = -LP(GP%TopLayer)%Z(i,j)
                    else
                        LP(GP%TopLayer)%H(1,i,j) = 0.0d0
                    endif
                endif
            elseif(iHMN.eq.2) then
                z1 = LP(GP%TopLayer)%Z(i,j);   z2 = LP(GP%TopLayer)%Z(i+1,j)
                h1 = LP(GP%TopLayer)%H(1,i,j); h2 = LP(GP%TopLayer)%H(1,i+1,j)
                call checkFluxDirection(GP%PermanentDryLimit,z1,z2,h1,h2,s)
                if(s.eq.999.or.s*LP(GP%TopLayer)%M(1,i,j).lt.0.0) LP(GP%TopLayer)%M(1,i,j) = 0.0d0
            elseif(iHMN.eq.3) then
                z1 = LP(GP%TopLayer)%Z(i,j);   z2 = LP(GP%TopLayer)%Z(i,j+1)
                h1 = LP(GP%TopLayer)%H(1,i,j); h2 = LP(GP%TopLayer)%H(1,i,j+1)
                call checkFluxDirection(GP%PermanentDryLimit,z1,z2,h1,h2,s)
                if(s.eq.999.or.s*LP(GP%TopLayer)%N(1,i,j).lt.0.0) LP(GP%TopLayer)%N(1,i,j) = 0.0d0
            endif
        enddo
        enddo    
        do iLayerLevel = 2, GP%NumLayerLevels
        do iLayer = 1, GP%NumLayers
        if(LP(iLayer)%Level.eq.iLayerLevel) then
            if(iHMN.eq.1) then
                iendx = LP(iLayer)%NX; iendy = LP(iLayer)%NY
            elseif(iHMN.eq.2) then
                iendx = LP(iLayer)%NX-1; iendy = LP(iLayer)%NY
            elseif(iHMN.eq.3) then
                iendx = LP(iLayer)%NX; iendy = LP(iLayer)%NY-1
            endif
            do j = 1,iendy
            do i = 1,iendx
                x = LP(iLayer)%X(i); y = LP(iLayer)%Y(j)
                if(iHMN.eq.2) x = x + 0.5*LP(iLayer)%dx
                if(iHMN.eq.3) y = y + 0.5*LP(iLayer)%dy
                call interpData(GP, LP, LP(iLayer)%Parent, -iHMN, x, y, val)
                if(iHMN.eq.1) LP(iLayer)%H(1,i,j) = val
                if(iHMN.eq.2) LP(iLayer)%M(1,i,j) = val
                if(iHMN.eq.3) LP(iLayer)%N(1,i,j) = val
                if(iHMN.eq.1) then
                    if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit) then 
                        LP(iLayer)%H(1,i,j) = 0.0d0
                    elseif(LP(iLayer)%Z(i,j)+LP(iLayer)%H(1,i,j).le.0.0) then
                        if(LP(iLayer)%Z(i,j).gt.0.0) then
                            LP(iLayer)%H(1,i,j) = -LP(iLayer)%Z(i,j)
                        else
                            LP(iLayer)%H(1,i,j) = 0.0d0
                        endif
                    endif
                elseif(iHMN.eq.2) then
                    z1 = LP(iLayer)%Z(i,j);   z2 = LP(iLayer)%Z(i+1,j)
                    h1 = LP(iLayer)%H(1,i,j); h2 = LP(iLayer)%H(1,i+1,j)
                    call checkFluxDirection(GP%PermanentDryLimit,z1,z2,h1,h2,s)
                    if(s.eq.999.or.s*LP(iLayer)%M(1,i,j).lt.0.0) LP(iLayer)%M(1,i,j) = 0.0d0
                elseif(iHMN.eq.3) then
                    z1 = LP(iLayer)%Z(i,j);   z2 = LP(iLayer)%Z(i,j+1)
                    h1 = LP(iLayer)%H(1,i,j); h2 = LP(iLayer)%H(1,i,j+1)
                    call checkFluxDirection(GP%PermanentDryLimit,z1,z2,h1,h2,s)
                    if(s.eq.999.or.s*LP(iLayer)%N(1,i,j).lt.0.0) LP(iLayer)%N(1,i,j) = 0.0d0
                endif
            enddo
            enddo
        endif
        enddo !loop for layers
        enddo !loop for layer levels
    endif
    !///broadcast initial H/M/N of all layers from master to other nodes///!
    if(irank.eq.master) then
        do iLayer = 1,GP%NumLayers
            do iNode = 0, LP(iLayer)%nsize-1
                if(iNode.ne.master) then
                    istartx = MAX(1,LP(iLayer)%PartitionInfo(iNode+1,1)-nRB)
                    istarty = MAX(1,LP(iLayer)%PartitionInfo(iNode+1,3)-nRB)
                    if(iHMN.eq.1) then
                        iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(iNode+1,2)+nRB)
                        iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(iNode+1,4)+nRB)
                    elseif(iHMN.eq.2) then
                        iendx   = MIN(LP(iLayer)%NX-1,LP(iLayer)%PartitionInfo(iNode+1,2)+nRB)
                        iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(iNode+1,4)+nRB)
                    elseif(iHMN.eq.3) then
                        iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(iNode+1,2)+nRB)
                        iendy   = MIN(LP(iLayer)%NY-1,LP(iLayer)%PartitionInfo(iNode+1,4)+nRB)
                    endif
                    do i = istartx,iendx
                    do j = istarty,iendy
                        if(iHMN.eq.1) then
                            LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1)) = LP(iLayer)%H(1,i,j)
                        elseif(iHMN.eq.2) then
                            LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1)) = LP(iLayer)%M(1,i,j)
                        elseif(iHMN.eq.3) then
                            LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1)) = LP(iLayer)%N(1,i,j)
                        endif
                    enddo
                    enddo
                    call MPI_SEND(LocalData,(iendx-istartx+1)*(iendy-istarty+1), &
                        MPI_DOUBLE_PRECISION,iNode,iLayer,MPI_COMM_WORLD,ierror)
                endif
            enddo
        enddo
    else
        do iLayer = 1, GP%NumLayers
        if(irank.lt.LP(iLayer)%nsize) then
            istartx = MAX(1,LP(iLayer)%PartitionInfo(irank+1,1)-nRB)
            istarty = MAX(1,LP(iLayer)%PartitionInfo(irank+1,3)-nRB)
            if(iHMN.eq.1) then
                iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(irank+1,2)+nRB)
                iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(irank+1,4)+nRB)
            elseif(iHMN.eq.2) then
                iendx   = MIN(LP(iLayer)%NX-1,LP(iLayer)%PartitionInfo(irank+1,2)+nRB)
                iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(irank+1,4)+nRB)
            elseif(iHMN.eq.3) then
                iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(irank+1,2)+nRB)
                iendy   = MIN(LP(iLayer)%NY-1,LP(iLayer)%PartitionInfo(irank+1,4)+nRB)
            endif
            call MPI_RECV(LocalData,(iendx-istartx+1)*(iendy-istarty+1), &
                MPI_DOUBLE_PRECISION,master,iLayer,MPI_COMM_WORLD,istatus,ierror)
            do i = istartx,iendx
            do j = istarty,iendy
                if(iHMN.eq.1) then
                    LP(iLayer)%H(1,i,j) = LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1))
                elseif(iHMN.eq.2) then
                    LP(iLayer)%M(1,i,j) = LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1))
                elseif(iHMN.eq.3) then
                    LP(iLayer)%N(1,i,j) = LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1))
                endif
                if(iHMN.eq.1) then
                    if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit) then 
                        LP(iLayer)%H(1,i,j) = 0.0d0
                    elseif(LP(iLayer)%Z(i,j)+LP(iLayer)%H(1,i,j).le.0.0) then
                        if(LP(iLayer)%Z(i,j).gt.0.0) then
                            LP(iLayer)%H(1,i,j) = -LP(iLayer)%Z(i,j)
                        else
                            LP(iLayer)%H(1,i,j) = 0.0d0
                        endif
                    endif
                elseif(iHMN.eq.2) then
                    z1 = LP(iLayer)%Z(i,j);   z2 = LP(iLayer)%Z(i+1,j)
                    h1 = LP(iLayer)%H(1,i,j); h2 = LP(iLayer)%H(1,i+1,j)
                    call checkFluxDirection(GP%PermanentDryLimit,z1,z2,h1,h2,s)
                    if(s.eq.999.or.s*LP(iLayer)%M(1,i,j).lt.0.0) LP(iLayer)%M(1,i,j) = 0.0d0
                elseif(iHMN.eq.3) then
                    z1 = LP(iLayer)%Z(i,j);   z2 = LP(iLayer)%Z(i,j+1)
                    h1 = LP(iLayer)%H(1,i,j); h2 = LP(iLayer)%H(1,i,j+1)
                    call checkFluxDirection(GP%PermanentDryLimit,z1,z2,h1,h2,s)
                    if(s.eq.999.or.s*LP(iLayer)%N(1,i,j).lt.0.0) LP(iLayer)%N(1,i,j) = 0.0d0
                endif
            enddo
            enddo
        endif
        enddo
    endif
enddo !loop for initial H/M/N
endif

if((GP%PurposeCalculation.eq.1.and.GP%InitialConditionType.eq.1.or.GP%PurposeCalculation.eq.2).and.&
    GP%KajiuraFilter.eq.0) then
    nFaultsNow = 0
    do iFault = 1,GP%NumFaults
        if((GP%PurposeCalculation.eq.1.and.FP(iFault)%NT.eq.iTimeStep.and.ABS(FP(iFault)%Slip).gt.1e-5).or. &
            GP%PurposeCalculation.eq.2.and.iTimeStep.eq.0.and.iFault.eq.iCal) then
                nFaultsNow = nFaultsNow + 1
            if(irank.eq.master) then
                if(iTimeStep.ne.0.and.nFaultsNow.eq.1) then
                    write(*,*)
                    write(*,*) '-------------------- New Tsunami Source -----------------------------'
                endif
                write(*,'(a,i5,a,i5)') 'calculating water elevation from Okada''s model: subfault ', &
                    iFault,'    of',GP%NumFaults
            endif
            do iLayer = 1,GP%NumLayers
            if(irank.lt.LP(iLayer)%nsize) then
                istartx = MAX(1,LP(iLayer)%PartitionInfo(irank+1,1)-nRowBoundary)
                iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(irank+1,2)+nRowBoundary)
                istarty = MAX(1,LP(iLayer)%PartitionInfo(irank+1,3)-nRowBoundary)
                iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(irank+1,4)+nRowBoundary)
                do i = istartx,iendx
                do j = istarty,iendy
                    if(LP(iLayer)%Z(i,j).gt.0.0.and.LP(iLayer)%H(1,i,j)+LP(iLayer)%Z(i,j).gt.GP%MinWaterDepth) then
                        x = LP(iLayer)%X(i); y = LP(iLayer)%Y(j)
                        depth = FP(iFault)%Depth; xlength = FP(iFault)%Length
                        xwidth = FP(iFault)%Width; strike = FP(iFault)%Strike
                        dip = FP(iFault)%Dip; slip = FP(iFault)%Slip
                        rake = FP(iFault)%Rake; x0 = FP(iFault)%X0; y0 = FP(iFault)%Y0
                        !///calculate water depth gradient///!
                        if(GP%CoordinatesType.eq.1) then
                            dx = LP(iLayer)%dx; dy = LP(iLayer)%dy
                        elseif(GP%CoordinatesType.eq.0) then
                            dx = GP%R_EARTH*LP(iLayer)%dx*GP%PI/180.0*COS(y*GP%PI/180.0)
                            dy = GP%R_EARTH*LP(iLayer)%dy*GP%PI/180.0
                        endif
                        if(i.eq.1.or.i.eq.LP(iLayer)%NX) then
                            zx = 0.0d0
                        else
                            zx = 0.5*(LP(iLayer)%Z(i+1,j)-LP(iLayer)%Z(i-1,j))/dx
                        endif
                        if(j.eq.1.or.j.eq.LP(iLayer)%NY) then
                            zy = 0.0d0
                        else
                            zy = 0.5*(LP(iLayer)%Z(i,j+1)-LP(iLayer)%Z(i,j-1))/dy
                        endif
                        z_slope = SQRT(zx*zx+zy*zy)
                        if(z_slope.gt.GP%MaxBottomSlope) then
                            zx = zx/z_slope*GP%MaxBottomSlope
                            zy = zy/z_slope*GP%MaxBottomSlope
                        endif
                        !///calculate bottom displacement using OKADA model///!
                        if(GP%CoordinatesType.eq.0) then
                            x = GP%R_Earth*COS(y0/180.0*GP%PI)*(x-x0)/180.0*GP%PI
                            y = GP%R_Earth*(y-y0)/180.0*GP%PI
                            x0 = 0.0
                            y0 = 0.0
                        endif
                        if(GP%PurposeCalculation.eq.2) slip = 1.0d0
                        call okada1985(x,y,u1,u2,h,depth,xlength,xwidth,slip,strike,dip,rake,x0,y0)
                        LP(iLayer)%H(1,i,j) = LP(iLayer)%H(1,i,j) + h
                        !///consider the contribution of horizontal motion///!
                        if(GP%HorizontalMotion.eq.1) then
                            h_horizontal = u1*zx + u2*zy
                            if(ABS(h_horizontal).gt.ABS(h)) h_horizontal = SIGN(h,h_horizontal)
                            LP(iLayer)%H(1,i,j) = LP(iLayer)%H(1,i,j) + h_horizontal
                        endif
                    endif
                enddo
                enddo
            endif
            enddo
        endif
    enddo
    if(irank.eq.master.and.iTimeStep.eq.0) write(*,*)
endif

if((GP%PurposeCalculation.eq.1.and.GP%InitialConditionType.eq.1.or.GP%PurposeCalculation.eq.2).and.&
    GP%KajiuraFilter.eq.1) then
    nFaultsNow = 0
    do iFault = 1,GP%NumFaults
        if((GP%PurposeCalculation.eq.1.and.FP(iFault)%NT.eq.iTimeStep.and.ABS(FP(iFault)%Slip).gt.1e-5).or. &
            GP%PurposeCalculation.eq.2.and.iTimeStep.eq.0.and.iFault.eq.iCal) then
            nFaultsNow = nFaultsNow + 1
        endif
    enddo
    if(nFaultsNow.eq.0) return
    LP(GP%TopLayer)%HK = 0.0d0; nFaultsNow = 0
    do iFault = 1,GP%NumFaults
        if((GP%PurposeCalculation.eq.1.and.FP(iFault)%NT.eq.iTimeStep.and.ABS(FP(iFault)%Slip).gt.1e-5).or. &
            GP%PurposeCalculation.eq.2.and.iTimeStep.eq.0.and.iFault.eq.iCal) then
            nFaultsNow = nFaultsNow + 1
            if(irank.eq.master) then
                if(iTimeStep.ne.0.and.nFaultsNow.eq.1) then
                    write(*,*)
                    write(*,*) '-------------------- New Tsunami Source -----------------------------'
                endif
                write(*,'(a,i5,a,i5)') 'calculating water elevation from Okada''s model: subfault ', &
                    iFault,'    of',GP%NumFaults
            endif
            !///calculate the ocean bottom deformation of the top layer///!
            if(irank.lt.LP(GP%TopLayer)%nsize) then
                istartx = LP(GP%TopLayer)%PartitionInfo(irank+1,1)
                iendx   = LP(GP%TopLayer)%PartitionInfo(irank+1,2)
                istarty = LP(GP%TopLayer)%PartitionInfo(irank+1,3)
                iendy   = LP(GP%TopLayer)%PartitionInfo(irank+1,4)
                do i = istartx,iendx
                do j = istarty,iendy
                    if(LP(GP%TopLayer)%Z(i,j).le.-GP%PermanentDryLimit) cycle
                    x = LP(GP%TopLayer)%X(i); y = LP(GP%TopLayer)%Y(j)
                    depth = FP(iFault)%Depth; xlength = FP(iFault)%Length
                    xwidth = FP(iFault)%Width; strike = FP(iFault)%Strike
                    dip = FP(iFault)%Dip; slip = FP(iFault)%Slip
                    rake = FP(iFault)%Rake; x0 = FP(iFault)%X0; y0 = FP(iFault)%Y0
                    !///calculate water depth gradient///!
                    if(GP%CoordinatesType.eq.1) then
                        dx = LP(GP%TopLayer)%dx; dy = LP(GP%TopLayer)%dy
                    elseif(GP%CoordinatesType.eq.0) then
                        dx = GP%R_EARTH*LP(GP%TopLayer)%dx*GP%PI/180.0*COS(y*GP%PI/180.0)
                        dy = GP%R_EARTH*LP(GP%TopLayer)%dy*GP%PI/180.0
                    endif
                    if(i.eq.1.or.i.eq.LP(GP%TopLayer)%NX) then
                        zx = 0.0d0
                    else
                        zx = 0.5*(LP(GP%TopLayer)%Z(i+1,j)-LP(GP%TopLayer)%Z(i-1,j))/dx
                    endif
                    if(j.eq.1.or.j.eq.LP(GP%TopLayer)%NY) then
                        zy = 0.0d0
                    else
                        zy = 0.5*(LP(GP%TopLayer)%Z(i,j+1)-LP(GP%TopLayer)%Z(i,j-1))/dy
                    endif
                    z_slope = SQRT(zx*zx+zy*zy)
                    if(z_slope.gt.GP%MaxBottomSlope) then
                        zx = zx/z_slope*GP%MaxBottomSlope
                        zy = zy/z_slope*GP%MaxBottomSlope
                    endif
                    !///calculate bottom displacement using OKADA model///!
                    if(GP%CoordinatesType.eq.0) then
                        x = GP%R_Earth*COS(y0/180.0*GP%PI)*(x-x0)/180.0*GP%PI
                        y = GP%R_Earth*(y-y0)/180.0*GP%PI
                        x0 = 0.0
                        y0 = 0.0
                    endif
                    if(GP%PurposeCalculation.eq.2) slip = 1.0d0
                    call okada1985(x,y,u1,u2,h,depth,xlength,xwidth,slip,strike,dip,rake,x0,y0)
                    LP(GP%TopLayer)%HK(i,j) = LP(GP%TopLayer)%HK(i,j) + h
                    !///consider the contribution of horizontal motion///!
                    if(GP%HorizontalMotion.eq.1) then
                        h_horizontal = u1*zx + u2*zy
                        if(ABS(h_horizontal).gt.ABS(h)) h_horizontal = SIGN(h,h_horizontal)
                        LP(GP%TopLayer)%HK(i,j) = LP(GP%TopLayer)%HK(i,j) + h_horizontal
                    endif
                enddo
                enddo
            endif
        endif
    enddo !loop for subfaults
    !///apply Kajiura filter to the bottom elevation///!
    call smoothElevation_Kajiura(GP, LP, GP%TopLayer, kajiuraRmax, LocalData, LocalDataLength)
    !///interpolate the smoothed bottom elevation to all child layers///!
    if(irank.eq.master) then
        do iLayerLevel = 2, GP%NumLayerLevels
        do iLayer = 1, GP%NumLayers
            if(LP(iLayer)%Level.eq.iLayerLevel) then
                parent = LP(iLayer)%Parent
                do j = 1,LP(iLayer)%NY
                do i = 1,LP(iLayer)%NX
                    if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit) then
                        LP(iLayer)%HK(i,j) = 0.0d0
                    else
                        x = LP(iLayer)%X(i); y = LP(iLayer)%Y(j)
                        call interpolate2D_linear(LP(parent)%HK, LP(parent)%X, LP(parent)%Y, &
                            LP(parent)%NX, LP(parent)%NY, x, y, val)
                        LP(iLayer)%HK(i,j) = val
                    endif
                enddo
                enddo
            endif
        enddo
        enddo
    endif !if master node

    !///broadcast smoothed bottom elevation of all layers///!
    do iLayer = 1,GP%NumLayers
    if(irank.eq.master) then
        do iNode = 0,LP(iLayer)%nsize-1
            istartx = MAX(1,LP(iLayer)%PartitionInfo(iNode+1,1)-nRowBoundary)
            iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(iNode+1,2)+nRowBoundary)
            istarty = MAX(1,LP(iLayer)%PartitionInfo(iNode+1,3)-nRowBoundary)
            iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(iNode+1,4)+nRowBoundary)
            if(iNode.ne.master) then
                do i = istartx,iendx
                do j = istarty,iendy
                    LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1)) = LP(iLayer)%HK(i,j)
                enddo
                enddo
                call MPI_SEND(LocalData,(iendx-istartx+1)*(iendy-istarty+1), &
                    MPI_DOUBLE_PRECISION,iNode,iLayer,MPI_COMM_WORLD,ierror)
            else
                do i = istartx,iendx
                do j = istarty,iendy
                    if(LP(iLayer)%Z(i,j).gt.0.0.and.LP(iLayer)%H(1,i,j)+LP(iLayer)%Z(i,j).gt.GP%MinWaterDepth) &
                    LP(iLayer)%H(1,i,j) = LP(iLayer)%H(1,i,j) + LP(iLayer)%HK(i,j)
                enddo
                enddo

            endif
        enddo
    elseif(irank.lt.LP(iLayer)%nsize) then
        istartx = MAX(1,LP(iLayer)%PartitionInfo(irank+1,1)-nRowBoundary)
        iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(irank+1,2)+nRowBoundary)
        istarty = MAX(1,LP(iLayer)%PartitionInfo(irank+1,3)-nRowBoundary)
        iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(irank+1,4)+nRowBoundary)
        call MPI_RECV(LocalData,(iendx-istartx+1)*(iendy-istarty+1), &
            MPI_DOUBLE_PRECISION,master,iLayer,MPI_COMM_WORLD,istatus,ierror)
        do i = istartx, iendx
        do j = istarty, iendy
            LP(iLayer)%HK(i,j) = LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1))
            if(LP(iLayer)%Z(i,j).gt.0.0.and.LP(iLayer)%H(1,i,j)+LP(iLayer)%Z(i,j).gt.GP%MinWaterDepth) &
                LP(iLayer)%H(1,i,j) = LP(iLayer)%H(1,i,j) + LP(iLayer)%HK(i,j)
        enddo
        enddo
    endif
    enddo
    if(irank.eq.master.and.iTimeStep.eq.0) write(*,*)
    
endif


if(GP%PurposeCalculation.eq.3.and.iTimeStep.eq.0) then
    ix = MOD(iCal-1, GP%SourceNX)+1
    iy = (iCal-1)/GP%SourceNX+1
    x0 = GP%SourceStartX+(ix-0.5)*GP%SourceDX
    y0 = GP%SourceStartY+(iy-0.5)*GP%SourceDY

    if(irank.eq.master) then
        if(GP%SourceBasisFunctionType.eq.1) then
            write(*,*) 'setting up initial water elevation of unit Gaussian shape ...'
            write(*,*)
        elseif(GP%SourceBasisFunctionType.eq.2) then
            write(*,*) 'setting up initial water elevation of approximate Sigmoid function 1/(1+exp(-x)) ...'
            write(*,*)
        endif
        fn = TRIM(ADJUSTL(GP%ResultPath))//TRIM(ADJUSTL(GP%GFParamFileName))
        if(iCal.eq.1) then
            open(23,file=TRIM(ADJUSTL(fn)),status='replace',form='formatted')
            write(23,'(a)') 'SourceBasisFunctionType: 1/Gaussian; 2/Sigmoid'
            write(23,'(a)') 'GaussianRatio, SigmoidCoefficient'
            write(23,'(a)') 'i  NX  NY  ix  iy  x0  y0  dx dy'
            write(23,'(i5)') GP%SourceBasisFunctionType
            write(23,'(2f20.10)') GP%GaussianRatio, GP%SigmoidCoefficient
        else
            open(23,file=TRIM(ADJUSTL(fn)),form='formatted',status='old',access='append')
        endif
        write(23,'(5i5,4f20.10)') iCal,GP%SourceNX,GP%SourceNY,ix,iy,x0,y0,GP%SourceDX,GP%SourceDY
        close(23)
    endif

    if(GP%SourceBasisFunctionType.eq.1) then
        lx2 = (GP%SourceDX*GP%GaussianRatio)**2
        ly2 = (GP%SourceDY*GP%GaussianRatio)**2
        do iLayer = 1,GP%NumLayers
        if(irank.lt.LP(iLayer)%nsize) then
            istartx = MAX(1,LP(iLayer)%PartitionInfo(irank+1,1)-nRowBoundary)
            iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(irank+1,2)+nRowBoundary)
            istarty = MAX(1,LP(iLayer)%PartitionInfo(irank+1,3)-nRowBoundary)
            iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(irank+1,4)+nRowBoundary)
            do i = istartx, iendx
            do j = istarty, iendy
                x = LP(iLayer)%X(i); y = LP(iLayer)%Y(j);
                h = EXP(-(x-x0)**2/lx2-(y-y0)**2/ly2)
                LP(iLayer)%H(1,i,j) = LP(iLayer)%H(1,i,j) + h
                if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit) then 
                    LP(iLayer)%H(1,i,j) = 0.0
                elseif(LP(iLayer)%Z(i,j)+LP(iLayer)%H(1,i,j).le.0.0) then
                    if(LP(iLayer)%Z(i,j).gt.0.0) then
                        LP(iLayer)%H(1,i,j) = -LP(iLayer)%Z(i,j)
                    else
                        LP(iLayer)%H(1,i,j) = 0.0
                    endif
                endif
            enddo
            enddo
        endif
        enddo

    elseif(GP%SourceBasisFunctionType.eq.2) then
        lx = GP%SourceDX*GP%SigmoidCoefficient
        ly = GP%SourceDY*GP%SigmoidCoefficient
        do iLayer = 1,GP%NumLayers
        if(irank.lt.LP(iLayer)%nsize) then
            istartx = MAX(1,LP(iLayer)%PartitionInfo(irank+1,1)-nRowBoundary)
            iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(irank+1,2)+nRowBoundary)
            istarty = MAX(1,LP(iLayer)%PartitionInfo(irank+1,3)-nRowBoundary)
            iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(irank+1,4)+nRowBoundary)
            do i = istartx, iendx
            do j = istarty, iendy
                x = LP(iLayer)%X(i); y = LP(iLayer)%Y(j);
                x_left  = x0 - 0.5*GP%SourceDX-lx
                x_right = x0 + 0.5*GP%SourceDX+lx
                y_left  = y0 - 0.5*GP%SourceDY-ly
                y_right = y0 + 0.5*GP%SourceDY+ly
                h = (1.0/(1.0+exp((-x+x_left)/lx))+1.0/(1.0+exp((x-x_right)/lx))-1.0)* &
                    (1.0/(1.0+exp((-y+y_left)/ly))+1.0/(1.0+exp((y-y_right)/ly))-1.0)
                LP(iLayer)%H(1,i,j) = LP(iLayer)%H(1,i,j) + h
                if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit) then 
                    LP(iLayer)%H(1,i,j) = 0.0
                elseif(LP(iLayer)%Z(i,j)+LP(iLayer)%H(1,i,j).le.0.0) then
                    if(LP(iLayer)%Z(i,j).gt.0.0) then
                        LP(iLayer)%H(1,i,j) = -LP(iLayer)%Z(i,j)
                    else
                        LP(iLayer)%H(1,i,j) = 0.0
                    endif
                endif
            enddo
            enddo
        endif
        enddo
    endif
endif

!/// store maximum and minimum water height in each layer ///!
do iLayer = 1,GP%NumLayers
    if(irank.lt.LP(iLayer)%nsize) then
        istartx = LP(iLayer)%PartitionInfo(irank+1,1)
        iendx   = LP(iLayer)%PartitionInfo(irank+1,2)
        istarty = LP(iLayer)%PartitionInfo(irank+1,3)
        iendy   = LP(iLayer)%PartitionInfo(irank+1,4)
        do j = istarty, iendy
        do i = istartx, iendx
            LP(iLayer)%Hmax(i,j) = MAX(LP(iLayer)%Hmax(i,j),LP(iLayer)%H(1,i,j))
            LP(iLayer)%Hmin(i,j) = MIN(LP(iLayer)%Hmin(i,j),LP(iLayer)%H(1,i,j))
        enddo
        enddo
    endif
enddo

end subroutine getInitialCondition



subroutine smoothElevation_Kajiura(GP, LP, iLayer, kajiuraRmax, LocalData, LocalDataLength)

use mpi
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4, intent(in)   ::  iLayer
real*8, intent(in)      ::  kajiuraRmax
integer*4, intent(in)   ::  LocalDataLength
real*8                  ::  LocalData(LocalDataLength)  
integer*4, parameter    ::  NInterpPoints = 601;
real*8      ::  kj_r(NInterpPoints), kj_g(NInterpPoints)
real*8      ::  smoothWeight(LP(iLayer)%NX,LP(iLayer)%NY), HS(LP(iLayer)%NX,LP(iLayer)%NY)
real*8      ::  smoothR, depth_local, dx, dy, lxy, gxy, GtermSum, GtermSum_Inv, val
real*8      ::  depth_average, ele_threshold, thresholdRatio, depth_average_part, ele_threshold_part
integer*4   ::  Ncells, Ncells_part
integer*4   ::  irank, master
integer*4   ::  ierror, istatus(MPI_STATUS_SIZE)
integer*4   ::  iPart, iNode, istartx, iendx, istarty, iendy
integer*4   ::  i, j, m, n, IR, JR, mstart, mend, nstart, nend
real*8      ::  CPUTime1, CPUTime2, t0, t1, t2, ETA
integer*4   ::  nh, nm, ns  
    
interface
    elemental function kajiura_G(r) result(output)
    implicit NONE
    real*8  ::  output
    real*8, intent(in)  ::  r
    end function kajiura_G
end interface
    
irank = GP%irank; master = GP%master
call CPU_TIME(CPUTime1)

if(irank.eq.master) then
    write(*,*)
    write(*,*) 'applying Kajiura filter to surface elevation ...'
endif

!///compute the average water depth of source area///!
if(GP%UseAverageDepth.eq.1) then
    if(irank.eq.master) then
        write(*,*) 'calculating average water depth of source area ...'
    endif
    !///estimate the threshold of deformation///!
    thresholdRatio = 0.2d0 !ratio of the threshold to the maximum elevation
    ele_threshold_part = 0.0d0
    if(irank.lt.LP(iLayer)%nsize) then
        istartx = LP(iLayer)%PartitionInfo(irank+1,1)
        iendx   = LP(iLayer)%PartitionInfo(irank+1,2)
        istarty = LP(iLayer)%PartitionInfo(irank+1,3)
        iendy   = LP(iLayer)%PartitionInfo(irank+1,4)
        do j = istarty,iendy
        do i = istartx,iendx
            if(LP(iLayer)%Z(i,j).le.0.0d0) cycle
            ele_threshold_part = MAX(ele_threshold_part, ABS(LP(iLayer)%HK(i,j)))
        enddo
        enddo
    endif
    ele_threshold_part = thresholdRatio * ele_threshold_part
    call MPI_ALLREDUCE(ele_threshold_part, ele_threshold, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)

    !///estimate the average water depth///!
    Ncells_part = 0; depth_average_part = 0.0d0
    if(irank.lt.LP(iLayer)%nsize) then
        istartx = LP(iLayer)%PartitionInfo(irank+1,1)
        iendx   = LP(iLayer)%PartitionInfo(irank+1,2)
        istarty = LP(iLayer)%PartitionInfo(irank+1,3)
        iendy   = LP(iLayer)%PartitionInfo(irank+1,4)
        do j = istarty,iendy
        do i = istartx,iendx
            if(LP(iLayer)%Z(i,j).le.0.0d0) cycle
            if(ABS(LP(iLayer)%HK(i,j)).ge.ele_threshold) then
                Ncells_part = Ncells_part + 1
                depth_average_part = depth_average_part + LP(iLayer)%Z(i,j)
            endif
        enddo
        enddo
        depth_average_part = depth_average_part/MAX(Ncells_part,1)
    endif

    call MPI_ALLREDUCE(Ncells_part, Ncells, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
    Ncells = MAX(Ncells,1)
    depth_average_part = depth_average_part*Ncells_part/Ncells
    call MPI_ALLREDUCE(depth_average_part, depth_average, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    depth_average = MAX(depth_average, GP%MinKajiuraDepth)
endif
    
!///compute Kajiura's Green function at regular points///!
do i = 1,NInterpPoints
    kj_r(i) = kajiuraRmax * (i-1.0d0)/(NInterpPoints-1.0d0)
enddo
kj_g = kajiura_G(kj_r)
    
!///convolve the Kajiura's smoothing kernel with ocean bottom deformation///!
if(irank.eq.master) then
    write(*,*) '==== progress of Kajiura filter ===='
    write(*,*) '    Percentage    ETA(hh:mm:ss)'
endif
HS = 0.0d0
call CPU_TIME(t0); t1 = t0
if(irank.lt.LP(iLayer)%nsize) then
    istartx = LP(iLayer)%PartitionInfo(irank+1,1)
    iendx   = LP(iLayer)%PartitionInfo(irank+1,2)
    istarty = LP(iLayer)%PartitionInfo(irank+1,3)
    iendy   = LP(iLayer)%PartitionInfo(irank+1,4)
    do j = istarty,iendy
        !///output the progress of Kajiura's filter///!
        if(irank.eq.master) then
            call CPU_TIME(t2)
            if(j.gt.istarty.and.(MOD(j-istarty,MAX(1,(iendy-istarty+1)/10)).eq.0.or.t2-t1.ge.30.0)) then
                ETA = (t2-t0)/(j-istarty)*(iendy-j+1)
                nh = FLOOR(ETA/3600)
                nm = FLOOR((ETA-nh*3600)/60)
                ns = INT(ETA-nh*3600-nm*60)
                write(*,'(f12.1,a1,7x,i3.2,a1,i2.2,a1,i2.2)') 1.0*(j-istarty)/(iendy-istarty+1)*100.0,'%', nh, ':', nm, ':', ns
                t1 = t2
            endif
        endif
    do i = istartx,iendx
        if(LP(iLayer)%Z(i,j).le.GP%MinKajiuraDepth) then
            HS(i,j) = HS(i,j) + LP(iLayer)%HK(i,j); cycle
        endif
        GtermSum = 0.0d0
        if(GP%UseAverageDepth.eq.1.and.LP(iLayer)%Z(i,j).gt.depth_average) then
            depth_local = depth_average
        else
            depth_local = LP(iLayer)%Z(i,j)
        endif
        smoothR = kajiuraRmax*depth_local
        if(GP%CoordinatesType.eq.1) then
            dx = LP(iLayer)%dx; dy = LP(iLayer)%dy
        elseif(GP%CoordinatesType.eq.0) then
            dx = GP%R_EARTH*COS(LP(iLayer)%Y(j)*GP%PI/180.0)*LP(iLayer)%dx*GP%PI/180.0
            dy = GP%R_EARTH*LP(iLayer)%dy*GP%PI/180.0
        endif
        IR = FLOOR(smoothR/dx); JR = FLOOR(smoothR/dy)
        mstart = MAX(i-IR,1); mend = MIN(i+IR,LP(iLayer)%NX)
        nstart = MAX(j-JR,1); nend = MIN(j+JR,LP(iLayer)%NY)
        do n = nstart,nend
        do m = mstart,mend
            dx = (m-i)*LP(iLayer)%dx; dy = (n-j)*LP(iLayer)%dy
            if(GP%CoordinatesType.eq.0) then
                dx = GP%R_EARTH*COS(LP(iLayer)%Y(j)*GP%PI/180.0)*dx*GP%PI/180.0
                dy = GP%R_EARTH*dy*GP%PI/180.0
            endif
            lxy = SQRT(dx*dx+dy*dy)/depth_local
            if(lxy.le.kajiuraRmax) then
                call interpolate1D_linear(kj_r, kj_g, NInterpPoints, NInterpPoints, lxy, gxy)
            else
                gxy = 0.0d0
            endif
            smoothWeight(m,n) = gxy; GtermSum = GtermSum + gxy
        enddo
        enddo
        GtermSum_Inv = 1.0d0/GtermSum
        do n = nstart,nend
        do m = mstart,mend
            HS(m,n) = HS(m,n) + GtermSum_Inv*smoothWeight(m,n)*LP(iLayer)%HK(i,j)
        enddo
        enddo
    enddo
    enddo
    LP(iLayer)%HK = HS
endif !if this node is used by this layer
    
!///suprimpose the smoothed elevation calculated on each node///!
do iPart = 1,LP(iLayer)%nsize
    istartx = LP(iLayer)%PartitionInfo(iPart,1)
    iendx   = LP(iLayer)%PartitionInfo(iPart,2)
    istarty = LP(iLayer)%PartitionInfo(iPart,3)
    iendy   = LP(iLayer)%PartitionInfo(iPart,4)
    if(irank.lt.LP(iLayer)%nsize.and.irank.ne.master) then
        do i = istartx,iendx
        do j = istarty,iendy
            LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1)) = LP(iLayer)%HK(i,j)
        enddo
        enddo
        call MPI_SEND(LocalData,(iendx-istartx+1)*(iendy-istarty+1), &
            MPI_DOUBLE_PRECISION,master,iPart,MPI_COMM_WORLD,ierror)
    endif
    if(irank.eq.master) then
        do iNode = 0,LP(iLayer)%nsize-1
            if(iNode.ne.master) then
                call MPI_RECV(LocalData,(iendx-istartx+1)*(iendy-istarty+1), &
                    MPI_DOUBLE_PRECISION,iNode,iPart,MPI_COMM_WORLD,istatus,ierror)
                do i = istartx,iendx
                do j = istarty,iendy
                    LP(iLayer)%HK(i,j) = LP(iLayer)%HK(i,j) + LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1))
                enddo
                enddo
            endif
        enddo
    endif
enddo !loop for all parts of the computation domain
if(irank.eq.master) then
    write(*,*) '==== Kajiura filtering finished! ===='
endif
    
call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1
    
end subroutine smoothElevation_Kajiura



elemental function kajiura_G(r) result(output)
!This function calculates Kajiura's G function using the infinite series representation.
!r is the distance normalized by water depth
implicit NONE
real*8  ::  output
real*8, intent(in)  ::  r
real*8, parameter       ::  RECURSIVE_STOP_FACTOR = 1.0e-5
integer*4, parameter    ::  NT = 1000
real*8, parameter       ::  PI_INV = 1.0d0/(atan(1.0d0)*4.0d0)
integer*4   ::  k, iter_counter
real*8      ::  g, g2
logical     ::  converged

converged = .FALSE.
!///compute the sum of the first NT terms in the series///!
g = 0.0d0
do k = 0,NT
    g = g + ((-1.0d0)**k)*(2.0d0*k+1.0d0)/((2.0d0*k+1.0d0)**2+r*r)**(1.5d0)
enddo
g = g * PI_INV
iter_counter = NT

!///calculate the trunction error///!
do while(.not.converged)
    g2 = 0.0d0
    do k = iter_counter+1,iter_counter+NT
        g2 = g2 + ((-1.0d0)**k)*(2.0d0*k+1.0d0)/((2.0d0*k+1.0d0)**2+r*r)**(1.5d0)
    enddo
    g2 = g2 * PI_INV
    g = g + g2
    iter_counter = iter_counter + NT
    if(ABS(g2).le.ABS(g)*RECURSIVE_STOP_FACTOR) converged = .TRUE.
enddo
output = g

end function kajiura_G



subroutine interpData(GP, LP, iLayer, iFlag, x, y, val)
!! iFlag = -1,-2,-3: interpolate H,M,N(1,i,j)
!! iFlag = 1,2,3: interpolate H,M,N(2,i,j)
!! iFlag = 0: interpolate Z(i,j)
!! iFlag = 4: interpolate Q(i,j)
!! iFlag = 5: interpolate Brkv(i,j)
!! iFlag = 6: interpolate BrkAge(i,j)

use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(LayerParameters)    ::  LP(100)
integer*4                ::  iLayer, iFlag, iTime
real*8                   ::  x, y, val
integer*4                ::  ii, jj, ii2, jj2
integer*4                ::  i1, i2, i3, i4
real*8                   ::  x1, y1, x2, y2, invx, invy
real*8                   ::  z1, z2, z3, z4

if(iFlag.lt.0) then
    iTime = 1
else
    iTime = 2
endif
if(iFlag.eq.0) then
    ii  = FLOOR((x-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
    jj  = FLOOR((y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
    ii2 = MIN(ii+1, LP(iLayer)%NX)
    jj2 = MIN(jj+1, LP(iLayer)%NY)
    z1  = LP(iLayer)%Z(ii,jj)
    z2  = LP(iLayer)%Z(ii2,jj)
    z3  = LP(iLayer)%Z(ii,jj2)
    z4  = LP(iLayer)%Z(ii2,jj2)
    x1  = LP(iLayer)%X(ii);  y1 = LP(iLayer)%Y(jj)
    x2  = LP(iLayer)%X(ii2); y2 = LP(iLayer)%Y(jj2)
elseif(iFlag.eq.-1.or.iFlag.eq.1) then
    ii  = FLOOR((x-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
    jj  = FLOOR((y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
    ii2 = MIN(ii+1, LP(iLayer)%NX)
    jj2 = MIN(jj+1, LP(iLayer)%NY)
    z1  = LP(iLayer)%H(iTime,ii,jj)
    z2  = LP(iLayer)%H(iTime,ii2,jj)
    z3  = LP(iLayer)%H(iTime,ii,jj2)
    z4  = LP(iLayer)%H(iTime,ii2,jj2)
    x1  = LP(iLayer)%X(ii);  y1 = LP(iLayer)%Y(jj)
    x2  = LP(iLayer)%X(ii2); y2 = LP(iLayer)%Y(jj2)
elseif(iFlag.eq.-2.or.iFlag.eq.2) then
    ii  = FLOOR((x-LP(iLayer)%X(1)-LP(iLayer)%dx*0.5)/LP(iLayer)%dx)+1
    jj  = FLOOR((y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
    ii2 = MIN(ii+1, LP(iLayer)%NX-1)
    jj2 = MIN(jj+1, LP(iLayer)%NY)
    if(ii.eq.0) then
        ii = 1; ii2 = 1
    endif
    z1 = LP(iLayer)%M(iTime,ii,jj)
    z2 = LP(iLayer)%M(iTime,ii2,jj)
    z3 = LP(iLayer)%M(iTime,ii,jj2)
    z4 = LP(iLayer)%M(iTime,ii2,jj2)
    x1 = LP(iLayer)%X(ii)+LP(iLayer)%dx*0.5;  y1 = LP(iLayer)%Y(jj)
    x2 = LP(iLayer)%X(ii2)+LP(iLayer)%dx*0.5; y2 = LP(iLayer)%Y(jj2)           
elseif(iFlag.eq.-3.or.iFlag.eq.3) then
    ii  = FLOOR((x-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
    jj  = FLOOR((y-LP(iLayer)%Y(1)-LP(iLayer)%dy*0.5)/LP(iLayer)%dy)+1
    ii2 = MIN(ii+1, LP(iLayer)%NX)
    jj2 = MIN(jj+1, LP(iLayer)%NY-1)
    if(jj.eq.0) then
        jj = 1; jj2 = 1
    endif
    z1 = LP(iLayer)%N(iTime,ii,jj)
    z2 = LP(iLayer)%N(iTime,ii2,jj)
    z3 = LP(iLayer)%N(iTime,ii,jj2)
    z4 = LP(iLayer)%N(iTime,ii2,jj2)
    x1 = LP(iLayer)%X(ii);  y1 = LP(iLayer)%Y(jj)+0.5*LP(iLayer)%dy
    x2 = LP(iLayer)%X(ii2); y2 = LP(iLayer)%Y(jj2)+0.5*LP(iLayer)%dy
else
    ii  = FLOOR((x-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
    jj  = FLOOR((y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
    ii2 = MIN(ii+1, LP(iLayer)%NX)
    jj2 = MIN(jj+1, LP(iLayer)%NY)
    x1 = LP(iLayer)%X(ii);  y1 = LP(iLayer)%Y(jj)
    x2 = LP(iLayer)%X(ii2); y2 = LP(iLayer)%Y(jj2)
    if(iFlag.eq.4) then
        z1 = LP(iLayer)%Q(ii,jj);  z2 = LP(iLayer)%Q(ii2,jj)
        z3 = LP(iLayer)%Q(ii,jj2); z4 = LP(iLayer)%Q(ii2,jj2)
    elseif(iFlag.eq.5) then
        z1 = LP(iLayer)%Brkv(ii,jj);  z2 = LP(iLayer)%Brkv(ii2,jj)
        z3 = LP(iLayer)%Brkv(ii,jj2); z4 = LP(iLayer)%Brkv(ii2,jj2)
    elseif(iFlag.eq.6) then
        z1 = LP(iLayer)%BrkAge(ii,jj);  z2 = LP(iLayer)%BrkAge(ii2,jj)
        z3 = LP(iLayer)%BrkAge(ii,jj2); z4 = LP(iLayer)%BrkAge(ii2,jj2)
    endif
endif
            
i1 = 1; i2 = 1; i3 = 1; i4 = 1
if(iFlag.eq.1.or.iFlag.eq.-1.or.(iFlag.ge.4.and.iFlag.le.6)) then
    if(LP(iLayer)%Z(ii,jj)+LP(iLayer)%H(iTime,ii,jj).le.0.0) i1 = 0
    if(LP(iLayer)%Z(ii2,jj)+LP(iLayer)%H(iTime,ii2,jj).le.0.0) i2 = 0
    if(LP(iLayer)%Z(ii,jj2)+LP(iLayer)%H(iTime,ii,jj2).le.0.0) i3 = 0
    if(LP(iLayer)%Z(ii2,jj2)+LP(iLayer)%H(iTime,ii2,jj2).le.0.0) i4 = 0
    if(i1.eq.0.and.i2.eq.0.and.i3.eq.0.and.i4.eq.0) then
        i1 = 1; z1 = 0.0d0
    endif
elseif(iFlag.ne.0) then
    if(LP(iLayer)%Z(ii,jj).le.-GP%PermanentDryLimit) i1 = 0
    if(LP(iLayer)%Z(ii2,jj).le.-GP%PermanentDryLimit) i2 = 0
    if(LP(iLayer)%Z(ii,jj2).le.-GP%PermanentDryLimit) i3 = 0
    if(LP(iLayer)%Z(ii2,jj2).le.-GP%PermanentDryLimit) i4 = 0
    if(i1.eq.0.and.i2.eq.0.and.i3.eq.0.and.i4.eq.0) then
        i1 = 1; z1 = 0.0d0
    endif
endif
if(i1.eq.1.and.i2.eq.0.and.i3.eq.0.and.i4.eq.0) then
    ii2 = ii; jj2 = jj
elseif(i1.eq.0.and.i2.eq.1.and.i3.eq.0.and.i4.eq.0) then
    ii = ii2; jj2 = jj; z1 = z2
elseif(i1.eq.0.and.i2.eq.0.and.i3.eq.1.and.i4.eq.0) then
    ii2 = ii; jj = jj2; z1 = z3
elseif(i1.eq.0.and.i2.eq.0.and.i3.eq.0.and.i4.eq.1) then
    ii = ii2; jj = jj2; z1 = z4
elseif(i1.eq.1.and.i2.eq.1.and.i3.eq.0.and.i4.eq.0) then
    jj2 = jj;
elseif(i1.eq.1.and.i2.eq.0.and.i3.eq.1.and.i4.eq.0) then
    ii2 = ii;
elseif(i1.eq.0.and.i2.eq.1.and.i3.eq.0.and.i4.eq.1) then
    ii = ii2; z1 = z2; z3 = z4
elseif(i1.eq.0.and.i2.eq.0.and.i3.eq.1.and.i4.eq.1) then
    jj = jj2; z1 = z3; z2 = z4
elseif(i1.eq.1.and.i2.eq.0.and.i3.eq.0.and.i4.eq.1) then
    jj2 = jj; z2 = z4
elseif(i1.eq.0.and.i2.eq.1.and.i3.eq.1.and.i4.eq.0) then
    jj2 = jj; z1 = z3
elseif((i1.eq.0.or.i2.eq.0).and.i3.eq.1.and.i4.eq.1) then
    jj = jj2; z1 = z3; z2 = z4
elseif(i1.eq.1.and.i2.eq.1.and.(i3.eq.0.or.i4.eq.0)) then
    jj2 = jj
endif

if(ii.ne.ii2) invx = 1.0/(x2-x1)
if(jj.ne.jj2) invy = 1.0/(y2-y1)
if(ii.eq.ii2.or.jj.eq.jj2) then
    if(ii.eq.ii2.and.jj.eq.jj2) then
        val = z1
    elseif(ii.eq.ii2) then
        val = -z1*(y-y2)*invy+z3*(y-y1)*invy
    else
        val = -z1*(x-x2)*invx+z2*(x-x1)*invx
    endif
else
    val = z1*(x2-x)*(y2-y)+z2*(x-x1)*(y2-y)+z3*(x2-x)*(y-y1)+z4*(x-x1)*(y-y1)
    val = val*invx*invy
endif

end subroutine interpData



subroutine interpDataByAverage(GP, LP, iLayer, iFlag, x, y, val)
!! iFlag = -1,-2,-3: interpolate H,M,N(1,i,j)
!! iFlag = 1,2,3: interpolate H,M,N(2,i,j)

use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(LayerParameters)    ::  LP(100)
integer*4                ::  iLayer, iFlag, iTime
real*8                   ::  x, y, val
integer*4                ::  ii, jj, ii2, jj2, ii0, jj0
real*8                   ::  z1, z2, z3, z4, z5, z6, z7, z8, z9, inv

if(iFlag.lt.0) then
    iTime = 1
else
    iTime = 2
endif
if(iFlag.eq.-1.or.iFlag.eq.1) then
    ii = NINT((x-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
    jj = NINT((y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
    ii0 = MAX(ii-1, 1); ii2 = MIN(ii+1, LP(iLayer)%NX)
    jj0 = MAX(jj-1, 1); jj2 = MIN(jj+1, LP(iLayer)%NY)
    z1 = LP(iLayer)%H(iTime,ii0,jj0)
    z2 = LP(iLayer)%H(iTime,ii0,jj)
    z3 = LP(iLayer)%H(iTime,ii0,jj2)
    z4 = LP(iLayer)%H(iTime,ii,jj0)
    z5 = LP(iLayer)%H(iTime,ii,jj)
    z6 = LP(iLayer)%H(iTime,ii,jj2)
    z7 = LP(iLayer)%H(iTime,ii2,jj0)
    z8 = LP(iLayer)%H(iTime,ii2,jj)
    z9 = LP(iLayer)%H(iTime,ii2,jj2)
elseif(iFlag.eq.-2.or.iFlag.eq.2) then
    ii = NINT((x-LP(iLayer)%X(1)-LP(iLayer)%dx*0.5)/LP(iLayer)%dx)+1
    jj = NINT((y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
    ii0 = MAX(ii-1, 1); ii2 = MIN(ii+1, LP(iLayer)%NX-1)
    jj0 = MAX(jj-1, 1); jj2 = MIN(jj+1, LP(iLayer)%NY)
    if(ii.eq.0) then
        ii0 = 1; ii = 1; ii2 = 1
    endif
    z1 = LP(iLayer)%M(iTime,ii0,jj0)
    z2 = LP(iLayer)%M(iTime,ii0,jj)
    z3 = LP(iLayer)%M(iTime,ii0,jj2)
    z4 = LP(iLayer)%M(iTime,ii,jj0)
    z5 = LP(iLayer)%M(iTime,ii,jj)
    z6 = LP(iLayer)%M(iTime,ii,jj2)
    z7 = LP(iLayer)%M(iTime,ii2,jj0)
    z8 = LP(iLayer)%M(iTime,ii2,jj)
    z9 = LP(iLayer)%M(iTime,ii2,jj2)
elseif(iFlag.eq.-3.or.iFlag.eq.3) then
    ii = FLOOR((x-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
    jj = FLOOR((y-LP(iLayer)%Y(1)-LP(iLayer)%dy*0.5)/LP(iLayer)%dy)+1
    ii0 = MAX(ii-1, 1); ii2 = MIN(ii+1, LP(iLayer)%NX)
    jj0 = MAX(jj-1, 1); jj2 = MIN(jj+1, LP(iLayer)%NY-1)
    if(jj.eq.0) then
        jj0 = 1; jj = 1; jj2 = 1
    endif
    z1 = LP(iLayer)%N(iTime,ii0,jj0)
    z2 = LP(iLayer)%N(iTime,ii0,jj)
    z3 = LP(iLayer)%N(iTime,ii0,jj2)
    z4 = LP(iLayer)%N(iTime,ii,jj0)
    z5 = LP(iLayer)%N(iTime,ii,jj)
    z6 = LP(iLayer)%N(iTime,ii,jj2)
    z7 = LP(iLayer)%N(iTime,ii2,jj0)
    z8 = LP(iLayer)%N(iTime,ii2,jj)
    z9 = LP(iLayer)%N(iTime,ii2,jj2)
endif
inv = 9.0
if(LP(iLayer)%Z(ii0,jj0).le.-GP%PermanentDryLimit) inv = inv-1.0
if(LP(iLayer)%Z(ii0,jj).le.-GP%PermanentDryLimit) inv = inv-1.0
if(LP(iLayer)%Z(ii0,jj2).le.-GP%PermanentDryLimit) inv = inv-1.0
if(LP(iLayer)%Z(ii,jj0).le.-GP%PermanentDryLimit) inv = inv-1.0
if(LP(iLayer)%Z(ii,jj).le.-GP%PermanentDryLimit) inv = inv-1.0
if(LP(iLayer)%Z(ii,jj2).le.-GP%PermanentDryLimit) inv = inv-1.0
if(LP(iLayer)%Z(ii2,jj0).le.-GP%PermanentDryLimit) inv = inv-1.0
if(LP(iLayer)%Z(ii2,jj).le.-GP%PermanentDryLimit) inv = inv-1.0
if(LP(iLayer)%Z(ii2,jj2).le.-GP%PermanentDryLimit) inv = inv-1.0
if(inv.eq.0) then
    val = 0.0
else
    val = (z1+z2+z3+z4+z5+z6+z7+z8+z9)/inv
endif

end subroutine interpDataByAverage



subroutine interpolate1D_linear(t,f,mlen,n,tq,fq)
!This subroutine finds the value of the function f=f(t) at the query point tq
!t: points at which the function is defined; f: function values at t
!tq: query point; fq: function value at tq
!mlen: maximum length of t and f; n: length of effective data in t and f
implicit NONE
real*8, intent(in)      ::  t(mlen), f(mlen), tq
real*8, intent(out)     ::  fq
integer*4, intent(in)   ::  mlen, n
integer*4   ::  j, jleft, jright
        
if(tq.le.t(1)) then
    fq = f(1)
elseif(tq.ge.t(n)) then
    fq = f(n)
else
    j = n/2; jleft = 1; jright = n;
    do
        if((tq.ge.t(j)).and.(tq.le.t(j+1))) exit
        if(tq.lt.t(j)) then
            jright = j; j = (j+jleft)/2
        else
            jleft = j; j = (j+jright)/2
        endif
    enddo
    fq = (f(j+1)-f(j))/(t(j+1)-t(j))*(tq-t(j))+f(j)
endif
    
end subroutine interpolate1D_linear



subroutine interpolate2D_linear(V,X,Y,xlen,ylen,xq,yq,val)
!This subroutine finds the value of the function V=V(x,y) at the query point (xq,yq)
!(X,Y): points at which the function is defined; V: function values at (X,Y)
!(xq,yq): query point; val: function value at (xq,yq)
!xlen, ylen: data length of X and Y
implicit NONE
real*8, intent(in)      ::  V(xlen,ylen), X(xlen), Y(ylen), xq, yq
real*8, intent(out)     ::  val
integer*4, intent(in)   ::  xlen, ylen
integer*4   ::  i, j, ileft, iright, jleft, jright
integer*4   ::  ii1, ii2, jj1, jj2
real*8      ::  x1, x2, y1, y2, invx, invy
real*8      ::  v1, v2, v3, v4 
     
if(xq.le.X(1)) then
    ii1 = 1; ii2 = ii1
elseif(xq.ge.X(xlen)) then
    ii1 = xlen; ii2 = ii1
else
    i = xlen/2; ileft = 1; iright = xlen
    do
        if((xq.ge.X(i)).and.xq.le.X(i+1)) exit
        if(xq.lt.X(i)) then
            iright = i; i = (ileft+iright)/2
        else
            ileft = i; i = (ileft+iright)/2
        endif
    enddo
    ii1 = i; ii2 = i+1
endif
    
if(yq.le.Y(1)) then
    jj1 = 1; jj2 = jj1
elseif(yq.ge.Y(ylen)) then
    jj1 = ylen; jj2 = jj1
else
    j = ylen/2; jleft = 1; jright = ylen
    do
        if((yq.ge.Y(j)).and.yq.le.Y(j+1)) exit
        if(yq.lt.Y(j)) then
            jright = j; j = (jleft+jright)/2
        else
            jleft = j; j = (jleft+jright)/2
        endif
    enddo
    jj1 = j; jj2 = j+1
endif
        
x1 = X(ii1); x2 = X(ii2)
y1 = Y(jj1); y2 = Y(jj2)
v1 = V(ii1, jj1); v2 = V(ii2,jj1)
v3 = V(ii1,jj2); v4 = V(ii2,jj2)
if(ii1.ne.ii2) invx = 1.0d0/(x2-x1)
if(jj1.ne.jj2) invy = 1.0d0/(y2-y1)
if(ii1.eq.ii2.or.jj1.eq.jj2) then
    if(ii1.eq.ii2.and.jj1.eq.jj2) then
        val = v1
    elseif(ii1.eq.ii2) then
        val = -v1*(yq-y2)*invy + v3*(yq-y1)*invy
    else
        val = -v1*(xq-x2)*invx + v2*(xq-x1)*invx
    endif
else
    val = v1*(x2-xq)*(y2-yq)+v2*(xq-x1)*(y2-yq)+v3*(x2-xq)*(yq-y1)+v4*(xq-x1)*(yq-y1)
    val = val*invx*invy
endif
        
end subroutine interpolate2D_linear



subroutine screenOutput(GP, LP, iCal, iTimeStep)

use mpi
use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(LayerParameters)    ::  LP(100)
integer*4  ::  iCal, iTimeStep
real*8     ::  CPUTime
real*8     ::  ETA
integer*4  ::  nh, nm, ns, i, j, iNode
integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
real*8     ::  LocalData(999*10), tt, calt, comt, writ

if(GP%irank.eq.GP%master) then
if(iTimeStep.eq.-1) then
    write(*,*)
    if(GP%PurposeCalculation.eq.1) then
        write(*,*) '==================== Start Computing : Regular Simulation ===================='
    elseif(GP%PurposeCalculation.eq.2) then
        write(*,'(a,i3,a,i5,a)') '=========== Start Computing : Green''s Functions (Fault ', &
            iCal,' of ',GP%nCalculations,') ==========='
    elseif(GP%PurposeCalculation.eq.3) then
        write(*,'(a,i3,a,i5,a)') '====== Start Computing : Green''s Functions (Initial Height ', &
            iCal,' of ',GP%nCalculations,') ======'
    endif
    write(*,*)
elseif(iTimeStep.gt.0.and.iTimeStep.le.GP%TotalTimeSteps) then
    call CPU_TIME(CPUTime)
    if((iTimeStep.eq.10).or. &
        (iTimeStep.gt.10.and.CPUTime-GP%CPUTimeInitial-GP%CPUTime(GP%master+1,1).ge.120.0).or.&
        (iTimeStep.gt.10.and.MOD(iTimeStep,MAX(1,GP%TotalTimeSteps/20)).eq.0)) then
        ETA = (GP%TotalTimeSteps-iTimeStep)*(CPUTime-GP%CPUTimeInitialCalculation)/iTimeStep
        nh = FLOOR(ETA/3600)
        nm = FLOOR((ETA-nh*3600)/60)
        ns = INT(ETA-nh*3600-nm*60)
        write(*,*)
        write(*,*) '===================== CALCULATION OUTPUT ========================'
        write(*,*) '    TimeStep       Seconds       Percentage       ETA(hh:mm:ss)'
        if(iTimeStep.eq.10) then
            write(*,'(i11,f16.1,9x,f4.1,a2,i13.2,a1,i2.2,a1,i2.2)') &
                0, 0.0d0, 0.0d0,'%',nh,':',nm,':',ns
        else
            write(*,'(i11,f16.1,8x,f5.1,a2,i13.2,a1,i2.2,a1,i2.2)') &
                iTimeStep,iTimeStep*GP%dt,iTimeStep*1.0/GP%TotalTimeSteps*100,'%',&
                nh,':',nm,':',ns
        endif
        GP%CPUTime(GP%master+1,1) = CPUTime-GP%CPUTimeInitial
    endif
endif
if(iTimeStep.eq.GP%TotalTimeSteps) then
    write(*,*) '==================================================='
endif
endif

if(iCal.eq.GP%nCalculations+100.and.iTimeStep.eq.GP%TotalTimeSteps+100) then
    if(GP%irank.ne.GP%master) then
        do j = 1,5
            LocalData(j) = GP%CPUTime(GP%irank+1,j)
        enddo
        call MPI_SEND(LocalData,5,MPI_DOUBLE_PRECISION, &
            GP%master,2015,MPI_COMM_WORLD,ierror)
    else
        do iNode=0, GP%nsizeTotal-1
        if(iNode.ne.GP%master) then
            call MPI_RECV(LocalData,5,MPI_DOUBLE_PRECISION, &
                iNode,2015,MPI_COMM_WORLD,istatus,ierror)
            do j = 1,5
                GP%CPUTime(iNode+1,j) = LocalData(j)
            enddo
        endif
        enddo
        tt   = MAXVAL(GP%CPUTime(1:GP%nsizeTotal,1)); calt = MAXVAL(GP%CPUTime(1:GP%nsizeTotal,3));
        comt = MAXVAL(GP%CPUTime(1:GP%nsizeTotal,4)); writ = MAXVAL(GP%CPUTime(1:GP%nsizeTotal,5));

        write(*,*)
        write(*,*) '========================================================================'
        write(*,*) 'pcomcot calculation complete.'
        write(*,*)
        write(*,'(a)')              ' Max Time On Node : '
        nh = FLOOR(tt/3600); nm = FLOOR((tt-nh*3600)/60); ns = MOD(NINT(tt),60)
        write(*,'(a,i2.2,a1,i2.2,a1,i2.2,f10.1,a)') &
            ' Total Time       : ',nh,':',nm,':',ns,100.0d0,' %'
        nh = FLOOR(calt/3600); nm = FLOOR((calt-nh*3600)/60); ns = MOD(NINT(calt),60)
        write(*,'(a,i2.2,a1,i2.2,a1,i2.2,f10.1,a)') &
            ' Calculation      : ',nh,':',nm,':',ns,calt/tt*100,' %'
        nh = FLOOR(comt/3600); nm = FLOOR((comt-nh*3600)/60); ns = MOD(NINT(comt),60)    
        write(*,'(a,i2.2,a1,i2.2,a1,i2.2,f10.1,a)') &
            ' Communication    : ',nh,':',nm,':',ns,comt/tt*100,' %'
        nh = FLOOR(writ/3600); nm = FLOOR((writ-nh*3600)/60); ns = MOD(NINT(writ),60)
        write(*,'(a,i2.2,a1,i2.2,a1,i2.2,f10.1,a)') &
            ' Write to File    : ',nh,':',nm,':',ns,writ/tt*100,' %'
        write(*,*)
    endif
endif

end subroutine screenOutput



subroutine sponge(GP, LP, iLayer)
    
use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(LayerParameters)    ::  LP(100)
integer*4  ::  iLayer
integer*4  ::  irank, nstartx, nendx, nstarty, nendy
integer*4  ::  i, j, Iwidth, Jwidth, iin, jin
real*8     ::  Cs, Csx, Csy

irank = GP%irank

if(irank.lt.LP(iLayer)%nsize) then
    nstartx = MAX(1,LP(iLayer)%PartitionInfo(irank+1,1)-GP%nRowBoundary)
    nendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(irank+1,2)+GP%nRowBoundary)
    nstarty = MAX(1,LP(iLayer)%PartitionInfo(irank+1,3)-GP%nRowBoundary)
    nendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(irank+1,4)+GP%nRowBoundary)

    if(GP%CoordinatesType.eq.1) then
        Iwidth = MAX(NINT(GP%SpongeWidthX/LP(iLayer)%dx), 0)
        Jwidth = MAX(NINT(GP%SpongeWidthY/LP(iLayer)%dy), 0)
    elseif(GP%CoordinatesType.eq.0) then
        Iwidth = NINT(GP%SpongeWidthX/GP%R_EARTH/LP(iLayer)%dx/GP%PI*180.0 &
            /COS(0.5*(LP(iLayer)%ymin+LP(iLayer)%ymax)*GP%PI/180.0))
        Jwidth = NINT(GP%SpongeWidthY/GP%R_EARTH/LP(iLayer)%dy/GP%PI*180.0)
        Iwidth = MAX(Iwidth, 0); Jwidth = MAX(Jwidth, 0)
    endif
    do j = nstarty, nendy
    do i = nstartx, nendx
        iin = MIN(i-1, LP(iLayer)%NX-i); jin = MIN(j-1, LP(iLayer)%NY-j)
        if(iin.le.Iwidth-1.or.jin.le.Jwidth-1) then !if this cell is in the sponge boundary
            if(iin.le.Iwidth-1) then
                Csx = GP%SpongeDampingR**iin; Csx = GP%SpongeDampingA**Csx
            else
                Csx = 1.0d0
            endif
            if(jin.le.Jwidth-1) then
                Csy = GP%SpongeDampingR**jin; Csy = GP%SpongeDampingA**Csy
            else
                Csy = 1.0d0
            endif
            Cs = MAX(Csx, Csy)
            LP(iLayer)%H(2,i,j) = LP(iLayer)%H(2,i,j)/Cs
            if(i.ne.LP(iLayer)%NX) LP(iLayer)%M(2,i,j) = LP(iLayer)%M(2,i,j)/Cs
            if(j.ne.LP(iLayer)%NY) LP(iLayer)%N(2,i,j) = LP(iLayer)%N(2,i,j)/Cs
        endif
    enddo
    enddo
endif !if: this node is used by this layer

    
end subroutine sponge



subroutine updateComputedResults(GP, LP, iLayer)

use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(LayerParameters)    ::  LP(100)
integer*4  ::  iLayer
integer*4  ::  irank, nstartx, nendx, nstarty, nendy
integer*4  ::  istartx, iendx, istarty, iendy
integer*4  ::  i, j

irank = GP%irank
if(irank.lt.LP(iLayer)%nsize) then
    istartx = MAX(1,LP(iLayer)%PartitionInfo(irank+1,1)-GP%nRowBoundary)
    iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(irank+1,2)+GP%nRowBoundary)
    istarty = MAX(1,LP(iLayer)%PartitionInfo(irank+1,3)-GP%nRowBoundary)
    iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(irank+1,4)+GP%nRowBoundary)

    do j = istarty, iendy
    do i = istartx, iendx
        LP(iLayer)%H(1,i,j) = LP(iLayer)%H(2,i,j)
        if(i.ne.LP(iLayer)%NX) LP(iLayer)%M(1,i,j) = LP(iLayer)%M(2,i,j)
        if(j.ne.LP(iLayer)%NY) LP(iLayer)%N(1,i,j) = LP(iLayer)%N(2,i,j)
    enddo
    enddo

    if(LP(iLayer)%Dispersion.eq.1) then
        nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
        nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
        nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
        nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
        do j = nstarty, nendy
        do i = nstartx, nendx
            LP(iLayer)%WS(1,i,j) = LP(iLayer)%WS(2,i,j)
            LP(iLayer)%WB(1,i,j) = LP(iLayer)%WB(2,i,j)
        enddo
        enddo
    endif

endif !if: irank is used

end subroutine updateComputedResults



subroutine calculateStationData(GP, LP, SP, iCal, iTimeStep)

use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(LayerParameters)    ::  LP(100)
type(StationParameters)  ::  SP(999)
integer*4  ::  iCal, iTimeStep, irank
integer*4  ::  iSta, iFlag
real*8     ::  val

irank = GP%irank
do iSta = 1,GP%NumStations
    if(irank.eq.SP(iSta)%nNode) then
        iFlag = -1
        call interpData(GP, LP, SP(iSta)%nLayer, iFlag, SP(iSta)%X, SP(iSta)%Y, val)
        SP(iSta)%H(iTimeStep+1) = val
        if(GP%SaveFlux.eq.1) then
            iFlag = -2
            call interpData(GP, LP, SP(iSta)%nLayer, iFlag, SP(iSta)%X, SP(iSta)%Y, val)
            SP(iSta)%M(iTimeStep+1) = val

            iFlag = -3
            call interpData(GP, LP, SP(iSta)%nLayer, iFlag, SP(iSta)%X, SP(iSta)%Y, val)
            SP(iSta)%N(iTimeStep+1) = val
        endif
        if(GP%SaveDynamicPressure.eq.1.and.LP(SP(iSta)%nLayer)%Dispersion.eq.1) then
            iFlag = 4
            call interpData(GP, LP, SP(iSta)%nLayer, iFlag, SP(iSta)%X, SP(iSta)%Y, val)
            SP(iSta)%Q(iTimeStep+1) = val
        endif
    endif
enddo

end subroutine calculateStationData



subroutine saveStationData(GP, LP, SP, iCal)

use mpi
use VariableDefination
implicit NONE
type(GlobalParameters)    ::  GP
type(LayerParameters)     ::  LP(100)
type(StationParameters)   ::  SP(999)
integer*4        ::  iCal
integer*4        ::  irank, nsize, master
integer*4        ::  ierror, istatus(MPI_STATUS_SIZE)
integer*4        ::  iSta, i, j
character(999)   ::  s

irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master

do iSta = 1,GP%NumStations
    if(irank.eq.SP(iSta)%nNode.and.irank.ne.master) then
        call MPI_SEND(SP(iSta)%H, GP%TotalTimeSteps+1, &
            MPI_DOUBLE_PRECISION, master, 2015, MPI_COMM_WORLD, ierror)
        if(GP%SaveFlux.eq.1) then
            call MPI_SEND(SP(iSta)%M, GP%TotalTimeSteps+1, &
                MPI_DOUBLE_PRECISION, master, 2015, MPI_COMM_WORLD, ierror)
            call MPI_SEND(SP(iSta)%N, GP%TotalTimeSteps+1, &
                MPI_DOUBLE_PRECISION, master, 2015, MPI_COMM_WORLD, ierror)
        endif
        if(GP%SaveDynamicPressure.eq.1.and.LP(SP(iSta)%nLayer)%Dispersion.eq.1) then
            call MPI_SEND(SP(iSta)%Q, GP%TotalTimeSteps+1, &
            MPI_DOUBLE_PRECISION, master, 2015, MPI_COMM_WORLD, ierror) 
        endif
    endif
enddo
do iSta = 1,GP%NumStations
    if(irank.eq.master.and.irank.ne.SP(iSta)%nNode) then
        call MPI_RECV(SP(iSta)%H, GP%TotalTimeSteps+1, &
            MPI_DOUBLE_PRECISION, SP(iSta)%nNode, 2015, MPI_COMM_WORLD, istatus, ierror)
        if(GP%SaveFlux.eq.1) then 
            call MPI_RECV(SP(iSta)%M, GP%TotalTimeSteps+1, &
                MPI_DOUBLE_PRECISION, SP(iSta)%nNode, 2015, MPI_COMM_WORLD, istatus, ierror)
            call MPI_RECV(SP(iSta)%N, GP%TotalTimeSteps+1, &
                MPI_DOUBLE_PRECISION, SP(iSta)%nNode, 2015, MPI_COMM_WORLD, istatus, ierror)
        endif
        if(GP%SaveDynamicPressure.eq.1.and.LP(SP(iSta)%nLayer)%Dispersion.eq.1) then
            call MPI_RECV(SP(iSta)%Q, GP%TotalTimeSteps+1, &
                MPI_DOUBLE_PRECISION, SP(iSta)%nNode, 2015, MPI_COMM_WORLD, istatus, ierror)
        endif
    endif
enddo

if(irank.eq.master) then
do i = 1,GP%NumStations
    write(s,'(a,i4.4,a4)') TRIM(ADJUSTL(GP%ResultPath))//'Station',i,'.dat'
    open(1000+i,file=s,status='old',form='unformatted',access='append')
    if(iCal.eq.1) then
        write(1000+i) (GP%t(j), j=1,GP%TotalTimeSteps+1)
    endif
    write(1000+i) (SP(i)%H(j), j=1,GP%TotalTimeSteps+1)
    close(1000+i)

    if(GP%SaveFlux.eq.1) then

    write(s,'(a,i4.4,a6)') TRIM(ADJUSTL(GP%ResultPath))//'Station',i,'_M.dat'
    open(1000+i,file=s,status ='old',form='unformatted',access='append')
    if(iCal.eq.1) then
        write(1000+i) (GP%t(j), j=1,GP%TotalTimeSteps+1)
    endif
    write(1000+i) (SP(i)%M(j), j=1,GP%TotalTimeSteps+1)
    close(1000+i)

    write(s,'(a,i4.4,a6)') TRIM(ADJUSTL(GP%ResultPath))//'Station',i,'_N.dat'
    open(1000+i,file=s,status='old',form='unformatted',access='append')
    if(iCal.eq.1) then
        write(1000+i) (GP%t(j), j=1,GP%TotalTimeSteps+1)
    endif
    write(1000+i) (SP(i)%N(j), j=1,GP%TotalTimeSteps+1)
    close(1000+i)

    endif

    if(GP%SaveDynamicPressure.eq.1.and.LP(SP(i)%nLayer)%Dispersion.eq.1) then

    write(s,'(a,i4.4,a6)') TRIM(ADJUSTL(GP%ResultPath))//'Station',i,'_Q.dat'
    open(1000+i,file=s,status ='old',form='unformatted',access='append')
    if(iCal.eq.1) then
        write(1000+i) (GP%t(j), j=1,GP%TotalTimeSteps+1)
    endif
    write(1000+i) (SP(i)%Q(j), j=1,GP%TotalTimeSteps+1)
    close(1000+i)

    endif

enddo
endif

end subroutine saveStationData



subroutine saveSnapshot(GP, LP, iCal, iTimeStep, LocalData, LocalDataLength)
    
use mpi
use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(LayerParameters)    ::  LP(100)
integer*4  ::  iCal, iTimeStep
integer*4  ::  LocalDataLength
real*8     ::  LocalData(LocalDataLength)
integer*4  ::  irank, nsize, master
integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
integer*4  ::  iLayer, iNode, i, j, iFlag
integer*4  ::  nstartx, nendx, nstarty, nendy
character(999)  ::  s
real*8     ::  CPUTime1, CPUTime2

irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
call CPU_TIME(CPUTime1)

do iLayer = 1,GP%NumLayers
    if(irank.lt.LP(iLayer)%nsize.and.irank.ne.master) then
        nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
        nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
        nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
        nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
        if(iTimeStep.eq.0) then
            do j=nstarty, nendy
            do i=nstartx, nendx
                LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%H(1,i,j)
            enddo
            enddo
        else
            do j=nstarty, nendy
            do i=nstartx, nendx
                LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%H(2,i,j)
            enddo
            enddo
        endif
        call MPI_SEND(LocalData, (nendx-nstartx+1)*(nendy-nstarty+1), &
            MPI_DOUBLE_PRECISION,master,2015,MPI_COMM_WORLD,ierror)
    endif
    if(irank.eq.master) then
        do iNode = 0, LP(iLayer)%nsize-1
            if(iNode.ne.master) then
                nstartx = LP(iLayer)%PartitionInfo(iNode+1,1)
                nendx   = LP(iLayer)%PartitionInfo(iNode+1,2)
                nstarty = LP(iLayer)%PartitionInfo(iNode+1,3)
                nendy   = LP(iLayer)%PartitionInfo(iNode+1,4)
                call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1), &
                    MPI_DOUBLE_PRECISION,iNode,2015,MPI_COMM_WORLD,istatus,ierror)
                if(iTimeStep.eq.0) then
                    do j=nstarty, nendy
                    do i=nstartx, nendx
                        LP(iLayer)%H(1,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                    enddo
                    enddo
                else
                    do j=nstarty, nendy
                    do i=nstartx, nendx
                        LP(iLayer)%H(2,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                    enddo
                    enddo
                endif
            endif
        enddo
        write(s,'(a,i2.2,a1,i6.6,a4)') TRIM(ADJUSTL(GP%ResultPath))//'z_',iLayer,'_',iTimeStep,'.dat'
        open(22,file=s,status='replace',form='unformatted')
        write(22) LP(iLayer)%NX, LP(iLayer)%NY
        if(iTimeStep.eq.0) then
            do j = 1,LP(iLayer)%NY
                write(22) (LP(iLayer)%H(1,i,j), i=1,LP(iLayer)%NX)
            enddo
        else
            do j = 1,LP(iLayer)%NY
                write(22) (LP(iLayer)%H(2,i,j), i=1,LP(iLayer)%NX)
            enddo
        endif
        close(22)
    endif
enddo


if(GP%SaveFlux.eq.1) then

do iLayer = 1,GP%NumLayers
    if(irank.lt.LP(iLayer)%nsize.and.irank.ne.master) then
        nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
        nendx   = MIN(LP(iLayer)%PartitionInfo(irank+1,2),LP(iLayer)%NX-1)
        nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
        nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
        if(iTimeStep.eq.0) then
            do j=nstarty, nendy
            do i=nstartx, nendx
                LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%M(1,i,j)
            enddo
            enddo
        else
            do j=nstarty, nendy
            do i=nstartx, nendx
                LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%M(2,i,j)
            enddo
            enddo
        endif
        call MPI_SEND(LocalData, (nendx-nstartx+1)*(nendy-nstarty+1),&
            MPI_DOUBLE_PRECISION,master,2015,MPI_COMM_WORLD,ierror)
    endif
    if(irank.eq.master) then
        do iNode=0, LP(iLayer)%nsize-1
            if(iNode.ne.master) then
                nstartx = LP(iLayer)%PartitionInfo(iNode+1,1)
                nendx   = MIN(LP(iLayer)%PartitionInfo(iNode+1,2),LP(iLayer)%NX-1)
                nstarty = LP(iLayer)%PartitionInfo(iNode+1,3)
                nendy   = LP(iLayer)%PartitionInfo(iNode+1,4)
                call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1), &
                    MPI_DOUBLE_PRECISION,iNode,2015,MPI_COMM_WORLD,istatus,ierror)
                if(iTimeStep.eq.0) then
                    do j=nstarty, nendy
                    do i=nstartx, nendx
                        LP(iLayer)%M(1,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                    enddo
                    enddo
                else
                    do j=nstarty, nendy
                    do i=nstartx, nendx
                        LP(iLayer)%M(2,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                    enddo
                    enddo
                endif
            endif
        enddo
        write(s,'(a,i2.2,a1,i6.6,a4)') TRIM(ADJUSTL(GP%ResultPath))//'M_',iLayer,'_',iTimeStep,'.dat'
        open(22,file=s,status='replace',form='unformatted')
        write(22) LP(iLayer)%NX-1, LP(iLayer)%NY
        if(iTimeStep.eq.0) then
            do j = 1,LP(iLayer)%NY
                write(22) (LP(iLayer)%M(1,i,j), i=1,LP(iLayer)%NX-1)
            enddo
        else
            do j = 1,LP(iLayer)%NY
                write(22) (LP(iLayer)%M(2,i,j), i=1,LP(iLayer)%NX-1)
            enddo
        endif
        close(22)
    endif
enddo

do iLayer = 1,GP%NumLayers
    if(irank.lt.LP(iLayer)%nsize.and.irank.ne.master) then
        nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
        nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
        nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
        nendy   = MIN(LP(iLayer)%PartitionInfo(irank+1,4),LP(iLayer)%NY-1)
        if(iTimeStep.eq.0) then
            do j=nstarty, nendy
            do i=nstartx, nendx
                LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%N(1,i,j)
            enddo
            enddo
        else
            do j=nstarty, nendy
            do i=nstartx, nendx
                LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%N(2,i,j)
            enddo
            enddo
        endif
        call MPI_SEND(LocalData, (nendx-nstartx+1)*(nendy-nstarty+1),&
            MPI_DOUBLE_PRECISION,master,2015,MPI_COMM_WORLD,ierror)
    endif
    if(irank.eq.master) then
        do iNode=0, LP(iLayer)%nsize-1
            if(iNode.ne.master) then
                nstartx = LP(iLayer)%PartitionInfo(iNode+1,1)
                nendx   = LP(iLayer)%PartitionInfo(iNode+1,2)
                nstarty = LP(iLayer)%PartitionInfo(iNode+1,3)
                nendy   = MIN(LP(iLayer)%PartitionInfo(iNode+1,4),LP(iLayer)%NY-1)
                call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1), &
                    MPI_DOUBLE_PRECISION,iNode,2015,MPI_COMM_WORLD,istatus,ierror)
                if(iTimeStep.eq.0) then
                    do j=nstarty, nendy
                    do i=nstartx, nendx
                        LP(iLayer)%N(1,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                    enddo
                    enddo
                else
                    do j=nstarty, nendy
                    do i=nstartx, nendx
                        LP(iLayer)%N(2,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                    enddo
                    enddo
                endif
            endif
        enddo
        write(s,'(a,i2.2,a1,i6.6,a4)') TRIM(ADJUSTL(GP%ResultPath))//'N_',iLayer,'_',iTimeStep,'.dat'
        open(22,file=s,status='replace',form='unformatted')
        write(22) LP(iLayer)%NX, LP(iLayer)%NY-1
        if(iTimeStep.eq.0) then
            do j = 1,LP(iLayer)%NY-1
                write(22) (LP(iLayer)%N(1,i,j), i=1,LP(iLayer)%NX)
            enddo
        else
            do j = 1,LP(iLayer)%NY-1
                write(22) (LP(iLayer)%N(2,i,j), i=1,LP(iLayer)%NX)
            enddo
        endif
        close(22)
    endif
enddo

endif

!/// save maximum water height in each layer ///!
if(iTimeStep.eq.GP%TotalTimeSteps) then
do iLayer = 1,GP%NumLayers
    if(irank.lt.LP(iLayer)%nsize.and.irank.ne.master) then
        nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
        nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
        nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
        nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
        do j = nstarty, nendy
        do i = nstartx, nendx
            LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%Hmax(i,j)
        enddo
        enddo
        call MPI_SEND(LocalData, (nendx-nstartx+1)*(nendy-nstarty+1), &
            MPI_DOUBLE_PRECISION,master,2015,MPI_COMM_WORLD,ierror)
    endif
    if(irank.eq.master) then
        do iNode = 0, LP(iLayer)%nsize-1
            if(iNode.ne.master) then
                nstartx = LP(iLayer)%PartitionInfo(iNode+1,1)
                nendx   = LP(iLayer)%PartitionInfo(iNode+1,2)
                nstarty = LP(iLayer)%PartitionInfo(iNode+1,3)
                nendy   = LP(iLayer)%PartitionInfo(iNode+1,4)
                call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1), &
                    MPI_DOUBLE_PRECISION,iNode,2015,MPI_COMM_WORLD,istatus,ierror)
                do j = nstarty, nendy
                do i = nstartx, nendx
                    LP(iLayer)%Hmax(i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                enddo
                enddo 
            endif
        enddo
        write(s,'(a,i2.2,a4)') TRIM(ADJUSTL(GP%ResultPath))//'zmax_',iLayer,'.dat'
        open(22,file=s,status='replace',form='unformatted')
        write(22) LP(iLayer)%NX, LP(iLayer)%NY
        do j = 1, LP(iLayer)%NY
            write(22) (LP(iLayer)%Hmax(i,j), i=1,LP(iLayer)%NX)
        enddo
        close(22)
    endif
enddo
endif !if: calculation is finished, and Hmax is obtained

!/// save minimum water height in each layer ///!
if(iTimeStep.eq.GP%TotalTimeSteps) then
do iLayer = 1,GP%NumLayers
    if(irank.lt.LP(iLayer)%nsize.and.irank.ne.master) then
        nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
        nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
        nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
        nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
        do j = nstarty, nendy
        do i = nstartx, nendx
            LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%Hmin(i,j)
        enddo
        enddo
        call MPI_SEND(LocalData, (nendx-nstartx+1)*(nendy-nstarty+1), &
            MPI_DOUBLE_PRECISION,master,2015,MPI_COMM_WORLD,ierror)
    endif
    if(irank.eq.master) then
        do iNode = 0, LP(iLayer)%nsize-1
            if(iNode.ne.master) then
                nstartx = LP(iLayer)%PartitionInfo(iNode+1,1)
                nendx   = LP(iLayer)%PartitionInfo(iNode+1,2)
                nstarty = LP(iLayer)%PartitionInfo(iNode+1,3)
                nendy   = LP(iLayer)%PartitionInfo(iNode+1,4)
                call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1), &
                    MPI_DOUBLE_PRECISION,iNode,2015,MPI_COMM_WORLD,istatus,ierror)
                do j = nstarty, nendy
                do i = nstartx, nendx
                    LP(iLayer)%Hmin(i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                enddo
                enddo 
            endif
        enddo
        write(s,'(a,i2.2,a4)') TRIM(ADJUSTL(GP%ResultPath))//'zmin_',iLayer,'.dat'
        open(22,file=s,status='replace',form='unformatted')
        write(22) LP(iLayer)%NX, LP(iLayer)%NY
        do j = 1, LP(iLayer)%NY
            write(22) (LP(iLayer)%Hmin(i,j), i=1,LP(iLayer)%NX)
        enddo
        close(22)
    endif
enddo
endif !if: calculation is finished, and Hmin is obtained

!///output the snapshot of dynamic pressure///!
do iLayer = 1,GP%NumLayers

if(GP%SaveDynamicPressure.eq.1.and.LP(iLayer)%Dispersion.eq.1) then

    if(irank.lt.LP(iLayer)%nsize.and.irank.ne.master) then
        nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
        nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
        nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
        nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
        do j = nstarty, nendy
        do i = nstartx, nendx
            LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%Q(i,j)
        enddo
        enddo
        call MPI_SEND(LocalData, (nendx-nstartx+1)*(nendy-nstarty+1), &
            MPI_DOUBLE_PRECISION,master,2015,MPI_COMM_WORLD,ierror)
    endif
    if(irank.eq.master) then
        do iNode = 0, LP(iLayer)%nsize-1
        if(iNode.ne.master) then
            nstartx = LP(iLayer)%PartitionInfo(iNode+1,1)
            nendx   = LP(iLayer)%PartitionInfo(iNode+1,2)
            nstarty = LP(iLayer)%PartitionInfo(iNode+1,3)
            nendy   = LP(iLayer)%PartitionInfo(iNode+1,4)
            call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1), &
                MPI_DOUBLE_PRECISION,iNode,2015,MPI_COMM_WORLD,istatus,ierror)
            do j = nstarty, nendy
            do i = nstartx, nendx
                LP(iLayer)%Q(i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
            enddo
            enddo 
        endif
        enddo
        write(s,'(a,i2.2,a1,i6.6,a4)') TRIM(ADJUSTL(GP%ResultPath))//'Q_',iLayer,'_',iTimeStep,'.dat'
        open(22,file=s,status='replace',form='unformatted')
        write(22) LP(iLayer)%NX, LP(iLayer)%NY
        do j = 1, LP(iLayer)%NY
            write(22) (LP(iLayer)%Q(i,j), i=1,LP(iLayer)%NX)
        enddo
        close(22)
    endif
endif !if: dispersion is calculated in this layer and Q needs saving
enddo 

call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,5) = GP%CPUTime(irank+1,5)+CPUTime2-CPUTime1

end subroutine saveSnapshot



subroutine deallocateVariables(GP, LP, SP)

use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
type(StationParameters) ::  SP(999)
integer*4               ::  iLayer, iSta, error

!/// Release the space assigned to GlobalPrameters ///!
DEALLOCATE(GP%t, stat=error)
    
!/// Release the space assigned to LayerParameters ///!
do iLayer = 1,GP%NumLayers
    DEALLOCATE(LP(iLayer)%X, LP(iLayer)%Y, LP(iLayer)%Z, stat=error)
    DEALLOCATE(LP(iLayer)%PermanentDryCells, LP(iLayer)%FloodingCells, stat = error)
    DEALLOCATE(LP(iLayer)%FloodingCellsWaterHeight, stat=error)
    DEALLOCATE(LP(iLayer)%CPX, LP(iLayer)%CPY, LP(iLayer)%CSY, stat=error)
    DEALLOCATE(LP(iLayer)%CS1, LP(iLayer)%CS2, LP(iLayer)%CS3, stat=error)
    DEALLOCATE(LP(iLayer)%CS6, LP(iLayer)%CS7, LP(iLayer)%CS8, stat=error)
    DEALLOCATE(LP(iLayer)%H, LP(iLayer)%M, LP(iLayer)%N, stat=error)
    DEALLOCATE(LP(iLayer)%D_M, LP(iLayer)%D_N, stat=error)
    DEALLOCATE(LP(iLayer)%Hmax, LP(iLayer)%Hmin, LP(iLayer)%H0, LP(iLayer)%M0, LP(iLayer)%N0, stat=error)
    DEALLOCATE(LP(iLayer)%HF, LP(iLayer)%MF, LP(iLayer)%NF, LP(iLayer)%HK, stat=error)
    DEALLOCATE(LP(iLayer)%Sign_M, LP(iLayer)%Sign_N, stat=error)
    DEALLOCATE(LP(iLayer)%CSD1, LP(iLayer)%CSD2, LP(iLayer)%CSD3, LP(iLayer)%CSD4, stat=error)
    DEALLOCATE(LP(iLayer)%Q, LP(iLayer)%Q0, LP(iLayer)%QF, LP(iLayer)%Has_Q, stat=error)
    DEALLOCATE( LP(iLayer)%WB, LP(iLayer)%WS, stat=error)
    DEALLOCATE(LP(iLayer)%PL, LP(iLayer)%PR, LP(iLayer)%PB, LP(iLayer)%PT, LP(iLayer)%PC, LP(iLayer)%PQ, stat=error)
    DEALLOCATE(LP(iLayer)%BrkAge, LP(iLayer)%BrkAge0, LP(iLayer)%BrkAgeF, stat=error)
    DEALLOCATE(LP(iLayer)%Brkv, LP(iLayer)%Brkv0, LP(iLayer)%BrkvF, stat=error)
    DEALLOCATE(LP(iLayer)%CSB1, LP(iLayer)%CSB2, LP(iLayer)%CSB3, LP(iLayer)%CSB4, stat=error)
enddo

!/// Release the space assigned to StationParameters ///!
do iSta = 1,GP%NumStations
    DEALLOCATE(SP(iSta)%H, SP(iSta)%M, SP(iSta)%N, SP(iSta)%Q, stat=error)
enddo

end subroutine deallocateVariables