!**************************************************************************************************!
! PCOMCOT is a parallel multi-grid dispersive tsunami package developed by Yifan Zhu & Chao An.    !                                      !
! PCOMCOT simulates nonlinear, dispersive and breaking tsunamis by combining shallow water         ! 
! equations and non-hydrostatic pressure model.                                                    !
! The first version for linear shallow water equations was develped by Chao An                     ! 
! at Cornell University, New York in 2012~2017.                                                    !
! The current version with nonlinearity, dispersion, breaking, and runup is developed by Yifan Zhu ! 
! at Shanghai Jiao Tong University, Shanghai in 2021~2023.                                         ! 
! Report bugs to: Yifan Zhu (zyftop@sjtu.edu.cn) or Chao An (anchao@sjtu.edu.cn)                   !                                                  !
! Last updated: April/20/2024                                                                      !
!**************************************************************************************************!

program pcomcot

use mpi
use VariableDefination
      
implicit NONE
integer*4   ::  irank, nsize, master
integer*4   ::  ierror, istatus(MPI_STATUS_SIZE)
      
type(GlobalParameters)      ::  GP
type(LayerParameters)       ::  LP(100)
type(StationParameters)     ::  SP(999)
type(FaultParameters)       ::  FP(4000)
integer*4           ::  iLayer, iLayerLevel, iSta, iCal, iTimeStep, iStep
real*8              ::  CPUTime1, CPUTime2
integer*4           ::  LocalDataLength
real*8, allocatable ::  LocalData(:)

!######### initialize MPI evironment #########!
call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nsize, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierror)
call CPU_TIME(GP%CPUTimeInitial); GP%CPUTime(irank+1,1) = 0.0d0
master = 0; GP%irank = irank; GP%nsizeTotal = nsize; GP%master = master

!########## read data on master node ##########!
if(irank.eq.master) then
    call readConfig(GP)
    call checkFiles(GP)
    call getBathymetryDataSize(GP, LP)
    call determineLayerDependency(GP, LP)
    call readLayerConfig(GP, LP)
    do iLayer = 1,GP%NumLayers
        ALLOCATE(LP(iLayer)%X(LP(iLayer)%NX))
        ALLOCATE(LP(iLayer)%Y(LP(iLayer)%NY))
        ALLOCATE(LP(iLayer)%Z(LP(iLayer)%NX,LP(iLayer)%NY))
    enddo
    call readBathymetry(GP, LP)
    call cflCheck(GP, LP)
    call readStations(GP, LP, SP)
    if((GP%PurposeCalculation.eq.1.and.GP%InitialConditionType.eq.1).or.(GP%PurposeCalculation.eq.2)) then
        call readFaultParameters(GP, FP)
    endif
    call partitionDomain(GP, LP, SP)
    call computeGlobalParameters(GP, LP, SP)
endif

!#### broadcast dimensions, allocate memory ####!
call bcastCommonInfo(GP, LP, SP, FP)
ALLOCATE(GP%t(GP%TotalTimeSteps+1))
call allocateLayerVariables(GP, LP, SP)
LocalDataLength = (GP%MaxNX+10)*(GP%MaxNY+10)
ALLOCATE(LocalData(LocalDataLength))

!#### broadcast, compute arrays after memories allocated ####!
call bcastBathymetry(GP, LP, SP, FP, LocalData, LocalDataLength)
call computeParameters(GP, LP, SP, FP)
call CPU_TIME(GP%CPUTime(irank+1,2))
GP%CPUTime(irank+1,2) = GP%CPUTime(irank+1,2)-GP%CPUTimeInitial

!############ start computing ############!
if(irank.eq.master) call removeStationFiles(GP, LP, SP)
      
do iCal = 1,GP%nCalculations
call screenOutput(GP, LP, iCal, -1)
!/// initialize variables ///!
do iLayer = 1,GP%NumLayers
    LP(iLayer)%H = 0.0d0; LP(iLayer)%M = 0.0d0; LP(iLayer)%N = 0.0d0 
    LP(iLayer)%D_M = 0.0d0; LP(iLayer)%D_N = 0.0d0
    LP(iLayer)%Hmax = 0.0d0; LP(iLayer)%Hmin = 0.0d0
    if(LP(iLayer)%Dispersion.eq.1) then
        LP(iLayer)%Q = 0.0d0; LP(iLayer)%WB = 0.0d0; LP(iLayer)%WS = 0.0d0
    endif
    if(LP(iLayer)%Breaking.eq.1) then
        LP(iLayer)%Brkv = 0.0d0; LP(iLayer)%BrkAge = 0.0d0
    endif
enddo
do iTimeStep = 0,GP%TotalTimeSteps
    GP%t(iTimeStep+1) = iTimeStep*GP%dt
    call getInitialCondition(GP, LP, SP, FP, iCal, iTimeStep, LocalData, LocalDataLength)
    if(iTimeStep.eq.1) call CPU_TIME(GP%CPUTimeInitialCalculation)
    call screenOutput(GP, LP, iCal, iTimeStep)

    if(iTimeStep.ne.0) then
    do iLayerLevel = 1,GP%NumLayerLevels
    do iLayer = 1,GP%NumLayers
    if(LP(iLayer)%Level.eq.iLayerLevel) then
        if(iTimeStep*GP%dt.lt.LP(iLayer)%ComputingStartTime) cycle            
        if(iLayerLevel.gt.1) then
            call getLayerBoundaryFromParent(GP, LP, iLayer, LocalData, LocalDataLength)
            if(LP(iLayer)%Dispersion.eq.1) then
                call getLayerBoundaryDynamicPressureFromParent(GP, LP, iLayer, Localdata, LocalDataLength)
            endif
            if(LP(iLayer)%Breaking.eq.1) then
                call getLayerBoundaryViscosityFromParent(GP, LP, iLayer, LocalData, LocalDataLength)
            endif
        endif
        do iStep = 1,LP(iLayer)%nStepsPerTimeStep
            if(iLayerLevel.gt.1) then
                call getLayerBoundaryAtFineTimeStep(GP, LP, iLayer, iStep)
                if(LP(iLayer)%Dispersion.eq.1) then
                    call getLayerBoundaryDynamicPressureAtFineTimeStep(GP, LP, iLayer, iStep)
                endif
                if(LP(iLayer)%Breaking.eq.1) then
                    call getLayerBoundaryViscosityAtFineTimeStep(GP, LP, iLayer, iStep)
                endif
            endif
            call mass(GP, LP, iLayer, iStep, iTimeStep, LocalData, LocalDataLength)
            call bcastComputeDomainBoundaryValue(GP, LP, iLayer, LocalData, LocalDataLength)
            call momentum(GP, LP, iLayer, iStep, iTimeStep, LocalData, LocalDataLength)
            call bcastComputeDomainBoundaryFlux(GP, LP, iLayer, LocalData, LocalDataLength)
            if(GP%BoundaryConditionType.eq.2.and.iLayer.eq.GP%TopLayer) then
                call sponge(GP, LP, iLayer)
            endif
            if(.not.(GP%FeedbackToParentLayer.eq.1.and.iStep.eq.LP(iLayer)%nStepsPerTimeStep)) then
                call updateComputedResults(GP, LP, iLayer)
            endif
        enddo
    endif
    enddo
    enddo
    endif

    if(iTimeStep.ne.0.and.GP%FeedbackToParentLayer.eq.1) then
    do iLayerLevel = GP%NumLayerLevels,2,-1
    do iLayer = 1,GP%NumLayers
    if(LP(iLayer)%Level.eq.iLayerLevel.and.iTimeStep*GP%dt.ge.LP(iLayer)%ComputingStartTime) then
        call feedbackToParentLayerByAverage(GP, LP, iLayer, LocalData, LocalDataLength)
        call momentum(GP, LP, LP(iLayer)%Parent, 0, iTimeStep, LocalData, LocalDataLength)
        call bcastComputeDomainBoundaryFlux(GP, LP, LP(iLayer)%Parent, LocalData, LocalDataLength)
    endif
    enddo
    enddo
    do iLayer = 1,GP%NumLayers
        call updateComputedResults(GP, LP, iLayer)
    enddo
    endif

    call calculateStationData(GP, LP, SP, iCal, iTimeStep)
    if(MOD(iTimeStep,GP%NDTSaveData).eq.0.or.iTimeStep.eq.GP%TotalTimeSteps) then
        call saveSnapshot(GP, LP, iCal, iTimeStep, LocalData, LocalDataLength)
    endif
enddo !do loop for iTimeStep
call screenOutput(GP, LP, iCal, GP%TotalTimeSteps+100)

call saveStationData(GP, LP, SP, iCal)
enddo !do loop for iCal

call CPU_TIME(GP%CPUTIME(irank+1,1))
GP%CPUTime(irank+1,1) = GP%CPUTime(irank+1,1)-GP%CPUTimeInitial
call screenOutput(GP, LP, GP%nCalculations+100, GP%TotalTimeSteps+100)

call deallocateVariables(GP, LP, SP)
call MPI_FINALIZE(ierror)

end program pcomcot

