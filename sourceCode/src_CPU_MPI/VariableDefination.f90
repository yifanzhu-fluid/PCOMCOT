
module VariableDefination

!// kind number of double-precidion floating-point
integer*4, parameter    ::  fp_kind = 8



type :: GlobalParameters

    !****************************** Read From Config ******************************!
    !/// Basic Control Parmeters ///!
    real*8      ::  Version = 2.0          ! pcomcot version
    integer*4   ::  PurposeCalculation     ! 1/Regular Simulation; 2,3/Compute Green Functions of (Subfaults/Initial Height)
    integer*4   ::  InitialConditionType   ! Initial condition  0/file; 1/fault
    integer*4   ::  CoordinatesType        ! 0/spherical; 1/cartisian
    real*8      ::  TotalTime              ! Total run time
    real*8      ::  dt                     ! Time step
    real*8      ::  DTSaveData             ! Time interval to save snapshots
    integer*4   ::  SaveFlux               ! 0/no; 1/save flux
    integer*4   ::  SaveDynamicPressure    ! 0/no; 1/save dynamic pressure Q
    integer*4   ::  MinGridsPerNode        ! Min grids on each computing node
    integer*4   ::  FeedbackToParentLayer  ! 0/no feedback; 1/yes

    !/// Parameters for Tsunami Source ///!
    integer*4   ::  HorizontalMotion        ! 0/only vertical displacement; 1/consider horizontal motion using Tanioka (1996)
    integer*4   ::  KajiuraFilter           ! 0/no Kajiura filtering; 1/determine surface elevation based on Kajiura' theory (1963)
    integer*4   ::  UseAverageDepth         ! 0/use local water depth for Kajiura filtering; 1/use average depth of source area
    real*8      ::  MinKajiuraDepth         ! Minimum water depth for applying Kajiura filter
    real*8      ::  MaxBottomSlope = 0.5d0  ! Maximum bottom slope allowed by the model

    !/// Parameters for Wave Physics and Numerics ///!
    integer*4   ::  Nonlinearity                ! 0/linear waves; 1/nonlinear waves
    integer*4   ::  Dispersion                  ! 0/non-dispersive waves; 1/dispersive waves
    integer*4   ::  DepthVariability            ! 0/smooth bottom; 1/relatively steep bottom
    real*8      ::  MinDispersionDepth          ! Minimum flow depth for calculating dynamic pressure
    integer*4   ::  Breaking                    ! 0/non-breaking waves; 1/breaking waves
    integer*4   ::  FluxCenter                  ! 0/TFSC scheme for LSWE; 1/use flux-centered scheme to improve stability
    real*8      ::  centerWeighting0 = 0.9995d0 ! center weighting for TFSC scheme
    real*8      ::  FrCap                       ! Upper bound of Froud number to avoid unreasonable velocity
    real*8      ::  DTFilterData                ! Time interval to filter out short wave components
    integer*4   ::  NDTFilterData               ! Time step interval to apply 9-point filter to wave field
    
    !/// Parameters for Boundary Condition ///!
    integer*4   ::  BoundaryConditionType  ! Boundary condition 0/open; 1/wall; 2/sponge **only wall and sponge are supported
    real*8      ::  SpongeWidthX           ! Width of sponge layer in x direction (unit: meter)
    real*8      ::  SpongeWidthY           ! Width of sponge layer in y direction (unit: meter)
    real*8      ::  MaxSpongeMANNING       ! Maximum Manning coefficient for friction-type sponge
    real*8      ::  SpongeDampingA         ! Damping coefficient A for L-D type sponge
    real*8      ::  SpongeDampingR         ! Damping coefficient R for L-D type sponge

    !/// Parameters for Inundation ///!
    real*8      ::  PermanentDryLimit      ! Considered to be permanently dry if ground elevation is no less than this
    real*8      ::  MinWaterDepth          ! Wet cells turn into dry, if total depth is no more than this
    real*8      ::  FrictionDepthLimit     ! Bottom friction is not calculated for too shallow area to avoid instability
    real*8      ::  MANNING                ! Manning coefficient for bottom friction

    !/// Parameters for Green Functions ///!
    real*8      ::  SourceStartX, SourceEndX, SourceDX  ! Starting/Ending/GridSize X for Source Area
    real*8      ::  SourceStartY, SourceEndY, SourceDY  ! Starting/Ending/GridSize Y for Source Area
    integer*4   ::  SourceBasisFunctionType             ! 1/Gaussian; 2/Sigmoid Function 1/(1+exp(-x+x_left))+1/(1+exp(x-x_right))-1
    real*8      ::  GaussianRatio                       ! Ratio of Gaussian Radius v.s. Source Grids Size
    real*8      ::  SigmoidCoefficient                  ! epsilon coefficient for sigmoid function to approximate step function
    
    !************************* Computed for Convinience **************************!
    integer*4         ::  NumLayers                  ! Value set in "checkFiles"
    character(999)    ::  InitialElevationFileName   ! Value set in "checkFiles"
    integer*4         ::  InitialElevationFileFormat ! Value set in "checkFiles" 1/nf; 2/xyz
    integer*4         ::  StartEastWest              ! Vaule set in "getBathymetryDataSize";
                                                     ! 0/Cartisian; 1/spherical long:160E~-160W(200E); 2/special long: 340W(-20E)~20E
    integer*4         ::  NumLayerLevels             ! Value set in "determineLayerDependency"
    integer*4         ::  TopLayer                   ! Value set in "determineLayerDependency"
    integer*4         ::  NumStations                ! Value set in "readStations"
    integer*4         ::  NumFaults                  ! Value set in "readFaultParameters"
    integer*4         ::  TotalTimeSteps             ! Value set in "computeGlobalParameters"
    integer*4         ::  NDTSaveData                ! Value set in "computeGlobalParameters"
    integer*4         ::  nCalculations              ! Value set in "computeGlobalParameters"
    integer*4         ::  SourceNX, SourceNY         ! Value set in "computeGlobalParaemters"
    real*8,dimension(:),pointer  ::  t               ! Value set when solving equations

    !************************* Default System Parameters *************************!
    character(999)    ::  COMCOTParametersFileName  = "pcomcot.ctl"
    character(999)    ::  COMCOTLayerCtlFileName  = "layers.ctl"
    character(999)    ::  BathymetryFilePrefix = "layer"                    ! Bathymetry: .nf/.xyz
    character(999)    ::  InitialElevationFilePrefix = "InitialElevation"   ! Initial elevation: .nf/.xyz
    character(999)    ::  InitialMFileName = 'InitialFluxM.xyz'             ! Initial flux M(optional)
    character(999)    ::  InitialNFileName = 'InitialFluxN.xyz'             ! Initial flux N(optional)
    character(999)    ::  FaultParametersFileName = "FaultParameters.ctl"   ! Fault parameters
    character(999)    ::  StationFileName = "Stations.ctl"                  ! File for stations coordinates
    character(999)    ::  GFParamFileName = "InitialElevationGFParameters.dat" ! File for initial height GF parameters
    character(999)    ::  GFParamFormatFileName = "InitialElevationGFParametersFormat.txt"
    character(999)    ::  ResultPath = 'PCOMCOToutput/'                     ! directory of all output files
    integer*4         ::  ComputeDivisionOpt = 1                            ! 1: divide domain averagely; 2: child layer on 1 node
    integer*4         ::  nRowBoundary = 2                                  ! Sync boundary rows between subdomains (<10)
    integer*4         ::  nRowBoundaryFlux = 2                              ! Sync boundary rows between subdomains for flux
    integer*4         ::  nGHOST = 4                                        ! Outermost 4 grid points in child layer are from parent
    integer*4         ::  nRowBathymetry = 3                                ! At least 3 cells are wet, smoothing bathymetry
    integer*4         ::  nRowBathymetryBoundary = 10                       ! At least 10 cells are wet from boundary
    real*8            ::  GRAV = 9.807
    real*8            ::  PI = 3.141593
    real*8            ::  R_EARTH = 6378000.0
    real*8            ::  OMEGA = 7.2921159E-5

    !******************************  MPI Variables  ******************************!
    integer*4         ::  irank, nsizeTotal, master
    integer*4         ::  MaxNX, MinNX, MaxNY, MinNY   ! Value set in "partitionDomain"; Max/Min grids on one node
    real*8            ::  CPUTimeInitial               ! Value set in "pcomcot"
    real*8            ::  CPUTimeInitialCalculation    ! Value set in "pcomcot"
    real*8            ::  CPUTime(999,5)               ! Value set in all subroutines; iNode;AllTime/Prepare/CalculationTime/Communication/WriteToDisk


end type GlobalParameters



type  ::  LayerParameters

    integer*4       ::  Level                           ! Value set in "determineLayerDependency"
    integer*4       ::  Parent                          ! Value set in "determineLayerDependency"
    character(999)  ::  BathymetryFileName              ! Value set in "getBathymetryDataSize"
    integer*4       ::  BathymetryFileFormat            ! Value set in "getBathymetryDataSize"
    real*8          ::  xmin, dx, xmax                  ! Value set in "getBathymetryDataSize"
    real*8          ::  ymin, dy, ymax                  ! Value set in "getBathymetryDataSize"
    integer*4       ::  NX, NY                          ! Value set in "getBathymetryDataSize"
    integer*4       ::  Nonlinearity                    ! Value set in "readLayerConfig" 0/linear; 1/nonlinear
    integer*4       ::  Dispersion                      ! Value set in "readLayerConfig" 0/non-dispersive; 1/dispersive
    integer*4       ::  DepthVariability                ! Value set in "readLayerConfig" 0/smooth bottom; 1/steep bathymetry
    integer*4       ::  Breaking                        ! Value set in "readLayerConfig" 0/non-breaking; 1/breaking
    integer*4       ::  FluxCenter                      ! Value set in "readLayerConfig" 0/FTCS scheme for LSWE; 1/flux-centered scheme
    integer*4       ::  NDTFilterData                   ! Value set in "readLayerConfig" time step interval to apply filter
    real*8          ::  DTFilterData                    ! Value set in "readLayerConfig" time interval to apply filter
    real*8          ::  ComputingStartTime              ! Value set in "readLayerConfig" skip this layer before this time
    real*8, dimension(:), pointer     ::  X, Y          ! Value set in "readBathymetry"
    real*8, dimension(:,:), pointer   ::  Z             ! Value set in "readBathymetry"
    real*8          ::  zmin, zmax                      ! Value set in "readBathymetry"
    real*8          ::  dt                              ! Value set in "cflCheck"
    integer*4       ::  nStepsPerTimeStep               ! Value set in "cflCheck"
    integer*4       ::  nsize                           ! Value set in "partitionDomain"; ComputeNodes for layer
    integer*4       ::  npartx, nparty                  ! Value set in "partitionDomain"; partition
    integer*4       ::  MaxNX, MinNX, MaxNY, MinNY      ! Value set in "partitionDomain"; Max/Min grids on one node
    integer*4       ::  PartitionInfo(999,4)            ! Value set in "partitionDomain"; ComputeNode; nstartx,nendx,nstarty,nendy
    integer*4       ::  BoundarySendRecvCount           ! Value set in "partitionDomain"; Count
    integer*4       ::  BoundarySendRecv(999,14)        ! Value set in "partitionDomain"; FromNode,ToNode,(nstartx,nendx...)*3 for H,M,N
    integer*4       ::  ParentToChildSendRecvCount      ! Value set in "partitionDomain"; Count
    integer*4       ::  ParentToChildSendRecv(9999,8)   ! Value set in "partitionDomain"; FromNode,ToNode,Boundary Flag,H/M/N Flag,nstartx,nendx,nstarty,nendy
                                                        ! nstart, nend: indices in child layer iLayer; coordinates (x_i, y_j)
                                                        ! check this position belongs to FromNode in pLayer; interpolate value in pLayer on FromNode
                                                        ! store the boundary value in H/M/N_left/right/bottom/top
    integer*4       ::  ChildToParentSendRecvCount      ! Value set in "partitionDomain"; Count
    integer*4       ::  ChildToParentSendRecv(9999,6)   ! Value set in "partitionDomain"; FromNode,ToNode,nstartx,nendx,nstarty,nendy
                                                        ! nstart, nend: indices in parent layer pLayer; coordinates (x_i, y_j)
                                                        ! check this position belongs to FromNode in cLayer; average value in cLayer on FromNode
                                                        ! store value in H(2,i,j) on FromNode; send to ToNode; save in H(2,i,j)

    !**** Values Computed as Parameters on each node  ****!
    !**** Values set in "computeParameters" ****!
    integer*4                          ::  PermanentDryCellsCount
    integer*4, dimension(:,:), pointer ::  PermanentDryCells               ! cells that are permanently dry on this computing node
    integer*4                          ::  FloodingCellsCount
    integer*4, dimension(:,:), pointer ::  FloodingCells                   ! cells that may be changed between dry and wet (-PDL<Z<PDL)
    real*8, dimension(:), pointer      ::  FloodingCellsWaterHeight        ! water height after being flooded
    real*8, dimension(:,:), pointer    ::  SpongeMANNING                   ! Manning coefficient of friction-type sponge layer
    real*8, dimension(:), pointer      ::  CPX, CPY                        ! Coriolis parameters CP=2*Omega*SIN(Y)
    real*8                             ::  CC1, CC2, CC3, CC4     ! coefficients for linear and nonlinear terms in Cartesian 
    real*8, dimension(:), pointer      ::  CSY, CS1, CS2, CS3     ! coefficients for linear terms in spherical 
    real*8                             ::  CS4
    real*8                             ::  CS5                    ! coefficients for nonlinear terms in spherical 
    real*8, dimension(:), pointer      ::  CS6, CS7, CS8          ! coefficients for nonlinear terms in spherical  

    !**** Values Computed by Solving Equations  ****!
    real*8, dimension(:,:,:), pointer   ::  H       ! water elevation, value in last/this step; x; y position
    real*8, dimension(:,:,:), pointer   ::  M       ! volume flux in x direction (h*u)
    real*8, dimension(:,:,:), pointer   ::  N       ! volume flux in y direction (h*v)
    real*8, dimension(:,:), pointer     ::  D_M, D_N! numerical Depth where flux is defined
    real*8, dimension(:,:), pointer     ::  Hmax, Hmin      ! Maximum/Minimum surface elevation in this layer
    real*8, dimension(:,:), pointer     ::  H0, M0, N0      ! Initial H/M/N in each iTimeStep
    real*8, dimension(:,:), pointer     ::  HF, MF, NF      ! Final H/M/N in each iTimeStep (for child layer)
    real*8, dimension(:,:), pointer     ::  HK              ! Surface elevation smoothed by Kajiura filter
    integer*4, dimension(:,:), pointer  ::  Sign_M, Sign_N  ! Allowed flow situations of M and N: 0, +1, -1, +2, -2, 999
                                                            ! 0: both directions are allowed; 999: cannot flow; 
                                                            ! +1: only positive and from higher to lower; +2: ~ from lower to higher
                                                            ! -1: only negative and from higher to lower; -2: ~ from lower to higher

    !/// Variables for Calculating Dispersion ///!
    !*** Values Computed as Parameters on each node ***!
    real*8                          ::  CCD1, CCD2, CCD3, CCD4  ! coefficients for constructing equations for q in cartesian coordinates
    real*8, dimension(:), pointer   ::  CSD1, CSD2, CSD3, CSD4  ! coefficients for constructing equations for q in spherical coordinates
    real*8                          ::  CSD5

    !*** Values Computed by Solving Equations for Dynamic Pressure 
    real*8, dimension(:,:), pointer     ::  Q, Q0, QF   ! Dynamic pressure, needs nesting and synchronization
    integer*4, dimension(:,:), pointer  ::  Has_Q       ! 1: has dynamic pressure; 0: has no dynamic pressure
    real*8, dimension(:,:,:), pointer   ::  WB, WS      ! Vertical velocity at bottom and surface, dont need nesting or sync
    real*8, dimension(:,:), pointer     ::  PL, PR, PB, PT, PC, PQ   ! coefficients and forcing vector of linear equations for dynamic pressure

    !/// Variables for eddy-visocity type of wave breaking ///!
    real*8, dimension(:,:), pointer     ::  BrkAge, BrkAge0, BrkAgeF
    real*8, dimension(:,:), pointer     ::  Brkv, Brkv0, BrkvF
    !*** Values Computed as Parameters on each node ***!
    real*8                          ::  CCB1, CCB2             ! coefficients for eddy viscosity terms in Cartesian
    real*8, dimension(:), pointer   ::  CSB1, CSB2, CSB3, CSB4 ! coefficients for eddy viscosity terms in Spherical
       

end type LayerParameters



type  ::  StationParameters

    real*8      ::  X,Y      ! Value set in "readStations";    Stations coordinates
    integer*4   ::  nLayer   ! Value set in "readStations";    In which layer is the station located
    integer*4   ::  nNode    ! Value set in "partitionDomain"; On which node is the station located

    real*8, dimension(:), pointer      ::  H, M, N, Q  !Values at this station
    
end type StationParameters



type FaultParameters

    real*8         ::  T0           ! rupture starting time
    integer*4      ::  NT           ! time step when rupture starts
    real*8         ::  Depth        ! focal depth
    real*8         ::  Length       ! length of fault plane
    real*8         ::  Width        ! width of fault plane
    real*8         ::  Slip         ! dislocation
    real*8         ::  Rake         ! rake
    real*8         ::  HSlip        ! horizontal dislocation     HSlip = Slip*cos(Rake)
    real*8         ::  PSlip        ! perpendicular dislocation  PSlip = Slip*sin(Rake)
    real*8         ::  Strike       ! strike
    real*8         ::  Dip          ! dip
    real*8         ::  Y0           ! epicenter(latitude)
    real*8         ::  X0           ! epicenter(longitude)

end type FaultParameters



end module VariableDefination
