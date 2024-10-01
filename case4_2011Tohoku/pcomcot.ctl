################################################################################
#                                                                              #
#                Control file for PCOMCOT program (v2.1)                       #
#                                                                              #
################################################################################
#===============================================:===============================
# Basic Control Parameters                      :     Values                   |
#===============================================:===============================
 Purpose of Calculation           (Choose 1,2,3):      1
 Initial Condition             (0:file, 1:fault):      1
 Coordinate System    (0:spherical, 1:cartesian):      0
 Total run time                         (second):   6000.0
 Time step                              (second):      2.0
 Time interval to Save Snapshots        (second):    200.0
 Save Flux                         (0:no, 1:yes):      0
 Save Non-hydrostatic Pressure     (0:no, 1:yes):      0
 Minimum grids on each computing node           :      2000
 Feedback to parent layer          (0:no, 1:yes):      0

#===============================================:===============================
# Parameters for Surface Deformation            :     Values                   |
#===============================================:===============================
 Consider Horizontal Motion        (0:no; 1:yes):      1
 Apply Kajiura filter              (0:no; 1:yes):      1
 Use average depth for Kajiura     (0:no; 1:yes):      1
 Water depth Limit for Kajiura           (meter):    100.0

#===============================================:===============================
# Parameters for Wave Physics and Numerics      :     Values                   |
#===============================================:===============================
 Nonlinearity     (0:no/linear; 1:yes/nonlinear):      1
 Dispersion                        (0:no; 1:yes):      1
 Depth change for Dispersion (0:smooth; 1:steep):      0
 Water depth Limit for Dispersion        (meter):    100.0
 Breaking                          (0:no; 1:yes):      0
 Scheme for LSWEs      (0:FTCS; 1:flux-centered):      1
 Froude number Cap to Limit velocity            :     10.0                    
 Time interval to Apply filter          (second):  50000.0

#===============================================:===============================
# Parameters for Boundary Condition             :     Values                   |
#===============================================:===============================
 Boundary Condition           (1:Wall; 2:Sponge):      2
 Width of Sponge (West-East)             (meter): 100000.0
 Width of Sponge (South-North)           (meter): 100000.0
 Maximum Manning coefficient in Sponge          :    100.0
 Damping coefficient A                          :      2.0
 Damping coefficient R                          :      0.9 

#===============================================:===============================
# Parameters for Inundation                     :     Values                   |
#===============================================:===============================
 Permanent dry Limit                     (meter):     50.0
 Water depth Limit for wet -> dry        (meter):     0.01
 Water depth limit for Bottom Friction   (meter):     0.05
 Manning coefficient for Bottom friction        :     0.03

#===============================================:===============================
# Parameters for Computing Green's Functions    :     Values                   |
#===============================================:===============================
 Source Area Starting Longitude                 :     00.00
 Source Area Ending   Longitude                 :     00.00
 Source Area Discretizatoin Grid Size X  (meter):     00.00
 Source Area Starting Latitude                  :     00.00
 Source Area Ending   Latitude                  :     00.00
 Source Area Discretizatoin Grid Size Y  (meter):     00.00
 Basis Function Type        (1/Guass; 2/Sigmoid):      1
 Ratio of Gaussian Raidus / Grid Size           :      1.0
 Sigmoid Coefficient 1/(1+exp(-x/(DX*coef)))    :      1.0
 
#===============================================:===============================
# System Parameters (DO NOT CHANGE)             :     Values                   |
#===============================================:===============================
 Purpose of Calculation                      (1):  Regular Simulation
 Purpose of Calculation                      (2):  Calculation of Green's Functions (Fault)
 Purpose of Calculation                      (3):  Calculation of Green's Functions (Initial Height)
 Initial Water Elevation File Name              :  InitialWaterElevation.nf/xyz
 File Name of Fault Parameters                  :  FaultParameters.ctl
 File Name of Bathymetry Data                   :  layerXX.nf/xyz
 File Name of Stations                          :  Stations.ctl
 File Name of Initial Height GF Parameters      :  InitialHeightGFParameters.dat
