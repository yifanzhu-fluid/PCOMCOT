
      subroutine getBathymetryDataSizeNetCDF(GP, LP, iLayer)

      use mpi
      use VariableDefination
      implicit NONE
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(GP%NumLayers)
      integer*4                ::  iLayer, errorcode, ierror

      write(*,*) 'ERROR: pcomcot compiled without netcdf lib. can''t read .nf files.'
      call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)

      end subroutine getBathymetryDataSizeNetCDF



      subroutine readBathymetryNetCDF(GP, LP, iLayer)

      use mpi
      use VariableDefination
      implicit NONE
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(GP%NumLayers)
      integer*4                ::  iLayer, errorcode, ierror

      write(*,*) 'ERROR: pcomcot compiled without netcdf lib. can''t read .nf files.'
      call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)

      end subroutine readBathymetryNetCDF



      subroutine readInitialElevationNetCDF(GP, LP)
      
      use mpi      
      use VariableDefination
      implicit NONE
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(GP%NumLayers)
      integer*4                ::  errorcode, ierror

      write(*,*) 'ERROR: pcomcot compiled without netcdf lib. can''t read .nf files.'
      call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)

      end subroutine readInitialElevationNetCDF
