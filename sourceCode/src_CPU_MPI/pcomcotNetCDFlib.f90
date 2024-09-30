
    subroutine getBathymetryDataSizeNetCDF(GP, LP, iLayer)

    use mpi  
    use netcdf
    use VariableDefination
    implicit NONE
    type(GlobalParameters)    ::  GP
    type(LayerParameters)     ::  LP(100)
    integer*4                 ::  iLayer, fid, vid
    real*8, allocatable       ::  x(:), y(:)
    character(NF90_MAX_NAME)  ::  chartmp

    call handleCDFError(NF90_OPEN(TRIM(ADJUSTL(LP(iLayer)%BathymetryFileName)), NF90_NOWRITE, fid))

    vid = 1; ! vid = 1: x coordinates, double (real*8)
    call handleCDFError(NF90_INQUIRE_DIMENSION(fid, vid, chartmp, LP(iLayer)%NX))
    ALLOCATE(x(LP(iLayer)%NX))
    call handleCDFError(NF90_GET_VAR(fid, vid, x))
    if(iLayer.eq.1) then
        GP%StartEastWest = 0
        if((GP%CoordinatesType.eq.0).and.(x(1).ge.0).and.(x(1).le.180.0)) then
            GP%StartEastWest = 1
        elseif((GP%CoordinatesType.eq.0).and.(x(1).gt.180.0.or.x(1).lt.0.0)) then
            GP%StartEastWest = 2
        endif
    endif
    if((GP%StartEastWest.eq.1).and.(x(1).lt.0)) x(1) = x(1)+360
    if((GP%StartEastWest.eq.2).and.(x(1).gt.180)) x(1) = x(1)-360
    if((GP%StartEastWest.eq.1).and.(x(2).lt.0)) x(2) = x(2)+360
    if((GP%StartEastWest.eq.2).and.(x(2).gt.180)) x(2) = x(2)-360
    if((GP%StartEastWest.eq.1).and.(x(LP(iLayer)%NX).lt.0)) x(LP(iLayer)%NX) = x(LP(iLayer)%NX)+360
    if((GP%StartEastWest.eq.2).and.(x(LP(iLayer)%NX).gt.180)) x(LP(iLayer)%NX) = x(LP(iLayer)%NX)-360
    LP(iLayer)%xmin = x(1);
    LP(iLayer)%xmax = x(LP(iLayer)%NX);
    LP(iLayer)%dx = x(2)-x(1)


    vid = 2; ! vid = 2: y coordinates, double (real*8)
    call handleCDFError(NF90_INQUIRE_DIMENSION(fid, vid, chartmp, LP(iLayer)%NY))
    ALLOCATE(y(LP(iLayer)%NY))
    call handleCDFError(NF90_GET_VAR(fid, vid, y))
    LP(iLayer)%ymin = y(1);
    LP(iLayer)%ymax = y(LP(iLayer)%NY);
    LP(iLayer)%dy = y(2)-y(1)

    call handleCDFError(NF90_CLOSE(fid))

    end subroutine getBathymetryDataSizeNetCDF



    subroutine readBathymetryNetCDF(GP, LP, iLayer)

    use mpi  
    use netcdf
    use VariableDefination
    implicit NONE
    type(GlobalParameters)    ::  GP
    type(LayerParameters)     ::  LP(100)
    integer*4                 ::  iLayer, fid, vid, i, j
    real*4, allocatable       ::  z(:,:)
    character(NF90_MAX_NAME)  ::  chartmp

    call handleCDFError(NF90_OPEN(TRIM(ADJUSTL(LP(iLayer)%BathymetryFileName)), NF90_NOWRITE, fid))

    vid = 1; ! vid = 1: x coordinates, double (real*8)
    call handleCDFError(NF90_GET_VAR(fid, vid, LP(iLayer)%X))
     do i = 1, LP(iLayer)%NX
        if(GP%StartEastWest.eq.1.and.LP(iLayer)%X(i).lt.0) LP(iLayer)%X(i) = LP(iLayer)%X(i)+360
        if(GP%StartEastWest.eq.2.and.LP(iLayer)%X(i).gt.180) LP(iLayer)%X(i) = LP(iLayer)%X(i)-360
    enddo

    vid = 2; ! vid = 2: y coordinates, double (real*8)
    call handleCDFError(NF90_GET_VAR(fid, vid, LP(iLayer)%Y))

    vid = 3; ! vid = 3: z values, float (real*4)
    ALLOCATE(z(LP(iLayer)%NX, LP(iLayer)%NY))
    call handleCDFError(NF90_GET_VAR(fid, vid, z))
    do i = 1, LP(iLayer)%NX
        do j = 1, LP(iLayer)%NY
            LP(iLayer)%Z(i,j) = DBLE(-z(i,j))
            LP(iLayer)%zmax = MAX(LP(iLayer)%zmax, LP(iLayer)%Z(i,j))
            LP(iLayer)%zmin = MIN(LP(iLayer)%zmin, LP(iLayer)%Z(i,j))
        enddo
    enddo

    call handleCDFError(NF90_CLOSE(fid))

    end subroutine readBathymetryNetCDF



    subroutine readInitialElevationNetCDF(GP, LP)
        
    use mpi  
    use netcdf
    use VariableDefination
    implicit NONE
    type(GlobalParameters)      ::  GP
    type(LayerParameters)       ::  LP(100)
    integer*4                   ::  fid, vid, i, j
    real*4, allocatable         ::  z(:,:)           
    character(NF90_MAX_NAME)  ::  chartmp

    call handleCDFError(NF90_OPEN(TRIM(ADJUSTL(GP%InitialElevationFileName)), NF90_NOWRITE, fid))

    vid = 3; ! vid = 3: z values, float (real*4)
    ALLOCATE(z(LP(GP%TopLayer)%NX, LP(GP%TopLayer)%NY))
    call handleCDFError(NF90_GET_VAR(fid, vid, z))

    do j = 1, LP(GP%TopLayer)%NY
    do i = 1, LP(GP%TopLayer)%NX
        LP(GP%TopLayer)%H(1,i,j) = DBLE(z(i,j))
    enddo
    enddo
    
    call handleCDFError(NF90_CLOSE(fid))

    end subroutine readInitialElevationNetCDF



    subroutine handleCDFError(ist)
      
    use netcdf
    implicit NONE
    integer*4  ::  ist

    if(ist.ne.NF90_NOERR) then
        write(*,*) TRIM(NF90_STRERROR(ist))
        stop "CDF ERROR"
    endif

    end subroutine handleCDFError

