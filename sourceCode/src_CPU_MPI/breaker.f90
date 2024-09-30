subroutine addBreaking(GP, LP, iLayer, iStep, iTimeStep, LocalData, LocalDataLength)

use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(LayerParameters)    ::  LP(100)
integer*4   ::  iLayer, iStep, iTimeStep, LocalDataLength
real*8      ::  LocalData(LocalDataLength)

if(iTimeStep.gt.1) then ! We assume that there is no breaking at the first major step
    call eddyViscosity(GP, LP, iLayer)
    call bcastComputeDomainBoundaryViscosity(GP, LP, iLayer, LocalData, LocalDataLength)
    call smoothEddyViscosity(GP, LP, iLayer)
    if(GP%CoordinatesType.eq.1) then
        call breakingCartesian(GP, LP, iLayer)
    elseif(GP%CoordinatesType.eq.0) then
        call breakingSpherical(GP, LP, iLayer)
    endif
else
    LP(iLayer)%BrkAge = 0.0d0; LP(iLayer)%Brkv = 0.0d0
endif

end subroutine addBreaking




subroutine breakingCartesian(GP, LP, iLayer)

use VariableDefination
implicit NONE
type(GlobalParameters)    ::  GP
type(LayerParameters)     ::  LP(100)
integer*4   ::  iLayer, irank
integer*4   ::  nstartx, nendx, nstarty, nendy
integer*4   ::  istart, iend, jstart, jend
integer*4   ::  i, j
real*8      ::  m0, m1, m2, m3, m5
real*8      ::  n0, n1, n2, n3, n5
real*8      ::  v0, v1, v2, v3, v4, v5, v6, v7
real*8      ::  CPUTime1, CPUTime2

irank = GP%irank
call CPU_TIME(CPUTime1)

if(irank.lt.LP(iLayer)%nsize) then
   
nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)

!///calculate eddy viscosity term for flux component M///!
if(LP(iLayer)%Level.eq.1) then !calculate all grids in top layer
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX-1,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY,nendy)
elseif(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids in child layer
    istart = MAX(GP%nGHOST+1,nstartx); iend = MIN(LP(iLayer)%NX-GP%nGHOST-1,nendx)
    jstart = MAX(GP%nGHOST+1,nstarty); jend = MIN(LP(iLayer)%NY-GP%nGHOST,nendy)
endif
do j = jstart,jend
do i = istart,iend
    !///exclude the cell where water cannot flow freely///!
    if(LP(iLayer)%Sign_M(i,j).ne.0) cycle
    !/// add eddy viscosity term to M ///!
    m0 = LP(iLayer)%M(1,i,j); v0 = LP(iLayer)%Brkv(i,j); v2 = LP(iLayer)%Brkv(i+1,j)
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
    if(j.eq.1) then
        m3 = 0.0d0; v3 = 0.0d0; v4 = 0.0d0
    else
        m3 = LP(iLayer)%M(1,i,j-1)
        v3 = LP(iLayer)%Brkv(i,j-1); v4 = LP(iLayer)%Brkv(i+1,j-1)
    endif
    if(j.eq.LP(iLayer)%NY) then
        m5 = 0.0d0; v5 = 0.0d0; v6 = 0.0d0
    else
        m5 = LP(iLayer)%M(1,i,j+1)
        v5 = LP(iLayer)%Brkv(i,j+1); v6 = LP(iLayer)%Brkv(i+1,j+1)
    endif
    LP(iLayer)%M(2,i,j) = LP(iLayer)%M(2,i,j) + LP(iLayer)%CCB1*(v2*(m2-m0)-v0*(m0-m1)) &
        + 0.25*LP(iLayer)%CCB2*((v0+v2+v5+v6)*(m5-m0)-(v0+v2+v3+v4)*(m0-m3))
enddo
enddo

!///calculate eddy viscosity term for flux component N///!
if(LP(iLayer)%Level.eq.1) then !calculate all grids in top layer 
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY-1,nendy)
elseif(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids on child layers
    istart = MAX(GP%nGHOST+1,nstartx); iend = MIN(LP(iLayer)%NX-GP%nGHOST,nendx)
    jstart = MAX(GP%nGHOST+1,nstarty); jend = MIN(LP(iLayer)%NY-GP%nGHOST-1,nendy)
endif
do j = jstart,jend
do i = istart,iend
    !///exclude the cell where water cannot flow freely///!
    if(LP(iLayer)%Sign_N(i,j).ne.0) cycle
    !/// add eddy viscosity term to N ///!
    n0 = LP(iLayer)%N(1,i,j); v0 = LP(iLayer)%Brkv(i,j); v5 = LP(iLayer)%Brkv(i,j+1)
    if(j.eq.1) then
        n3 = 0.0d0
    else
        n3 = LP(iLayer)%N(1,i,j-1)
    endif
    if(j.eq.LP(iLayer)%NY-1) then
        n5 = 0.0d0
    else
        n5 = LP(iLayer)%N(1,i,j+1)
    endif
    if(i.eq.1) then
        n1 = 0.0d0; v1 = 0.0d0; v7 = 0.0d0
    else
        n1 = LP(iLayer)%N(1,i-1,j)
        v1 = LP(iLayer)%Brkv(i-1,j); v7 = LP(iLayer)%Brkv(i-1,j+1)
    endif
    if(i.eq.LP(iLayer)%NX) then
        n2 = 0.0d0; v2 = 0.0d0; v6 = 0.0d0
    else
        n2 = LP(iLayer)%N(1,i+1,j)
        v2 = LP(iLayer)%Brkv(i+1,j); v6 = LP(iLayer)%Brkv(i+1,j+1)
    endif
    LP(iLayer)%N(2,i,j) = LP(iLayer)%N(2,i,j) + LP(iLayer)%CCB2*(v5*(n5-n0)-v0*(n0-n3)) &
        + 0.25*LP(iLayer)%CCB1*((v0+v2+v5+v6)*(n2-n0)-(v0+v1+v5+v7)*(n0-n1))
enddo
enddo

endif !if:this node is used by this layer

call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1

end subroutine breakingCartesian




subroutine breakingSpherical(GP, LP, iLayer)

use VariableDefination
implicit NONE
type(GlobalParameters)    ::  GP
type(LayerParameters)     ::  LP(100)
integer*4   ::  iLayer, irank
integer*4   ::  nstartx, nendx, nstarty, nendy
integer*4   ::  istart, iend, jstart, jend
integer*4   ::  i, j
real*8      ::  m0, m1, m2, m3, m5
real*8      ::  n0, n1, n2, n3, n5
real*8      ::  csy1, csy2
real*8      ::  v0, v1, v2, v3, v4, v5, v6, v7
real*8      ::  CPUTime1, CPUTime2
    
irank = GP%irank
call CPU_TIME(CPUTime1)
    
if(irank.lt.LP(iLayer)%nsize) then
       
nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
    
!///calculate eddy viscosity term for flux component M///!
if(LP(iLayer)%Level.eq.1) then !calculate all grids in top layer
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX-1,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY,nendy)
elseif(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids in child layer
    istart = MAX(GP%nGHOST+1,nstartx); iend = MIN(LP(iLayer)%NX-GP%nGHOST-1,nendx)
    jstart = MAX(GP%nGHOST+1,nstarty); jend = MIN(LP(iLayer)%NY-GP%nGHOST,nendy)
endif
do j = jstart,jend
do i = istart,iend
    !///exclude the cell where water cannot flow freely///!
    if(LP(iLayer)%Sign_M(i,j).ne.0) cycle
    !/// add eddy viscosity term to M ///!
    m0 = LP(iLayer)%M(1,i,j); v0 = LP(iLayer)%Brkv(i,j); v2 = LP(iLayer)%Brkv(i+1,j)
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
    if(j.eq.1) then
        m3 = 0.0d0; v3 = 0.0d0; v4 = 0.0d0
        csy1 = COS((LP(iLayer)%ymin-0.5*LP(iLayer)%dy)*GP%PI/180.0)
    else
        m3 = LP(iLayer)%M(1,i,j-1)
        v3 = LP(iLayer)%Brkv(i,j-1); v4 = LP(iLayer)%Brkv(i+1,j-1)
        csy1 = LP(iLayer)%CSY(j-1)
    endif
    if(j.eq.LP(iLayer)%NY) then
        m5 = 0.0d0; v5 = 0.0d0; v6 = 0.0d0
        csy2 = COS((LP(iLayer)%ymax+0.5*LP(iLayer)%dy)*GP%PI/180.0)
    else
        m5 = LP(iLayer)%M(1,i,j+1)
        v5 = LP(iLayer)%Brkv(i,j+1); v6 = LP(iLayer)%Brkv(i+1,j+1)
        csy2 = LP(iLayer)%CSY(j)
    endif

    LP(iLayer)%M(2,i,j) = LP(iLayer)%M(2,i,j) + LP(iLayer)%CSB1(j)*(v2*(m2-m0)-v0*(m0-m1)) &
        + 0.25*LP(iLayer)%CSB2(j)*((v0+v2+v5+v6)*(m5-m0)*csy2-(v0+v2+v3+v4)*(m0-m3)*csy1)
enddo
enddo
    
!///calculate eddy viscosity term for flux component N///!
if(LP(iLayer)%Level.eq.1) then !calculate all grids in top layer 
    istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX,nendx)
    jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY-1,nendy)
elseif(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids on child layers
    istart = MAX(GP%nGHOST+1,nstartx); iend = MIN(LP(iLayer)%NX-GP%nGHOST,nendx)
    jstart = MAX(GP%nGHOST+1,nstarty); jend = MIN(LP(iLayer)%NY-GP%nGHOST-1,nendy)
endif
do j = jstart,jend
do i = istart,iend
    !///exclude the cell where water cannot flow freely///!
    if(LP(iLayer)%Sign_N(i,j).ne.0) cycle
    !/// add eddy viscosity term to N ///!
    n0 = LP(iLayer)%N(1,i,j); v0 = LP(iLayer)%Brkv(i,j); v5 = LP(iLayer)%Brkv(i,j+1)
    csy1 = COS(LP(iLayer)%Y(j)*GP%PI/180.0); csy2 = COS(LP(iLayer)%Y(j+1)*GP%PI/180.0)
    if(j.eq.1) then
        n3 = 0.0d0
    else
        n3 = LP(iLayer)%N(1,i,j-1)
    endif
    if(j.eq.LP(iLayer)%NY-1) then
        n5 = 0.0d0
    else
        n5 = LP(iLayer)%N(1,i,j+1)
    endif
    if(i.eq.1) then
        n1 = 0.0d0; v1 = 0.0d0; v7 = 0.0d0
    else
        n1 = LP(iLayer)%N(1,i-1,j)
        v1 = LP(iLayer)%Brkv(i-1,j); v7 = LP(iLayer)%Brkv(i-1,j+1)
    endif
    if(i.eq.LP(iLayer)%NX) then
        n2 = 0.0d0; v2 = 0.0d0; v6 = 0.0d0
    else
        n2 = LP(iLayer)%N(1,i+1,j)
        v2 = LP(iLayer)%Brkv(i+1,j); v6 = LP(iLayer)%Brkv(i+1,j+1)
    endif
    LP(iLayer)%N(2,i,j) = LP(iLayer)%N(2,i,j) + LP(iLayer)%CSB4(j)*(v5*(n5-n0)*csy2-v0*(n0-n3)*csy1) &
        + 0.25*LP(iLayer)%CSB3(j)*((v0+v2+v5+v6)*(n2-n0)-(v0+v1+v5+v7)*(n0-n1))
enddo
enddo
    
endif !if:this node is used by this layer
    
call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1
    
end subroutine breakingSpherical
    



subroutine eddyViscosity(GP, LP, iLayer)

use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer, irank
integer*4   ::  nstartx, nendx, nstarty, nendy, i, j
real*8      ::  Cbrk1, Cbrk2, T_brk
real*8      ::  etat_star, T_star, etat_I, etat_F, B
real*8      ::  ETAt, ETAx, ETAy, D, age, tmp
real*8      ::  h0, h1, h2, h3, h4, z1, z2, z3, z4, lx, ly
real*8      ::  C, angle, age1, age2, age3, propx, propy, propxy, small  
real*8      ::  CPUTime1, CPUTime2

irank = GP%irank
call CPU_TIME(CPUTime1)

Cbrk1 = 0.65; Cbrk2 = 0.15; small = 1.0e-6

if(irank.lt.LP(iLayer)%nsize) then

nstartx = LP(iLayer)%PartitionInfo(irank+1,1); nendx = LP(iLayer)%PartitionInfo(irank+1,2)
nstarty = LP(iLayer)%PartitionInfo(irank+1,3); nendy = LP(iLayer)%PartitionInfo(irank+1,4)
if(LP(iLayer)%Level.gt.1) then !do not calculate GHOST grids in child layers
    nstartx = MAX(nstartx,GP%nGHOST+1); nendx = MIN(nendx,LP(iLayer)%NX-GP%nGHOST)
    nstarty = MAX(nstarty,GP%nGHOST+1); nendy = MIN(nendy,LP(iLayer)%NY-GP%nGHOST)
endif

!/// calculate edd-viscosity coefficients at wet cells ///!
do j = nstarty, nendy
do i = nstartx, nendx
    !/// exclude outermost cells in the top layer ///!
    if(i.eq.1.or.i.eq.LP(iLayer)%NX.or.j.eq.1.or.j.eq.LP(iLayer)%NY) then
        LP(iLayer)%BrkAge(i,j) = 0.0d0
        LP(iLayer)%Brkv(i,j) = 0.0d0
        cycle
    endif
    !/// exclude dry cells ///!
    if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit.or. &
        LP(iLayer)%Z(i,j)+LP(iLayer)%H(2,i,j).le.0.0) then
        LP(iLayer)%BrkAge(i,j) = 0.0d0
        LP(iLayer)%Brkv(i,j) = 0.0d0
        cycle
    endif

    D = MAX(LP(iLayer)%Z(i,j)+LP(iLayer)%H(2,i,j), GP%FrictionDepthLimit)
    T_brk = 10.0*SQRT(MAX(LP(iLayer)%Z(i,j),GP%FrictionDepthLimit)/GP%GRAV)
    tmp = SQRT(GP%GRAV*D); etat_I = Cbrk1*tmp; etat_F = Cbrk2*tmp

    !/// estimate the age of breaking event ///!
    age = LP(iLayer)%BrkAge(i,j)
    ETAt = LP(iLayer)%H(2,i,j) - MAX(LP(iLayer)%H(1,i,j),-LP(iLayer)%Z(i,j))
    ETAt = ETAt/LP(iLayer)%dt

    if(ETAt.ge.etat_I.and.(age.eq.0.0.or.age.gt.T_brk)) then
        LP(iLayer)%BrkAge(i,j) = LP(iLayer)%dt
    elseif(age.gt.0.0) then
        LP(iLayer)%BrkAge(i,j) = LP(iLayer)%BrkAge(i,j) + LP(iLayer)%dt
    else
        h0 = LP(iLayer)%H(2,i,j)
        z1 = LP(iLayer)%Z(i-1,j); h1 = LP(iLayer)%H(2,i-1,j)
        z2 = LP(iLayer)%Z(i+1,j); h2 = LP(iLayer)%H(2,i+1,j)
        z3 = LP(iLayer)%Z(i,j-1); h3 = LP(iLayer)%H(2,i,j-1)
        z4 = LP(iLayer)%Z(i,j+1); h4 = LP(iLayer)%H(2,i,j+1)
        if(GP%CoordinatesType.eq.1) then
            lx = LP(iLayer)%dx; ly = LP(iLayer)%dy
        elseif(GP%CoordinatesType.eq.0) then
            lx = GP%R_EARTH*COS(LP(iLayer)%Y(j)*GP%PI/180.0)*LP(iLayer)%dx*GP%PI/180.0
            ly = GP%R_EARTH*LP(iLayer)%dy*GP%PI/180.0
        endif
        if(z1+h1.le.0.0) h1 = h0
        if(z2+h2.le.0.0) h2 = h0
        if(z3+h3.le.0.0) h3 = h0
        if(z4+h4.le.0.0) h4 = h0
        ETAx = 0.5*(h2-h1)/lx; ETAy = 0.5*(h4-h3)/ly

        tmp = MAX(SQRT(ETAx*ETAx+ETAy*ETAy), small)
        C = MIN(ABS(ETAt)/tmp, SQRT(GP%GRAV*D))
        propx = lx/MAX(C,small)
        propy = ly/MAX(C,small)
        propxy = SQRT(lx*lx+ly*ly)/MAX(C,small)

        if(ETAt.ge.etat_F.and.(ABS(ETAx).gt.small.or.ABS(ETAy).gt.small)) then
            angle = ATAN2(-ETAy,-ETAx)*180.0/GP%PI
            !/// quadrant 1, there are 4 quadrants in total ///!
            if(angle.ge.0.0.and.angle.lt.90.0) then
                age1 = LP(iLayer)%BrkAge(i-1,j)
                age2 = LP(iLayer)%BrkAge(i-1,j-1)
                age3 = LP(iLayer)%BrkAge(i,j-1)
                if((age1.ge.LP(iLayer)%dt.and.age1.gt.propx).or. &
                   (age2.ge.LP(iLayer)%dt.and.age2.gt.propxy).or. &
                   (age3.ge.LP(iLayer)%dt.and.age3.gt.propy)) then
                    LP(iLayer)%BrkAge(i,j) = LP(iLayer)%dt
                endif
            endif
            !/// quadrant 2 ///!
            if(angle.ge.90.0.and.angle.lt.180.0) then
                age1 = LP(iLayer)%BrkAge(i+1,j)
                age2 = LP(iLayer)%BrkAge(i+1,j-1)
                age3 = LP(iLayer)%BrkAge(i,j-1)
                if((age1.ge.LP(iLayer)%dt.and.age1.gt.propx).or. &
                   (age2.ge.LP(iLayer)%dt.and.age2.gt.propxy).or. &
                   (age3.ge.LP(iLayer)%dt.and.age3.gt.propy)) then
                    LP(iLayer)%BrkAge(i,j) = LP(iLayer)%dt
                endif
            endif
            !/// quadrant 3 ///!
            if(angle.ge.-180.0.and.angle.lt.-90.0) then
                age1 = LP(iLayer)%BrkAge(i+1,j)
                age2 = LP(iLayer)%BrkAge(i+1,j+1)
                age3 = LP(iLayer)%BrkAge(i,j+1)
                if((age1.ge.LP(iLayer)%dt.and.age1.gt.propx).or. &
                   (age2.ge.LP(iLayer)%dt.and.age2.gt.propxy).or. &
                   (age3.ge.LP(iLayer)%dt.and.age3.gt.propy)) then
                    LP(iLayer)%BrkAge(i,j) = LP(iLayer)%dt
                endif
            endif
            !/// quadrant 4 ///!
            if(angle.ge.-90.0.and.angle.lt.0.0) then
                age1 = LP(iLayer)%BrkAge(i-1,j)
                age2 = LP(iLayer)%BrkAge(i-1,j+1)
                age3 = LP(iLayer)%BrkAge(i,j+1)
                if((age1.ge.LP(iLayer)%dt.and.age1.gt.propx).or. &
                   (age2.ge.LP(iLayer)%dt.and.age2.gt.propxy).or. &
                   (age3.ge.LP(iLayer)%dt.and.age3.gt.propy)) then
                    LP(iLayer)%BrkAge(i,j) = LP(iLayer)%dt
                endif
            endif
        endif
    endif

    !/// estimate eddy viscosity based on the age of breaking event ///!
    age = LP(iLayer)%BrkAge(i,j)
    if(age.gt.0.0.and.age.lt.T_brk.and.ETAt.gt.etat_F) then
        etat_star = etat_F
        T_star = 5.0*SQRT(MAX(LP(iLayer)%Z(i,j),GP%FrictionDepthLimit)/GP%GRAV)
        if(age.lt.T_star) etat_star = etat_I + age/T_star*(etat_F-etat_I)

        B = 0.0d0
        if(ETAt.gt.etat_star.and.ETAt.le.2.0*etat_star) then
            B = ETAt/etat_star-1
        elseif(ETAt.gt.2.0*etat_star) then
            B = 1.0d0
        endif
        LP(iLayer)%Brkv(i,j) = B*D*ETAt
    else
        LP(iLayer)%Brkv(i,j) = 0.0d0
    endif

enddo
enddo

endif !if: this node is used by this layer

call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1


end subroutine eddyViscosity




subroutine smoothEddyViscosity(GP, LP, iLayer)

use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer, irank
integer*4   ::  nstartx, nendx, nstarty, nendy, i, j
real*8      ::  Brkv_smooth(LP(iLayer)%NX,LP(iLayer)%NY)
real*8      ::  v0, v1, v2, v3, v4, v5, v6, v7, v8
real*8      ::  CPUTime1, CPUTime2

irank = GP%irank
call CPU_TIME(CPUTime1)

if(irank.lt.LP(iLayer)%nsize) then

nstartx = MAX(LP(iLayer)%PartitionInfo(irank+1,1)-1,1)
nendx   = MIN(LP(iLayer)%PartitionInfo(irank+1,2)+1, LP(iLayer)%NX)
nstarty = MAX(LP(iLayer)%PartitionInfo(irank+1,3)-1, 1)
nendy   = MIN(LP(iLayer)%PartitionInfo(irank+1,4)+1, LP(iLayer)%NY)

!/// smooth eddy viscosity with a two-dimensional 9-point filter ///!
do j = nstarty, nendy
do i = nstartx, nendx
    !/// leave the value on the outermost cells unchanged ///!
    if(i.eq.1.or.i.eq.LP(iLayer)%NX.or.j.eq.1.or.j.eq.LP(iLayer)%NY) then
        Brkv_smooth(i,j) = LP(iLayer)%Brkv(i,j); cycle
    endif
    v0 = LP(iLayer)%Brkv(i,j)
    v1 = LP(iLayer)%Brkv(i-1,j);   v2 = LP(iLayer)%Brkv(i+1,j)
    v3 = LP(iLayer)%Brkv(i,j-1);   v4 = LP(iLayer)%Brkv(i,j+1)
    v5 = LP(iLayer)%Brkv(i-1,j-1); v6 = LP(iLayer)%Brkv(i+1,j-1)
    v7 = LP(iLayer)%Brkv(i-1,j+1); v8 = LP(iLayer)%Brkv(i+1,j+1)
    Brkv_smooth(i,j) = v0 + (v1+v2+v3+v4-4.0*v0)/8.0 + (v5+v6+v7+v8-4.0*v0)/16.0
enddo
enddo
do j = nstarty, nendy
do i = nstartx, nendx
    LP(iLayer)%Brkv(i,j) = Brkv_smooth(i,j)
enddo
enddo

endif !if:this node is used by this layer

call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1


end subroutine smoothEddyViscosity



subroutine bcastComputeDomainBoundaryViscosity(GP, LP, iLayer, LocalData, LocalDataLength)

use mpi
use VariableDefination
implicit NONE
type(GlobalParameters)   ::  GP
type(LayerParameters)    ::  LP(100)
integer*4  ::  iLayer, LocalDataLength
real*8     ::  LocalData(LocalDataLength)
integer*4  ::  irank, nsize, master
integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
integer*4  ::  iBC, iVA, i, j
integer*4  ::  nstartx, nendx, nstarty, nendy
real*8     ::  CPUTime1, CPUTime2

irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
call CPU_TIME(CPUTime1)

do iBC = 1,LP(iLayer)%BoundarySendRecvCount
    do iVA = 1,2 !1:Brkv; 2:BrkAge
    if(irank.eq.LP(iLayer)%BoundarySendRecv(iBC, 1)) then
        nstartx = LP(iLayer)%BoundarySendRecv(iBC,3)
        nendx   = LP(iLayer)%BoundarySendRecv(iBC,4)
        nstarty = LP(iLayer)%BoundarySendRecv(iBC,5)
        nendy   = LP(iLayer)%BoundarySendRecv(iBC,6)
        if(nstartx.le.nendx.and.nstarty.le.nendy) then
        do j = nstarty, nendy
        do i = nstartx, nendx
            if(iVA.eq.1) then
                LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%Brkv(i,j)
            else
                LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%BrkAge(i,j)
            endif
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
            if(iVA.eq.1) then
                LP(iLayer)%Brkv(i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
            else
                LP(iLayer)%BrkAge(i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
            endif
        enddo
        enddo
        endif
    endif
    enddo !loop for Brkv and BrkAge
enddo !loop for iBC

call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,4) = GP%CPUTime(irank+1,4)+CPUTime2-CPUTime1


end subroutine bcastComputeDomainBoundaryViscosity



subroutine getLayerBoundaryViscosityFromParent(GP, LP, iLayer, Localdata, LocalDataLength)
    
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
integer*4   ::  istart, iend, jstart, jend, i, j, iVA
real*8      ::  x, y, val
real*8      ::  CPUTime1, CPUTime2

irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master

!/// For each major step iTimeStep, store initial value of viscosity in Brkv0, and final value in BrkvF
!/// For each major step iTimeStep, store initial value of breaking age in BrkAge0, and final value in BrkAgeF
pLayer = LP(iLayer)%Parent
LP(iLayer)%Brkv0 = LP(iLayer)%Brkv; LP(iLayer)%BrkvF = 0.0d0
LP(iLayer)%BrkAge0 = LP(iLayer)%BrkAge; LP(iLayer)%BrkAgeF = 0.0d0
if(LP(pLayer)%Breaking.eq.0) return

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

    if(iHMN.eq.1) then ! use only the table of H for Brkv and BrkAge
        do iVA = 1,2 ! 1:Brkv; 2:BrkAge
        !/// interpolate final value of Brkv or BrkAge on iNodeFrom ///!
        if(irank.eq.iNodeFrom) then
            do j = jstart,jend
            do i = istart,iend
                x = LP(iLayer)%X(i); y = LP(iLayer)%Y(j)
                if(iVA.eq.1) then
                    call interpData(GP,LP,pLayer,5,x,y,val)
                    LP(iLayer)%BrkvF(i,j) = val
                else
                    call interpData(GP,LP,pLayer,6,x,y,val)
                    LP(iLayer)%BrkAgeF(i,j) = val
                endif

                if(iNodeFrom.eq.iNodeTo) then ! modify the value of Brkv or BrkAge
                    if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit.or. &
                        LP(iLayer)%Z(i,j)+LP(iLayer)%HF(i,j).le.0.0) then
                        if(iVA.eq.1) then
                            LP(iLayer)%BrkvF(i,j) = 0.0d0
                        else
                            LP(iLayer)%BrkAgeF(i,j) = 0.0d0
                        endif
                    endif
                endif
            enddo
            enddo
        endif

        !/// send interpolated Brkv or BrkAge from iNodeFrom to iNodeTo ///!
        if(irank.eq.iNodeFrom.and.iNodeTo.ne.iNodeFrom) then
            do j = jstart,jend
            do i = istart,iend
                if(iVA.eq.1) then
                    LocalData((i-istart)*(jend-jstart+1)+j-jstart+1) = LP(iLayer)%BrkvF(i,j)
                else
                    LocalData((i-istart)*(jend-jstart+1)+j-jstart+1) = LP(iLayer)%BrkAgeF(i,j)
                endif
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
                if(iVA.eq.1) then
                    LP(iLayer)%BrkvF(i,j) = LocalData((i-istart)*(jend-jstart+1)+j-jstart+1)
                else
                    LP(iLayer)%BrkAgeF(i,j) = LocalData((i-istart)*(jend-jstart+1)+j-jstart+1)
                endif
                if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit.or. &
                    LP(iLayer)%Z(i,j)+LP(iLayer)%HF(i,j).le.0.0) then
                    if(iVA.eq.1) then
                        LP(iLayer)%BrkvF(i,j) = 0.0d0
                    else
                        LP(iLayer)%BrkAgeF(i,j) = 0.0d0
                    endif
                endif
            enddo
            enddo
        endif

    enddo !loop for Brkv and BrkAge
    endif !if this table is for H
enddo !loop for all ParenToChildSendRecv

call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,4) = GP%CPUTime(irank+1,4)+CPUTime2-CPUTime1


end subroutine getLayerBoundaryViscosityFromParent




subroutine getLayerBoundaryViscosityAtFineTimeStep(GP, LP, iLayer, iStep)
    
use VariableDefination
implicit NONE
type(GlobalParameters)  ::  GP
type(LayerParameters)   ::  LP(100)
integer*4   ::  iLayer, iStep, irank
integer*4   ::  iPC, iNodeFrom, iNodeTo, iBoundary, iHMN
integer*4   ::  istart, iend, jstart, jend, i, j
real*8      ::  t1, t2, t, deltat, age0, agef

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

    if(iHMN.eq.1.and.irank.eq.iNodeTo) then ! only the table of H is used for Brkv and  BrkAge
    do j = jstart,jend
    do i = istart,iend
        LP(iLayer)%Brkv(i,j) = LP(iLayer)%Brkv0(i,j) + &
        (LP(iLayer)%BrkvF(i,j)-LP(iLayer)%Brkv0(i,j))*deltat*(t-t1)
        age0 = LP(iLayer)%BrkAge0(i,j); agef = LP(iLayer)%BrkAgeF(i,j)
        if(age0.gt.agef) age0 = 0.0d0
        LP(iLayer)%BrkAge(i,j) = age0 + (agef-age0)*deltat*(t-t1)
        if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit.or. &
            LP(iLayer)%Z(i,j)+LP(iLayer)%H(2,i,j).le.0.0) then
            LP(iLayer)%Brkv(i,j) = 0.0d0
            LP(iLayer)%BrkAge(i,j) = 0.0d0
        endif
    enddo
    enddo
    endif !if: this table is for H

enddo ! loop for all ParentToChildSendRecv
endif ! if this node is used by this layer


end subroutine getLayerBoundaryViscosityAtFineTimeStep