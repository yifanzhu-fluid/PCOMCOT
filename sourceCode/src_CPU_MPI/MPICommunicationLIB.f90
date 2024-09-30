
    subroutine partitionDomain(GP, LP, SP)
      
    use VariableDefination
    implicit NONE 
    type(GlobalParameters)  ::  GP
    type(LayerParameters)   ::  LP(100)
    type(StationParameters) ::  SP(999)
    integer*4  ::  nsizeTotal, nsize
    integer*4  ::  iLayerLevel, iLayer, pLayer, iSta, i, j, ii, jj
    integer*4  ::  ii1, ii2, jj1, jj2
    integer*4  ::  npartx, nparty, ndiff, ntmp, nlength, nCBoundary
    integer*4  ::  nstartx, nendx, nstarty, nendy
    integer*4  ::  istart, iend, jstart, jend
    integer*4  ::  istartx, iendx, jstarty, jendy
    integer*4  ::  ixy, nBlocks, iBlock, AllBlocks(99,2), IsNewBlock, IsMergeBlock
    integer*4  ::  xpartInfo(999,2), ypartInfo(999,2), partInfo(999,2), partInfotmp(999,2)
    integer*4  ::  npart, istartpart, iendpart 
    integer*4  ::  iMarginx, iMarginy, MaxGrids, MaxGridsOld, MinGrids, AvgGrids
    integer*4  ::  MaxNX, MinNX, MaxNY, MinNY
    integer*4  ::  nCnt, iNode, iCnt, iBoundary, iHMN
    integer*4  ::  iNodeFrom, iNodeTo, HasBoundary, InDomain, IsExist
    real*8     ::  x, y

    write(*,*) 'dividing each computational domain into subdomains ...'
    nsizeTotal = GP%nsizeTotal
    do iLayerLevel = 1,GP%NumLayerLevels
    do iLayer = 1,GP%NumLayers
    if(LP(iLayer)%Level.eq.iLayerLevel) then

        if(GP%MinGridsPerNode.le.0) then
            nsize = nsizeTotal
        else
            nsize = MIN(nsizeTotal,LP(iLayer)%NX*LP(iLayer)%NY/GP%MinGridsPerNode)
            nsize = MAX(1,nsize)
        endif
        npartx = 1; nparty = nsize
        ndiff = ABS(LP(iLayer)%NX/npartx-LP(iLayer)%NY/nparty)
        do i = 1,nsize
            if(MOD(nsize,i).eq.0) then
                ntmp = ABS(LP(iLayer)%NX/i-LP(iLayer)%NY/(nsize/i))
                if(ntmp.lt.ndiff) then
                    npartx = i; nparty = nsize/i; ndiff = ntmp;
                endif
            endif
        enddo
        LP(iLayer)%nsize = nsize
        LP(iLayer)%npartx = npartx; LP(iLayer)%nparty = nparty

        LP(iLayer)%PartitionInfo = -1
        do i=0,nsize-1 ! the layer is divided averagely into all compute nodes
            nlength = LP(iLayer)%NX/npartx
            if(MOD(i,npartx)+1.le.(MOD(LP(iLayer)%NX,npartx))) then
                LP(iLayer)%PartitionInfo(i+1,1) = MOD(i,npartx)*(nlength+1)+1
                LP(iLayer)%PartitionInfo(i+1,2) = LP(iLayer)%PartitionInfo(i+1,1) + nlength
            else
                LP(iLayer)%PartitionInfo(i+1,1) = MOD(LP(iLayer)%NX,npartx)*(nlength+1)+ &
                    (MOD(i,npartx)-MOD(LP(iLayer)%NX,npartx))*nlength + 1
                LP(iLayer)%PartitionInfo(i+1,2) = LP(iLayer)%PartitionInfo(i+1,1)+nlength-1
            endif
            nlength = LP(iLayer)%NY/nparty
            if(i/npartx+1.le.MOD(LP(iLayer)%NY,nparty)) then
                LP(iLayer)%PartitionInfo(i+1,3) = (i/npartx)*(nlength+1) + 1
                LP(iLayer)%PartitionInfo(i+1,4) = LP(iLayer)%PartitionInfo(i+1,3) + nlength
            else
                LP(iLayer)%PartitionInfo(i+1,3) = MOD(LP(iLayer)%NY,nparty)*(nlength+1)+ &
                    (i/npartx-MOD(LP(iLayer)%NY,nparty))*nlength+1
                LP(iLayer)%PartitionInfo(i+1,4) = LP(iLayer)%PartitionInfo(i+1,3)+nlength-1
            endif
        enddo

        ! To save communication time, a child layer is contained in one single node of parent layer!
        ! Partition of parent layers have to be optimized !
        if(GP%ComputeDivisionOpt.eq.2) then !if: child layer on ONE node
	    if(LP(iLayer)%Level.lt.GP%NumLayerLevels) then !if: this layer is a parent layer
            do i = 1,npartx
                xpartInfo(i,1) = LP(iLayer)%PartitionInfo(i,1)
                xpartInfo(i,2) = LP(iLayer)%PartitionInfo(i,2)
            enddo
            do i = 1,nparty
                ypartInfo(i,1) = LP(iLayer)%PartitionInfo((i-1)*npartx+1,3)
                ypartInfo(i,2) = LP(iLayer)%PartitionInfo((i-1)*npartx+1,4)
            enddo
            do ixy = 1,2
                nBlocks = 0; AllBlocks = -1;
                do i = 1,GP%NumLayers
                if(LP(i)%Parent.eq.iLayer) then
                    IsNewBlock = 1
                    if(ixy.eq.1) then
                        istart = NINT((LP(i)%xmin-LP(iLayer)%xmin)/LP(iLayer)%dx)+1
                        iend   = NINT((LP(i)%xmax-LP(iLayer)%xmin)/LP(iLayer)%dx)+1
                    else
                        istart = NINT((LP(i)%ymin-LP(iLayer)%ymin)/LP(iLayer)%dy)+1
                        iend   = NINT((LP(i)%ymax-LP(iLayer)%ymin)/LP(iLayer)%dy)+1
                    endif
                    do iBlock = 1,nBlocks
                        if(istart.ge.AllBlocks(iBlock,1).and.istart.le.AllBlocks(iBlock,2)) then
                            IsNewBlock = 0; AllBlocks(iBlock,2) = MAX(AllBlocks(iBlock,2),iend); exit
                        elseif(iend.ge.AllBlocks(iBlock,1).and.iend.le.AllBlocks(iBlock,2)) then
                            IsNewBlock = 0; AllBlocks(iBlock,1) = MIN(AllBlocks(iBlock,1),istart); exit
                        elseif(istart.lt.AllBlocks(iBlock,1).and.iend.gt.AllBlocks(iBlock,2)) then
                            IsNewBlock = 0; AllBlocks(iBlock,1) = istart; AllBlocks(iBlock,2) = iend; exit
                        endif
                    enddo
                    if(IsNewBlock.eq.1) then
                        nBlocks = nBlocks+1; AllBlocks(nBlocks,1) = istart; ALlBlocks(nBlocks,2) = iend
                    endif
                endif
                enddo
                do i = 1,nBlocks-1
                do j = i+1, nBlocks
                    if(AllBlocks(i,1).gt.AllBlocks(j,1)) then
                        ntmp = AllBlocks(i,1); AllBlocks(i,1) = AllBlocks(j,1); AllBlocks(j,1) = ntmp
                        ntmp = AllBlocks(i,2); AllBlocks(i,2) = AllBlocks(j,2); AllBlocks(j,2) = ntmp
                    endif
                enddo
                enddo
                do
                    IsMergeBlock = 0
                    do iBlock = 1,nBlocks-1
                        if(AllBlocks(iBlock,2).ge.AllBlocks(iBlock+1,1)) then
                            AllBlocks(iBlock,2) = MAX(AllBlocks(iBlock,2),AllBlocks(iBlock+1,2))
                            do i = iBlock+1, nBlocks-1
                                AllBlocks(i,:) = AllBlocks(i+1,:)
                            enddo
                            IsMergeBlock = 1; nBlocks = nBlocks-1; exit
                        endif
                    enddo
                    if(IsMergeBlock.eq.0) exit
                enddo
                if(ixy.eq.1) then
                    npart = npartx; partInfo = xpartInfo;
                else
                    npart = nparty; partInfo = ypartInfo;
                endif
                do iBlock = 1,nBlocks
                    do i = 1,npart
                        if(AllBlocks(iBlock,1)-3.ge.partInfo(i,1).and. &
                            AllBlocks(iBlock,1)-3.le.partInfo(i,2)) istartpart = i
                        if(AllBlocks(iBlock,2)+3.ge.partInfo(i,1).and. &
                                AllBlocks(iBlock,2)+3.le.partInfo(i,2)) iendpart = i
                    enddo
                    if(iendpart.ne.istartpart) then !if this block covers more than one partition, merge these partitions
                        if(AllBlocks(iBlock,1)-partInfo(istartpart,1)-3.le.&
                            partInfo(istartpart,2)-AllBlocks(iBlock,1)+3) then
                            if(AllBlocks(iBlock,2)-partInfo(iendpart,1)+3.le.&
                                partInfo(iendpart,2)-AllBlocks(iBlock,2)-3) then
                                partInfo(istartpart,2) = AllBlocks(iBlock,2)+3
                                partInfo(istartpart+1,1) = partInfo(istartpart,2)+1
                                PartInfo(istartpart+1,2) = partInfo(iendpart,2)
                                do i = iendpart+1,npart
                                    partInfo(istartpart+1+i-iendpart,:) = partInfo(i,:)
                                enddo
                                npart = npart-(iendpart-istartpart-1)
                            else
                                partInfo(istartpart,2) = partInfo(iendpart,2)
                                do i = iendpart+1,npart
                                    partInfo(istartpart+i-iendpart,:) = partInfo(i,:)
                                enddo
                                npart = npart-(iendpart-istartpart)
                            endif
                        else
                            partInfotmp = partInfo
                            partInfo(istartpart,2) = AllBlocks(iBlock,1)-4
                            partInfo(istartpart+1,1) = partInfo(istartpart,2)+1
                            if(AllBlocks(iBlock,2)-partInfotmp(iendpart,1)+3.le.&
                                partInfotmp(iendpart,2)-AllBlocks(iBlock,2)-3) then
                                partInfo(istartpart+1,2) = AllBlocks(iBlock,2)+3
                                partInfo(istartpart+2,1) = partInfo(istartpart+1,2)+1
                                partInfo(istartpart+2,2) = partInfotmp(iendpart,2)
                                do i = iendpart+1, npart
                                    partInfo(istartpart+2+i-iendpart,:) = partInfotmp(i,:)
                                enddo
                                npart = npart-(iendpart-istartpart-2)
                            else
                                partInfo(istartpart+1,2) = partInfotmp(iendpart,2)
                                do i = iendpart+1,npart
                                    partInfo(istartpart+1+i-iendpart,:) = partInfotmp(i,:)
                                enddo
                                npart = npart-(iendpart-istartpart-1)
                            endif
                        endif
                    endif !if this block covers more than one partition
                enddo
                if(ixy.eq.1) then
                    npartx = npart; xpartInfo = partInfo
                else
                    nparty = npart; ypartInfo = partInfo
                endif
            enddo !loop for ixy
                  
            do !loop for removing too small partitions
                MinNX = LP(iLayer)%NX; MinNY = LP(iLayer)%NY
                MaxNX = 0; MaxNY = 0
                do i = 1,npartx
                    MinNX = MIN(MinNX,xpartInfo(i,2)-xpartInfo(i,1)+1)
                    MaxNX = MAX(MaxNX,xpartInfo(i,2)-xpartInfo(i,1)+1)
                enddo
                do i = 1,nparty
                    MinNY = MIN(MinNY,ypartInfo(i,2)-ypartInfo(i,1)+1)
                    MaxNY = MAX(MaxNY,ypartInfo(i,2)-ypartInfo(i,1)+1)
                enddo
                if((MinNX*MinNY.lt.GP%MinGridsPerNode.or.MinNX*MinNY.lt.0.2*MaxNX*MaxNY) &
                    .and.npartx*nparty.gt.1) then !if there are partitions need to and can be merged
                    if((MinNX.le.MinNY.or.nparty.eq.1).and.npartx.gt.1) then
                        do i = 1,npartx
                            if(xpartInfo(i,2)-xpartInfo(i,1)+1.eq.MinNX) then
                                if(i.eq.1) then
                                    xpartInfo(i,2) = xpartInfo(i+1,2)
                                    do j = i+1, npartx-1
                                        xpartInfo(j,:) = xpartInfo(j+1,:)
                                    enddo
                                    npartx = npartx-1; exit
                                elseif(i.eq.LP(iLayer)%NX) then
                                    xpartInfo(i-1,2) = xpartInfo(i,2)
                                    npartx = npartx-1; exit
                                elseif(xpartInfo(i-1,2)-xpartInfo(i-1,1).le. &
                                    xpartInfo(i+1,2)-xpartInfo(i+1,1)) then
                                    xpartInfo(i-1,2) = xpartInfo(i,2)
                                    do j = i,npartx-1
                                        xpartInfo(j,:) = xpartInfo(j+1,:)
                                    enddo
                                    npartx = npartx-1; exit
                                else
                                    xpartInfo(i,2) = xpartInfo(i+1,2)
                                    do j = i+1, npartx-1
                                        xpartInfo(j,:) = xpartInfo(j+1,:)
                                    enddo
                                    npartx = npartx-1; exit
                                endif
                            endif
                        enddo
                    else
                        do i = 1,nparty
                            if(ypartInfo(i,2)-ypartInfo(i,1)+1.eq.MinNY) then
                                if(i.eq.1) then
                                    ypartInfo(i,2) = ypartInfo(i+1,2)
                                    do j = i+1, nparty-1
                                        ypartInfo(j,:) = ypartInfo(j+1,:)
                                    enddo
                                    nparty = nparty-1; exit
                                elseif(i.eq.LP(iLayer)%NY) then
                                    ypartInfo(i-1,2) = ypartInfo(i,2)
                                    nparty = nparty-1; exit
                                elseif(ypartInfo(i-1,2)-ypartInfo(i-1,1).le.&
                                    ypartInfo(i+1,2)-ypartInfo(i+1,1)) then
                                    ypartInfo(i-1,2) = ypartInfo(i,2)
                                    do j = i,nparty-1
                                        ypartInfo(j,:) = ypartInfo(j+1,:)
                                    enddo
                                    nparty = nparty-1; exit
                                else
                                    ypartInfo(i,2) = ypartInfo(i+1,2)
                                    do j = i+1, nparty-1
                                        ypartInfo(j,:) = ypartInfo(j+1,:)
                                    enddo
                                    nparty = nparty-1; exit
                                endif
                            endif
                        enddo
                    endif
                else !if no partitions need to or can be merged
                    exit
                endif
            enddo !loop for removing too small partitions
            nsize = npartx*nparty; LP(iLayer)%nsize = nsize
            LP(iLayer)%npartx = npartx; LP(iLayer)%nparty = nparty
            LP(iLayer)%PartitionInfo = -1
            do j = 1, nparty
            do i = 1, npartx
                LP(iLayer)%PartitionInfo((j-1)*npartx+i,1) = xpartInfo(i,1)
                LP(iLayer)%PartitionInfo((j-1)*npartx+i,2) = xpartInfo(i,2)
                LP(iLayer)%PartitionInfo((j-1)*npartx+i,3) = ypartInfo(j,1)
                LP(iLayer)%PartitionInfo((j-1)*npartx+i,4) = ypartInfo(j,2)
            enddo
            enddo
        endif !if this layer is a parent layer
        endif !if child layer on ONE parent node
    endif
    enddo
    enddo

    !////// calculate MaxNX, MinNX, MaxNY, MinNY //////!
    do iLayer = 1,GP%NumLayers
        LP(iLayer)%MaxNX = 0; LP(iLayer)%MinNX = LP(iLayer)%NX;
        LP(iLayer)%MaxNY = 0; LP(iLayer)%MinNY = LP(iLayer)%NY;
        do i = 1,LP(iLayer)%nsize
            LP(iLayer)%MaxNX = MAX(LP(iLayer)%MaxNX,LP(iLayer)%PartitionInfo(i,2)-LP(iLayer)%PartitionInfo(i,1)+1)
            LP(iLayer)%MinNX = MIN(LP(iLayer)%MinNX,LP(iLayer)%PartitionInfo(i,2)-LP(iLayer)%PartitionInfo(i,1)+1)
            LP(iLayer)%MaxNY = MAX(LP(iLayer)%MaxNY,LP(iLayer)%PartitionInfo(i,4)-LP(iLayer)%PartitionInfo(i,3)+1)
            LP(iLayer)%MinNY = MIN(LP(iLayer)%MinNY,LP(iLayer)%PartitionInfo(i,4)-LP(iLayer)%PartitionInfo(i,3)+1)
        enddo
    enddo
    GP%MaxNX = LP(1)%MaxNX; GP%MinNX = LP(1)%MinNX;
    GP%MaxNY = LP(1)%MaxNY; GP%MinNY = LP(1)%MinNY;
    do iLayer = 2, GP%NumLayers
        GP%MaxNX = MAX(GP%MaxNX, LP(iLayer)%MaxNX); GP%MinNX = MIN(GP%MinNX, LP(iLayer)%MinNX);
        GP%MaxNY = MAX(GP%MaxNY, LP(iLayer)%MaxNY); GP%MinNY = MIN(GP%MinNY, LP(iLayer)%MinNY);
    enddo

    !****** Calculate Boundary Send Receive Table for Each Layer ******!
    do iLayer = 1,GP%NumLayers
        nsize = LP(iLayer)%nsize; npartx = LP(iLayer)%npartx; nparty = LP(iLayer)%nparty
        LP(iLayer)%BoundarySendRecvCount = 0; LP(iLayer)%BoundarySendRecv = -1; nCnt = 0
        do i = 0, nsize-1
            nstartx = LP(iLayer)%PartitionInfo(i+1,1); nendx = LP(iLayer)%PartitionInfo(i+1,2)
            nstarty = LP(iLayer)%PartitionInfo(i+1,3); nendy = LP(iLayer)%PartitionInfo(i+1,4)
            !////// left/right <--> bottom/top must in order to correctly sync coners //////!
            !////// send H/M/N on left boundary of this compute domain //////!
            if(MOD(i,npartx).ne.0) then
                nCnt = nCnt+1
                LP(iLayer)%BoundarySendRecv(nCnt, 1) = i
                LP(iLayer)%BoundarySendRecv(nCnt, 2) = i-1
                LP(iLayer)%BoundarySendRecv(nCnt, 3) = nstartx
                LP(iLayer)%BoundarySendRecv(nCnt, 4) = nstartx+GP%nRowBoundary-1
                LP(iLayer)%BoundarySendRecv(nCnt, 5) = MAX(nstarty-GP%nRowBoundary,1)
                LP(iLayer)%BoundarySendRecv(nCnt, 6) = MIN(nendy+GP%nRowBoundary,LP(iLayer)%NY)
                LP(iLayer)%BoundarySendRecv(nCnt, 7) = nstartx
                LP(iLayer)%BoundarySendRecv(nCnt, 8) = nstartx+GP%nRowBoundaryFlux-1
                LP(iLayer)%BoundarySendRecv(nCnt, 9) = MAX(nstarty-GP%nRowBoundaryFlux,1)
                LP(iLayer)%BoundarySendRecv(nCnt,10) = MIN(nendy+GP%nRowBoundaryFlux,LP(iLayer)%NY)
                LP(iLayer)%BoundarySendRecv(nCnt,11) = nstartx
                LP(iLayer)%BoundarySendRecv(nCnt,12) = nstartx+GP%nRowBoundaryFlux-1
                LP(iLayer)%BoundarySendRecv(nCnt,13) = MAX(nstarty-GP%nRowBoundaryFlux,1)
                LP(iLayer)%BoundarySendRecv(nCnt,14) = MIN(nendy+GP%nRowBoundaryFlux,LP(iLayer)%NY-1)
            endif
        enddo
        do i=0, nsize-1
            nstartx = LP(iLayer)%PartitionInfo(i+1,1); nendx = LP(iLayer)%PartitionInfo(i+1,2)
            nstarty = LP(iLayer)%PartitionInfo(i+1,3); nendy = LP(iLayer)%PartitionInfo(i+1,4)
            !////// send H/M/N on right boundary of this compute domain //////!
            if(MOD(i,npartx).ne.npartx-1) then
                nCnt = nCnt+1
                LP(iLayer)%BoundarySendRecv(nCnt, 1) = i
                LP(iLayer)%BoundarySendRecv(nCnt, 2) = i+1
                LP(iLayer)%BoundarySendRecv(nCnt, 3) = nendx-GP%nRowBoundary+1
                LP(iLayer)%BoundarySendRecv(nCnt, 4) = nendx
                LP(iLayer)%BoundarySendRecv(nCnt, 5) = MAX(nstarty-GP%nRowBoundary,1)
                LP(iLayer)%BoundarySendRecv(nCnt, 6) = MIN(nendy+GP%nRowBoundary,LP(iLayer)%NY)
                LP(iLayer)%BoundarySendRecv(nCnt, 7) = nendx-GP%nRowBoundaryFlux+1
                LP(iLayer)%BoundarySendRecv(nCnt, 8) = nendx
                LP(iLayer)%BoundarySendRecv(nCnt, 9) = MAX(nstarty-GP%nRowBoundaryFlux,1)
                LP(iLayer)%BoundarySendRecv(nCnt,10) = MIN(nendy+GP%nRowBoundaryFlux,LP(iLayer)%NY)
                LP(iLayer)%BoundarySendRecv(nCnt,11) = nendx-GP%nRowBoundaryFlux+1
                LP(iLayer)%BoundarySendRecv(nCnt,12) = nendx
                LP(iLayer)%BoundarySendRecv(nCnt,13) = MAX(nstarty-GP%nRowBoundaryFlux,1)
                LP(iLayer)%BoundarySendRecv(nCnt,14) = MIN(nendy+GP%nRowBoundaryFlux,LP(iLayer)%NY-1)
            endif
        enddo
        do i=0, nsize-1
            nstartx = LP(iLayer)%PartitionInfo(i+1,1); nendx = LP(iLayer)%PartitionInfo(i+1,2)
            nstarty = LP(iLayer)%PartitionInfo(i+1,3); nendy = LP(iLayer)%PartitionInfo(i+1,4)
            !////// send H/M/N on bottom boundary of this compute domain //////!
            if(i/npartx.ne.0) then
                nCnt = nCnt+1
                LP(iLayer)%BoundarySendRecv(nCnt, 1) = i
                LP(iLayer)%BoundarySendRecv(nCnt, 2) = i-npartx
                LP(iLayer)%BoundarySendRecv(nCnt, 3) = MAX(nstartx-GP%nRowBoundary,1)
                LP(iLayer)%BoundarySendRecv(nCnt, 4) = MIN(nendx+GP%nRowBoundary,LP(iLayer)%NX)
                LP(iLayer)%BoundarySendRecv(nCnt, 5) = nstarty
                LP(iLayer)%BoundarySendRecv(nCnt, 6) = nstarty+GP%nRowBoundary-1
                LP(iLayer)%BoundarySendRecv(nCnt, 7) = MAX(nstartx-GP%nRowBoundaryFlux,1)
                LP(iLayer)%BoundarySendRecv(nCnt, 8) = MIN(nendx+GP%nRowBoundaryFlux,LP(iLayer)%NX-1)
                LP(iLayer)%BoundarySendRecv(nCnt, 9) = nstarty
                LP(iLayer)%BoundarySendRecv(nCnt,10) = nstarty+GP%nRowBoundaryFlux-1
                LP(iLayer)%BoundarySendRecv(nCnt,11) = MAX(nstartx-GP%nRowBoundaryFlux,1)
                LP(iLayer)%BoundarySendRecv(nCnt,12) = MIN(nendx+GP%nRowBoundaryFlux,LP(iLayer)%NX)
                LP(iLayer)%BoundarySendRecv(nCnt,13) = nstarty
                LP(iLayer)%BoundarySendRecv(nCnt,14) = nstarty+GP%nRowBoundaryFlux-1
            endif
        enddo
        do i=0, nsize-1
            nstartx = LP(iLayer)%PartitionInfo(i+1,1); nendx = LP(iLayer)%PartitionInfo(i+1,2)
            nstarty = LP(iLayer)%PartitionInfo(i+1,3); nendy = LP(iLayer)%PartitionInfo(i+1,4)
            !////// send H/M/N on top boundary of this compute domain //////!
            if(i/npartx.ne.nparty-1) then
                nCnt = nCnt+1
                LP(iLayer)%BoundarySendRecv(nCnt, 1) = i
                LP(iLayer)%BoundarySendRecv(nCnt, 2) = i+npartx
                LP(iLayer)%BoundarySendRecv(nCnt, 3) = MAX(nstartx-GP%nRowBoundary,1)
                LP(iLayer)%BoundarySendRecv(nCnt, 4) = MIN(nendx+GP%nRowBoundary,LP(iLayer)%NX)
                LP(iLayer)%BoundarySendRecv(nCnt, 5) = nendy-GP%nRowBoundary+1
                LP(iLayer)%BoundarySendRecv(nCnt, 6) = nendy
                LP(iLayer)%BoundarySendRecv(nCnt, 7) = MAX(nstartx-GP%nRowBoundaryFlux,1)
                LP(iLayer)%BoundarySendRecv(nCnt, 8) = MIN(nendx+GP%nRowBoundaryFLux,LP(iLayer)%NX-1)
                LP(iLayer)%BoundarySendRecv(nCnt, 9) = nendy-GP%nRowBoundaryFlux+1
                LP(iLayer)%BoundarySendRecv(nCnt,10) = nendy
                LP(iLayer)%BoundarySendRecv(nCnt,11) = MAX(nstartx-GP%nRowBoundaryFlux,1)
                LP(iLayer)%BoundarySendRecv(nCnt,12) = MIN(nendx+GP%nRowBoundaryFLux,LP(iLayer)%NX)
                LP(iLayer)%BoundarySendRecv(nCnt,13) = nendy-GP%nRowBoundaryFlux+1
                LP(iLayer)%BoundarySendRecv(nCnt,14) = nendy
            endif
        enddo
        LP(iLayer)%BoundarySendRecvCount = nCnt
    enddo

    !****** Calculate Parent to Child Send Receive Table ******!
    nCBoundary = MAX(GP%nGHOST+1, GP%nRowBoundary, GP%nRowBoundaryFlux)
    do iLayer = 1,GP%NumLayers
    if(LP(iLayer)%Level.gt.1) then ! this layer needs to get boundary value from parent
        pLayer = LP(iLayer)%Parent
        LP(iLayer)%ParentToChildSendRecvCount = 0; LP(iLayer)%ParentToChildSendRecv = -1; nCnt = 0
        do iNodeTo = 0,LP(iLayer)%nsize-1
            HasBoundary = 0
            nstartx = LP(iLayer)%PartitionInfo(iNodeTo+1,1)
            nendx   = LP(iLayer)%PartitionInfo(iNodeTo+1,2)
            nstarty = LP(iLayer)%PartitionInfo(iNodeTo+1,3)
            nendy   = LP(iLayer)%PartitionInfo(iNodeTo+1,4)
            if(nstartx.eq.1.or.nendx.eq.LP(iLayer)%NX.or.&
               nstarty.eq.1.or.nendy.eq.LP(iLayer)%NY) HasBoundary = 1

            if(HasBoundary.eq.1) then
            do iBoundary = 1,4   !1/2/3/4 -->left/right/bottom/top
                istart = 0; iend = -1; jstart = 0; jend = -1
                if(nstartx.eq.1.and.iBoundary.eq.1) then
                    istart = 1 
                    iend   = nCBoundary
                    jstart = MAX(nstarty-nCBoundary,1)
                    jend   = MIN(nendy+nCBoundary,LP(iLayer)%NY)
                elseif(nendx.eq.LP(iLayer)%NX.and.iBoundary.eq.2) then
                    istart = LP(iLayer)%NX-nCBoundary+1
                    iend   = LP(iLayer)%NX
                    jstart = MAX(nstarty-nCBoundary,1)
                    jend   = MIN(nendy+nCBoundary,LP(iLayer)%NY)
                elseif(nstarty.eq.1.and.iBoundary.eq.3) then
                    istart = MAX(nstartx-nCBoundary,1)
                    iend   = MIN(nendx+nCBoundary,LP(iLayer)%NX)
                    jstart = 1
                    jend   = nCBoundary
                elseif(nendy.eq.LP(iLayer)%NY.and.iBoundary.eq.4) then
                    istart = MAX(nstartx-nCBoundary,1)
                    iend   = MIN(nendx+nCBoundary,LP(iLayer)%NX)
                    jstart = LP(iLayer)%NY-nCBoundary+1
                    jend   = LP(iLayer)%NY
                endif
                do iHMN = 1,3   !1/2/3-->H/M/N
                    istartx = istart; iendx = iend
                    jstarty = jstart; jendy = jend
                    if(iHMN.eq.2) then
                        iendx = iend - 1
                    elseif(iHMN.eq.3) then
                        jendy = jend - 1
                    endif
                    do j = jstarty, jendy
                    do i = istartx, iendx
                        iNodeFrom = -1; x = LP(iLayer)%X(i); y = LP(iLayer)%Y(j)
                        if(iHMN.eq.2) x = x + 0.5*LP(iLayer)%dx
                        if(iHMN.eq.3) y = y + 0.5*LP(iLayer)%dy
                        ii = FLOOR((x-LP(pLayer)%X(1))/LP(pLayer)%dx)+1
                        jj = FLOOR((y-LP(pLayer)%Y(1))/LP(pLayer)%dy)+1 
                        do iNode = 0, LP(pLayer)%nsize-1
                            ii1 = LP(pLayer)%PartitionInfo(iNode+1,1)
                            ii2 = LP(pLayer)%PartitionInfo(iNode+1,2)
                            jj1 = LP(pLayer)%PartitionInfo(iNode+1,3)
                            jj2 = LP(pLayer)%PartitionInfo(iNode+1,4)
                            if(ii.ge.ii1.and.ii.le.ii2.and.jj.ge.jj1.and.jj.le.jj2) then
                                iNodeFrom = iNode; exit
                            endif
                        enddo
                        IsExist = 0
                        do iCnt = 1,nCnt
                            if (iNodeFrom.eq.LP(iLayer)%ParentToChildSendRecv(iCnt,1).and.&
                                iNodeTo.eq.LP(iLayer)%ParentToChildSendRecv(iCnt,2).and.&
                                iBoundary.eq.LP(iLayer)%ParentToChildSendRecv(iCnt,3).and.&
                                iHMN.eq.LP(iLayer)%ParentToChildSendRecv(iCnt,4)) then
                                IsExist = iCnt; exit
                            endif
                        enddo
                        if (IsExist.eq.0) then
                            nCnt = nCnt + 1
                            LP(iLayer)%ParentToChildSendRecv(nCnt,1) = iNodeFrom
                            LP(iLayer)%ParentToChildSendRecv(nCnt,2) = iNodeTo
                            LP(iLayer)%ParentToChildSendRecv(nCnt,3) = iBoundary
                            LP(iLayer)%ParentToChildSendRecv(nCnt,4) = iHMN
                            LP(iLayer)%ParentToChildSendRecv(nCnt,5) = i
                            LP(iLayer)%ParentToChildSendRecv(nCnt,6) = i
                            LP(iLayer)%ParentToChildSendRecv(nCnt,7) = j
                            LP(iLayer)%ParentToChildSendRecv(nCnt,8) = j
                        else
                            LP(iLayer)%ParentToChildSendRecv(IsExist,6) = &
                                MAX(LP(iLayer)%ParentToChildSendRecv(IsExist,6),i)
                            LP(iLayer)%ParentToChildSendRecv(IsExist,8) = &
                                MAX(LP(iLayer)%ParentToChildSendRecv(IsExist,8),j)
                        endif
                    enddo
                    enddo
                enddo !iHMN = 1,3
            enddo !iBoundary = 1,4
            endif ! if has boundary
        enddo ! iNodeTo of child layer
        LP(iLayer)%ParentToChildSendRecvCount = nCnt
    endif ! if is a child layer
    enddo


    !****** Calculate Child to Parent Send Receive Table ******!
    do iLayer = 1,GP%NumLayers
    if(LP(iLayer)%Level.gt.1) then !this layer needs to send H to its parent
        pLayer = LP(iLayer)%Parent
        LP(iLayer)%ChildToParentSendRecvCount = 0; LP(iLayer)%ChildToParentSendRecv = -1; nCnt = 0
        do iNodeTo = 0, LP(pLayer)%nsize-1
            nstartx = LP(pLayer)%PartitionInfo(iNodeTo+1,1)
            nendx   = LP(pLayer)%PartitionInfo(iNodeTo+1,2)
            nstarty = LP(pLayer)%PartitionInfo(iNodeTo+1,3)
            nendy   = LP(pLayer)%PartitionInfo(iNodeTo+1,4)
            istart = MAX(nstartx-GP%nRowBoundary,1); iend = MIN(nendx+GP%nRowBoundary,LP(pLayer)%NX)
            jstart = MAX(nstarty-GP%nRowBoundary,1); jend = MIN(nendy+GP%nRowBoundary,LP(pLayer)%NY)
            do j = jstart,jend
            do i = istart,iend
                iNodeFrom = -1; InDomain = 0
                x = LP(pLayer)%X(i); y = LP(pLayer)%Y(j)
                if(x.ge.LP(iLayer)%xmin+2.0*LP(iLayer)%dx.and.&
                   x.le.LP(iLayer)%xmax-2.0*LP(iLayer)%dx.and.&
                   y.ge.LP(iLayer)%ymin+2.0*LP(iLayer)%dy.and.&
                   y.le.LP(iLayer)%ymax-2.0*LP(iLayer)%dy) then
                    InDomain = 1
                endif
                if(InDomain.eq.1) then
                ii = FLOOR((x-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
                jj = FLOOR((y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1  
                do iNode = 0,LP(iLayer)%nsize-1
                    if(ii.ge.LP(iLayer)%PartitionInfo(iNode+1,1).and.ii.le.LP(iLayer)%PartitionInfo(iNode+1,2).and.&
                        jj.ge.LP(iLayer)%PartitionInfo(iNode+1,3).and.jj.le.LP(iLayer)%PartitionInfo(iNode+1,4)) then
                        iNodeFrom = iNode; exit
                    endif
                enddo
                IsExist = 0
                do iCnt = 1,nCnt
                    if(iNodeFrom.eq.LP(iLayer)%ChildToParentSendRecv(iCnt,1).and.&
                        iNodeTo.eq.LP(iLayer)%ChildToParentSendRecv(iCnt,2)) then
                        IsExist = iCnt; exit
                    endif
                enddo
                if(IsExist.eq.0) then
                    nCnt = nCnt+1
                    LP(iLayer)%ChildToParentSendRecv(nCnt,1) = iNodeFrom
                    LP(iLayer)%ChildToParentSendRecv(nCnt,2) = iNodeTo
                    LP(iLayer)%ChildToParentSendRecv(nCnt,3) = i
                    LP(iLayer)%ChildToParentSendRecv(nCnt,4) = i 
                    LP(iLayer)%ChildToParentSendRecv(nCnt,5) = j
                    LP(iLayer)%ChildToParentSendRecv(nCnt,6) = j
                else
                    LP(iLayer)%ChildToParentSendRecv(IsExist,4) = MAX(LP(iLayer)%ChildToParentSendRecv(IsExist,4),i)
                    LP(iLayer)%ChildToParentSendRecv(IsExist,6) = MAX(LP(iLayer)%ChildToParentSendRecv(IsExist,6),j)
                endif
                endif !if parent grid in child layer
            enddo
            enddo
        enddo !loop for all parent nodes
        LP(iLayer)%ChildToParentSendRecvCount = nCnt
    endif
    enddo   
                                                          
    !////// Calculate on which node a station is located //////!
    do iSta = 1,GP%NumStations
        iLayer  = SP(iSta)%nLayer; SP(iSta)%nNode = -1
        ii = FLOOR((SP(iSta)%X-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
        jj = FLOOR((SP(iSta)%Y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
        do i=0, LP(iLayer)%nsize-1
            nstartx = LP(iLayer)%PartitionInfo(i+1,1)
            nendx   = LP(iLayer)%PartitionInfo(i+1,2)
            nstarty = LP(iLayer)%PartitionInfo(i+1,3)
            nendy   = LP(iLayer)%PartitionInfo(i+1,4)
            if(ii.ge.nstartx.and.ii.le.nendx.and.jj.ge.nstarty.and.jj.le.nendy) then
                SP(iSta)%nNode = i; exit
            endif
        enddo
    enddo

    !////// Display partition infomation //////!
    do iLayer = 1,GP%NumLayers
        write(*,*)
        write(*,'(a,a11,a)',advance='no') '        ',ADJUSTL(LP(iLayer)%BathymetryFileName),': '
        write(*,'(a,i5,a,i5,a,i10,a,i4,a)') 'NX, NY, NX*NY = ', LP(iLayer)%NX, ',', &
            LP(iLayer)%NY, ',', LP(iLayer)%NX*LP(iLayer)%NY, ',  use', LP(iLayer)%nsize,' nodes...'
        MaxGrids = LP(iLayer)%MaxNX*LP(iLayer)%MaxNY
        MinGrids = LP(iLayer)%MinNX*LP(iLayer)%MinNY
        AvgGrids = NINT((LP(iLayer)%NX*LP(iLayer)%NY*1.0)/(LP(iLayer)%npartx*LP(iLayer)%nparty))
        write(*,'(a,i7,a,i7,a,i7,a,f6.2,a,f6.2,a,f6.2)')           &
            '        Max, Min, Average grids on compute node = ',  &
            MaxGrids, ',', MinGrids, ',', AvgGrids, ',  ',         &
            MaxGrids*1.0/AvgGrids, ',', MinGrids*1.0/AvgGrids, ',', 1.0d0
        if(LP(iLayer)%nsize.ne.nsizeTotal) then
            write(*,'(a,a,a,i4,a,i4,a)') &
                '        WARNING: ',TRIM(ADJUSTL(LP(iLayer)%BathymetryFileName)), &
                ' uses ',LP(iLayer)%nsize,'  out of ',nsizeTotal,'  computing nodes.'
        endif
    enddo
    write(*,*)
    open(99,file=TRIM(ADJUSTL(GP%ResultPath))//'PartitionInfo.dat',form='formatted')
    do iLayer = 1,GP%NumLayers
        do i=0, LP(iLayer)%nsize-1
            write(99,'(9i7)') iLayer,LP(iLayer)%nsize,LP(iLayer)%npartx,LP(iLayer)%nparty, &
                i,(LP(iLayer)%PartitionInfo(i+1,j), j=1,4)
        enddo
    enddo
    close(99)            
    
    end subroutine partitionDomain



    subroutine bcastCommonInfo(GP, LP, SP, FP)

    use mpi 
    use VariableDefination
    implicit NONE
    type(GlobalParameters)   ::  GP
    type(LayerParameters)    ::  LP(100)
    type(StationParameters)  ::  SP(999)
    type(FaultParameters)    ::  FP(4000)
    integer*4  ::  irank, nsize, master
    integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
    integer*4  ::  iLayer, iSta, iFault, LocalData(999*10*14), i, j

    irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
    if(irank.eq.master) then
        write(*,*) 'broadcasting common info from master node ...'
        write(*,*)
    endif

    !/// Broadcast basic control parameters ///!
    call MPI_BCAST(GP%Version,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%PurposeCalculation,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%InitialConditionType,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%CoordinatesType,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%TotalTime,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%dt,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%DTSaveData,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SaveFlux,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SaveDynamicPressure,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%MinGridsPerNode,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%FeedbackToParentLayer,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    !/// Broadcast parameters for Tsunami Source ///!
    call MPI_BCAST(GP%HorizontalMotion,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%KajiuraFilter,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%UseAverageDepth,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%MinKajiuraDepth,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    !/// Broadcast parameters for Wave Physics & Numerics ///!
    call MPI_BCAST(GP%Nonlinearity,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%Dispersion,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%DepthVariability,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%MinDispersionDepth,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%Breaking,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%FluxCenter,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%FrCap,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%DTFilterData,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%NDTFilterData,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    !/// Broadcast Parameters for Boundary  Condition ///!
    call MPI_BCAST(GP%BoundaryConditionType,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SpongeWidthX,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SpongeWidthY,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%MaxSpongeMANNING,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SpongeDampingA,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SpongeDampingR,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    !/// Broadcast parameters for inundation ///!
    call MPI_BCAST(GP%PermanentDryLimit,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%MinWaterDepth,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%FrictionDepthLimit,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%MANNING,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    !/// Broadcast parameters for Green functions ///!    
    call MPI_BCAST(GP%SourceStartX,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SourceEndX,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SourceDX,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SourceStartY,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SourceEndY,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SourceDY,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SourceBasisFunctionType,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%GaussianRatio,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SigmoidCoefficient,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)

    call MPI_BCAST(GP%NumLayers,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      !**** InitialElevationFileName   not broadcasted
      !**** InitialElevationFileFormat not broadcasted
    call MPI_BCAST(GP%StartEastWest,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%NumLayerLevels,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%TopLayer,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%NumStations,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%NumFaults,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%TotalTimeSteps,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%NDTSaveData,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%nCalculations,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SourceNX,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%SourceNY,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)

    call MPI_BCAST(GP%MaxNX,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%MinNX,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%MaxNY,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(GP%MinNY,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)

    do iLayer = 1,GP%NumLayers
        call MPI_BCAST(LP(iLayer)%Level,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%Parent,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%Nonlinearity,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%Dispersion,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror) 
        call MPI_BCAST(LP(iLayer)%DepthVariability,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror) 
        call MPI_BCAST(LP(iLayer)%Breaking,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%FluxCenter,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%ComputingStartTime,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%DTFilterData,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%NDTFilterData,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          !**** BathymetryFileName    not broadcasted
          !**** BathymetryFileFormat  not broadcasted
        call MPI_BCAST(LP(iLayer)%xmin,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%dx,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%xmax,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%ymin,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%dy,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%ymax,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%NX,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%NY,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          !**** X, Y, Z broadcasted after memory allocated
        call MPI_BCAST(LP(iLayer)%zmin,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%zmax,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%dt,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%nStepsPerTimeStep,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%nsize,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%npartx,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%nparty,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%MaxNX,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%MinNX,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%MaxNY,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%MinNY,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        if(irank.eq.master) then
            do i = 1,LP(iLayer)%nsize
                do j = 1,4
                    LocalData((i-1)*4+j) = LP(iLayer)%PartitionInfo(i,j)
                enddo
            enddo
            do i = 0, nsize-1 ! common info is sent to all nodes including unused
                if(i.ne.master)  then
                    call MPI_SEND(LocalData,LP(iLayer)%nsize*4,MPI_INTEGER,i,2015,MPI_COMM_WORLD,ierror)
                endif
            enddo
        else
            call MPI_RECV(LocalData,LP(iLayer)%nsize*4,MPI_INTEGER,master,2015,MPI_COMM_WORLD,istatus,ierror)
            do i = 1,LP(iLayer)%nsize
                do j = 1,4
                    LP(iLayer)%PartitionInfo(i,j) = LocalData((i-1)*4+j)
                enddo
            enddo
        endif
        call MPI_BCAST(LP(iLayer)%BoundarySendRecvCount,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        if(irank.eq.master) then
            do i = 1,LP(iLayer)%BoundarySendRecvCount
                do j = 1,14
                    LocalData((i-1)*14+j) = LP(iLayer)%BoundarySendRecv(i,j)
                enddo
            enddo
            do i = 0, nsize-1
                if(i.ne.master) call MPI_SEND(LocalData, LP(iLayer)%BoundarySendRecvCount*14, &
                    MPI_INTEGER, i, 2015, MPI_COMM_WORLD,ierror)
            enddo
        else
            call MPI_RECV(LocalData, LP(iLayer)%BoundarySendRecvCount*14, &
                MPI_INTEGER, master, 2015, MPI_COMM_WORLD, istatus, ierror)
            do i = 1,LP(iLayer)%BoundarySendRecvCount
                do j = 1,14
                    LP(iLayer)%BoundarySendRecv(i,j) = LocalData((i-1)*14+j)
                enddo
            enddo
        endif
        call MPI_BCAST(LP(iLayer)%ParentToChildSendRecvCount,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        if(irank.eq.master) then
            do i = 1,LP(iLayer)%ParentToChildSendRecvCount
                do j = 1,8
                    LocalData((i-1)*8+j) = LP(iLayer)%ParentToChildSendRecv(i,j)
                enddo
            enddo
            do i = 0, nsize-1
                if(i.ne.master) call MPI_SEND(LocalData, LP(iLayer)%ParentToChildSendRecvCount*8, &
                    MPI_INTEGER, i, 2015, MPI_COMM_WORLD,ierror)
            enddo
        else
            call MPI_RECV(LocalData, LP(iLayer)%ParentToChildSendRecvCount*8, &
                MPI_INTEGER, master, 2015, MPI_COMM_WORLD, istatus, ierror)
            do i = 1,LP(iLayer)%ParentToChildSendRecvCount
                do j = 1,8
                    LP(iLayer)%ParentToChildSendRecv(i,j) = LocalData((i-1)*8+j)
                enddo
            enddo
        endif
        call MPI_BCAST(LP(iLayer)%ChildToParentSendRecvCount,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        if(irank.eq.master) then
            do i = 1,LP(iLayer)%ChildToParentSendRecvCount
                do j = 1,6
                    LocalData((i-1)*6+j) = LP(iLayer)%ChildToParentSendRecv(i,j)
                enddo
            enddo
            do i = 0, nsize-1
                if(i.ne.master) call MPI_SEND(LocalData, LP(iLayer)%ChildToParentSendRecvCount*6, &
                    MPI_INTEGER, i, 2015, MPI_COMM_WORLD,ierror)
            enddo
        else
            call MPI_RECV(LocalData, LP(iLayer)%ChildToParentSendRecvCount*6, &
                MPI_INTEGER, master, 2015, MPI_COMM_WORLD, istatus, ierror)
            do i = 1,LP(iLayer)%ChildToParentSendRecvCount
                do j = 1,6
                    LP(iLayer)%ChildToParentSendRecv(i,j) = LocalData((i-1)*6+j)
                enddo
            enddo
        endif
    enddo

    do iSta = 1,GP%NumStations
        call MPI_BCAST(SP(iSta)%X,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(SP(iSta)%Y,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(SP(iSta)%nLayer,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(SP(iSta)%nNode,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    enddo

    if(GP%InitialConditionType.eq.1) then
        do iFault = 1,GP%NumFaults
            call MPI_BCAST(FP(iFault)%T0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
            call MPI_BCAST(FP(iFault)%NT,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
            call MPI_BCAST(FP(iFault)%Depth,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
            call MPI_BCAST(FP(iFault)%Length,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
            call MPI_BCAST(FP(iFault)%Width,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
            call MPI_BCAST(FP(iFault)%Slip,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
            call MPI_BCAST(FP(iFault)%Rake,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
            call MPI_BCAST(FP(iFault)%HSlip,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
            call MPI_BCAST(FP(iFault)%PSlip,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
            call MPI_BCAST(FP(iFault)%Strike,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
            call MPI_BCAST(FP(iFault)%Dip,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
            call MPI_BCAST(FP(iFault)%X0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
            call MPI_BCAST(FP(iFault)%Y0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        enddo
    endif

    end subroutine bcastCommonInfo



    subroutine bcastBathymetry(GP, LP, SP, FP, LocalData, LocalDataLength)

    use mpi
    use VariableDefination
    implicit NONE
    type(GlobalParameters)   ::  GP
    type(LayerParameters)    ::  LP(100)
    type(StationParameters)  ::  SP(999)
    type(FaultParameters)    ::  FP(4000)
    integer*4  ::  LocalDataLength
    real*8     ::  LocalData(LocalDataLength)
    integer*4  ::  irank, nsize, master
    integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
    integer*4  ::  iLayer, iNode, i, j
    integer*4  ::  istartx, iendx, istarty, iendy, nCBoundary

    irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
    if(irank.eq.master) then
        write(*,*) 'broadcasting bathymetry data from master node ...'
        write(*,*)
    endif

    do iLayer = 1,GP%NumLayers
        call MPI_BCAST(LP(iLayer)%X,LP(iLayer)%NX,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(LP(iLayer)%Y,LP(iLayer)%NY,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    enddo

    nCBoundary = MAX(GP%nGHOST, GP%nRowBoundary+1) ! one extra row/column of water depth for calculating gradient
    if(irank.eq.master) then
        do iLayer = 1,GP%NumLayers
        do iNode = 0, LP(iLayer)%nsize-1
            if(iNode.ne.master) then
                istartx = MAX(1,LP(iLayer)%PartitionInfo(iNode+1,1)-nCBoundary)
                iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(iNode+1,2)+nCBoundary)
                istarty = MAX(1,LP(iLayer)%PartitionInfo(iNode+1,3)-nCBoundary)
                iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(iNode+1,4)+nCBoundary)
                do i = istartx,iendx
                do j = istarty,iendy
                    LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1)) = LP(iLayer)%Z(i,j)
                enddo
                enddo
                call MPI_SEND(LocalData,(iendx-istartx+1)*(iendy-istarty+1), &
                    MPI_DOUBLE_PRECISION,iNode,iLayer,MPI_COMM_WORLD,ierror)
            endif
        enddo
        enddo
    else
        do iLayer = 1,GP%NumLayers
            if(irank.lt.LP(iLayer)%nsize) then
                istartx = MAX(1,LP(iLayer)%PartitionInfo(irank+1,1)-nCBoundary)
                iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(irank+1,2)+nCBoundary)
                istarty = MAX(1,LP(iLayer)%PartitionInfo(irank+1,3)-nCBoundary)
                iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(irank+1,4)+nCBoundary)
                call MPI_RECV(LocalData,(iendx-istartx+1)*(iendy-istarty+1),  &
                    MPI_DOUBLE_PRECISION,master,iLayer,MPI_COMM_WORLD,istatus,ierror)
                do i = istartx,iendx
                do j = istarty,iendy
                    LP(iLayer)%Z(i,j) = LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1))
                enddo
                enddo
            endif
        enddo
    endif
      
    end subroutine bcastBathymetry



    subroutine bcastComputeDomainBoundaryValue(GP, LP, iLayer, LocalData, LocalDataLength)

    use mpi
    use VariableDefination
    implicit NONE
    type(GlobalParameters)   ::  GP
    type(LayerParameters)    ::  LP(100)
    integer*4  ::  iLayer
    integer*4  ::  LocalDataLength
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
                do j=nstarty, nendy
                do i=nstartx, nendx
                    LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%H(2,i,j)
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
                do j=nstarty, nendy
                do i=nstartx, nendx
                    LP(iLayer)%H(2,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                enddo
                enddo
            endif
        endif
    enddo
    call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,4) = GP%CPUTime(irank+1,4)+CPUTime2-CPUTime1

    end subroutine bcastComputeDomainBoundaryValue



    subroutine bcastComputeDomainBoundaryFlux(GP, LP, iLayer, LocalData, LocalDataLength)

    use mpi
    use VariableDefination
    implicit NONE
    type(GlobalParameters)   ::  GP
    type(LayerParameters)    ::  LP(100)
    integer*4  ::  iLayer
    integer*4  ::  LocalDataLength
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
            nstartx = LP(iLayer)%BoundarySendRecv(iBC, 7)
            nendx   = LP(iLayer)%BoundarySendRecv(iBC, 8)
            nstarty = LP(iLayer)%BoundarySendRecv(iBC, 9)
            nendy   = LP(iLayer)%BoundarySendRecv(iBC,10)
            if(nstartx.le.nendx.and.nstarty.le.nendy) then
                do j=nstarty, nendy
                do i=nstartx, nendx
                    LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%M(2,i,j)
                enddo
                enddo
                call MPI_SEND(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1),MPI_DOUBLE_PRECISION, &
                    LP(iLayer)%BoundarySendRecv(iBC,2),iLayer,MPI_COMM_WORLD,ierror)
            endif
            nstartx = LP(iLayer)%BoundarySendRecv(iBC,11)
            nendx   = LP(iLayer)%BoundarySendRecv(iBC,12)
            nstarty = LP(iLayer)%BoundarySendRecv(iBC,13)
            nendy   = LP(iLayer)%BoundarySendRecv(iBC,14)
            if(nstartx.le.nendx.and.nstarty.le.nendy) then
                do j=nstarty, nendy
                do i=nstartx, nendx
                    LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%N(2,i,j)
                enddo
                enddo
                call MPI_SEND(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1),MPI_DOUBLE_PRECISION, &
                    LP(iLayer)%BoundarySendRecv(iBC,2),iLayer,MPI_COMM_WORLD,ierror)
            endif
        elseif(irank.eq.LP(iLayer)%BoundarySendRecv(iBC,2)) then
            nstartx = LP(iLayer)%BoundarySendRecv(iBC, 7)
            nendx   = LP(iLayer)%BoundarySendRecv(iBC, 8)
            nstarty = LP(iLayer)%BoundarySendRecv(iBC, 9)
            nendy   = LP(iLayer)%BoundarySendRecv(iBC,10)
            if(nstartx.le.nendx.and.nstarty.le.nendy) then
                call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1),MPI_DOUBLE_PRECISION, &
                    LP(iLayer)%BoundarySendRecv(iBC,1),iLayer,MPI_COMM_WORLD,istatus,ierror)
                do j=nstarty, nendy
                do i=nstartx, nendx
                    LP(iLayer)%M(2,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                enddo
                enddo
            endif
            nstartx = LP(iLayer)%BoundarySendRecv(iBC,11)
            nendx   = LP(iLayer)%BoundarySendRecv(iBC,12)
            nstarty = LP(iLayer)%BoundarySendRecv(iBC,13)
            nendy   = LP(iLayer)%BoundarySendRecv(iBC,14)
            if(nstartx.le.nendx.and.nstarty.le.nendy) then
                call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1),MPI_DOUBLE_PRECISION, &
                    LP(iLayer)%BoundarySendRecv(iBC,1),iLayer,MPI_COMM_WORLD,istatus,ierror)
                do j=nstarty, nendy
                do i=nstartx, nendx
                    LP(iLayer)%N(2,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                enddo
                enddo
            endif
        endif
    enddo
    call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,4) = GP%CPUTime(irank+1,4)+CPUTime2-CPUTime1

    end subroutine bcastComputeDomainBoundaryFlux



    subroutine getLayerBoundaryFromParent(GP, LP, iLayer, LocalData, LocalDataLength)

    use mpi
    use VariableDefination
    implicit NONE
    type(GlobalParameters)   ::  GP
    type(LayerParameters)    ::  LP(100)
    integer*4  ::  iLayer
    integer*4  ::  LocalDataLength
    real*8     ::  LocalData(LocalDataLength)
    integer*4  ::  irank, nsize, master
    integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
    integer*4  ::  pLayer
    integer*4  ::  iPC, iNodeFrom, iNodeTo, iBoundary, iHMN
    integer*4  ::  istart, iend, jstart, jend, i, j, s
    real*8     ::  x, y, val, z1, z2, h1, h2
    real*8     ::  CPUTime1, CPUTime2

    irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master

    !/// for each major step iTimestep store initial value of H/M/N in H0/M0/N0, and final value in HF/MF/NF
    pLayer = LP(iLayer)%Parent   
    LP(iLayer)%H0 = LP(iLayer)%H(1,:,:); LP(iLayer)%M0 = LP(iLayer)%M(1,:,:); LP(iLayer)%N0 = LP(iLayer)%N(1,:,:)
    LP(iLayer)%HF = 0.0d0; LP(iLayer)%MF = 0.0d0; LP(iLayer)%NF = 0.0d0

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
     
        !///interpolate child layer boundary height and fluxes on iNodeFrom///!
        if(irank.eq.iNodeFrom) then
            do j = jstart,jend
            do i = istart,iend
                x = LP(iLayer)%X(i); y = LP(iLayer)%Y(j)
                if(iHMN.eq.2) x = x + 0.5*LP(iLayer)%dx
                if(iHMN.eq.3) y = y + 0.5*LP(iLayer)%dy
                call interpData(GP,LP,pLayer,iHMN,x,y,val)
                if(iHMN.eq.1) then
                    LP(iLayer)%HF(i,j) = val
                elseif(iHMN.eq.2) then
                    LP(iLayer)%MF(i,j) = val
                elseif(iHMN.eq.3) then
                    LP(iLayer)%NF(i,j) = val
                endif
                if(iNodeFrom.eq.iNodeTo) then
                    if(iHMN.eq.1) then
                        if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit) then
                            LP(iLayer)%HF(i,j) = 0.0d0
                        elseif(LP(iLayer)%Z(i,j)+LP(iLayer)%HF(i,j).le.0.0) then
                            if(LP(iLayer)%Z(i,j).gt.0.0) then
                                LP(iLayer)%HF(i,j) = -LP(iLayer)%Z(i,j)
                            else
                                LP(iLayer)%HF(i,j) = 0.0d0
                            endif
                        endif
                    elseif(iHMN.eq.2) then
                        z1 = LP(iLayer)%Z(i,j);  z2 = LP(iLayer)%Z(i+1,j)
                        h1 = LP(iLayer)%HF(i,j); h2 = LP(iLayer)%HF(i+1,j)
                        call checkFluxDirection(GP%PermanentDryLimit,z1,z2,h1,h2,s)
                        if(s.eq.999.or.s*LP(iLayer)%MF(i,j).lt.0.0) LP(iLayer)%MF(i,j) = 0.0d0
                    elseif(iHMN.eq.3) then
                        z1 = LP(iLayer)%Z(i,j);  z2 = LP(iLayer)%Z(i,j+1)
                        h1 = LP(iLayer)%HF(i,j); h2 = LP(iLayer)%HF(i,j+1)
                        call checkFluxDirection(GP%PermanentDryLimit,z1,z2,h1,h2,s)
                        if(s.eq.999.or.s*LP(iLayer)%NF(i,j).lt.0.0) LP(iLayer)%NF(i,j) = 0.0d0
                    endif
                endif
            enddo
            enddo
        endif

        !///send interpolated value from iNodeFrom to iNodeTo///!
        if(irank.eq.iNodeFrom.and.iNodeTo.ne.iNodeFrom) then
            do j = jstart,jend
            do i = istart,iend
                if(iHMN.eq.1) then
                    LocalData((i-istart)*(jend-jstart+1)+j-jstart+1) = LP(iLayer)%HF(i,j)
                elseif(iHMN.eq.2) then
                    LocalData((i-istart)*(jend-jstart+1)+j-jstart+1) = LP(iLayer)%MF(i,j)
                elseif(iHMN.eq.3) then
                    LocalData((i-istart)*(jend-jstart+1)+j-jstart+1) = LP(iLayer)%NF(i,j)
                endif
            enddo
            enddo
            call MPI_SEND(LocalData,(iend-istart+1)*(jend-jstart+1),MPI_DOUBLE_PRECISION,&
                iNodeTo,2015,MPI_COMM_WORLD,ierror)
        endif

        if(irank.eq.iNodeTo.and.iNodeTo.ne.iNodeFrom) then
            call MPI_RECV(LocalData,(iend-istart+1)*(jend-jstart+1),MPI_DOUBLE_PRECISION,&
                iNodeFrom,2015,MPI_COMM_WORLD,istatus,ierror)
            do j = jstart,jend
            do i = istart,iend
                if(iHMN.eq.1) then
                    LP(iLayer)%HF(i,j) = LocalData((i-istart)*(jend-jstart+1)+j-jstart+1)
                elseif(iHMN.eq.2) then
                    LP(iLayer)%MF(i,j) = LocalData((i-istart)*(jend-jstart+1)+j-jstart+1)
                elseif(iHMN.eq.3) then
                    LP(iLayer)%NF(i,j) = LocalData((i-istart)*(jend-jstart+1)+j-jstart+1)
                endif
                if(iHMN.eq.1) then
                    if(LP(iLayer)%Z(i,j).le.-GP%PermanentDryLimit) then
                        LP(iLayer)%HF(i,j) = 0.0d0
                    elseif(LP(iLayer)%Z(i,j)+LP(iLayer)%HF(i,j).le.0.0) then
                        if(LP(iLayer)%Z(i,j).gt.0.0) then
                            LP(iLayer)%HF(i,j) = -LP(iLayer)%Z(i,j)
                        else
                            LP(iLayer)%HF(i,j) = 0.0d0
                        endif
                    endif
                elseif(iHMN.eq.2) then
                    z1 = LP(iLayer)%Z(i,j);  z2 = LP(iLayer)%Z(i+1,j)
                    h1 = LP(iLayer)%HF(i,j); h2 = LP(iLayer)%HF(i+1,j)
                    call checkFluxDirection(GP%PermanentDryLimit,z1,z2,h1,h2,s)
                    if(s.eq.999.or.s*LP(iLayer)%MF(i,j).lt.0.0) LP(iLayer)%MF(i,j) = 0.0d0
                elseif(iHMN.eq.3) then
                    z1 = LP(iLayer)%Z(i,j);  z2 = LP(iLayer)%Z(i,j+1)
                    h1 = LP(iLayer)%HF(i,j); h2 = LP(iLayer)%HF(i,j+1)
                    call checkFluxDirection(GP%PermanentDryLimit,z1,z2,h1,h2,s)
                    if(s.eq.999.or.s*LP(iLayer)%NF(i,j).lt.0.0) LP(iLayer)%NF(i,j) = 0.0d0
                endif
            enddo
            enddo
        endif

    enddo !loop for all ParentToChildSendRecv
    call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,4) = GP%CPUTime(irank+1,4)+CPUTime2-CPUTime1

    end subroutine getLayerBoundaryFromParent



    subroutine getLayerBoundaryAtFineTimeStep(GP, LP, iLayer, iStep)

    use VariableDefination
    implicit NONE
    type(GlobalParameters)  ::  GP
    type(LayerParameters)   ::  LP(100)
    integer*4   ::  iLayer, iStep, irank
    integer*4   ::  iPC, iNodeFrom, iNodeTo, iBoundary, iHMN
    integer*4   ::  istart, iend, jstart, jend, i, j, s
    real*8      ::  t1, t2, t, deltat
    real*8      ::  z1, z2, h1, h2

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
        
        if(irank.eq.iNodeTo) then
        do j = jstart,jend
        do i = istart,iend
            if(iHMN.eq.1) then
                LP(iLayer)%H(2,i,j) = LP(iLayer)%H0(i,j) + &
                    (LP(iLayer)%HF(i,j)-LP(iLayer)%H0(i,j))*deltat*(t-t1)
                if(LP(iLayer)%H(2,i,j)+LP(iLayer)%Z(i,j).le.0.0) then
                    if(LP(iLayer)%Z(i,j).gt.0.0) then
                        LP(iLayer)%H(2,i,j) = -LP(iLayer)%Z(i,j)
                    else
                        LP(iLayer)%H(2,i,j) = 0.0d0
                    endif
                endif
            elseif(iHMN.eq.2) then
                LP(iLayer)%M(2,i,j) = LP(iLayer)%M0(i,j) + &
                    (LP(iLayer)%MF(i,j)-LP(iLayer)%M0(i,j))*deltat*(t-t1)
                z1 = LP(iLayer)%Z(i,j);   h1 = LP(iLayer)%H(2,i,j)
                z2 = LP(iLayer)%Z(i+1,j); h2 = LP(iLayer)%H(2,i+1,j)
                call checkFluxDirection(GP%PermanentDryLimit, z1, z2, h1, h2, s)
                if(s.eq.999.or.s*LP(iLayer)%M(2,i,j).lt.0.0) LP(iLayer)%M(2,i,j) = 0.0d0
            elseif(iHMN.eq.3) then
                LP(iLayer)%N(2,i,j) = LP(iLayer)%N0(i,j) + &
                    (LP(iLayer)%NF(i,j)-LP(iLayer)%N0(i,j))*deltat*(t-t1)
                z1 = LP(iLayer)%Z(i,j);   h1 = LP(iLayer)%H(2,i,j)
                z2 = LP(iLayer)%Z(i,j+1); h2 = LP(iLayer)%H(2,i,j+1)
                call checkFluxDirection(GP%PermanentDryLimit, z1, z2, h1, h2, s)
                if(s.eq.999.or.s*LP(iLayer)%N(2,i,j).lt.0.0) LP(iLayer)%N(2,i,j) = 0.0d0
            endif
        enddo
        enddo
        endif !if this node has boundary

    enddo !loop for all ParentToChildSendRecv
    endif !if this node is used by this layer

    end subroutine getLayerBoundaryAtFineTimeStep



    subroutine feedbackToParentLayer(GP, LP, iLayer, LocalData, LocalDataLength)

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
    integer*4   ::  iCP, iNodeFrom, iNodeTo, istart, iend, jstart, jend
    integer*4   ::  i, j
    real*8      ::  x, y, h
    real*8      ::  CPUTime1, CPUTime2

    irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
    pLayer = LP(iLayer)%Parent   ! interpolate data at iLayer, feedback to pLayer
    call CPU_TIME(CPUTime1)

    do iCP = 1,LP(iLayer)%ChildToParentSendRecvCount
        iNodeFrom  = LP(iLayer)%ChildToParentSendRecv(iCP,1)
        iNodeTo    = LP(iLayer)%ChildToParentSendRecv(iCP,2)
        istart     = LP(iLayer)%ChildToParentSendRecv(iCP,3)
        iend       = LP(iLayer)%ChildToParentSendRecv(iCP,4)
        jstart     = LP(iLayer)%ChildToParentSendRecv(iCP,5)
        jend       = LP(iLayer)%ChildToParentSendRecv(iCP,6)

        if(irank.eq.iNodeFrom) then
            do j = jstart,jend
            do i = istart,iend
                x = LP(pLayer)%X(i); y = LP(pLayer)%Y(j)
                call interpData(GP, LP, iLayer, 1, x, y, h)
                LP(pLayer)%H(2,i,j) = h
                if(iNodeFrom.eq.iNodeTo) then
                    if(LP(pLayer)%Z(i,j).le.-GP%PermanentDryLimit) then
                        LP(pLayer)%H(2,i,j) = 0.0d0
                    elseif(LP(pLayer)%Z(i,j)+LP(pLayer)%H(2,i,j).le.0.0) then
                        if(LP(pLayer)%Z(i,j).gt.0.0) then
                            LP(pLayer)%H(2,i,j) = -LP(pLayer)%Z(i,j)
                        else
                            LP(pLayer)%H(2,i,j) = 0.0d0
                        endif
                    endif
                endif
            enddo
            enddo
        endif

        if(irank.eq.iNodeFrom.and.iNodeFrom.ne.iNodeTo) then
            do j = jstart,jend
            do i = istart,iend
                LocalData((i-istart)*(jend-jstart+1)+j-jstart+1) = LP(pLayer)%H(2,i,j)
            enddo
            enddo
            call MPI_SEND(LocalData,(iend-istart+1)*(jend-jstart+1), &
                MPI_DOUBLE_PRECISION,iNodeTo,2015,MPI_COMM_WORLD,ierror)
        endif

        if(irank.eq.iNodeTo.and.iNodeFrom.ne.iNodeTo) then
            call MPI_RECV(LocalData,(iend-istart+1)*(jend-jstart+1), &
                MPI_DOUBLE_PRECISION,iNodeFrom,2015,MPI_COMM_WORLD,istatus,ierror)
            do j = jstart,jend
            do i = istart,iend
                LP(pLayer)%H(2,i,j) = LocalData((i-istart)*(jend-jstart+1)+j-jstart+1)
                if(LP(pLayer)%Z(i,j).le.-GP%PermanentDryLimit) then
                    LP(pLayer)%H(2,i,j) = 0.0d0
                elseif(LP(pLayer)%Z(i,j)+LP(pLayer)%H(2,i,j).le.0.0) then
                    if(LP(pLayer)%Z(i,j).gt.0.0) then
                        LP(pLayer)%H(2,i,j) = -LP(pLayer)%Z(i,j)
                    else
                        LP(pLayer)%H(2,i,j) = 0.0d0
                    endif
                endif
            enddo
            enddo
        endif
    enddo !loop for all ChildToParent Table
    call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,4) = GP%CPUTime(irank+1,4)+CPUTime2-CPUTime1

    end subroutine feedbackToParentLayer



    subroutine feedbackToParentLayerByAverage(GP, LP, iLayer, LocalData, LocalDataLength)
    
    use mpi
    use VariableDefination
    implicit NONE
    type(GlobalParameters)   ::  GP
    type(LayerParameters)    ::  LP(100)
    integer*4  ::  iLayer
    integer*4  ::  LocalDataLength
    real*8     ::  LocalData(LocalDataLength)
    integer*4  ::  irank, nsize, master
    integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
    integer*4  ::  pLayer
    integer*4  ::  iCP, iNodeFrom, iNodeTo, istart, iend, jstart, jend
    integer*4  ::  i, j
    real*8     ::  x, y, h
    real*8     ::  CPUTime1, CPUTime2

    irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
    pLayer = LP(iLayer)%Parent   ! interpolate data at iLayer, feedback to pLayer
    call CPU_TIME(CPUTime1)

    do iCP = 1,LP(iLayer)%ChildToParentSendRecvCount
        iNodeFrom  = LP(iLayer)%ChildToParentSendRecv(iCP,1)
        iNodeTo    = LP(iLayer)%ChildToParentSendRecv(iCP,2)
        istart     = LP(iLayer)%ChildToParentSendRecv(iCP,3)
        iend       = LP(iLayer)%ChildToParentSendRecv(iCP,4)
        jstart     = LP(iLayer)%ChildToParentSendRecv(iCP,5)
        jend       = LP(iLayer)%ChildToParentSendRecv(iCP,6)

        if(irank.eq.iNodeFrom) then
            do j = jstart,jend
            do i = istart,iend
                x = LP(pLayer)%X(i); y = LP(pLayer)%Y(j)
                call interpDataByAverage(GP, LP, iLayer, 1, x, y, h)
                LP(pLayer)%H(2,i,j) = h
                if(iNodeFrom.eq.iNodeTo) then
                    if(LP(pLayer)%Z(i,j).le.-GP%PermanentDryLimit) then
                        LP(pLayer)%H(2,i,j) = 0.0d0
                    elseif(LP(pLayer)%Z(i,j)+LP(pLayer)%H(2,i,j).le.0.0) then
                        if(LP(pLayer)%Z(i,j).gt.0.0) then
                            LP(pLayer)%H(2,i,j) = -LP(pLayer)%Z(i,j)
                        else
                            LP(pLayer)%H(2,i,j) = 0.0d0
                        endif
                    endif
                endif
            enddo
            enddo
        endif

        if(irank.eq.iNodeFrom.and.iNodeFrom.ne.iNodeTo) then
            do j = jstart,jend
            do i = istart,iend
                LocalData((i-istart)*(jend-jstart+1)+j-jstart+1) = LP(pLayer)%H(2,i,j)
            enddo
            enddo
            call MPI_SEND(LocalData,(iend-istart+1)*(jend-jstart+1), &
                MPI_DOUBLE_PRECISION,iNodeTo,2015,MPI_COMM_WORLD,ierror)
        endif

        if(irank.eq.iNodeTo.and.iNodeFrom.ne.iNodeTo) then
            call MPI_RECV(LocalData,(iend-istart+1)*(jend-jstart+1), &
                MPI_DOUBLE_PRECISION,iNodeFrom,2015,MPI_COMM_WORLD,istatus,ierror)
            do j = jstart,jend
            do i = istart,iend
                LP(pLayer)%H(2,i,j) = LocalData((i-istart)*(jend-jstart+1)+j-jstart+1)
                if(LP(pLayer)%Z(i,j).le.-GP%PermanentDryLimit) then
                    LP(pLayer)%H(2,i,j) = 0.0d0
                elseif(LP(pLayer)%Z(i,j)+LP(pLayer)%H(2,i,j).le.0.0) then
                    if(LP(pLayer)%Z(i,j).gt.0.0) then
                        LP(pLayer)%H(2,i,j) = -LP(pLayer)%Z(i,j)
                    else
                        LP(pLayer)%H(2,i,j) = 0.0d0
                    endif
                endif
            enddo
            enddo
        endif
    enddo !loop for all ChildToParent Table
    call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,4) = GP%CPUTime(irank+1,4)+CPUTime2-CPUTime1

    end subroutine feedbackToParentLayerByAverage
