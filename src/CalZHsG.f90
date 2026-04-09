subroutine CalZHG(nx, ny, nz, vsf, goxd, gozd, dvxd, dvzd, &
                  nZH, zhx, zhz, izhx, izhz, depz, minthk, &
                  nperZH, perZH, dallZH, adaptive_radius_min, &
                  dsyn, maxnar, &
                  nar, row, col, Gvar)

   implicit none

   integer, parameter :: maxnperZH = 400
   real,    parameter :: ftol      = 1.0e-4

   integer, intent(in) :: nx, ny, nz, nZH, maxnar, dallZH
   real,    intent(in) :: goxd, gozd, dvxd, dvzd, minthk

   real,    intent(in) :: vsf(nx, ny, nz)
   real,    intent(in) :: depz(nz)
   integer, intent(in) :: nperZH(nZH)
   real,    intent(in) :: perZH(nZH, maxnperZH)
   real,    intent(in) :: zhx(nZH), zhz(nZH)
   integer, intent(in) :: izhx(nZH), izhz(nZH)
   real,    intent(in) :: adaptive_radius_min
   real,    intent(out) :: dsyn(dallZH)
   integer, intent(out) :: nar
   integer, intent(out) :: row(maxnar), col(maxnar)
   real,    intent(out) :: Gvar(maxnar)

   integer :: nvx, nvy, nvz, nfull
   integer :: iZH, iDat, kk, k
   integer :: iz, ix, iy, iwt
   integer :: ix_full, iy_full, iz_full, k_full
   integer :: nWeight
   integer :: ip
   integer :: izh_dat
   integer :: ix_in, iy_in, idx_in

   real    :: grid_spacing, Gij
   real    :: radii(nZH)
   real    :: vels(nz-1)
   real    :: vs_loc, vp_loc
   real    :: coef_vp, coef_rho_vp, coef_rho

   real    :: sen_vpZH(dallZH, nz-1)
   real    :: sen_vsZH(dallZH, nz-1)
   real    :: sen_rhoZH(dallZH, nz-1)
   real    :: pvZH_tmp(dallZH)

   integer, allocatable :: xUse(:), zUse(:)
   real,    allocatable :: Wcoef(:)

   integer :: idx2izh(dallZH)
   real, allocatable :: Tmatrix(:,:,:)
   real, allocatable :: vels_sta(:,:)

   real    :: sumW, w
   integer :: x0, y0

   nvx   = nx - 2
   nvy   = ny - 2
   nvz   = nz - 1
   nfull = nx * ny * nz

   allocate(Tmatrix(nZH, nvz, nvx*nvy))
   allocate(vels_sta(nZH, nvz))
   Tmatrix  = 0.0
   vels_sta = 0.0

   allocate(xUse(nx*ny))
   allocate(zUse(nx*ny))
   allocate(Wcoef(nx*ny))

   row  = 0
   col  = 0
   Gvar = 0.0
   dsyn = 0.0

   grid_spacing = 0.5 * (dvxd + dvzd)

   ip = 0
   do kk = 1, nZH
      do k = 1, nperZH(kk)
         ip = ip + 1
         if (ip > dallZH) then
            print *, 'ERROR: ip > dallZH in index mapping'
            stop
         end if
         idx2izh(ip) = kk
      end do
   end do

   if (ip /= dallZH) then
      print *, 'ERROR: sum(nperZH) != dallZH', ip, dallZH
      stop
   end if

   call calculate_adaptive_radius(zhx, zhz, nZH, grid_spacing,adaptive_radius_min, radii)
   call depthkernelZH(nx, ny, nz, vsf, goxd, gozd, dvxd, dvzd, &
                      nZH, nperZH, perZH, dallZH, zhx, zhz, izhx, izhz, &
                      depz, minthk, radii, &
                      pvZH_tmp, sen_vsZH, sen_vpZH, sen_rhoZH)
   dsyn(1:dallZH) = pvZH_tmp(1:dallZH)

   do iZH = 1, nZH
      call weight(radii(iZH), nZH, zhx, zhz, dvxd, dvzd, goxd, gozd, &
                  nx, ny, nWeight, Wcoef, xUse, zUse, iZH)

      if (nWeight < 1 .or. nWeight > nx*ny) then
         print *, 'ERROR: nWeight out of range: ', nWeight, ' station=', iZH
         stop
      end if

      do iz = 1, nvz

         vels(iz) = 0.0
         sumW     = 0.0

         do iwt = 1, nWeight
            x0 = xUse(iwt)
            y0 = zUse(iwt)

            ix_full = x0 + 1
            iy_full = y0 + 1

            if (ix_full >= 2 .and. ix_full <= nx-1 .and. &
                iy_full >= 2 .and. iy_full <= ny-1) then
               sumW = sumW + Wcoef(iwt)
            end if
         end do

         if (sumW <= 0.0) then
            print *, 'ERROR: sumW <= 0 (no interior weights). station=', iZH, ' iz=', iz
            stop
         end if

         do iwt = 1, nWeight
            x0 = xUse(iwt)
            y0 = zUse(iwt)

            ix_full = x0 + 1
            iy_full = y0 + 1

            if (ix_full >= 2 .and. ix_full <= nx-1 .and. &
                iy_full >= 2 .and. iy_full <= ny-1) then

               w = Wcoef(iwt) / sumW

               ix_in  = ix_full - 1
               iy_in  = iy_full - 1
               idx_in = (iy_in-1)*nvx + ix_in

               if (idx_in < 1 .or. idx_in > nvx*nvy) then
                  print *, 'ERROR: idx_in out of range: ', idx_in, ' station=', iZH
                  stop
               end if

               vels(iz) = vels(iz) + w * vsf(ix_full, iy_full, iz)
               Tmatrix(iZH, iz, idx_in) = Tmatrix(iZH, iz, idx_in) + w

            end if
         end do

         vels_sta(iZH, iz) = vels(iz)

      end do

   end do

   nar = 0

   do iDat = 1, dallZH

      izh_dat = idx2izh(iDat)
      if (izh_dat < 1 .or. izh_dat > nZH) then
         print *, 'ERROR: invalid izh_dat: ', izh_dat, ' iDat=', iDat
         stop
      end if

      do iz = 1, nvz

         vs_loc = vels_sta(izh_dat, iz)

         vp_loc = 0.9409 + 2.0947*vs_loc - 0.8206*vs_loc**2 + &
                  0.2683*vs_loc**3 - 0.0251*vs_loc**4

         coef_vp = 2.0947 - 1.6412*vs_loc + &
                   0.8049*vs_loc**2 - 0.1004*vs_loc**3

         coef_rho_vp = 1.6612 - 0.9442*vp_loc + &
                       0.2013*vp_loc**2 - 0.0172*vp_loc**3 + &
                       0.00053*vp_loc**4

         coef_rho = coef_rho_vp * coef_vp

         do iy = 1, nvy
            do ix = 1, nvx

               ix_full = ix + 1
               iy_full = iy + 1
               iz_full = iz

               k_full  = ix_full + (iy_full-1)*nx + (iz_full-1)*nx*ny
               if (k_full < 1 .or. k_full > nfull) cycle

               idx_in = (iy-1)*nvx + ix

               Gij = ( sen_vsZH(iDat, iz) + &
                       coef_vp  * sen_vpZH(iDat, iz) + &
                       coef_rho * sen_rhoZH(iDat, iz) ) * Tmatrix(izh_dat, iz, idx_in)

               if (abs(Gij) > ftol) then
                  nar = nar + 1
                  if (nar > maxnar) then
                     print *, 'ERROR: nar exceeds maxnar'
                     stop
                  end if
                  row(nar)  = iDat
                  col(nar)  = k_full
                  Gvar(nar) = Gij
               end if

            end do
         end do

      end do

   end do

   deallocate(xUse, zUse, Wcoef)
   deallocate(Tmatrix, vels_sta)

end subroutine CalZHG 

subroutine depthkernelZH(nx, ny, nz, vel, goxdf, goydf, dvxdf, dvydf, &
                         nZH, nperZH, perZH, dallZH, zhx, zhy, izhx, izhy, &
                         depz, minthk, radii, &
                         pvZH, sen_vsZH, sen_vpZH, sen_rhoZH)

   use omp_lib
   use utils
   implicit none

   !------------------------------
   ! 参数定义
   !------------------------------
   integer, parameter :: maxnperZH = 400
   integer, parameter :: NL_base = 500
   integer, parameter :: NPERIOD = 400

   integer :: nx, ny, nz, nZH, dallZH
   real    :: vel(nx, ny, nz)
   real    :: goxdf, goydf, dvxdf, dvydf
   real    :: zhx(nZH), zhy(nZH), depz(nz)
   real    :: minthk, wavelength
   integer :: izhx(nZH), izhy(nZH)
   integer nperZH(nZH)
   real    :: perZH(nZH, maxnperZH)
   real    :: radii(nZH)

   ! 输出
   real    :: pvZH(dallZH)
   real    :: sen_vsZH(dallZH, nz-1), sen_vpZH(dallZH, nz-1), sen_rhoZH(dallZH, nz-1)

   ! 内部常量/临时
   integer :: nvx, nvy, nvz, mmax
   integer :: ii, iz, iwt, k, kk, i, nn, j
   integer :: nWeight, start_index, rmax
   real    :: dlnVs, dlnVp, dlnrho

   ! 注意：下列数组在并行区域每个线程都要有独立拷贝
   real    :: vsz(nz), vpz(nz), rhoz(nz)
   real    :: vsz_temp(nz), vpz_temp(nz), rhoz_temp(nz)
   real    :: depm(NL_base), vsm(NL_base), vpm(NL_base), rhom(NL_base), thkm(NL_base)
   real    :: rdep(NL_base), rvp(NL_base), rvs(NL_base), rrho(NL_base), rthk(NL_base)
   real    :: dfdvs(NPERIOD, nz-1), dfdvp(NPERIOD, nz-1), dfdr(NPERIOD, nz-1)
   real*8  :: ZH(maxnperZH), ZH2(maxnperZH)
   integer :: xUse(nx*ny), zUse(nx*ny)
   real    :: Wcoef(nx*ny)

   ! 新增：内网格过滤后的权重与索引
   integer :: xUse_in(nx*ny), zUse_in(nx*ny)
   real    :: Wcoef_in(nx*ny)
   integer :: nUse
   integer :: ix_full, iy_full
   real    :: sumW, w
   integer :: nper0(nZH)
   integer :: sum0
   logical :: dbg
   dbg = .false.

   ! 初始化
   nvx = nx - 2
   nvy = ny - 2
   nvz = nz - 1
   mmax = nz - 1
   dlnVs = 0.01
   dlnVp = 0.01
   dlnrho = 0.01

   ! 清零输出
   pvZH = 0.0
   sen_vsZH = 0.0
   sen_vpZH = 0.0
   sen_rhoZH = 0.0
   nper0 = nperZH


   ! 并行化：对每个台站并行
   !$omp parallel do default(shared) &
   !$omp& private(ii,iz,iwt,k,kk,i,nn,j,nWeight,start_index,rmax, &
   !$omp&        vsz,vpz,rhoz,vsz_temp,vpz_temp,rhoz_temp, &
   !$omp&        depm,vsm,vpm,rhom,thkm, &
   !$omp&        rdep,rvp,rvs,rrho,rthk, &
   !$omp&        dfdvs,dfdvp,dfdr,ZH,ZH2, &
   !$omp&        xUse,zUse,Wcoef, &
   !$omp&        xUse_in,zUse_in,Wcoef_in,nUse,ix_full,iy_full,sumW,w) &
   !$omp& shared(nx,ny,nz,nZH,vel,goxdf,goydf,dvxdf,dvydf,zhx,zhy,izhx,izhy, &
   !$omp&        nperZH,perZH,depz,minthk,radii,pvZH,sen_vsZH,sen_vpZH,sen_rhoZH,dbg)
   do ii = 1, nZH

      ! 计算在全局 pv/sen 中的起始位置
      if (ii == 1) then
         start_index = 1
      else
         start_index = sum(nperZH(1:ii-1)) + 1
      end if

      !----------------------------------------
      ! 1) 原始权重
      !----------------------------------------
      call weight(radii(ii), nZH, zhx, zhy, dvxdf, dvydf, goxdf, goydf, &
                  nx, ny, nWeight, Wcoef, xUse, zUse, ii)


      if (nWeight < 1 .or. nWeight > nx*ny) then
         print *, 'ERROR: nWeight out of range in depthkernelZH:', nWeight, ' station=', ii
         stop
      end if

      !----------------------------------------
      ! 2) 过滤边界点，并对内点重新归一化
      !    与 CalZHG 完全一致
      !----------------------------------------
      nUse = 0
      sumW = 0.0

      do iwt = 1, nWeight
         ix_full = xUse(iwt) + 1
         iy_full = zUse(iwt) + 1

         if (ix_full >= 2 .and. ix_full <= nx-1 .and. &
             iy_full >= 2 .and. iy_full <= ny-1) then
            nUse = nUse + 1
            xUse_in(nUse) = xUse(iwt)
            zUse_in(nUse) = zUse(iwt)
            Wcoef_in(nUse) = Wcoef(iwt)
            sumW = sumW + Wcoef(iwt)
         end if
      end do

      if (nUse < 1 .or. sumW <= 0.0) then
         print *, 'ERROR: no interior weights in depthkernelZH. station=', ii
         stop
      end if

      do iwt = 1, nUse
         Wcoef_in(iwt) = Wcoef_in(iwt) / sumW
      end do

      !----------------------------------------
      ! 3) 用“内网格重归一化后的权重”构造局部 1D 模型
      !----------------------------------------
      vsz = 0.0
      do iz = 1, nz
         do iwt = 1, nUse
            vsz(iz) = vsz(iz) + Wcoef_in(iwt) * vel(xUse_in(iwt)+1, zUse_in(iwt)+1, iz)
         end do
      end do

      !----------------------------------------
      ! 4) 经验关系计算 Vp, rho
      !----------------------------------------
      do k = 1, nz
         vpz(k)  = 0.9409 + 2.0947*vsz(k) - 0.8206*vsz(k)**2 + 0.2683*vsz(k)**3 - 0.0251*vsz(k)**4
         rhoz(k) = 1.6612*vpz(k) - 0.4721*vpz(k)**2 + 0.0671*vpz(k)**3 - 0.0043*vpz(k)**4 + 0.000106*vpz(k)**5
      end do

      ! 分层模型
      call refineGrid2LayerMdl(minthk, nz, depz, vpz, vsz, rhoz, rmax, rdep, rvp, rvs, rrho, rthk)
      if (any(nperZH /= nper0)) then
         print *, 'ERROR: nperZH corrupted after refine(base), station=', ii
         print *, 'old/new sums =', sum0, sum(nperZH)
         stop
      end if

      call zhFwd(rthk, rvp, rvs, rrho, rmax, 1, nperZH(ii), perZH(ii,1:nperZH(ii)), ZH)
  

      if (any(nperZH /= nper0)) then
         print *, 'ERROR: nperZH corrupted after zhFwd(base), station=', ii
         print *, 'old/new sums =', sum0, sum(nperZH)
         stop
      end if

      ! 存储预测值
      do j = 1, nperZH(ii)
         pvZH(start_index + j - 1) = ZH(j)
      end do

      !----------------------------------------
      ! 5) 构建用于灵敏度计算的分层模型
      !----------------------------------------
      do kk = 1, mmax
         depm(kk) = depz(kk)
         vsm(kk)  = vsz(kk)
         vpm(kk)  = vpz(kk)
         rhom(kk) = rhoz(kk)
         thkm(kk) = depz(kk+1) - depz(kk)
      end do
      depm(mmax+1) = depz(mmax+1)

      !----------------------------------------
      ! 6) 灵敏度计算
      !----------------------------------------
      do i = 1, mmax
         ! Vs 灵敏度
         vsz_temp = vsz
         vsz_temp(i) = vsz(i) * (1.0 + dlnVs)
         call refineGrid2LayerMdl(minthk, nz, depz, vpz, vsz_temp, rhoz, rmax, rdep, rvp, rvs, rrho, rthk)
         ! Vs
         call zhFwd(rthk, rvp, rvs, rrho, rmax, 1, nperZH(ii), perZH(ii,1:nperZH(ii)), ZH2)
   
         do nn = 1, nperZH(ii)
            dfdvs(nn, i) = (ZH2(nn) - ZH(nn)) / (dlnVs * vsz(i))
         end do

         ! Vp 灵敏度
         vpz_temp = vpz
         vpz_temp(i) = vpz(i) * (1.0 + dlnVp)
         call refineGrid2LayerMdl(minthk, nz, depz, vpz_temp, vsz, rhoz, rmax, rdep, rvp, rvs, rrho, rthk)
         call zhFwd(rthk, rvp, rvs, rrho, rmax, 1, nperZH(ii), perZH(ii,1:nperZH(ii)), ZH2)
         do nn = 1, nperZH(ii)
            dfdvp(nn, i) = (ZH2(nn) - ZH(nn)) / (dlnVp * vpz(i))
         end do

         ! Rho 灵敏度
         rhoz_temp = rhoz
         rhoz_temp(i) = rhoz(i) * (1.0 + dlnrho)
         call refineGrid2LayerMdl(minthk, nz, depz, vpz, vsz, rhoz_temp, rmax, rdep, rvp, rvs, rrho, rthk)
         call zhFwd(rthk, rvp, rvs, rrho, rmax, 1, nperZH(ii), perZH(ii,1:nperZH(ii)), ZH2)
         do nn = 1, nperZH(ii)
            dfdr(nn, i) = (ZH2(nn) - ZH(nn)) / (dlnrho * rhoz(i))
         end do
      end do

      ! 存储灵敏度
      do nn = 1, nperZH(ii)
         sen_vsZH(start_index + nn - 1, 1:mmax) = dfdvs(nn, 1:mmax)
         sen_vpZH(start_index + nn - 1, 1:mmax) = dfdvp(nn, 1:mmax)
         sen_rhoZH(start_index + nn - 1, 1:mmax) = dfdr(nn, 1:mmax)
      end do

   end do
   !$omp end parallel do

end subroutine depthkernelZH

subroutine calculate_adaptive_radius(station_lons, station_lats, nstations, &
                                    grid_spacing, min_physical_radius, radii)
    implicit none
    
    ! 输入参数
    integer, intent(in) :: nstations
    real, intent(in) :: station_lons(nstations), station_lats(nstations)
    real, intent(in) :: grid_spacing, min_physical_radius
    real, intent(out) :: radii(nstations)
    
    ! 局部变量
    integer :: i, j, k
    integer :: ix, iy
    integer :: count_nonzero
    integer :: nearest_sta
    integer :: nlon, nlat
    logical :: covered, all_covered
    
    real :: distances(nstations), temp_dist
    real :: sorted_distances(nstations)
    real :: lon_min, lon_max, lat_min, lat_max
    real :: lon_range, lat_range, max_radius
    real :: gx, gy, dist, min_dist
    
    !-----------------------------
    ! 计算台站范围
    !-----------------------------
    lon_min = minval(station_lons)
    lon_max = maxval(station_lons)
    lat_min = minval(station_lats)
    lat_max = maxval(station_lats)

    lon_range = lon_max - lon_min
    lat_range = lat_max - lat_min
    max_radius = min(lon_range, lat_range) / 2.0

    ! 网格数（覆盖整个台站范围）
    nlon = int((lon_max - lon_min) / grid_spacing + 0.5) + 1
    nlat = int((lat_max - lat_min) / grid_spacing + 0.5) + 1

    !-----------------------------
    ! 1) 初始半径：最近邻非零距离
    !-----------------------------
    do i = 1, nstations

        do j = 1, nstations
            distances(j) = sqrt((station_lons(i) - station_lons(j))**2 + &
                                (station_lats(i) - station_lats(j))**2)
        end do

        sorted_distances = distances

        ! 排序（从小到大）
        do j = 1, nstations-1
            do k = j+1, nstations
                if (sorted_distances(j) > sorted_distances(k)) then
                    temp_dist = sorted_distances(j)
                    sorted_distances(j) = sorted_distances(k)
                    sorted_distances(k) = temp_dist
                end if
            end do
        end do

        count_nonzero = 0
        do j = 1, nstations
            if (sorted_distances(j) > 1.0e-8) then
                count_nonzero = count_nonzero + 1
            end if
        end do

        if (count_nonzero > 0) then
            radii(i) = min_physical_radius
            do j = 2, nstations
                if (sorted_distances(j) > 1.0e-8) then
                    radii(i) = sorted_distances(j)
                    exit
                end if
            end do
        else
            radii(i) = min_physical_radius
        end if

        if (radii(i) < min_physical_radius) radii(i) = min_physical_radius
        if (radii(i) > max_radius) radii(i) = max_radius

    end do

    !-----------------------------
    ! 2) 扩大半径直到规则网格全覆盖
    !    网格定义：台站范围内，步长为 grid_spacing
    !-----------------------------
    do
        all_covered = .true.

        do iy = 0, nlat-1
            gy = lat_min + real(iy) * grid_spacing

            do ix = 0, nlon-1
                gx = lon_min + real(ix) * grid_spacing

                covered = .false.
                min_dist = 1.0e30
                nearest_sta = 1

                ! 找最近台站，并判断该点是否已被覆盖
                do i = 1, nstations
                    dist = sqrt((gx - station_lons(i))**2 + &
                                (gy - station_lats(i))**2)

                    if (dist < min_dist) then
                        min_dist = dist
                        nearest_sta = i
                    end if

                    if (dist <= radii(i)) then
                        covered = .true.
                    end if
                end do

                ! 若未覆盖，则扩大最近台站半径到该网格点
                if (.not. covered) then
                    all_covered = .false.
                    radii(nearest_sta) = max(radii(nearest_sta), min_dist)

                    if (radii(nearest_sta) < min_physical_radius) then
                        radii(nearest_sta) = min_physical_radius
                    end if

                    if (radii(nearest_sta) > max_radius) then
                        radii(nearest_sta) = max_radius
                    end if
                end if

            end do
        end do

        if (all_covered) exit
    end do

end subroutine calculate_adaptive_radius

subroutine weight(radius, nsta, zhx, zhz, dvx, dvz, gox, goz, &
                  nx, nz, nWeight, Wcoef, xUse, zUse, ista)

    implicit none
    !------ input ------
    real, intent(in) :: radius         ! 影响半径
    integer, intent(in) :: nsta        ! 总台站数
    real, intent(in) :: zhx(nsta)      ! 台站 x
    real, intent(in) :: zhz(nsta)      ! 台站 z
    real, intent(in) :: dvx, dvz       ! 网格间距
    real, intent(in) :: gox, goz       ! 网格左上角坐标
    integer, intent(in) :: nx, nz      ! FULL-GRID 尺寸
    integer, intent(in) :: ista        ! 当前台站编号

    !------ output ------
    integer, intent(out) :: nWeight
    real,    intent(out) :: Wcoef(nx*nz)
    integer, intent(out) :: xUse(nx*nz), zUse(nx*nz)

    !------ local ------
    integer :: ix, iz, idx
    real :: dx, dz, dist, sigma2
    real :: sumw
    real :: xr, zr
   
    ! 台站坐标
    xr = zhx(ista)
    zr = zhz(ista)
    sigma2 = (radius * radius) / 3.0
    nWeight = 0
    !----------------------------------------
    !  FULL-GRID 搜索所有点
    !----------------------------------------
    do iz = 0, nz-1
        do ix = 0, nx-1

            dx = gox - ix*dvx - xr
            dz = goz + iz*dvz - zr
            dist = sqrt(dx*dx + dz*dz)
            if (dist <= radius) then
                nWeight = nWeight + 1
                idx = nWeight
                Wcoef(idx) = exp(- dist*dist / (2.0*sigma2))
                xUse(idx) = ix
                zUse(idx) = iz
            end if

        end do
    end do

    !----------------------------------------
    !  归一化权重
    !----------------------------------------
    if (nWeight > 0) then
        sumw = sum(Wcoef(1:nWeight))
        Wcoef(1:nWeight) = Wcoef(1:nWeight) / sumw
    end if

end subroutine weight

subroutine refineGrid2LayerMdl(minthk0,mmax,dep,vp,vs,rho,&
    rmax,rdep,rvp,rvs,rrho,rthk)
  implicit none
  integer, parameter :: NL=500
  integer mmax,rmax
  real minthk0
  real minthk
  real dep(*),vp(*),vs(*),rho(*)
  real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
  integer nsublay(NL)
  real thk,newthk,initdep
  integer i,j,k

  k = 0
  initdep = 0.0

  do i = 1, mmax-1
    thk = dep(i+1)-dep(i)
    minthk = thk/minthk0
    nsublay(i) = int((thk+1.0e-4)/minthk) + 1
    newthk = thk/nsublay(i)

    do j = 1, nsublay(i)
      k = k + 1
      if (k > NL) then
        print *, 'ERROR: refineGrid2LayerMdl overflow: k=', k, ' NL=', NL, ' i=', i
        stop
      end if
      rthk(k) = newthk
      rdep(k) = initdep + rthk(k)
      initdep = rdep(k)
      rvp(k) = vp(i)+(2*j-1)*(vp(i+1)-vp(i))/(2*nsublay(i))
      rvs(k) = vs(i)+(2*j-1)*(vs(i+1)-vs(i))/(2*nsublay(i))
      rrho(k) = rho(i)+(2*j-1)*(rho(i+1)-rho(i))/(2*nsublay(i))
    end do
  end do

  k = k + 1
  if (k > NL) then
    print *, 'ERROR: refineGrid2LayerMdl halfspace overflow: k=', k, ' NL=', NL
    stop
  end if

  rthk(k) = 0.0
  rvp(k)  = vp(mmax)
  rvs(k)  = vs(mmax)
  rrho(k) = rho(mmax)
  rdep(k) = dep(mmax)

  rmax = k
end subroutine

subroutine syntheticZH(nx, ny, nz, vels, goxdf, goydf, dvxdf, dvydf, &
                       depz, minthk, nZH, zhx, zhy, izhx, izhy, noiselevel, &
                       nperZH, perZH, adaptive_radius_min, obsZH)

   implicit none

   integer, parameter :: NL = 500
   integer, parameter :: maxnperZH = 400

   integer, intent(in) :: nx, ny, nz, nZH
   real,    intent(in) :: vels(nx, ny, nz)
   real,    intent(in) :: goxdf, goydf, dvxdf, dvydf
   real,    intent(in) :: depz(nz), minthk
   real,    intent(in) :: zhx(nZH), zhy(nZH)
   integer, intent(in) :: izhx(nZH), izhy(nZH)
   real,    intent(in) :: noiselevel
   integer, intent(in) :: nperZH(nZH)
   real,    intent(in) :: perZH(nZH, maxnperZH)
   real,    intent(in) :: adaptive_radius_min

   real, intent(out) :: obsZH(sum(nperZH))

   integer :: i, j, k, iz, iwt, rmax
   integer :: nWeight, nUse
   integer :: ix_full, iy_full
   integer :: start_index

   real :: vpz(nz), rhoz(nz), vsz(nz)
   real :: rdep(NL), rrho(NL), rthk(NL), rvp(NL), rvs(NL)
   real*8 :: ZH(maxnperZH)
   real :: radii(nZH)
   real :: grid_spacing
   real :: sumW, randv

   integer, allocatable :: xUse(:), zUse(:)
   real,    allocatable :: Wcoef(:)
   integer, allocatable :: xUse_in(:), zUse_in(:)
   real,    allocatable :: Wcoef_in(:)

   allocate(xUse(nx*ny))
   allocate(zUse(nx*ny))
   allocate(Wcoef(nx*ny))
   allocate(xUse_in(nx*ny))
   allocate(zUse_in(nx*ny))
   allocate(Wcoef_in(nx*ny))

   obsZH = 0.0

   grid_spacing = 0.5 * (dvxdf + dvydf)
   call calculate_adaptive_radius(zhx, zhy, nZH, grid_spacing, adaptive_radius_min, radii)

   do i = 1, nZH

      if (i == 1) then
         start_index = 1
      else
         start_index = sum(nperZH(1:i-1)) + 1
      end if

      call weight(radii(i), nZH, zhx, zhy, dvxdf, dvydf, goxdf, goydf, &
                  nx, ny, nWeight, Wcoef, xUse, zUse, i)

      if (nWeight < 1 .or. nWeight > nx*ny) then
         print *, 'ERROR: nWeight out of range in syntheticZH: ', nWeight, ' station=', i
         stop
      end if

      nUse = 0
      sumW = 0.0
      do iwt = 1, nWeight
         ix_full = xUse(iwt) + 1
         iy_full = zUse(iwt) + 1

         if (ix_full >= 2 .and. ix_full <= nx-1 .and. &
             iy_full >= 2 .and. iy_full <= ny-1) then
            nUse = nUse + 1
            xUse_in(nUse) = xUse(iwt)
            zUse_in(nUse) = zUse(iwt)
            Wcoef_in(nUse) = Wcoef(iwt)
            sumW = sumW + Wcoef(iwt)
         end if
      end do

      if (nUse < 1 .or. sumW <= 0.0) then
         print *, 'ERROR: no interior weights in syntheticZH. station=', i
         stop
      end if

      do iwt = 1, nUse
         Wcoef_in(iwt) = Wcoef_in(iwt) / sumW
      end do

      vsz = 0.0
      do iz = 1, nz
         do iwt = 1, nUse
            vsz(iz) = vsz(iz) + Wcoef_in(iwt) * vels(xUse_in(iwt)+1, zUse_in(iwt)+1, iz)
         end do
      end do

      do k = 1, nz
         vpz(k) = 0.9409 + 2.0947*vsz(k) - 0.8206*vsz(k)**2 + &
                  0.2683*vsz(k)**3 - 0.0251*vsz(k)**4
         rhoz(k) = 1.6612*vpz(k) - 0.4721*vpz(k)**2 + &
                   0.0671*vpz(k)**3 - 0.0043*vpz(k)**4 + &
                   0.000106*vpz(k)**5
      end do

      call refineGrid2LayerMdl(minthk, nz, depz, vpz, vsz, rhoz, &
                               rmax, rdep, rvp, rvs, rrho, rthk)

      call zhFwd(rthk, rvp, rvs, rrho, rmax, 1, nperZH(i), perZH(i,1:nperZH(i)), ZH)

      do j = 1, nperZH(i)
         obsZH(start_index + j - 1) = real(ZH(j))
         if (noiselevel > 0.0) then
            call random_number(randv)
            randv = 2.0*randv - 1.0
            obsZH(start_index + j - 1) = obsZH(start_index + j - 1) + &
                                         noiselevel * abs(obsZH(start_index + j - 1)) * randv
         end if
      end do

   end do

   deallocate(xUse, zUse, Wcoef)
   deallocate(xUse_in, zUse_in, Wcoef_in)

end subroutine syntheticZH
