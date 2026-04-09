program ZH_inversion_3d
    use Def_zh
    use Def_variable
    use Read_Data
    use utils
    use lsmrModule, only: lsmr
    use lsmrblasInterface, only: dnrm2
    implicit none

    integer :: synthetic_flag
    real :: syn_noise_level
    real :: adaptive_radius_min
    real, allocatable :: vsf_true(:,:,:)

    character(len=32)  :: out_prefix
    character(len=256) :: log_file
    integer            :: log_unit
    character(len=8)   :: cdate
    character(len=10)  :: ctime

    ! Initialize grid and model parameters
    !--------------------------------------------------------------------------------
    call read_input_file("3DZHTomo.in", params, ierr)
    if (ierr /= 0) stop

    zh_path = params%files%zh_data_file
    nx = params%grid%nx
    ny = params%grid%ny
    nz = params%grid%nz
    goxd = params%grid%goxd
    gozd = params%grid%gozd
    dvxd = params%grid%dvxd
    dvzd = params%grid%dvzd

    maxiter = params%inversion%max_iter
    vs_min = params%inversion%vmin
    vs_max = params%inversion%vmax
    Tik = params%inversion%reg_order
    reg_weight = params%inversion%reg_weight
    axis = params%inversion%direction_flag
    adaptive_radius_min = params%inversion%adaptive_radius_min
    synthetic_flag = params%inversion%synthetic_flag
    syn_noise_level = params%inversion%syn_noise_level

    if (synthetic_flag == 1) then
        out_prefix = 'syn_'
        log_file   = 'syn_run.log'
    else
        out_prefix = 'real_'
        log_file   = 'real_run.log'
    endif

    !--------------------------------------------------------------------------------
    ! Open log file
    !--------------------------------------------------------------------------------
    call date_and_time(cdate, ctime)
    open(newunit=log_unit, file=trim(log_file), status='replace', action='write')

    write(log_unit, '(A)') '======================================'
    write(log_unit, '(A)') 'ZH inversion run log'
    write(log_unit, '(A,A,A)') 'Start time: ', cdate, ' ' // ctime
    write(log_unit, '(A)') '======================================'

    if (synthetic_flag == 1) then
        print *, '======================================'
        print *, ' Running SYNTHETIC TEST mode'
        print *, ' Initial model  : MOD.true'
        print *, ' Min radius     : ', adaptive_radius_min
        print *, ' Noise level    : ', syn_noise_level
        print *, '======================================'

        write(log_unit, '(A)') 'Mode                : SYNTHETIC TEST'
        write(log_unit, '(A)') 'True model file     : MOD.true'
        write(log_unit, '(A,F10.4)') 'Synthetic noise     : ', syn_noise_level
    else
        print *, '======================================'
        print *, ' Running REAL DATA inversion mode'
        print *, ' Min radius     : ', adaptive_radius_min
        print *, '======================================'

        write(log_unit, '(A)') 'Mode                : REAL DATA'
    endif

    write(log_unit, '(A,A)') 'ZH data path        : ', trim(zh_path)
    write(log_unit, '(A,I6,1X,I6,1X,I6)') 'Grid (nx,ny,nz)     : ', nx, ny, nz
    write(log_unit, '(A,F10.4,1X,F10.4)') 'Origin (goxd,gozd)  : ', goxd, gozd
    write(log_unit, '(A,F10.4,1X,F10.4)') 'Spacing (dvxd,dvzd) : ', dvxd, dvzd
    write(log_unit, '(A,F10.4,1X,F10.4)') 'Vs range            : ', vs_min, vs_max
    write(log_unit, '(A,I6)') 'Max iteration       : ', maxiter
    write(log_unit, '(A,I6)') 'Tik order           : ', Tik
    write(log_unit, '(A,F10.4)') 'Reg weight          : ', reg_weight
    write(log_unit, '(A,I6)') 'Direction flag      : ', axis
    write(log_unit, '(A,F10.4)') 'Adaptive radius min : ', adaptive_radius_min
    write(log_unit, '(A,A)') 'Output prefix       : ', trim(out_prefix)
    write(log_unit, '(A)') '======================================'

    call ReadZHDirectory(zh_path, zh)

    allocate(depz(nz))
    allocate(vsf(nx,ny,nz))
    call ReadModel("MOD", nx, ny, nz, depz, vsf)

    print *, "Grid origin (x, y) = ", goxd, gozd
    print *, "Grid spacing (dx, dy) = ", dvxd, dvzd

    ! Initialize ZH data
    !--------------------------------------------------------------------------------
    nZH = zh%nstations
    allocate(zhx(nZH), zhz(nZH), izhx(nZH), izhz(nZH), nperZH(nZH))

    do i = 1, nZH
        zhx(i) = zh%stations(i)%lat
        zhz(i) = zh%stations(i)%lon
        call NearestGridIndex(zhx(i), zhz(i), goxd, gozd, dvxd, dvzd, nx, ny, izhx(i), izhz(i))
        nperZH(i) = zh%stations(i)%ndata
    end do

    dallZH = sum(nperZH)
    print *, 'Total number of ZH data points:', dallZH
    print *, 'Number of ZH stations:', nZH

    write(log_unit, '(A,I10)') 'Number of ZH stations: ', nZH
    write(log_unit, '(A,I10)') 'Total ZH data points : ', dallZH

    allocate(perZH(nZH, maxnperZH), obsZH(dallZH), dsynZH(dallZH), cbsZH(dallZH))
    allocate(rowZH(maxnar), colZH(maxnar), GvarZH(maxnar))
    allocate(row(maxnar), col(maxnar), Gvar(maxnar))

    ! 周期数组
    do i = 1, nZH
        nper = min(nperZH(i), maxnperZH)
        do j = 1, nper
            per = zh%stations(i)%points(j)%period
            perZH(i, j) = per
        end do
    end do

    ! obsZH 来源：真实数据 or synthetic
    if (synthetic_flag == 1) then
        allocate(vsf_true(nx, ny, nz))
        call ReadModel('MOD.true', nx, ny, nz, depz, vsf_true)

        call syntheticZH(nx, ny, nz, vsf_true, goxd, gozd, dvxd, dvzd, &
                         depz, minthk, nZH, zhx, zhz, izhx, izhz, syn_noise_level, &
                         nperZH, perZH, adaptive_radius_min, obsZH)

        deallocate(vsf_true)

        write(log_unit, '(A)') 'obsZH source         : syntheticZH from MOD.true'
    else
        k = 0
        do i = 1, nZH
            do j = 1, nperZH(i)
                k = k + 1
                obsZH(k) = zh%stations(i)%points(j)%zh_ratio
            end do
        end do

        write(log_unit, '(A)') 'obsZH source         : real ZH observations'
    endif
    allocate(radii(nZH))
    grid_spacing = 0.5 * (dvxd + dvzd)
    call calculate_adaptive_radius(zhx, zhz, nZH, grid_spacing,adaptive_radius_min, radii)
    write(log_unit, '(A)') 'Adaptive radii by station'
    write(log_unit, '(A)') 'Format: RADIUS i=... lat=... lon=... radius=... nper=...'

    do i = 1, nZH
        write(log_unit,'(A,I6,A,F12.4,A,F12.4,A,F12.6,A,I6)') &
            'RADIUS i=', i, &
            ' lat=', zhx(i), &
            ' lon=', zhz(i), &
            ' radius=', radii(i), &
            ' nper=', nperZH(i)
    end do

    write(log_unit, '(A)') '--------------------------------------'

    write(log_unit, '(A)') '======================================'
    write(log_unit, '(A)') 'Iteration log'
    write(log_unit, '(A)') 'iter    RMS_ZH      m_total      nar'
    write(log_unit, '(A)') '--------------------------------------'

    ! Main Iteration
    !--------------------------------------------------------------------------------
    do iter = 1, maxiter
        print *, 'Before iter=', iter, ' sum(nperZH)=', sum(nperZH), ' dallZH=', dallZH
        print *, 'min/max nperZH=', minval(nperZH), maxval(nperZH)
        print *, "Computing ZHG...iter", iter

        call CalZHG(nx, ny, nz, vsf, goxd, gozd, dvxd, dvzd, &
                    nZH, zhx, zhz, izhx, izhz, depz, minthk, &
                    nperZH, perZH, dallZH, adaptive_radius_min, &
                    dsynZH, maxnar, &
                    narZH, rowZH, colZH, GvarZH)

        print *, "=== End of CalZHG ==="
        print *, 'after iter=', iter, ' sum(nperZH)=', sum(nperZH), ' dallZH=', dallZH
        print *, 'min/max nperZH=', minval(nperZH), maxval(nperZH)

        do i = 1, dallZH
            cbsZH(i) = obsZH(i) - dsynZH(i)
        end do

        rmsZH = sqrt(sum(cbsZH**2) / dallZH)
        print *, "RMS ZH:", rmsZH

        nar = narZH
        m = dallZH
        n = nx * ny * nz

        row(1:nar) = rowZH(1:nar)
        col(1:nar) = colZH(1:nar)
        Gvar(1:nar) = GvarZH(1:nar)

        call calc_reg_nnz(Tik, axis, nx, ny, nz, nnz_reg)
        print*, 'Regularization added:', nnz_reg, 'non-zero elements'

        storage_start = nar + 1
        reg_row_start = m + 1

        call generate_reg_matrix(Tik, reg_weight, axis, nx, ny, nz, &
            row(storage_start), col(storage_start), Gvar(storage_start), reg_row_start)

        nar = nar + nnz_reg

        nreg = 0
        select case(Tik)
        case(0)
            nreg = nx*ny*nz
        case(1)
            if(axis == 0 .or. axis == 1) nreg = nreg + (nx-1)*ny*nz
            if(axis == 0 .or. axis == 2) nreg = nreg + nx*(ny-1)*nz
            if(axis == 0 .or. axis == 3) nreg = nreg + nx*ny*(nz-1)
        case(2)
            nreg = nx*ny*nz
        end select

        m_total = m + nreg
        print*, 'Total equations:', m_total, '(data:', m, ', reg:', nreg, ')'

        write(log_unit, '(I4,2X,F12.6,2X,I10,2X,I10)') iter, rmsZH, m_total, nar

        allocate(iw_lsmr(2*nar+1), rw_lsmr(nar))
        iw_lsmr(1) = nar

        do i = 1, nar
            iw_lsmr(1+i) = row(i)
            iw_lsmr(1+nar+i) = col(i)
            rw_lsmr(i) = Gvar(i)
        end do

        allocate(b_lsmr(m_total))
        b_lsmr(1:m) = cbsZH(1:m)
        b_lsmr(m+1:m_total) = 0.0d0

        damp = 0.0d0
        atol = 1.0d-6
        btol = 1.0d-6
        conlim = 1.0d12
        itnlim = 1000
        localSize = 0
        nout = 0

        allocate(x_lsmr(n))
        istop = 0
        itn = 0
        normA = 0.0d0
        condA = 0.0d0
        normr = 0.0d0
        normAr = 0.0d0
        normx = 0.0d0

        call LSMR(m_total, n, 2*nar+1, nar, iw_lsmr, rw_lsmr, b_lsmr, damp, &
                  atol, btol, conlim, itnlim, localSize, nout, &
                  x_lsmr, istop, itn, normA, condA, normr, normAr, normx)

        deallocate(iw_lsmr, rw_lsmr, b_lsmr)

        allocate(delta_vsf(nx, ny, nz))
        k = 0
        do iz = 1, nz
            do iy = 1, ny
                do ix = 1, nx
                    k = k + 1
                    delta_vsf(ix, iy, iz) = x_lsmr(k)
                end do
            end do
        end do

        do iz = 1, nz
            do iy = 1, ny
                do ix = 1, nx
                    if (delta_vsf(ix, iy, iz) >  max_update) delta_vsf(ix, iy, iz) =  max_update
                    if (delta_vsf(ix, iy, iz) < -max_update) delta_vsf(ix, iy, iz) = -max_update

                    vsf_new = vsf(ix, iy, iz) + delta_vsf(ix, iy, iz)
                    vsf_new = max(vs_min, min(vs_max, vsf_new))
                    vsf(ix, iy, iz) = vsf_new
                end do
            end do
        end do

        write(iter_str, '(I0)') iter
        model_file = trim(out_prefix) // 'model_iter' // trim(iter_str) // '.dat'
        call save_model_to_file(vsf, depz, model_file, nx, ny, nz, goxd, gozd, dvxd, dvzd)

        write(log_unit, '(A,A)') 'Saved model file     : ', trim(model_file)
        write(log_unit, '(A)') '--------------------------------------'

        deallocate(delta_vsf)
        deallocate(x_lsmr)
    end do

    call date_and_time(cdate, ctime)
    write(log_unit, '(A)') '======================================'
    write(log_unit, '(A,A,A)') 'End time: ', cdate, ' ' // ctime
    write(log_unit, '(A)') 'Run finished successfully.'
    close(log_unit)

end program ZH_inversion_3d