module utils
implicit none
contains

subroutine NearestGridIndex(lat, lon, lat0, lon0, dlat, dlon, nx, ny, ix, iy)
    implicit none
    ! 输入
    real, intent(in) :: lat, lon
    real, intent(in) :: lat0, lon0
    real, intent(in) :: dlat, dlon
    integer, intent(in) :: nx, ny

    ! 输出
    integer, intent(out) :: ix, iy

    ! 局部变量
    real :: rx, ry

    ! 纬度方向: 起点在上, 向南纬度减小 => 用 (lat0 - lat)
    rx = (lat0 - lat) / dlat + 1.0
    ! 经度方向: 起点在左, 向东递增
    ry = (lon - lon0) / dlon + 1.0

    ! 四舍五入到最近网格点
    ix = nint(rx)
    iy = nint(ry)

    ! 边界控制
    if (ix < 1) ix = 1
    if (ix > nx) ix = nx
    if (iy < 1) iy = 1
    if (iy > ny) iy = ny
end subroutine NearestGridIndex

subroutine calc_reg_nnz(Tik, axis, nx, ny, nz, nnz)
    implicit none
    integer, intent(in) :: Tik, axis, nx, ny, nz
    integer, intent(out) :: nnz
    
    nnz = 0
    select case(Tik)
    case(0)  ! 0阶
        nnz = nx * ny * nz
        
    case(1)  ! 1阶
        if(axis == 0 .or. axis == 1) nnz = nnz + 2*(nx-1)*ny*nz
        if(axis == 0 .or. axis == 2) nnz = nnz + 2*nx*(ny-1)*nz
        if(axis == 0 .or. axis == 3) nnz = nnz + 2*nx*ny*(nz-1)
        
    case(2)  ! 2阶
        nnz = 7*(nx*ny*nz)  ! 中心点+6个邻居
        ! 边界修正
        nnz = nnz - 2*((nx-1)*ny*nz + nx*(ny-1)*nz + nx*ny*(nz-1))
    end select
end subroutine calc_reg_nnz
subroutine generate_reg_matrix(Tik, weight, axis, nx, ny, nz, &
                               row, col, Gvar, row_start)
    implicit none
    integer, intent(in) :: Tik, axis, nx, ny, nz, row_start
    real, intent(in) :: weight
    integer, intent(out) :: row(*), col(*)
    real, intent(out) :: Gvar(*)

    integer :: i, j, k, idx, count, reg_eq
    real :: w(3)

    w = [weight, weight, weight]

    ! 这里必须从1开始，因为传进来的row/col/Gvar已经是切片起点
    count = 1
    reg_eq = 0

    select case(Tik)

    case(0)  ! 0阶
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    reg_eq = reg_eq + 1
                    idx = (k-1)*ny*nx + (j-1)*nx + i

                    row(count) = row_start + reg_eq - 1
                    col(count) = idx
                    Gvar(count) = weight
                    count = count + 1
                end do
            end do
        end do

    case(1)  ! 1阶
        ! x方向
        if(axis == 0 .or. axis == 1) then
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx-1
                        reg_eq = reg_eq + 1
                        idx = (k-1)*ny*nx + (j-1)*nx + i

                        row(count) = row_start + reg_eq - 1
                        col(count) = idx
                        Gvar(count) = -w(1)
                        count = count + 1

                        row(count) = row_start + reg_eq - 1
                        col(count) = idx + 1
                        Gvar(count) =  w(1)
                        count = count + 1
                    end do
                end do
            end do
        end if

        ! y方向
        if(axis == 0 .or. axis == 2) then
            do k = 1, nz
                do j = 1, ny-1
                    do i = 1, nx
                        reg_eq = reg_eq + 1
                        idx = (k-1)*ny*nx + (j-1)*nx + i

                        row(count) = row_start + reg_eq - 1
                        col(count) = idx
                        Gvar(count) = -w(2)
                        count = count + 1

                        row(count) = row_start + reg_eq - 1
                        col(count) = idx + nx
                        Gvar(count) =  w(2)
                        count = count + 1
                    end do
                end do
            end do
        end if

        ! z方向
        if(axis == 0 .or. axis == 3) then
            do k = 1, nz-1
                do j = 1, ny
                    do i = 1, nx
                        reg_eq = reg_eq + 1
                        idx = (k-1)*ny*nx + (j-1)*nx + i

                        row(count) = row_start + reg_eq - 1
                        col(count) = idx
                        Gvar(count) = -w(3)
                        count = count + 1

                        row(count) = row_start + reg_eq - 1
                        col(count) = idx + nx*ny
                        Gvar(count) =  w(3)
                        count = count + 1
                    end do
                end do
            end do
        end if

    case(2)  ! 2阶
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    reg_eq = reg_eq + 1
                    idx = (k-1)*ny*nx + (j-1)*nx + i

                    row(count) = row_start + reg_eq - 1
                    col(count) = idx
                    Gvar(count) = -6.0 * weight
                    count = count + 1

                    if(i > 1) then
                        row(count) = row_start + reg_eq - 1
                        col(count) = idx - 1
                        Gvar(count) = weight
                        count = count + 1
                    end if
                    if(i < nx) then
                        row(count) = row_start + reg_eq - 1
                        col(count) = idx + 1
                        Gvar(count) = weight
                        count = count + 1
                    end if

                    if(j > 1) then
                        row(count) = row_start + reg_eq - 1
                        col(count) = idx - nx
                        Gvar(count) = weight
                        count = count + 1
                    end if
                    if(j < ny) then
                        row(count) = row_start + reg_eq - 1
                        col(count) = idx + nx
                        Gvar(count) = weight
                        count = count + 1
                    end if

                    if(k > 1) then
                        row(count) = row_start + reg_eq - 1
                        col(count) = idx - nx*ny
                        Gvar(count) = weight
                        count = count + 1
                    end if
                    if(k < nz) then
                        row(count) = row_start + reg_eq - 1
                        col(count) = idx + nx*ny
                        Gvar(count) = weight
                        count = count + 1
                    end if
                end do
            end do
        end do

    end select
end subroutine generate_reg_matrix

subroutine save_model_to_file(vsf, depz, filename, nx, ny, nz, goxd, gozd, dvxd, dvzd)
    implicit none

    integer, intent(in) :: nx, ny, nz
    real, intent(in) :: vsf(nx, ny, nz)
    real, intent(in) :: depz(nz)
    real, intent(in) :: goxd, gozd, dvxd, dvzd
    character(len=*), intent(in) :: filename

    integer :: ix, iy, iz, unit, ios
    real :: lon, lat, dep, vel

    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening output file: ', trim(filename)
        stop
    endif

    ! 输出格式：经度 纬度 深度 速度
    ! 与参考程序一致：
    ! lon = gozd + (iy-1)*dvzd
    ! lat = goxd - (ix-1)*dvxd
    do iz = 1, nz
        dep = depz(iz)
        do iy = 1, ny
            lon = gozd + (iy - 1) * dvzd
            do ix = 1, nx
                lat = goxd - (ix - 1) * dvxd
                vel = vsf(ix, iy, iz)
                write(unit, '(4f11.5)') lon, lat, dep, vel
            end do
        end do
    end do

    close(unit)
    print *, 'Model saved successfully: ', trim(filename)

end subroutine save_model_to_file

end module utils