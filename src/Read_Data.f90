
!==============================================================
! Read input parameter file
!==============================================================
module Read_Data
    use Def_model
    implicit none
contains

subroutine read_input_file(filename, params, ierr)
    implicit none
    character(len=*), intent(in)  :: filename
    type(InputParameters), intent(out) :: params
    integer, intent(out) :: ierr

    integer :: iu, ios
    character(len=512) :: line

    ierr = 0

    open(newunit=iu, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) "ERROR: cannot open input file ", trim(filename)
        ierr = 1
        return
    endif

    !--------------------------------------------------
    ! 1. ZH data list
    !--------------------------------------------------
    call read_data_line(iu, line, ios)
    read(line, '(A)') params%files%zh_data_file

    !--------------------------------------------------
    ! 2. nx ny nz
    !--------------------------------------------------
    call read_data_line(iu, line, ios)
    read(line, *) params%grid%nx, params%grid%ny, params%grid%nz

    !--------------------------------------------------
    ! 3. goxd gozd
    !--------------------------------------------------
    call read_data_line(iu, line, ios)
    read(line, *) params%grid%goxd, params%grid%gozd

    !--------------------------------------------------
    ! 4. dvxd dvzd
    !--------------------------------------------------
    call read_data_line(iu, line, ios)
    read(line, *) params%grid%dvxd, params%grid%dvzd

    !--------------------------------------------------
    ! 5. velocity range
    !--------------------------------------------------
    call read_data_line(iu, line, ios)
    read(line, *) params%inversion%vmin, params%inversion%vmax

    !--------------------------------------------------
    ! 6. max iteration
    !--------------------------------------------------
    call read_data_line(iu, line, ios)
    read(line, *) params%inversion%max_iter

    !--------------------------------------------------
    ! 7. regularization order
    !--------------------------------------------------
    call read_data_line(iu, line, ios)
    read(line, *) params%inversion%reg_order

    !--------------------------------------------------
    ! 8. regularization weight
    !--------------------------------------------------
    call read_data_line(iu, line, ios)
    read(line, *) params%inversion%reg_weight

    !--------------------------------------------------
    ! 9. regularization direction flag
    !--------------------------------------------------
    call read_data_line(iu, line, ios)
    read(line, *) params%inversion%direction_flag

    !--------------------------------------------------
    ! 10. synthetic flag
    !--------------------------------------------------
    call read_data_line(iu, line, ios)
    read(line, *) params%inversion%synthetic_flag

    !--------------------------------------------------
    ! 11. adaptive radius minimum
    !--------------------------------------------------
    call read_data_line(iu, line, ios)
    read(line, *) params%inversion%adaptive_radius_min

    !--------------------------------------------------
    ! 12. adaptive radius minimum
    !--------------------------------------------------
    call read_data_line(iu, line, ios)
    read(line, *) params%inversion%syn_noise_level

    close(iu)
end subroutine read_input_file

!==============================================================
! Read next non-comment, non-empty line
!==============================================================
subroutine read_data_line(iu, line, ios)
    implicit none
    integer, intent(in) :: iu
    character(len=*), intent(out) :: line
    integer, intent(out) :: ios

    do
        read(iu, '(A)', iostat=ios) line
        if (ios /= 0) return

        line = adjustl(line)
        if (len_trim(line) == 0) cycle
        if (line(1:1) == '@') cycle

        if (index(line, '@') > 0) then
            line = line(1:index(line,'@')-1)
        endif

        line = trim(line)
        exit
    end do
end subroutine read_data_line

!=======================================================================
! Subroutine: ReadModel
! Description: Reads a 3D shear wave velocity model from a MOD file
! Input:
!   filename - Path to MOD file
!   nx, ny, nz - Grid dimensions
! Output:
!   depth - Depth array (km)
!   vel - 3D shear wave velocity model (nx, ny, nz)
!
! MOD file format:
!   Line 1: depth(1:nz)
!   Then for each depth k = 1..nz:
!       for each latitude j = 1..ny:
!           one row containing vel(i,j,k), i = 1..nx
!
!   That is: depth -> latitude -> longitude
!=======================================================================

subroutine ReadModel(filename, nx, ny, nz, depth, vel)
    implicit none

    ! 输入参数
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nx, ny, nz

    ! 输出参数
    real, allocatable, intent(out) :: depth(:)
    real, allocatable, intent(out) :: vel(:,:,:)

    ! 局部变量
    integer :: i, j, k, ios
    integer :: unit

    ! 分配数组
    allocate(depth(nz))
    allocate(vel(nx, ny, nz))

    ! 打开文件
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening file: ', trim(filename)
        stop
    endif

    !---------------------------------
    ! 第一行：深度（nz个点）
    !---------------------------------
    read(unit, *, iostat=ios) (depth(k), k=1, nz)
    if (ios /= 0) then
        print *, 'Error reading depth array from file: ', trim(filename)
        close(unit)
        stop
    endif

    !---------------------------------
    ! 后续行：速度值
    ! 顺序：先深度(k)，再纬度(j)，行内为经度(i)
    ! 即每一行读 vel(i,j,k), i=1,nx
    !---------------------------------
    do k = 1, nz
        do j = 1, ny
            read(unit, *, iostat=ios) (vel(i, j, k), i=1, nx)
            if (ios /= 0) then
                print *, 'Error reading velocity at depth index = ', k, &
                         ', latitude index = ', j
                close(unit)
                stop
            endif
        enddo
    enddo

    close(unit)

    print *, 'Model file read successfully: ', trim(filename)
    print *, 'Grid size: ', nx, ny, nz

end subroutine ReadModel

end module Read_Data