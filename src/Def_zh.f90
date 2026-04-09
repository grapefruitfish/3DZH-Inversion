module Def_zh
    implicit none

    !-----------------------------------------------------------------------
    ! Type: StationDataPoint
    ! Description: Represents a single data point for a station
    ! Fields:
    !   period - Period value (seconds)
    !   zh_ratio - ZH ratio value
    !-----------------------------------------------------------------------
    type :: StationDataPoint
        real :: period
        real :: zh_ratio
    end type StationDataPoint

    !-----------------------------------------------------------------------
    ! Type: StationData
    ! Description: Represents data for a single station
    ! Fields:
    !   lon - Station longitude
    !   lat - Station latitude
    !   ndata - Number of data points for this station
    !   points - Array of data points (StationDataPoint)
    !-----------------------------------------------------------------------
    type :: StationData
        real :: lon
        real :: lat
        integer :: ndata
        type(StationDataPoint), allocatable :: points(:)
    end type StationData

    !-----------------------------------------------------------------------
    ! Type: ZHData
    ! Description: Represents ZH ratio observation data
    ! Fields:
    !   data_path - Path to ZH data directory
    !   nstations - Number of stations
    !   stations - Array of station data
    !   data_vector - Flattened vector of all ZH ratio observations
    !-----------------------------------------------------------------------
    type :: ZH_Data
        character(len=256) :: data_path
        integer :: nstations
        type(StationData), allocatable :: stations(:)
        real, allocatable :: data_vector(:)
    end type ZH_Data

contains

    !-----------------------------------------------------------------------
    ! Subroutine: ReadZHDirectory
    ! Description: Reads all txt files in a directory and fills ZHData
    ! Input:
    !   dirpath - directory containing txt files
    ! Output:
    !   zhdata - populated ZHData structure
    !-----------------------------------------------------------------------
    subroutine ReadZHDirectory(dirpath, zhdata)
        implicit none
        character(len=*), intent(in) :: dirpath
        type(ZH_Data), intent(out) :: zhdata

        !---- 局部变量 ----
        character(len=1024) :: tmpfile, cmd, line
        character(len=1024), allocatable :: filenames(:)
        integer, allocatable :: counts(:)
        integer :: unit_tmp, unit_f, ios, iret
        integer :: nfiles, i, pcount, total_obs, k
        real :: lon, lat, zh , per

        ! ---- 准备临时文件名 ----
        tmpfile = trim(dirpath)
        if (len_trim(tmpfile) == 0) then
            zhdata%data_path = ''
            zhdata%nstations = 0
            allocate(zhdata%stations(0))
            allocate(zhdata%data_vector(0))
            return
        end if
        if (tmpfile(len_trim(tmpfile):len_trim(tmpfile)) == '/') then
            tmpfile = tmpfile(1:len_trim(tmpfile)-1)
        else
            tmpfile = tmpfile(1:len_trim(tmpfile))
        end if
        tmpfile = trim(tmpfile)//"/.zh_file_list.tmp"

        ! ---- 列出目录下的 .txt 文件到临时文件 ----
        cmd = "find '"//trim(dirpath)//"' -maxdepth 1 -type f -name '*.txt' -print > '"//trim(tmpfile)//"' 2>/dev/null"
        call execute_command_line(trim(cmd), wait=.true., exitstat=iret)

        ! ---- 打开临时文件并统计文件数 ----
        unit_tmp = 99
        open(unit=unit_tmp, file=trim(tmpfile), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            zhdata%data_path = dirpath
            zhdata%nstations = 0
            allocate(zhdata%stations(0))
            allocate(zhdata%data_vector(0))
            return
        end if

        nfiles = 0
        do
            read(unit_tmp,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) > 0) nfiles = nfiles + 1
        end do

        if (nfiles == 0) then
            close(unit_tmp)
            cmd = "rm -f '"//trim(tmpfile)//"'"
            call execute_command_line(trim(cmd), wait=.true., exitstat=iret)
            zhdata%data_path = dirpath
            zhdata%nstations = 0
            allocate(zhdata%stations(0))
            allocate(zhdata%data_vector(0))
            return
        end if

        rewind(unit_tmp)
        allocate(filenames(nfiles))
        i = 1
        do
            read(unit_tmp,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) > 0) then
                filenames(i) = trim(line)
                i = i + 1
            end if
        end do
        close(unit_tmp)

        ! ---- 第一遍统计每个文件的数据行数 ----
        allocate(counts(nfiles))
        counts = 0
        total_obs = 0
        unit_f = 98
        do i = 1, nfiles
            open(unit=unit_f, file=trim(filenames(i)), status='old', action='read', iostat=ios)
            if (ios /= 0) then
                counts(i) = 0
                cycle
            end if
            ! 读首行 lon lat
            read(unit_f,*, iostat=ios) lat, lon
            if (ios /= 0) then
                close(unit_f)
                counts(i) = 0
                cycle
            end if
            do
                read(unit_f,'(A)', iostat=ios) line
                if (ios /= 0) exit
                if (len_trim(line) == 0) cycle
                if (line(1:1) == '#') cycle
                counts(i) = counts(i) + 1
                total_obs = total_obs + 1
            end do
            close(unit_f)
        end do

        ! ---- 分配输出结构 ----
        zhdata%data_path = dirpath
        zhdata%nstations = nfiles
        allocate(zhdata%stations(nfiles))
        if (total_obs > 0) then
            allocate(zhdata%data_vector(total_obs))
        else
            allocate(zhdata%data_vector(0))
        end if

        ! ---- 第二遍读取数据并填充 ----
        k = 0
        do i = 1, nfiles
            open(unit=unit_f, file=trim(filenames(i)), status='old', action='read', iostat=ios)
            if (ios /= 0) then
                zhdata%stations(i)%lon = 0.0
                zhdata%stations(i)%lat = 0.0
                zhdata%stations(i)%ndata = 0
                cycle
            end if
            ! 读首行 lon lat
            read(unit_f,*, iostat=ios) lat, lon
            if (ios /= 0) then
                close(unit_f)
                zhdata%stations(i)%lon = 0.0
                zhdata%stations(i)%lat = 0.0
                zhdata%stations(i)%ndata = 0
                cycle
            end if
            zhdata%stations(i)%lon = lon
            zhdata%stations(i)%lat = lat
            zhdata%stations(i)%ndata = counts(i)
            if (counts(i) > 0) allocate(zhdata%stations(i)%points(counts(i)))

            pcount = 0
            do
                read(unit_f,'(A)', iostat=ios) line
                if (ios /= 0) exit
                if (len_trim(line) == 0) cycle
                if (line(1:1) == '#') cycle
                read(line,*, iostat=ios) per, zh
                if (ios /= 0) cycle
                pcount = pcount + 1
                zhdata%stations(i)%points(pcount)%period = per
                zhdata%stations(i)%points(pcount)%zh_ratio = zh
                k = k + 1
                zhdata%data_vector(k) = zh
            end do
            close(unit_f)
        end do

        ! ---- 清理临时文件 ----
        cmd = "rm -f '"//trim(tmpfile)//"'"
        call execute_command_line(trim(cmd), wait=.true., exitstat=iret)

        if (allocated(filenames)) deallocate(filenames)
        if (allocated(counts)) deallocate(counts)

    end subroutine ReadZHDirectory

end module Def_zh