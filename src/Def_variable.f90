module Def_variable
    use Def_model
    use Def_zh
    use lsmrDataModule, only : dp
    implicit none

    !========================================================
    ! 全局常量
    !========================================================
    integer, parameter :: maxnperZH = 100
    real,    parameter :: minthk = 0.5
    integer, parameter :: maxnar = 10000000
    real,    parameter :: pi = 3.1415926

    !========================================================
    ! 网格与模型
    !========================================================
    integer :: nx, ny, nz
    integer :: nvx, nvy, nvz
    integer :: count_reg
    integer :: idx


    real :: goxd, gozd, dvxd, dvzd
    real :: grid_spacing
    real, allocatable :: vsf(:,:,:)
    real, allocatable :: depz(:)

    !========================================================
    ! ZH 数据
    !========================================================
    integer :: nZH
    integer :: dallZH
    integer, allocatable :: nperZH(:)
    real,    allocatable :: zhx(:), zhz(:)
    integer, allocatable :: izhx(:), izhz(:)
    real,    allocatable :: perZH(:,:)
    real,    allocatable :: wavelengthZH(:)
    real,    allocatable :: obsZH(:), dsynZH(:)
    real,    allocatable :: cbsZH(:), alphaZH(:)

    !========================================================
    ! 稀疏矩阵 G
    !========================================================
    integer :: narZH
    integer, allocatable :: rowZH(:)
    integer, allocatable :: colZH(:)
    real,    allocatable :: GvarZH(:)

    !========================================================
    ! 迭代与统计
    !========================================================
    integer :: maxiter, iter, i, j, k, nper, nar, m, n, ii, jj
    real    :: meanZH, rmsZH, std_devsZH
    real    :: weight
    real    :: per

    !========================================================
    ! 反演存储
    !========================================================
    integer, allocatable :: row(:)
    integer, allocatable :: col(:)
    real,    allocatable :: Gvar(:)
    real,    allocatable :: cbst(:)

    !========================================================
    ! 正则化
    !========================================================
    integer :: reg_start, Tik, axis, nnz_reg, nreg
    integer :: m_total
    real    :: reg_weight

    !========================================================
    ! LSMR 输入/输出
    !========================================================
    integer,  allocatable :: iw_lsmr(:)
    real,     allocatable :: rw_lsmr(:)
    real(dp), allocatable :: b_lsmr(:)
    real(dp), allocatable :: x_lsmr(:)

    real(dp) :: damp, atol, btol, conlim
    integer  :: itnlim, localSize, nout
    integer  :: istop, itn
    real(dp) :: normA, condA, normr, normAr, normx

    !========================================================
    ! 模型更新
    !========================================================
    real, allocatable :: delta_vsf(:,:,:)
    real, allocatable :: vsf_smoothed(:,:,:)

    real(8) :: vsf_new
    real(8), parameter :: max_update = 0.2d0
    real(8) :: vs_min
    real(8) :: vs_max

    integer :: ix, iy, iz
    logical, parameter :: apply_smoothing = .true.
    real(8), parameter :: smooth_factor = 0.2d0
    character(len=10) :: iter_str
    real(8) :: max_delta, min_delta, avg_delta, neighbor_avg

    !========================================================
    ! 输入参数与其他全局变量
    !========================================================
    type(InputParameters) :: params
    integer :: ierr
    real, allocatable :: radii(:)
    character(len=512) :: zh_path
    type(ZH_Data) :: zh
    character(len=256) :: model_file
    integer :: storage_start, reg_row_start

end module Def_variable