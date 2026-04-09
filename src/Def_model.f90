module Def_model
    implicit none

    !=============================================================
    ! Grid parameters
    !=============================================================
    type :: GridParameters
        integer :: nx, ny, nz
        real    :: goxd, gozd      ! upper-left corner (lat, lon)
        real    :: dvxd, dvzd      ! grid spacing (lat, lon)
    end type GridParameters

    !=============================================================
    ! Inversion parameters (ZH inversion)
    !=============================================================
    type :: InversionParameters
        real    :: vmin, vmax      ! a priori velocity bounds
        integer :: max_iter        ! maximum iteration

        integer :: reg_order       ! Tikhonov order (0/1/2)
        real    :: reg_weight      ! regularization weight

        integer :: direction_flag  ! 0=all, 1=x, 2=y, 3=z
        real    :: adaptive_radius_min ! minimum radius for adaptive smoothing
        
        integer :: synthetic_flag  ! 0=real data, 1=synthetic test
        real    :: syn_noise_level  ! noise level for synthetic data
    end type InversionParameters

    !=============================================================
    ! File paths
    !=============================================================
    type :: FilePaths
        character(len=512) :: zh_data_file
    end type FilePaths

    !=============================================================
    ! Main input container
    !=============================================================
    type :: InputParameters
        type(FilePaths)           :: files
        type(GridParameters)      :: grid
        type(InversionParameters) :: inversion
    end type InputParameters

end module Def_model