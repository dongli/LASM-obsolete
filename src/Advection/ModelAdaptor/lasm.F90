!
! Fortran interface for adapting LASM to GAMIL.
!
! Note: Currently, the parallel computing is not implemented yet in LASM.
!

module lasm

    implicit none

    private

    public lasm_init
    public lasm_advance
    public lasm_final

    real(8), parameter :: re = 6371000.0d0 ! Earth radius (m)
    real(8), parameter :: pi = 4.0d0*atan(1.0d0)

contains

    subroutine lasm_init(dt, num_lon, num_lat, num_lev, lat)

        ! Note: The lat grids is from North Pole to South Pole.

        real(8), intent(in) :: dt ! time step size (s)
        integer, intent(in) :: num_lon, num_lat, num_lev
        real(8), intent(in) :: lat(num_lat)

        real(8) lon(num_lon), dlon
        real(8) lat_(num_lat) ! lat grid coordinates from South Pole to North Pole.
        integer i, j

        dlon = 2.0d0*pi/num_lon
        do i = 1, num_lon
            lon(i) = (i-1)*dlon
        end do

        lat_(1) = -pi*0.5d0
        do j = 2, num_lat-1
            lat_(j) = lat(j+1)-pi*0.5d0
        end do
        lat_(num_lat) = pi*0.5d0

        call lasm_init_cpp(re, dt, num_lon, lon, num_lat, lat_, num_lev)

    end subroutine lasm_init

    subroutine lasm_advance(u, v, q)

        ! Note: u, v and q contain one level halos.

        real(8), intent(in) :: u(:,:,:), v(:,:,:), q(:,:,:)

        call lasm_advance_cpp(u, v, q)

    end subroutine lasm_advance

    subroutine lasm_final()

        call lasm_final_cpp()

    end subroutine lasm_final

end module lasm
