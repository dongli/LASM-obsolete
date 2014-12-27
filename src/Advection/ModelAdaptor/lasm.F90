module lasm

    implicit none

    private

    public lasm_init
    public lasm_advance
    public lasm_final

contains

    subroutine lasm_init(re, dt, num_lon, lon, num_lat, lat, num_lev)

        real(8), intent(in) :: re, dt
        integer, intent(in) :: num_lon, num_lat, num_lev
        real(8), intent(in) :: lon(num_lon), lat(num_lat)

        real(8) lon_half(num_lon), lat_half(num_lat-1)

        call lasm_init_cpp(re, dt, num_lon, lon, lon_half, num_lat, lat, lat_half, num_lev)

    end subroutine lasm_init

    subroutine lasm_advance(u, v, q)

        real(8), intent(in) :: u(:,:,:), v(:,:,:), q(:,:,:)

        ! u, v and q should not contain halos.

        call lasm_advance_cpp(u, v, q)

    end subroutine lasm_advance

    subroutine lasm_final()

        call lasm_final_cpp()

    end subroutine lasm_final

end module lasm
