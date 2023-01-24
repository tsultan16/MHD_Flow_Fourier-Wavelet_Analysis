MODULE grid_arrays_mod

USE constants_mod
IMPLICIT NONE



REAL*4, ALLOCATABLE :: fx(:,:,:,:), buffer(:,:,:), vp(:,:,:,:)
REAL*4, ALLOCATABLE :: vxk_re(:,:,:), vxk_im(:,:,:),vyk_re(:,:,:), vyk_im(:,:,:),vzk_re(:,:,:), vzk_im(:,:,:), &
                       bxk_re(:,:,:), bxk_im(:,:,:),byk_re(:,:,:), byk_im(:,:,:),bzk_re(:,:,:), bzk_im(:,:,:), &
                       rhok_re(:,:,:), rhok_im(:,:,:), va_k_re(:,:,:), va_k_im(:,:,:), vf_k_re(:,:,:), vf_k_im(:,:,:), &
                       vs_k_re(:,:,:), vs_k_im(:,:,:), &
                       vwx(:,:,:), vwy(:,:,:), vwz(:,:,:), fwbk_re(:,:,:), fwbk_im(:,:,:)

INTEGER :: nmax, mem_bytes


CONTAINS


SUBROUTINE create_grid_arrays()

    nmax = MIN(nranks_x*nx,nranks_y*ny,nranks_z*nz)

    ALLOCATE(vxk_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vxk_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vyk_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vyk_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vzk_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vzk_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(bxk_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(bxk_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(byk_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(byk_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(bzk_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(bzk_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(rhok_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(rhok_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(va_k_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(va_k_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vf_k_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vf_k_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vs_k_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vs_k_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(fwbk_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(fwbk_im(0:nmax-1,0:nmax-1,0:nmax-1))


    ALLOCATE(fx(1:nx*nranks_x,1:ny*nranks_y,1:nz*nranks_z,7))
    ALLOCATE(vp(1:nx*nranks_x,1:ny*nranks_y,1:nz*nranks_z,3))
    ALLOCATE(buffer(1:nx*nranks_x,1:ny*nranks_y,1:nz*nranks_z))
    ALLOCATE(vwx(1:nx*nranks_x,1:ny*nranks_y,1:nz*nranks_z))
    ALLOCATE(vwy(1:nx*nranks_x,1:ny*nranks_y,1:nz*nranks_z))
    ALLOCATE(vwz(1:nx*nranks_x,1:ny*nranks_y,1:nz*nranks_z))
    
    mem_bytes = 22*SIZEOF(vxk_re) + SIZEOF(fx) + 4*SIZEOF(buffer) + SIZEOF(vp)
     
    PRINT*,''
    PRINT*,'Memory allocated for work arrays (Mb) = ', mem_bytes*1.e-6
    PRINT*,''

END SUBROUTINE create_grid_arrays


SUBROUTINE destroy_grid_arrays()

    DEALLOCATE(fx, vp, buffer)
    DEALLOCATE(vxk_re, vxk_im, vyk_re, vyk_im, vzk_re, vzk_im)
    DEALLOCATE(bxk_re, bxk_im, byk_re, byk_im, bzk_re, bzk_im)
    DEALLOCATE(rhok_re, rhok_im)
    DEALLOCATE(va_k_re, va_k_im, vf_k_re, vf_k_im, vs_k_re, vs_k_im)
    DEALLOCATE(vwx, vwy, vwz, fwbk_re, fwbk_im)

END SUBROUTINE destroy_grid_arrays


END MODULE grid_arrays_mod