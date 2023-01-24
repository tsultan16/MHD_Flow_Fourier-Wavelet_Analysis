!#################################################
! 3D Fast Wavelet Transform  #
! Reference: Numerical Recipes, Press et' al'    #
!#################################################

MODULE fwt_mod

USE OMP_LIB

USE constants_mod

IMPLICIT NONE


INTEGER, PARAMETER :: daub_num = 20 ! set this to either 4 or 12 or 20
INTEGER :: nx_tot 


CONTAINS


SUBROUTINE fwt_init(nx_in)

    INTEGER, INTENT(IN) :: nx_in

    nx_tot = nx_in

END SUBROUTINE fwt_init



SUBROUTINE generate_wavelet_basis(x)

    REAL*4, INTENT(INOUT) :: x(nx_tot, nx_tot, nx_tot,3)
    INTEGER :: i, j, k, ix, iy, iz
    REAL*4 :: xmin(3)
    REAL*4, ALLOCATABLE :: fx(:,:,:)


    ALLOCATE(fx(nx_tot,nx_tot,nx_tot))


    !$OMP PARALLEL PRIVATE(fx, xmin, i, j, k) SHARED(x)
    !$OMP DO
    DO k = 1, nx_tot 
    DO j = 1, nx_tot 
    DO i = 1, nx_tot 

        !PRINT*,'Wavelet # ',i,j,k

        ! compute basis wavelet funciton
        CALL compute_basis_wavelet(i,j,k,fx)


        ! find wavelet center
        xmin = MINLOC(fx) 
        x(i,j,k,:) = xmin(:)

        
    END DO
    END DO
    END DO   
    !$OMP END PARALLEL

    DEALLOCATE(fx)


END SUBROUTINE generate_wavelet_basis


! computes basis wavelet
SUBROUTINE compute_basis_wavelet(i,j,k,fx)

    INTEGER, INTENT(IN) :: i,j,k
    REAL, INTENT(INOUT) :: fx(nx,nx,nx)
    REAL, ALLOCATABLE :: fw(:,:,:)
    INTEGER :: ix, iy, iz
    
    
    ALLOCATE(fw(nx_tot,nx_tot,nx_tot))
    
    
    fw = 0.0
    fw(i,j,k) = 1.0

    ! inverse dwt
    CALL fwt_3d(fw,fx,-1)
    
    
    DEALLOCATE(FW)

END SUBROUTINE compute_basis_wavelet


SUBROUTINE fwt_3d(fx, fw, sgn)

    REAL*4, INTENT(IN) :: fx(1:nx_tot,1:nx_tot,1:nx_tot) 
    REAL*4, INTENT(INOUT) :: fw(1:nx_tot,1:nx_tot,1:nx_tot) 
    INTEGER, INTENT(IN) :: sgn ! 1: forward dwt, -1: inverse dwt
    REAL*4 :: buffer(1:nx_tot)   
    INTEGER :: ix, iy, iz
  
    
    ! clear output array    
    fw = 0.d0
        
    PRINT*,' Doing x pass..'    
   
    ! DWT in x-direction
    !$OMP PARALLEL PRIVATE(buffer,ix,iy,iz) SHARED(fx,fw)
    !$OMP DO
    DO iz = 1, nx_tot
        DO iy = 1, nx_tot
                                
            ! copy strips-along x into 1d buffer    
            DO ix = 1, nx_tot
                buffer(ix) = fx(ix,iy,iz)
            END DO
        
            ! perform 1D dwt 
            CALL fwt_1d(buffer,nx_tot,sgn)
        
            ! copy back
            DO ix = 1, nx_tot
                fw(ix,iy,iz) = buffer(ix)  
            END DO

        END DO       
    END DO       
    !$OMP END PARALLEL    
       
    PRINT*,' Doing y pass..'    

    ! DWT in y-direction
    !$OMP PARALLEL PRIVATE(buffer,ix,iy,iz) SHARED(fw)
    !$OMP DO
    DO iz = 1, nx_tot 
        DO ix = 1, nx_tot 
                
            ! copy strips-along y into 1d buffer    
            DO iy = 1, nx_tot
                buffer(iy) = fw(ix,iy,iz)
            END DO
        
            ! perform 1D dwt 
            CALL fwt_1d(buffer,nx_tot,sgn)
                    
            ! copy back
            DO iy = 1, nx_tot
                fw(ix,iy,iz) = buffer(iy)  
            END DO            
                
        END DO
    END DO
    !$OMP END PARALLEL    
    
    PRINT*,' Doing z pass..'    

    ! DWT in z-direction
    !$OMP PARALLEL PRIVATE(buffer,ix,iy,iz) SHARED(fw)
    !$OMP DO
    DO iy = 1, nx_tot 
        DO ix = 1, nx_tot 
                
            ! copy strips-along y into 1d buffer    
            DO iz = 1, nx_tot
                buffer(iz) = fw(ix,iy,iz)
            END DO
        
            ! perform 1D dwt 
            CALL fwt_1d(buffer,nx_tot,sgn)
                    
            ! copy back
            DO iz = 1, nx_tot
                fw(ix,iy,iz) = buffer(iz)  
            END DO            
                
        END DO
    END DO
    !$OMP END PARALLEL    


END SUBROUTINE fwt_3d



! This subroutine computes the (inverse) wavelet transform of an input vector of length 'n' for sgn = (-1) 1 
! n has to be power of 2
SUBROUTINE fwt_1d(a,n,sgn)

    REAL*4, INTENT(INOUT) :: a(:)    ! input vector
    INTEGER, INTENT(IN) :: sgn, n 
    INTEGER :: nn
    
    IF(n .LT. 4) RETURN
      
    ! compute the wavelet transform  
    IF(sgn .GE. 0) THEN
        
        nn = n ! start at largest hierarchy
        
        DO WHILE(nn .GE. 4) 
            CALL pwt(a,nn,daub_num,sgn)
            nn = nn/2
        END DO
        
    ELSE ! inverse transform

        nn = 4 ! start at smallest

        DO WHILE(nn .LE. n)
            CALL pwt(a,nn,daub_num,sgn)
            nn = nn*2        
        END DO
        
    END IF    


END SUBROUTINE fwt_1d



! partial wavelet transform, i.e. multiplication by the wavelet coefficient matrix followed by a permutation that rearranges the resulting vector
! so that all smooth components (low-pass) occupy the first half followed by the detail coefficients (high-pass) in the second half
SUBROUTINE pwt(a,n,filter,sgn)

    REAL*4, INTENT(INOUT) :: a(:)    ! input data vector
    INTEGER, INTENT(IN) :: sgn, n, filter
    INTEGER, PARAMETER :: nmax = 1024 ! maximum allowed value of n    
    INTEGER, PARAMETER :: ncmax = 50
    REAL*4 :: wksp(nmax), cc(ncmax), cr(ncmax) 
    INTEGER :: ncof, ioff, joff
    INTEGER :: i, ii, j, jf, jr, k, n1, ni, nj, nh, nmod
    REAL*4 :: ai, ai1

    IF(n .LT. 4) RETURN

    IF(n .GT. nmax) THEN 
        PRINT*,'nmax too small in daub4...'
        STOP
    END IF

    ! set filter co-efficients
    
    ncof = filter
    CALL wtset(ncof, ioff, joff, cc, cr)

    nmod = ncof*n        ! A positive constant equal to zero mod n.
    n1 = n-1             ! Mask of all bits, since n a power of 2.
    nh = n/2
    DO j=1,n
        wksp(j) = 0.
    END DO

    ! apply filter
    IF(sgn .GT. 0) THEN
        
        ii = 1
        
        DO i = 1, n, 2
        
            ni = i + nmod + ioff ! Pointer to be incremented and wrapped-around.
            nj = i + nmod + joff
           
            DO k = 1, ncof
                jf = IAND(n1,ni+k) ! We use bitwise and to wrap-around the pointers.
                jr = IAND(n1,nj+k)
                wksp(ii) = wksp(ii) + cc(k) * a(jf+1)        ! these are the smooth coefficients (stored in first half of array)
                wksp(ii+nh) = wksp(ii+nh) + cr(k) * a(jr+1)  ! these are the detail coefficients (stored in the second half of array)
            END DO

            ii = ii + 1
        
        END DO

    ELSE ! inverse transform

        ii = 1
        
        DO i = 1, n, 2
        
            ai = a(ii)
            ai1 = a(ii+nh)
            ni = i + nmod + ioff ! See comments above.
            nj = i + nmod + joff
            
            DO k = 1, ncof
                jf = IAND(n1,ni+k) + 1
                jr = IAND(n1,nj+k) + 1
                wksp(jf) = wksp(jf) + cc(k) * ai
                wksp(jr) = wksp(jr) + cr(k) * ai1
            END DO

            ii=ii+1

        END DO
                
    END IF    

    ! copy from buffer array into input array
    a(1:n) = wksp(1:n)
 

END SUBROUTINE pwt



! initilaization for wavelet filter co-efficients
SUBROUTINE wtset(ncof, ioff, joff, cc, cr)

    INTEGER, INTENT(IN) :: ncof
    INTEGER, INTENT(INOUT) :: ioff, joff
    REAL*4, INTENT(INOUT) :: cc(:), cr(:)
    REAL*4 :: c4(4), c12(12), c20(20)
    INTEGER :: sig, k
    
    ! DAUB4 filter co-efficients
    c4 = (/ 0.4829629131445341, 0.8365163037378079, 0.2241438680420134,-0.1294095225512604 /)    
        
    ! DAUB12 filter co-efficients
    c12 = (/.111540743350, .494623890398, .751133908021, .315250351709,-.226264693965,-.129766867567, &
            .097501605587, .027522865530,-.031582039318, .000553842201, .004777257511,-.001077301085/)
    
    ! DAUB20 filter co-efficients
    c20 = (/ .026670057901, .188176800078, .527201188932, .688459039454, .281172343661,-.249846424327, &
            -.195946274377, .127369340336, .093057364604, -.071394147166,-.029457536822, .033212674059, &
             .003606553567,-.010733175483, .001395351747, .001992405295,-.000685856695,-.000116466855, &
             .000093588670,-.000013264203 /)

    sig = -1
    
    DO k = 1, ncof
    
        IF(ncof .EQ. 4)THEN
            cc(k) = c4(k)
        ELSE IF(ncof .EQ. 12) THEN
            cc(k) = c12(k)
        ELSE IF(ncof .EQ. 20) THEN
            cc(k) = c20(k)
        ELSE
            PRINT*,'Unimplemented value ncof in pwtset. Need to choose from 4, 12 and 20.'
            STOP
        END IF
        
        cr(ncof+1-k) = sig * cc(k)
        sig = -sig

    END DO
    
    joff = -ncof/2 ! center for wavelet function support
    ioff = -ncof/2


END SUBROUTINE wtset


END MODULE fwt_mod