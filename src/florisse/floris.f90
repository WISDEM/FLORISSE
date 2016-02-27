! Flow field calculations have been intentionally left out to save development time.
! The flow field can be calculated using the pure python version of floris 

! This implementation is fully smooth and differentiable with the exception of a 
! discontinuity at the hub of each turbine. The discontinuity only presents issues if
! turbines are place within 1E-15 * rotor diameter of one another, which should never
! happen in a standard optimization problem with a separation constraint.

    
subroutine Hermite_Spline(x, x0, x1, y0, dy0, y1, dy1, y)
    !    This function produces the y and dy values for a hermite cubic spline
    !    interpolating between two end points with known slopes
    !
    !    :param x: x position of output y
    !    :param x0: x position of upwind endpoint of spline
    !    :param x1: x position of downwind endpoint of spline
    !    :param y0: y position of upwind endpoint of spline
    !    :param dy0: slope at upwind endpoint of spline
    !    :param y1: y position of downwind endpoint of spline
    !    :param dy1: slope at downwind endpoint of spline
    !
    !    :return: y: y value of spline at location x
    
    implicit none
        
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)
    
    ! in
    real(dp), intent(in) :: x, x0, x1, y0, dy0, y1, dy1
    
    ! out
    real(dp), intent(out) :: y !, dy_dx
    
    ! local
    real(dp) :: c3, c2, c1, c0

    ! initialize coefficients for parametric cubic spline
    c3 = (2.0_dp*(y1))/(x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3) - &
         (2.0_dp*(y0))/(x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3) + &
         (dy0)/(x0**2 - 2.0_dp*x0*x1 + x1**2) + &
         (dy1)/(x0**2 - 2.0_dp*x0*x1 + x1**2)
         
    c2 = (3.0_dp*(y0)*(x0 + x1))/(x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3) - &
         ((dy1)*(2.0_dp*x0 + x1))/(x0**2 - 2.0_dp*x0*x1 + x1**2) - ((dy0)*(x0 + &
         2.0_dp*x1))/(x0**2 - 2.0_dp*x0*x1 + x1**2) - (3.0_dp*(y1)*(x0 + x1))/(x0**3 - &
         3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3)
         
    c1 = ((dy0)*(x1**2 + 2.0_dp*x0*x1))/(x0**2 - 2.0_dp*x0*x1 + x1**2) + ((dy1)*(x0**2 + &
         2.0_dp*x1*x0))/(x0**2 - 2.0_dp*x0*x1 + x1**2) - (6.0_dp*x0*x1*(y0))/(x0**3 - &
         3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3) + (6.0_dp*x0*x1*(y1))/(x0**3 - &
         3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3)
         
    c0 = ((y0)*(- x1**3 + 3.0_dp*x0*x1**2))/(x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - &
         x1**3) - ((y1)*(- x0**3 + 3.0_dp*x1*x0**2))/(x0**3 - 3.0_dp*x0**2*x1 + &
         3.0_dp*x0*x1**2 - x1**3) - (x0*x1**2*(dy0))/(x0**2 - 2.0_dp*x0*x1 + x1**2) - &
         (x0**2*x1*(dy1))/(x0**2 - 2.0_dp*x0*x1 + x1**2)
!    print *, 'c3 = ', c3
!    print *, 'c2 = ', c2
!    print *, 'c1 = ', c1
!    print *, 'c0 = ', c0
    ! Solve for y and dy values at the given point
    y = c3*x**3 + c2*x**2 + c1*x + c0
    !dy_dx = c3*3*x**2 + c2*2*x + c1

end subroutine Hermite_Spline


subroutine calcOverlapAreas(nTurbines, turbineX, turbineY, rotorDiameter, wakeDiameters, &
                            wakeCenters, wakeOverlapTRel_m)
!    calculate overlap of rotors and wake zones (wake zone location defined by wake 
!    center and wake diameter)
!   turbineX,turbineY is x,y-location of center of rotor
!
!    wakeOverlap(TURBI,TURB,ZONEI) = overlap area of zone ZONEI of wake of turbine 
!     TURB with rotor of downstream turbine
!    TURBI

    implicit none
        
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)
    
    ! in
    integer, intent(in) :: nTurbines
    real(dp), dimension(nTurbines), intent(in) :: turbineX, turbineY, rotorDiameter
    real(dp), dimension(nTurbines, nTurbines, 3), intent(in) :: wakeDiameters
    real(dp), dimension(nTurbines, nTurbines), intent(in) :: wakeCenters
    
    ! out    
    real(dp), dimension(nTurbines, nTurbines, 3), intent(out) :: wakeOverlapTRel_m
    
    ! local
    integer :: turb, turbI, zone
    real(dp), parameter :: pi = 3.141592653589793_dp, tol = 0.000001_dp
    real(dp) :: OVdYd, OVr, OVRR, OVL, OVz
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeOverlap
        
    wakeOverlapTRel_m = 0.0_dp
    wakeOverlap = 0.0_dp
    
    do turb = 1, nTurbines
        do turbI = 1, nTurbines
            if (turbineX(turbI) > turbineX(turb)) then
                OVdYd = wakeCenters(turbI, turb)-turbineY(turbI)    ! distance between wake center and rotor center
                OVr = rotorDiameter(turbI)/2                        ! rotor diameter
                do zone = 1, 3
                    OVRR = wakeDiameters(turbI, turb, zone)/2.0_dp        ! wake diameter
                    OVdYd = abs(OVdYd)
                    if (OVdYd >= 0.0_dp + tol) then
                        ! calculate the distance from the wake center to the vertical line between
                        ! the two circle intersection points
                        OVL = (-OVr*OVr+OVRR*OVRR+OVdYd*OVdYd)/(2.0_dp*OVdYd)
                    else
                        OVL = 0.0_dp
                    end if

                    OVz = OVRR*OVRR-OVL*OVL

                    ! Finish calculating the distance from the intersection line to the outer edge of the wake zone
                    if (OVz > 0.0_dp + tol) then
                        OVz = sqrt(OVz)
                    else
                        OVz = 0.0_dp
                    end if

                    if (OVdYd < (OVr+OVRR)) then ! if the rotor overlaps the wake zone

                        if (OVL < OVRR .and. (OVdYd-OVL) < OVr) then
                            wakeOverlap(turbI, turb, zone) = OVRR*OVRR*dacos(OVL/OVRR) + OVr*OVr*dacos((OVdYd-OVL)/OVr) - OVdYd*OVz
                        else if (OVRR > OVr) then
                            wakeOverlap(turbI, turb, zone) = pi*OVr*OVr
                        else
                            wakeOverlap(turbI, turb, zone) = pi*OVRR*OVRR
                        end if
                    else
                        wakeOverlap(turbI, turb, zone) = 0.0_dp
                    end if
                    
                end do
                
            end if
            
        end do
        
    end do


    do turb = 1, nTurbines
    
        do turbI = 1, nTurbines
    
            wakeOverlap(turbI, turb, 3) = wakeOverlap(turbI, turb, 3)-wakeOverlap(turbI, turb, 2)
            wakeOverlap(turbI, turb, 2) = wakeOverlap(turbI, turb, 2)-wakeOverlap(turbI, turb, 1)
    
        end do
    
    end do
    
    wakeOverlapTRel_m = wakeOverlap

    do turbI = 1, nTurbines
            wakeOverlapTRel_m(turbI, :, :) = wakeOverlapTRel_m(turbI, :, &
                                                         :)/((pi*rotorDiameter(turbI) &
                                                       *rotorDiameter(turbI))/4.0_dp)
    end do
   
                                    
end subroutine calcOverlapAreas    


subroutine CTtoAxialInd(CT, nTurbines, axial_induction)
    
    implicit none
    
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: nTurbines
    real(dp), dimension(nTurbines), intent(in) :: CT

    ! local
    integer :: i

    ! out
    real(dp), dimension(nTurbines), intent(out) :: axial_induction

    axial_induction = 0.0_dp

    ! execute
    do i = 1, nTurbines
        if (CT(i) > 0.96) then  ! Glauert condition
            axial_induction(i) = 0.143_dp + sqrt(0.0203_dp-0.6427_dp*(0.889_dp - CT(i)))
        else
            axial_induction(i) = 0.5_dp*(1.0_dp-sqrt(1.0_dp-CT(i)))
        end if
    end do
    
end subroutine CTtoAxialInd
    

subroutine floris_wcent_wdiam(nTurbines, kd, initialWakeDisplacement, &
                              initialWakeAngle, ke_in, keCorrCT, Region2CT, yaw_deg, Ct, &
                              turbineXw, turbineYw, rotorDiameter, me, bd, useWakeAngle, &
                              adjustInitialWakeDiamToYaw, wakeCentersYT_vec, &
                              wakeDiametersT_vec)
                              
    ! independent variables: yaw_deg Ct turbineXw turbineYw rotorDiameter
    ! dependent variables: wakeCentersYT_vec wakeDiametersT_vec
    
    implicit none
    
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)
    
    ! in
    integer, intent(in) :: nTurbines
    real(dp), intent(in) :: kd, initialWakeDisplacement, initialWakeAngle, ke_in
    real(dp), intent(in) :: keCorrCT, Region2CT, bd
    real(dp), dimension(nTurbines), intent(in) :: yaw_deg, Ct, turbineXw, turbineYw
    real(dp), dimension(nTurbines), intent(in) :: rotorDiameter
    real(dp), dimension(3), intent(in) :: me
    logical, intent(in) :: useWakeAngle, adjustInitialWakeDiamToYaw
    
    ! out
    real(dp), dimension(nTurbines*nTurbines), intent(out) :: wakeCentersYT_vec
    real(dp), dimension(3*nTurbines*nTurbines), intent(out) :: wakeDiametersT_vec
    
    ! local
    real(dp) :: spline_bound ! in rotor diameters
    real(dp), dimension(nTurbines) :: ke, yaw
    Integer :: turb, turbI, zone
    real(dp) :: wakeAngleInit, zeroloc
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp) :: deltax, factor, displacement, x, x1, x2, y1, y2, dy1, dy2
    real(dp) :: wakeDiameter0
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeDiametersT_mat
    real(dp), dimension(nTurbines, nTurbines) :: wakeCentersYT_mat
    
    spline_bound = 1.0_dp
       
    ! execute
    !    if (CTcorrected) then
    !        Ct = Ct_in
    !    else
    !        Ct = Ct_in*cos(yaw)*cos(yaw)
    !    end if

    ! convert yaw from degrees to radians
    yaw = yaw_deg*pi/180.0_dp
        
    ! calculate y-locations of wake centers in wind ref. frame
    wakeCentersYT_mat = 0.0_dp
    
    do turb = 1, nTurbines
        wakeAngleInit = 0.5_dp*sin(yaw(turb))*Ct(turb)
        
        if (useWakeAngle) then
            wakeAngleInit = wakeAngleInit + initialWakeAngle*pi/180.0_dp
        end if
        
        do turbI = 1, nTurbines
            deltax = turbineXw(turbI) - turbineXw(turb)
            factor = (2.0_dp*kd*deltax/rotorDiameter(turb)) + 1.0_dp
            
            if (turbineXw(turb) < turbineXw(turbI)) then
                wakeCentersYT_mat(turbI, turb) = turbineYw(turb) - initialWakeDisplacement
                ! yaw-induced wake center displacement   
                displacement = (wakeAngleInit*(15.0_dp*(factor*factor*factor*factor) &
                               +(wakeAngleInit*wakeAngleInit))/((30.0_dp*kd* &
                               (factor*factor*factor*factor*factor))/ &
                               rotorDiameter(turb))) - (wakeAngleInit* &
                               rotorDiameter(turb)*(15.0_dp + &
                               (wakeAngleInit*wakeAngleInit))/(30.0_dp*kd))
                               
                if (useWakeAngle .eqv. .false.) then
                    displacement = displacement + bd*deltax
                end if
                
                wakeCentersYT_mat(turbI, turb) = wakeCentersYT_mat(turbI, turb) + displacement
                    
            else
                wakeCentersYT_mat(turbI, turb) = turbineYw(turb) - initialWakeDisplacement
            end if

        end do
        
    end do
        
    do turbI = 1, nTurbines
        wakeCentersYT_vec(nTurbines*(turbI-1)+1:nTurbines*(turbI-1)+nTurbines) &
                                 = wakeCentersYT_mat(turbI, :)
    end do
    
    !adjust k_e to C_T, adjusted to yaw
    ke = ke_in + keCorrCT*(Ct-Region2CT)
    
    wakeDiametersT_mat = 0.0_dp

    do turb = 1, nTurbines
        
        if (adjustInitialWakeDiamToYaw) then
            wakeDiameter0 = rotorDiameter(turb)*cos(yaw(turb))
        else
            wakeDiameter0 = rotorDiameter(turb)        
        end if
        
        do turbI = 1, nTurbines
        
            ! turbine separation
            deltax = turbineXw(turbI) - turbineXw(turb)            
            
            ! x position of interest
            x = turbineXw(turbI)                         
              
            zone = 1
            ! define centerpoint of spline
            zeroloc = turbineXw(turb) - wakeDiameter0/(2.0_dp*ke(turb)*me(zone))
            
            if (zeroloc + spline_bound*rotorDiameter(turb) < turbineXw(turbI)) then
                wakeDiametersT_mat(turbI, turb, zone) = 0.0_dp
            
            else if (zeroloc - spline_bound*rotorDiameter(turb) < turbineXw(turbI)) then
                               
                !!!!!!!!!!!!!!!!!!!!!! calculate spline values !!!!!!!!!!!!!!!!!!!!!!!!!!
                
                ! position of upwind point
                x1 = zeroloc - spline_bound*rotorDiameter(turb)
                
                ! diameter of upwind point
                y1 = wakeDiameter0+2.0_dp*ke(turb)*me(zone)*(x1 - turbineXw(turb))
                                
                ! slope at upwind point
                dy1 = 2.0_dp*ke(turb)*me(zone)
                
                ! position of downwind point
                x2 = zeroloc+spline_bound*rotorDiameter(turb)             

                ! diameter at downwind point
                y2 = 0.0_dp
                
                ! slope at downwind point
                dy2 = 0.0_dp
                
                ! solve for the wake zone diameter and its derivative w.r.t. the downwind
                ! location at the point of interest
                call Hermite_Spline(x, x1, x2, y1, dy1, y2, dy2, wakeDiametersT_mat(turbI, turb, zone))
            
            else if (turbineXw(turb) < turbineXw(turbI)) then
                wakeDiametersT_mat(turbI, turb, zone) = wakeDiameter0+2.0_dp*ke(turb)*me(zone)*deltax            
            end if
            
                        
            if (turbineXw(turb) < turbineXw(turbI)) then
                zone = 2
                wakeDiametersT_mat(turbI, turb, zone) = wakeDiameter0 + 2.0_dp*ke(turb)*me(zone)*deltax                   
                zone = 3
                wakeDiametersT_mat(turbI, turb, zone) = wakeDiameter0+2.0_dp*ke(turb)*me(zone)*deltax
            end if      
            
        end do
        
    end do
    
    do turbI = 1, nTurbines
        wakeDiametersT_vec(3*nTurbines*(turbI-1)+1:3*nTurbines*(turbI-1)+nTurbines) &
                                 = wakeDiametersT_mat(turbI, :, 1)
        wakeDiametersT_vec(3*nTurbines*(turbI-1)+nTurbines+1:3*nTurbines*(turbI-1) &
                                   +2*nTurbines) = wakeDiametersT_mat(turbI, :, 2)
        wakeDiametersT_vec(3*nTurbines*(turbI-1)+2*nTurbines+1:nTurbines*(turbI-1) &
                                   +3*nTurbines) = wakeDiametersT_mat(turbI, :, 3) 
    end do
    
    
end subroutine floris_wcent_wdiam


subroutine floris_overlap(nTurbines, turbineXw, turbineYw, rotorDiameter, &
                          wakeDiametersT_vec, wakeCentersYT_vec, cos_spread, &
                          wakeOverlapTRel_vec, cosFac_vec)
                          
    ! dependent variables: wakeOverlapTRel_vec cosFac_vec
    ! independent variables: turbineYw rotorDiameter wakeDiametersT_vec wakeCentersYT_vec
                              
    implicit none
        
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)
    
    ! in
    integer, intent(in) :: nTurbines
    real(dp), intent(in) :: cos_spread
    real(dp), dimension(nTurbines), intent(in) :: turbineXw, turbineYw, rotorDiameter
    real(dp), dimension(3*nTurbines*nTurbines), intent(in) :: wakeDiametersT_vec
    real(dp), dimension(nTurbines*nTurbines), intent(in) :: wakeCentersYT_vec
    
    ! out
    real(dp), dimension(3*nTurbines*nTurbines), intent(out) :: wakeOverlapTRel_vec
    real(dp), dimension(3*nTurbines*nTurbines), intent(out) :: cosFac_vec
    
    ! local
    real(dp) :: rmax
    integer :: turbI, turb, zone
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeOverlapTRel_mat, wakeDiametersT_mat
    real(dp), dimension(nTurbines, nTurbines, 3) :: cosFac_mat
    real(dp), dimension(nTurbines, nTurbines) :: wakeCentersYT_mat

    ! execute    
    ! pack wakeDiametersT_v vector into a matrix
    do turbI = 1, nTurbines
        wakeDiametersT_mat(turbI, :, 1) = wakeDiametersT_vec((3*nTurbines*(turbI-1)+1): &
                                   (3*nTurbines*(turbI-1)+nTurbines))
        wakeDiametersT_mat(turbI, :, 2) = wakeDiametersT_vec(3*nTurbines*(turbI-1)+nTurbines+1:&
                                   3*nTurbines*(turbI-1)+2*nTurbines)
        wakeDiametersT_mat(turbI, :, 3) = wakeDiametersT_vec(3*nTurbines*(turbI-1)+2*nTurbines+1:&
                                   3*nTurbines*(turbI-1)+3*nTurbines)
    end do
    
    ! pack wakeCentersYT_v vector into a matrix
    do turbI = 1, nTurbines
        wakeCentersYT_mat(turbI, :) = wakeCentersYT_vec((nTurbines*(turbI-1)+1): &
                                   (nTurbines*(turbI-1)+nTurbines))
    end do
    
    ! calculate relative overlap
    call calcOverlapAreas(nTurbines, turbineXw, turbineYw, rotorDiameter, wakeDiametersT_mat, &
                            wakeCentersYT_mat, wakeOverlapTRel_mat)
    
    ! calculate cosine factor
    do turbI = 1, nTurbines
        do turb = 1, nTurbines
            do zone = 1, 3
                rmax = cos_spread*0.5_dp*(wakeDiametersT_mat(turbI, turb, 3) + rotorDiameter(turbI))
                cosFac_mat(turbI, turb, zone) = 0.5_dp*(1.0_dp + cos(pi*dabs( &
                    wakeCentersYT_mat(turbI, turb)-turbineYw(turbI))/rmax))
                !cosFac_mat(turbI, turb, zone) = 1.0_dp
            end do
        end do
    end do
    
    wakeOverlapTRel_vec = 0.0_dp
    
    ! pack relative wake overlap and cosine factor into vectors
    do turbI = 1, nTurbines
            wakeOverlapTRel_vec(3*nTurbines*(turbI-1)+1:3*nTurbines*(turbI-1)+nTurbines) &
                                 = wakeOverlapTRel_mat(turbI, :, 1)
            wakeOverlapTRel_vec(3*nTurbines*(turbI-1)+nTurbines+1:3*nTurbines*(turbI-1) &
                                   +2*nTurbines) = wakeOverlapTRel_mat(turbI, :, 2)
            wakeOverlapTRel_vec(3*nTurbines*(turbI-1)+2*nTurbines+1:3*nTurbines*(turbI-1) &
                                   +3*nTurbines) = wakeOverlapTRel_mat(turbI, :, 3)
            
            cosFac_vec(3*nTurbines*(turbI-1)+1:3*nTurbines*(turbI-1)+nTurbines) &
                                 = cosFac_mat(turbI, :, 1)
            cosFac_vec(3*nTurbines*(turbI-1)+nTurbines+1:3*nTurbines*(turbI-1) &
                                   +2*nTurbines) = cosFac_mat(turbI, :, 2)
            cosFac_vec(3*nTurbines*(turbI-1)+2*nTurbines+1:3*nTurbines*(turbI-1) &
                                   +3*nTurbines) = cosFac_mat(turbI, :, 3)            
    end do
    
end subroutine floris_overlap


subroutine floris_velocity(nTurbines, wakeOverlapTRel_v, CosFac_v, Ct, a_in, &
                 axialIndProvided, useaUbU, keCorrCT, Region2CT, ke_in, Vinf, & 
                 keCorrArray, turbineXw, yaw_deg, rotorDiameter, MU, aU, bU, &
                 velocitiesTurbines)

    ! dependent variables: velocitiesTurbines, wt_power, power
    ! independent variables: wakeOverlapTRel_v, cosFac_v, Ct, a_in, turbineXw, yaw_deg, rotorDiameter
    

    implicit none
    
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)


    ! in
    integer, intent(in) :: nTurbines                                            ! Number of turbines
    real(dp), dimension(3*nTurbines*nTurbines), intent(in) :: wakeOverlapTRel_v, cosFac_v ! relative overlap and cosine factor (vector of length 3*nTurbines**2)
    real(dp), dimension(nTurbines), intent(in) :: Ct, a_in             ! thrust coeff, axial induction, turbine yaw
    logical,  intent(in) :: axialIndProvided, useaUbU                   ! option logicals
    real(dp), intent(in) :: keCorrCT, Region2CT, ke_in, Vinf, keCorrArray
    real(dp), dimension(nTurbines), intent(in) :: turbineXw, yaw_deg
    real(dp), dimension(nTurbines), intent(in) :: rotorDiameter
    real(dp), dimension(3), intent(in) :: MU
    real(dp), intent(in) :: aU, bU
        
    ! local
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeOverlapTRel, cosFac    ! relative overlap and Cosine factor (NxNx3)
    real(dp), dimension(nTurbines) :: a, ke, keArray, yaw
    real(dp), dimension(3) :: mmU
    real(dp) :: s, wakeEffCoeff, wakeEffCoeffPerZone, deltax
    real(dp), parameter :: pi = 3.141592653589793_dp
    integer :: turb, turbI, zone
    
    ! out
    real(dp), dimension(nTurbines), intent(out) :: velocitiesTurbines
    
    intrinsic cos

    ! convert yaw from degrees to radians
    yaw = yaw_deg*pi/180.0_dp

    ! pack overlap and cosine factor vectors into matrices
    do turbI = 1, nTurbines
        wakeOverlapTRel(turbI, :, 1) = wakeOverlapTRel_v(3*nTurbines*(turbI-1)+1: &
                                   3*nTurbines*(turbI-1)+nTurbines)
        wakeOverlapTRel(turbI, :, 2) = wakeOverlapTRel_v(3*nTurbines*(turbI-1)+nTurbines+1:&
                                   3*nTurbines*(turbI-1)+2*nTurbines)
        wakeOverlapTRel(turbI, :, 3) = wakeOverlapTRel_v(3*nTurbines*(turbI-1)+2*nTurbines+1:&
                                   3*nTurbines*(turbI-1)+3*nTurbines)
        cosFac(turbI, :, 1) = cosFac_v(3*nTurbines*(turbI-1)+1: &
                                   3*nTurbines*(turbI-1)+nTurbines)
        cosFac(turbI, :, 2) = cosFac_v(3*nTurbines*(turbI-1)+nTurbines+1:&
                                   3*nTurbines*(turbI-1)+2*nTurbines)
        cosFac(turbI, :, 3) = cosFac_v(3*nTurbines*(turbI-1)+2*nTurbines+1:&
                                   3*nTurbines*(turbI-1)+3*nTurbines)
    end do

    if (axialIndProvided) then
        a = a_in
    else
        call CTtoAxialInd(Ct, nTurbines, a)
    end if
    
    ! adjust ke to Ct as adjusted to yaw
    ke = ke_in + keCorrCT*(Ct-Region2CT)

    do turb = 1, nTurbines
        s = sum(wakeOverlapTrel(turb, :, 1) + wakeOverlapTRel(turb, :, 2))
        keArray(turb) = ke(turb)*(1+s*keCorrArray) 
    end do
    
    ! find effective wind speeds at downstream turbines, then predict the power of 
    ! downstream turbine    
    velocitiesTurbines = Vinf
    do turbI = 1, nTurbines
        wakeEffCoeff = 0.0_dp
        
        ! find overlap-area weighted effect of each wake zone
        do turb = 1, nTurbines
            wakeEffCoeffPerZone = 0.0_dp
            deltax = turbineXw(turbI) - turbineXw(turb)
            
            if (useaUbU) then
                mmU = MU/cos(aU*pi/180.0_dp + bU*yaw(turb))
            end if
            
            if (deltax > 0.0_dp) then
                do zone = 1, 3
                
                    if (useaUbU) then
                        wakeEffCoeffPerZone = wakeEffCoeffPerZone + &
                        ((((rotorDiameter(turb))/(rotorDiameter(turb)+2.0_dp*keArray(turb) &
                        *mmU(zone)*deltax))*(cosFac(turbI, turb, zone)))**2)*wakeOverlapTRel(turbI, turb, zone)   
                    else
                        wakeEffCoeffPerZone = wakeEffCoeffPerZone + &
                        (((rotorDiameter(turb)/(rotorDiameter(turb)+2.0_dp*keArray(turb) &
                        *MU(zone)*deltax))*(cosFac(turbI, turb, zone)))**2)*wakeOverlapTRel(turbI, turb, zone)   
                    end if                     
                            
                end do
                wakeEffCoeff = wakeEffCoeff + (a(turb)*wakeEffCoeffPerZone)**2
            end if
        end do
        wakeEffCoeff = 1.0_dp - 2.0_dp*sqrt(wakeEffCoeff)
        
        ! multiply the inflow speed with the wake coefficients to find effective wind 
        ! speed at turbine
        velocitiesTurbines(turbI) = velocitiesTurbines(turbI)*wakeEffCoeff
    end do
    
end subroutine floris_velocity


!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5902M) - 15 Dec 2015 09:00
!
!  Differentiation of hermite_spline in reverse (adjoint) mode:
!   gradient     of useful results: y
!   with respect to varying inputs: x x0 x1 dy0 y0
SUBROUTINE HERMITE_SPLINE_BV(x, xb, x0, x0b, x1, x1b, y0, y0b, dy0, dy0b&
& , y1, dy1, y, yb, nbdirs)
!  USE DIFFSIZES
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
!dy_dx = c3*3*x**2 + c2*2*x + c1
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: x, x0, x1, y0, dy0, y1, dy1
  REAL(dp), DIMENSION(nbdirs) :: xb, x0b, x1b, y0b, dy0b
! out
!, dy_dx
  REAL(dp) :: y
  REAL(dp), DIMENSION(nbdirs) :: yb
! local
  REAL(dp) :: c3, c2, c1, c0
  REAL(dp), DIMENSION(nbdirs) :: c3b, c2b, c1b, c0b
  INTRINSIC KIND
  INTEGER :: nd
  INTEGER :: nbdirs
  REAL(dp) :: temp3
  REAL(dp) :: temp2
  REAL(dp) :: temp1
  REAL(dp) :: temp0
  REAL(dp) :: tempb9(nbdirs)
  REAL(dp) :: tempb8(nbdirs)
  REAL(dp) :: tempb7(nbdirs)
  REAL(dp) :: tempb6(nbdirs)
  REAL(dp) :: tempb5(nbdirs)
  REAL(dp) :: tempb4(nbdirs)
  REAL(dp) :: tempb19(nbdirs)
  REAL(dp) :: tempb3(nbdirs)
  REAL(dp) :: tempb18(nbdirs)
  REAL(dp) :: tempb2(nbdirs)
  REAL(dp) :: tempb17(nbdirs)
  REAL(dp) :: tempb1(nbdirs)
  REAL(dp) :: tempb16(nbdirs)
  REAL(dp) :: tempb0(nbdirs)
  REAL(dp) :: tempb15(nbdirs)
  REAL(dp) :: tempb14(nbdirs)
  REAL(dp) :: tempb13(nbdirs)
  REAL(dp) :: tempb12(nbdirs)
  REAL(dp) :: tempb11(nbdirs)
  REAL(dp) :: tempb10(nbdirs)
  REAL(dp) :: temp16
  REAL(dp) :: temp15
  REAL(dp) :: temp14
  REAL(dp) :: temp13
  REAL(dp) :: temp12
  REAL(dp) :: temp11
  REAL(dp) :: temp10
  REAL(dp) :: tempb(nbdirs)
  REAL(dp) :: tempb30(nbdirs)
  REAL(dp) :: tempb29(nbdirs)
  REAL(dp) :: tempb28(nbdirs)
  REAL(dp) :: tempb27(nbdirs)
  REAL(dp) :: tempb26(nbdirs)
  REAL(dp) :: tempb25(nbdirs)
  REAL(dp) :: temp
  REAL(dp) :: tempb24(nbdirs)
  REAL(dp) :: tempb23(nbdirs)
  REAL(dp) :: tempb22(nbdirs)
  REAL(dp) :: temp9
  REAL(dp) :: tempb21(nbdirs)
  REAL(dp) :: temp8
  REAL(dp) :: tempb20(nbdirs)
  REAL(dp) :: temp7
  REAL(dp) :: temp6
  REAL(dp) :: temp5
  REAL(dp) :: temp4
! initialize coefficients for parametric cubic spline
  c3 = 2.0_dp*y1/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) - 2.0_dp*&
&   y0/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) + dy0/(x0**2-2.0_dp&
&   *x0*x1+x1**2) + dy1/(x0**2-2.0_dp*x0*x1+x1**2)
  c2 = 3.0_dp*y0*(x0+x1)/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) -&
&   dy1*(2.0_dp*x0+x1)/(x0**2-2.0_dp*x0*x1+x1**2) - dy0*(x0+2.0_dp*x1)/(&
&   x0**2-2.0_dp*x0*x1+x1**2) - 3.0_dp*y1*(x0+x1)/(x0**3-3.0_dp*x0**2*x1&
&   +3.0_dp*x0*x1**2-x1**3)
  c1 = dy0*(x1**2+2.0_dp*x0*x1)/(x0**2-2.0_dp*x0*x1+x1**2) + dy1*(x0**2+&
&   2.0_dp*x1*x0)/(x0**2-2.0_dp*x0*x1+x1**2) - 6.0_dp*x0*x1*y0/(x0**3-&
&   3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) + 6.0_dp*x0*x1*y1/(x0**3-&
&   3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3)
!    print *, 'c3 = ', c3
!    print *, 'c2 = ', c2
!    print *, 'c1 = ', c1
!    print *, 'c0 = ', c0
! Solve for y and dy values at the given point
  temp13 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp12 = 3.0_dp*x0*x1**2 - x1**3
  temp14 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp15 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp16 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp8 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp7 = x1**2 + 2.0_dp*x0*x1
  temp9 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp10 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp11 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp3 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp4 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp5 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp6 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp0 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp1 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp2 = x0**2 - 2.0_dp*x0*x1 + x1**2
  DO nd=1,nbdirs
    c3b(nd) = x**3*yb(nd)
    xb(nd) = (c1+c2*2*x+c3*3*x**2)*yb(nd)
    c2b(nd) = x**2*yb(nd)
    c1b(nd) = x*yb(nd)
    c0b(nd) = yb(nd)
    tempb(nd) = c0b(nd)/temp13
    tempb0(nd) = y0*tempb(nd)
    tempb1(nd) = -(y0*temp12*tempb(nd)/temp13)
    tempb2(nd) = -(y1*c0b(nd)/temp14)
    tempb3(nd) = -((3.0_dp*(x1*x0**2)-x0**3)*tempb2(nd)/temp14)
    tempb4(nd) = -(c0b(nd)/temp15)
    tempb5(nd) = x1**2*tempb4(nd)
    tempb6(nd) = -(x1**2*x0*dy0*tempb4(nd)/temp15)
    tempb7(nd) = -(dy1*c0b(nd)/temp16)
    tempb8(nd) = -(x0**2*x1*tempb7(nd)/temp16)
    tempb30(nd) = c1b(nd)/temp8
    tempb12(nd) = dy0*tempb30(nd)
    tempb13(nd) = -(dy0*temp7*tempb30(nd)/temp8)
    tempb14(nd) = dy1*c1b(nd)/temp9
    tempb15(nd) = -((x0**2+2.0_dp*(x1*x0))*tempb14(nd)/temp9)
    tempb16(nd) = y1*6.0_dp*c1b(nd)/temp10
    tempb17(nd) = -(x0*x1*tempb16(nd)/temp10)
    tempb9(nd) = -(6.0_dp*c1b(nd)/temp11)
    tempb18(nd) = -(x0*x1*y0*tempb9(nd)/temp11)
    tempb11(nd) = 3.0_dp*c2b(nd)/temp3
    tempb29(nd) = -(y0*(x0+x1)*tempb11(nd)/temp3)
    tempb28(nd) = -(dy1*c2b(nd)/temp4)
    tempb27(nd) = -((2.0_dp*x0+x1)*tempb28(nd)/temp4)
    tempb26(nd) = -(c2b(nd)/temp5)
    dy0b(nd) = temp7*tempb30(nd) + c3b(nd)/temp1 + (x0+2.0_dp*x1)*&
&     tempb26(nd) + x0*tempb5(nd)
    tempb25(nd) = -(dy0*(x0+2.0_dp*x1)*tempb26(nd)/temp5)
    tempb24(nd) = -(y1*3.0_dp*c2b(nd)/temp6)
    tempb23(nd) = -((x0+x1)*tempb24(nd)/temp6)
    tempb19(nd) = -(y1*2.0_dp*c3b(nd)/temp**2)
    tempb10(nd) = -(2.0_dp*c3b(nd)/temp0)
    y0b(nd) = x0*x1*tempb9(nd) + tempb10(nd) + (x0+x1)*tempb11(nd) + &
&     temp12*tempb(nd)
    tempb20(nd) = -(y0*tempb10(nd)/temp0)
    tempb21(nd) = -(dy0*c3b(nd)/temp1**2)
    tempb22(nd) = -(dy1*c3b(nd)/temp2**2)
    x0b(nd) = 2.0_dp*x1*tempb12(nd) + (2*x0-2.0_dp*x1)*tempb13(nd) + (&
&     2.0_dp*x1+2*x0)*tempb14(nd) + (2*x0-2.0_dp*x1)*tempb15(nd) + x1*&
&     tempb16(nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb17(nd) + &
&     y0*x1*tempb9(nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb18(&
&     nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb19(nd) + (3.0_dp*&
&     x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb20(nd) + (2*x0-2.0_dp*x1)*&
&     tempb21(nd) + (2*x0-2.0_dp*x1)*tempb22(nd) + (3.0_dp*x1**2-x1*&
&     3.0_dp*2*x0+3*x0**2)*tempb23(nd) + tempb24(nd) + (2*x0-2.0_dp*x1)*&
&     tempb25(nd) + dy0*tempb26(nd) + (2*x0-2.0_dp*x1)*tempb27(nd) + &
&     2.0_dp*tempb28(nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb29&
&     (nd) + y0*tempb11(nd) + (2*x0-2.0_dp*x1)*tempb8(nd) + x1*2*x0*&
&     tempb7(nd) + (2*x0-2.0_dp*x1)*tempb6(nd) + dy0*tempb5(nd) + (&
&     3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb3(nd) + (x1*3.0_dp*2*x0-&
&     3*x0**2)*tempb2(nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb1&
&     (nd) + 3.0_dp*x1**2*tempb0(nd)
    x1b(nd) = (2.0_dp*x0+2*x1)*tempb12(nd) + (2*x1-2.0_dp*x0)*tempb13(nd&
&     ) + 2.0_dp*x0*tempb14(nd) + (2*x1-2.0_dp*x0)*tempb15(nd) + x0*&
&     tempb16(nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb17(nd) + &
&     y0*x0*tempb9(nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb18(&
&     nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb19(nd) + (x0*&
&     3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb20(nd) + (2*x1-2.0_dp*x0)*&
&     tempb21(nd) + (2*x1-2.0_dp*x0)*tempb22(nd) + (x0*3.0_dp*2*x1-3*x1&
&     **2-3.0_dp*x0**2)*tempb23(nd) + tempb24(nd) + (2*x1-2.0_dp*x0)*&
&     tempb25(nd) + dy0*2.0_dp*tempb26(nd) + (2*x1-2.0_dp*x0)*tempb27(nd&
&     ) + tempb28(nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb29(nd&
&     ) + y0*tempb11(nd) + (2*x1-2.0_dp*x0)*tempb8(nd) + x0**2*tempb7(nd&
&     ) + (2*x1-2.0_dp*x0)*tempb6(nd) + x0*dy0*2*x1*tempb4(nd) + (x0*&
&     3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb3(nd) + 3.0_dp*x0**2*tempb2&
&     (nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb1(nd) + (x0*&
&     3.0_dp*2*x1-3*x1**2)*tempb0(nd)
  END DO
END SUBROUTINE HERMITE_SPLINE_BV


!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5902M) - 15 Dec 2015 09:00
!
!  Differentiation of calcoverlapareas in reverse (adjoint) mode:
!   gradient     of useful results: rotordiameter turbiney wakeoverlaptrel_m
!                wakediameters wakecenters
!   with respect to varying inputs: rotordiameter turbiney wakediameters
!                wakecenters
SUBROUTINE CALCOVERLAPAREAS_BV(nturbines, turbinex, turbiney, turbineyb&
& , rotordiameter, rotordiameterb, wakediameters, wakediametersb, &
& wakecenters, wakecentersb, wakeoverlaptrel_m, wakeoverlaptrel_mb, &
& nbdirs)
!  USE DIFFSIZES
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: turbinex, turbiney, &
& rotordiameter
  REAL(dp), DIMENSION(nbdirs, nturbines) :: turbineyb, rotordiameterb
  REAL(dp), DIMENSION(nturbines, nturbines, 3), INTENT(IN) :: &
& wakediameters
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakediametersb
  REAL(dp), DIMENSION(nturbines, nturbines), INTENT(IN) :: wakecenters
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines) :: wakecentersb
! out    
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlaptrel_m
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlaptrel_mb
! local
  INTEGER :: turb, turbi, zone
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp, tol=0.000001_dp
  REAL(dp) :: ovdyd, ovr, ovrr, ovl, ovz
  REAL(dp), DIMENSION(nbdirs) :: ovdydb, ovrb, ovrrb, ovlb, ovzb
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlap
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlapb
  INTRINSIC KIND
  INTRINSIC ABS
  INTRINSIC SQRT
  INTRINSIC DACOS
  INTEGER :: nd
  INTEGER :: branch
  INTEGER :: nbdirs
  REAL(dp) :: temp1
  REAL(dp) :: temp0
  REAL(dp) :: tempb2(nbdirs, nturbines, 3)
  REAL(dp) :: tempb1(nbdirs)
  REAL(dp) :: tempb0(nbdirs)
  REAL(dp) :: tempb(nbdirs)
  REAL(dp) :: temp
  wakeoverlap = 0.0_dp
  DO turb=1,nturbines
    DO turbi=1,nturbines
      IF (turbinex(turbi) .GT. turbinex(turb)) THEN
! distance between wake center and rotor center
        CALL PUSHREAL4ARRAY(ovdyd, dp/4)
        ovdyd = wakecenters(turbi, turb) - turbiney(turbi)
! rotor diameter
        CALL PUSHREAL4ARRAY(ovr, dp/4)
        ovr = rotordiameter(turbi)/2
        DO zone=1,3
! wake diameter
          ovrr = wakediameters(turbi, turb, zone)/2.0_dp
          IF (ovdyd .GE. 0.) THEN
            CALL PUSHREAL4ARRAY(ovdyd, dp/4)
            ovdyd = ovdyd
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREAL4ARRAY(ovdyd, dp/4)
            ovdyd = -ovdyd
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ovdyd .GE. 0.0_dp + tol) THEN
! calculate the distance from the wake center to the vertical line between
! the two circle intersection points
            CALL PUSHREAL4ARRAY(ovl, dp/4)
            ovl = (-(ovr*ovr)+ovrr*ovrr+ovdyd*ovdyd)/(2.0_dp*ovdyd)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREAL4ARRAY(ovl, dp/4)
            ovl = 0.0_dp
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL4ARRAY(ovz, dp/4)
          ovz = ovrr*ovrr - ovl*ovl
! Finish calculating the distance from the intersection line to the outer edge of the wake zone
          IF (ovz .GT. 0.0_dp + tol) THEN
            CALL PUSHREAL4ARRAY(ovz, dp/4)
            ovz = SQRT(ovz)
            CALL PUSHCONTROL1B(0)
          ELSE
            ovz = 0.0_dp
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ovdyd .LT. ovr + ovrr) THEN
! if the rotor overlaps the wake zone
            IF (ovl .LT. ovrr .AND. ovdyd - ovl .LT. ovr) THEN
              wakeoverlap(turbi, turb, zone) = ovrr*ovrr*DACOS(ovl/ovrr)&
&               + ovr*ovr*DACOS((ovdyd-ovl)/ovr) - ovdyd*ovz
              CALL PUSHCONTROL2B(3)
            ELSE IF (ovrr .GT. ovr) THEN
              wakeoverlap(turbi, turb, zone) = pi*ovr*ovr
              CALL PUSHCONTROL2B(2)
            ELSE
              wakeoverlap(turbi, turb, zone) = pi*ovrr*ovrr
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE
            wakeoverlap(turbi, turb, zone) = 0.0_dp
            CALL PUSHCONTROL2B(0)
          END IF
        END DO
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
  END DO
  DO turb=1,nturbines
    DO turbi=1,nturbines
      wakeoverlap(turbi, turb, 3) = wakeoverlap(turbi, turb, 3) - &
&       wakeoverlap(turbi, turb, 2)
      wakeoverlap(turbi, turb, 2) = wakeoverlap(turbi, turb, 2) - &
&       wakeoverlap(turbi, turb, 1)
    END DO
  END DO
  wakeoverlaptrel_m = wakeoverlap
  DO turbi=nturbines,1,-1
    temp1 = pi*rotordiameter(turbi)**2
    DO nd=1,nbdirs
      tempb2(nd, :, :) = 4.0_dp*wakeoverlaptrel_mb(nd, turbi, :, :)/&
&       temp1
      rotordiameterb(nd, turbi) = rotordiameterb(nd, turbi) + pi*2*&
&       rotordiameter(turbi)*SUM(-(wakeoverlaptrel_m(turbi, :, :)*tempb2&
&       (nd, :, :)/temp1))
      wakeoverlaptrel_mb(nd, turbi, :, :) = tempb2(nd, :, :)
    END DO
  END DO
  DO nd=1,nbdirs
    wakeoverlapb(nd, :, :, :) = 0.0
    wakeoverlapb(nd, :, :, :) = wakeoverlaptrel_mb(nd, :, :, :)
  END DO
  DO turb=nturbines,1,-1
    DO turbi=nturbines,1,-1
      DO nd=1,nbdirs
        wakeoverlapb(nd, turbi, turb, 1) = wakeoverlapb(nd, turbi, turb&
&         , 1) - wakeoverlapb(nd, turbi, turb, 2)
        wakeoverlapb(nd, turbi, turb, 2) = wakeoverlapb(nd, turbi, turb&
&         , 2) - wakeoverlapb(nd, turbi, turb, 3)
      END DO
    END DO
  END DO
  DO turb=nturbines,1,-1
    DO turbi=nturbines,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        DO nd=1,nbdirs
          ovdydb(nd) = 0.0
          ovrb(nd) = 0.0
        END DO
        DO zone=3,1,-1
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              DO nd=1,nbdirs
                wakeoverlapb(nd, turbi, turb, zone) = 0.0
              END DO
              ovrr = wakediameters(turbi, turb, zone)/2.0_dp
              DO nd=1,nbdirs
                ovlb(nd) = 0.0
                ovrrb(nd) = 0.0
                ovzb(nd) = 0.0
              END DO
              GOTO 100
            ELSE
              ovrr = wakediameters(turbi, turb, zone)/2.0_dp
              DO nd=1,nbdirs
                ovrrb(nd) = pi*2*ovrr*wakeoverlapb(nd, turbi, turb, zone&
&                 )
                wakeoverlapb(nd, turbi, turb, zone) = 0.0
              END DO
            END IF
          ELSE IF (branch .EQ. 2) THEN
            DO nd=1,nbdirs
              ovrb(nd) = ovrb(nd) + pi*2*ovr*wakeoverlapb(nd, turbi, &
&               turb, zone)
              wakeoverlapb(nd, turbi, turb, zone) = 0.0
            END DO
            ovrr = wakediameters(turbi, turb, zone)/2.0_dp
            DO nd=1,nbdirs
              ovrrb(nd) = 0.0
            END DO
          ELSE
            ovrr = wakediameters(turbi, turb, zone)/2.0_dp
            temp = ovl/ovrr
            temp0 = (ovdyd-ovl)/ovr
            DO nd=1,nbdirs
              IF (temp .EQ. 1.0 .OR. temp .EQ. (-1.0)) THEN
                tempb0(nd) = 0.0
              ELSE
                tempb0(nd) = -(ovrr*wakeoverlapb(nd, turbi, turb, zone)/&
&                 SQRT(1.D0-temp**2))
              END IF
              IF (temp0 .EQ. 1.0 .OR. temp0 .EQ. (-1.0)) THEN
                tempb1(nd) = 0.0
              ELSE
                tempb1(nd) = -(ovr*wakeoverlapb(nd, turbi, turb, zone)/&
&                 SQRT(1.D0-temp0**2))
              END IF
              ovrrb(nd) = DACOS(temp)*2*ovrr*wakeoverlapb(nd, turbi, &
&               turb, zone) - temp*tempb0(nd)
              ovlb(nd) = tempb0(nd) - tempb1(nd)
              ovrb(nd) = ovrb(nd) + DACOS(temp0)*2*ovr*wakeoverlapb(nd, &
&               turbi, turb, zone) - temp0*tempb1(nd)
              ovdydb(nd) = ovdydb(nd) + tempb1(nd) - ovz*wakeoverlapb(nd&
&               , turbi, turb, zone)
              ovzb(nd) = -(ovdyd*wakeoverlapb(nd, turbi, turb, zone))
              wakeoverlapb(nd, turbi, turb, zone) = 0.0
            END DO
            GOTO 100
          END IF
          DO nd=1,nbdirs
            ovlb(nd) = 0.0
            ovzb(nd) = 0.0
          END DO
 100      CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4ARRAY(ovz, dp/4)
            DO nd=1,nbdirs
              IF (ovz .EQ. 0.0) THEN
                ovzb(nd) = 0.0
              ELSE
                ovzb(nd) = ovzb(nd)/(2.0*SQRT(ovz))
              END IF
            END DO
          ELSE
            DO nd=1,nbdirs
              ovzb(nd) = 0.0
            END DO
          END IF
          CALL POPREAL4ARRAY(ovz, dp/4)
          DO nd=1,nbdirs
            ovrrb(nd) = ovrrb(nd) + 2*ovrr*ovzb(nd)
            ovlb(nd) = ovlb(nd) - 2*ovl*ovzb(nd)
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4ARRAY(ovl, dp/4)
            DO nd=1,nbdirs
              tempb(nd) = ovlb(nd)/(2.0_dp*ovdyd)
              ovrrb(nd) = ovrrb(nd) + 2*ovrr*tempb(nd)
              ovrb(nd) = ovrb(nd) - 2*ovr*tempb(nd)
              ovdydb(nd) = ovdydb(nd) + (2*ovdyd-(ovrr**2-ovr**2+ovdyd**&
&               2)/ovdyd)*tempb(nd)
            END DO
          ELSE
            CALL POPREAL4ARRAY(ovl, dp/4)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4ARRAY(ovdyd, dp/4)
          ELSE
            CALL POPREAL4ARRAY(ovdyd, dp/4)
            DO nd=1,nbdirs
              ovdydb(nd) = -ovdydb(nd)
            END DO
          END IF
          DO nd=1,nbdirs
            wakediametersb(nd, turbi, turb, zone) = wakediametersb(nd, &
&             turbi, turb, zone) + ovrrb(nd)/2.0_dp
          END DO
        END DO
        CALL POPREAL4ARRAY(ovr, dp/4)
        CALL POPREAL4ARRAY(ovdyd, dp/4)
        DO nd=1,nbdirs
          rotordiameterb(nd, turbi) = rotordiameterb(nd, turbi) + ovrb(&
&           nd)/2
          wakecentersb(nd, turbi, turb) = wakecentersb(nd, turbi, turb) &
&           + ovdydb(nd)
          turbineyb(nd, turbi) = turbineyb(nd, turbi) - ovdydb(nd)
        END DO
      END IF
    END DO
  END DO
END SUBROUTINE CALCOVERLAPAREAS_BV


!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5902M) - 15 Dec 2015 09:00
!
!  Differentiation of cttoaxialind in reverse (adjoint) mode:
!   gradient     of useful results: axial_induction ct
!   with respect to varying inputs: ct
! Flow field calculations have been intentionally left out to save development time.
! The flow field can be calculated using the pure python version of floris 
! This implementation is fully smooth and differentiable with the exception of a 
! discontinuity at the hub of each turbine. The discontinuity only presents issues if
! turbines are place within 1E-15 * rotor diameter of one another, which should never
! happen in a standard optimization problem with a separation constraint.
SUBROUTINE CTTOAXIALIND_BV(ct, ctb, nturbines, axial_induction, &
& axial_inductionb, nbdirs)
!  USE DIFFSIZES
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: ct
  REAL(dp), DIMENSION(nbdirs, nturbines) :: ctb
! local
  INTEGER :: i
! out
  REAL(dp), DIMENSION(nturbines) :: axial_induction
  REAL(dp), DIMENSION(nbdirs, nturbines) :: axial_inductionb
  INTRINSIC KIND
  INTRINSIC SQRT
  INTEGER :: nd
  INTEGER :: branch
  INTEGER :: nbdirs
! execute
  DO i=1,nturbines
    IF (ct(i) .GT. 0.96) THEN
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
  DO i=nturbines,1,-1
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO nd=1,nbdirs
        IF (.NOT.1.0_dp - ct(i) .EQ. 0.0) ctb(nd, i) = ctb(nd, i) + &
&           0.5_dp*axial_inductionb(nd, i)/(2.0*SQRT(1.0_dp-ct(i)))
        axial_inductionb(nd, i) = 0.0
      END DO
    ELSE
      DO nd=1,nbdirs
        IF (.NOT.0.0203_dp - 0.6427_dp*(0.889_dp-ct(i)) .EQ. 0.0) ctb(nd&
&         , i) = ctb(nd, i) + 0.6427_dp*axial_inductionb(nd, i)/(2.0*&
&           SQRT(0.0203_dp-0.6427_dp*(0.889_dp-ct(i))))
        axial_inductionb(nd, i) = 0.0
      END DO
    END IF
  END DO
END SUBROUTINE CTTOAXIALIND_BV


!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5902M) - 15 Dec 2015 09:00
!
!  Differentiation of floris_wcent_wdiam in reverse (adjoint) mode:
!   gradient     of useful results: wakediameterst_vec wakecentersyt_vec
!   with respect to varying inputs: rotordiameter turbinexw yaw_deg
!                turbineyw wakediameterst_vec ct wakecentersyt_vec
!   RW status of diff variables: rotordiameter:out turbinexw:out
!                yaw_deg:out turbineyw:out wakediameterst_vec:in-out
!                ct:out wakecentersyt_vec:in-out
SUBROUTINE FLORIS_WCENT_WDIAM_BV(nturbines, kd, initialwakedisplacement&
& , initialwakeangle, ke_in, kecorrct, region2ct, yaw_deg, yaw_degb, ct&
& , ctb, turbinexw, turbinexwb, turbineyw, turbineywb, rotordiameter, &
& rotordiameterb, me, bd, usewakeangle, adjustinitialwakediamtoyaw, &
& wakecentersyt_vecb, &
& wakediameterst_vecb, nbdirs)
!  USE DIFFSIZES
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), INTENT(IN) :: kd, initialwakedisplacement, initialwakeangle&
& , ke_in
  REAL(dp), INTENT(IN) :: kecorrct, region2ct, bd
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: yaw_deg, ct, turbinexw, &
& turbineyw
  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: yaw_degb, ctb, turbinexwb&
& , turbineywb
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: rotordiameter
  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: rotordiameterb
  REAL(dp), DIMENSION(3), INTENT(IN) :: me
  LOGICAL, INTENT(IN) :: usewakeangle, adjustinitialwakediamtoyaw
! out
  REAL(dp), DIMENSION(nturbines*nturbines) :: wakecentersyt_vec
  REAL(dp), DIMENSION(nbdirs, nturbines*nturbines) :: &
& wakecentersyt_vecb
  REAL(dp), DIMENSION(3*nturbines*nturbines) :: wakediameterst_vec
  REAL(dp), DIMENSION(nbdirs, 3*nturbines*nturbines) :: &
& wakediameterst_vecb
! local
! in rotor diameters
  REAL(dp) :: spline_bound
  REAL(dp), DIMENSION(nturbines) :: ke, yaw
  REAL(dp), DIMENSION(nbdirs, nturbines) :: keb, yawb
  INTEGER :: turb, turbi, zone
  REAL(dp) :: wakeangleinit, zeroloc
  REAL(dp), DIMENSION(nbdirs) :: wakeangleinitb, zerolocb
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp
  REAL(dp) :: deltax, factor, displacement, x, x1, x2, y1, y2, dy1, dy2
  REAL(dp), DIMENSION(nbdirs) :: deltaxb, factorb, displacementb, xb&
& , x1b, x2b, y1b, dy1b
  REAL(dp) :: wakediameter0
  REAL(dp), DIMENSION(nbdirs) :: wakediameter0b
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakediameterst_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakediameterst_matb
  REAL(dp), DIMENSION(nturbines, nturbines) :: wakecentersyt_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines) :: &
& wakecentersyt_matb
  INTRINSIC KIND
  INTRINSIC SIN
  INTRINSIC COS
  INTEGER :: nd
  INTEGER :: branch
  INTEGER :: nbdirs
  REAL(dp) :: temp2
  REAL(dp) :: temp1
  REAL(dp) :: temp0
  REAL(dp) :: tempb5(nbdirs)
  REAL(dp) :: tempb4(nbdirs)
  REAL(dp) :: tempb3(nbdirs)
  REAL(dp) :: tempb2(nbdirs)
  REAL(dp) :: tempb1(nbdirs)
  REAL(dp) :: tempb0(nbdirs)
  REAL(dp) :: tempb(nbdirs)
  REAL(dp) :: temp
  spline_bound = 1.0_dp
! execute
!    if (CTcorrected) then
!        Ct = Ct_in
!    else
!        Ct = Ct_in*cos(yaw)*cos(yaw)
!    end if
! convert yaw from degrees to radians
  yaw = yaw_deg*pi/180.0_dp
! calculate y-locations of wake centers in wind ref. frame
  DO turb=1,nturbines
    CALL PUSHREAL4ARRAY(wakeangleinit, dp/4)
    wakeangleinit = 0.5_dp*SIN(yaw(turb))*ct(turb)
    IF (usewakeangle) wakeangleinit = wakeangleinit + initialwakeangle*&
&       pi/180.0_dp
    DO turbi=1,nturbines
      CALL PUSHREAL4ARRAY(deltax, dp/4)
      deltax = turbinexw(turbi) - turbinexw(turb)
      factor = 2.0_dp*kd*deltax/rotordiameter(turb) + 1.0_dp
      IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
! yaw-induced wake center displacement   
        IF (usewakeangle .EQV. .false.) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
  END DO
!adjust k_e to C_T, adjusted to yaw
  ke = ke_in + kecorrct*(ct-region2ct)
  DO turb=1,nturbines
    IF (adjustinitialwakediamtoyaw) THEN
      CALL PUSHREAL4ARRAY(wakediameter0, dp/4)
      wakediameter0 = rotordiameter(turb)*COS(yaw(turb))
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHREAL4ARRAY(wakediameter0, dp/4)
      wakediameter0 = rotordiameter(turb)
      CALL PUSHCONTROL1B(0)
    END IF
    DO turbi=1,nturbines
! turbine separation
      CALL PUSHREAL4ARRAY(deltax, dp/4)
      deltax = turbinexw(turbi) - turbinexw(turb)
! x position of interest
      x = turbinexw(turbi)
      CALL PUSHINTEGER4(zone)
      zone = 1
! define centerpoint of spline
      zeroloc = turbinexw(turb) - wakediameter0/(2.0_dp*ke(turb)*me(zone&
&       ))
      IF (zeroloc + spline_bound*rotordiameter(turb) .LT. turbinexw(&
&         turbi)) THEN
        CALL PUSHCONTROL2B(0)
      ELSE IF (zeroloc - spline_bound*rotordiameter(turb) .LT. turbinexw&
&         (turbi)) THEN
        CALL PUSHCONTROL2B(1)
      ELSE IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        CALL PUSHCONTROL2B(2)
      ELSE
        CALL PUSHCONTROL2B(3)
      END IF
      IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        CALL PUSHINTEGER4(zone)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
  END DO
  DO nd=1,nbdirs
    wakediameterst_matb(nd, :, :, :) = 0.0
  END DO
  DO turbi=nturbines,1,-1
    DO nd=1,nbdirs
      wakediameterst_matb(nd, turbi, :, 3) = wakediameterst_matb(nd, &
&       turbi, :, 3) + wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+2*&
&       nturbines+1:nturbines*(turbi-1)+3*nturbines)
      wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+2*nturbines+1:&
&     nturbines*(turbi-1)+3*nturbines) = 0.0
      wakediameterst_matb(nd, turbi, :, 2) = wakediameterst_matb(nd, &
&       turbi, :, 2) + wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+&
&       nturbines+1:3*nturbines*(turbi-1)+2*nturbines)
      wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+nturbines+1:3*&
&     nturbines*(turbi-1)+2*nturbines) = 0.0
      wakediameterst_matb(nd, turbi, :, 1) = wakediameterst_matb(nd, &
&       turbi, :, 1) + wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+1:3&
&       *nturbines*(turbi-1)+nturbines)
      wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+1:3*nturbines*(turbi&
&     -1)+nturbines) = 0.0
    END DO
  END DO
  DO nd=1,nbdirs
    rotordiameterb(nd, :) = 0.0
    turbinexwb(nd, :) = 0.0
    yawb(nd, :) = 0.0
    keb(nd, :) = 0.0
  END DO
  DO turb=nturbines,1,-1
    DO nd=1,nbdirs
      wakediameter0b(nd) = 0.0
    END DO
    DO turbi=nturbines,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO nd=1,nbdirs
          deltaxb(nd) = 0.0
        END DO
      ELSE
        deltax = turbinexw(turbi) - turbinexw(turb)
        zone = 3
        DO nd=1,nbdirs
          tempb4(nd) = me(zone)*2.0_dp*wakediameterst_matb(nd, turbi, &
&           turb, zone)
          wakediameter0b(nd) = wakediameter0b(nd) + wakediameterst_matb(&
&           nd, turbi, turb, zone)
          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
        END DO
        zone = 2
        DO nd=1,nbdirs
          tempb5(nd) = me(zone)*2.0_dp*wakediameterst_matb(nd, turbi, &
&           turb, zone)
          keb(nd, turb) = keb(nd, turb) + deltax*tempb5(nd) + deltax*&
&           tempb4(nd)
          deltaxb(nd) = ke(turb)*tempb5(nd) + ke(turb)*tempb4(nd)
          wakediameter0b(nd) = wakediameter0b(nd) + wakediameterst_matb(&
&           nd, turbi, turb, zone)
          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
        END DO
        CALL POPINTEGER4(zone)
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          DO nd=1,nbdirs
            wakediameterst_matb(nd, turbi, turb, zone) = 0.0
          END DO
          DO nd=1,nbdirs
            xb(nd) = 0.0
            zerolocb(nd) = 0.0
          END DO
        ELSE
          zone = 1
          dy1 = 2.0_dp*ke(turb)*me(zone)
          zeroloc = turbinexw(turb) - wakediameter0/(2.0_dp*ke(turb)*me(&
&           zone))
          x1 = zeroloc - spline_bound*rotordiameter(turb)
          y1 = wakediameter0 + 2.0_dp*ke(turb)*me(zone)*(x1-turbinexw(&
&           turb))
          dy2 = 0.0_dp
          y2 = 0.0_dp
          x = turbinexw(turbi)
          x2 = zeroloc + spline_bound*rotordiameter(turb)
          CALL HERMITE_SPLINE_BV(x, xb, x1, x1b, x2, x2b, y1, y1b, dy1, &
&                          dy1b, y2, dy2, wakediameterst_mat(turbi, turb&
&                          , zone), wakediameterst_matb(1, turbi, turb, &
&                          zone), nbdirs)
          DO nd=1,nbdirs
            tempb2(nd) = me(zone)*2.0_dp*y1b(nd)
            x1b(nd) = x1b(nd) + ke(turb)*tempb2(nd)
            wakediameterst_matb(nd, turbi, turb, zone) = 0.0
            zerolocb(nd) = x1b(nd) + x2b(nd)
            rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&             spline_bound*x2b(nd) - spline_bound*x1b(nd)
            keb(nd, turb) = keb(nd, turb) + (x1-turbinexw(turb))*tempb2(&
&             nd) + me(zone)*2.0_dp*dy1b(nd)
            wakediameter0b(nd) = wakediameter0b(nd) + y1b(nd)
            turbinexwb(nd, turb) = turbinexwb(nd, turb) - ke(turb)*&
&             tempb2(nd)
          END DO
        END IF
      ELSE
        IF (branch .EQ. 2) THEN
          DO nd=1,nbdirs
            tempb3(nd) = me(zone)*2.0_dp*wakediameterst_matb(nd, turbi, &
&             turb, zone)
            wakediameter0b(nd) = wakediameter0b(nd) + &
&             wakediameterst_matb(nd, turbi, turb, zone)
            keb(nd, turb) = keb(nd, turb) + deltax*tempb3(nd)
            deltaxb(nd) = deltaxb(nd) + ke(turb)*tempb3(nd)
            wakediameterst_matb(nd, turbi, turb, zone) = 0.0
          END DO
        END IF
        DO nd=1,nbdirs
          xb(nd) = 0.0
          zerolocb(nd) = 0.0
        END DO
      END IF
      temp2 = 2.0_dp*me(zone)*ke(turb)
      DO nd=1,nbdirs
        turbinexwb(nd, turb) = turbinexwb(nd, turb) + zerolocb(nd)
        wakediameter0b(nd) = wakediameter0b(nd) - zerolocb(nd)/temp2
        keb(nd, turb) = keb(nd, turb) + wakediameter0*2.0_dp*me(zone)*&
&         zerolocb(nd)/temp2**2
        turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + deltaxb(nd) + xb&
&         (nd)
        turbinexwb(nd, turb) = turbinexwb(nd, turb) - deltaxb(nd)
      END DO
      CALL POPINTEGER4(zone)
      CALL POPREAL4ARRAY(deltax, dp/4)
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4ARRAY(wakediameter0, dp/4)
      DO nd=1,nbdirs
        rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&         wakediameter0b(nd)
      END DO
    ELSE
      CALL POPREAL4ARRAY(wakediameter0, dp/4)
      DO nd=1,nbdirs
        rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + COS(yaw(&
&         turb))*wakediameter0b(nd)
        yawb(nd, turb) = yawb(nd, turb) - rotordiameter(turb)*SIN(yaw(&
&         turb))*wakediameter0b(nd)
      END DO
    END IF
  END DO
  DO nd=1,nbdirs
    ctb(nd, :) = 0.0
    ctb(nd, :) = kecorrct*keb(nd, :)
  END DO
  DO nd=1,nbdirs
    wakecentersyt_matb(nd, :, :) = 0.0
  END DO
  DO turbi=nturbines,1,-1
    DO nd=1,nbdirs
      wakecentersyt_matb(nd, turbi, :) = wakecentersyt_matb(nd, turbi, :&
&       ) + wakecentersyt_vecb(nd, nturbines*(turbi-1)+1:nturbines*(&
&       turbi-1)+nturbines)
      wakecentersyt_vecb(nd, nturbines*(turbi-1)+1:nturbines*(turbi-1)+&
&     nturbines) = 0.0
    END DO
  END DO
  DO nd=1,nbdirs
    turbineywb(nd, :) = 0.0
  END DO
  DO turb=nturbines,1,-1
    DO nd=1,nbdirs
      wakeangleinitb(nd) = 0.0
    END DO
    DO turbi=nturbines,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO nd=1,nbdirs
          turbineywb(nd, turb) = turbineywb(nd, turb) + &
&           wakecentersyt_matb(nd, turbi, turb)
          wakecentersyt_matb(nd, turbi, turb) = 0.0
        END DO
        DO nd=1,nbdirs
          deltaxb(nd) = 0.0
          factorb(nd) = 0.0
        END DO
      ELSE
        DO nd=1,nbdirs
          displacementb(nd) = wakecentersyt_matb(nd, turbi, turb)
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO nd=1,nbdirs
            deltaxb(nd) = bd*displacementb(nd)
          END DO
        ELSE
          DO nd=1,nbdirs
            deltaxb(nd) = 0.0
          END DO
        END IF
        deltax = turbinexw(turbi) - turbinexw(turb)
        factor = 2.0_dp*kd*deltax/rotordiameter(turb) + 1.0_dp
        temp1 = 30.0_dp*kd*factor**5
        temp0 = wakeangleinit*rotordiameter(turb)
        temp = 15.0_dp*factor**4 + wakeangleinit**2
        DO nd=1,nbdirs
          tempb0(nd) = displacementb(nd)/temp1
          tempb1(nd) = -((wakeangleinit**2+15.0_dp)*displacementb(nd)/(&
&           30.0_dp*kd))
          factorb(nd) = (15.0_dp*temp0*4*factor**3-30.0_dp*kd*temp*temp0&
&           *5*factor**4/temp1)*tempb0(nd)
          wakeangleinitb(nd) = wakeangleinitb(nd) + rotordiameter(turb)*&
&           tempb1(nd) - wakeangleinit**2*rotordiameter(turb)*2*&
&           displacementb(nd)/(30.0_dp*kd) + (temp*rotordiameter(turb)+&
&           temp0*2*wakeangleinit)*tempb0(nd)
          rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&           wakeangleinit*tempb1(nd) + temp*wakeangleinit*tempb0(nd)
          turbineywb(nd, turb) = turbineywb(nd, turb) + &
&           wakecentersyt_matb(nd, turbi, turb)
          wakecentersyt_matb(nd, turbi, turb) = 0.0
        END DO
      END IF
      DO nd=1,nbdirs
        tempb(nd) = kd*2.0_dp*factorb(nd)/rotordiameter(turb)
        deltaxb(nd) = deltaxb(nd) + tempb(nd)
        rotordiameterb(nd, turb) = rotordiameterb(nd, turb) - deltax*&
&         tempb(nd)/rotordiameter(turb)
        turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + deltaxb(nd)
        turbinexwb(nd, turb) = turbinexwb(nd, turb) - deltaxb(nd)
      END DO
      CALL POPREAL4ARRAY(deltax, dp/4)
    END DO
    CALL POPREAL4ARRAY(wakeangleinit, dp/4)
    DO nd=1,nbdirs
      yawb(nd, turb) = yawb(nd, turb) + ct(turb)*0.5_dp*COS(yaw(turb))*&
&       wakeangleinitb(nd)
      ctb(nd, turb) = ctb(nd, turb) + 0.5_dp*SIN(yaw(turb))*&
&       wakeangleinitb(nd)
    END DO
  END DO
  DO nd=1,nbdirs
    yaw_degb(nd, :) = 0.0
    yaw_degb(nd, :) = pi*yawb(nd, :)/180.0_dp
  END DO
END SUBROUTINE FLORIS_WCENT_WDIAM_BV


!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5902M) - 15 Dec 2015 09:00
!
!  Differentiation of floris_overlap in reverse (adjoint) mode:
!   gradient     of useful results: cosfac_vec wakeoverlaptrel_vec
!   with respect to varying inputs: rotordiameter turbineyw cosfac_vec
!                wakeoverlaptrel_vec wakediameterst_vec wakecentersyt_vec
!   RW status of diff variables: rotordiameter:out turbineyw:out
!                cosfac_vec:in-out wakeoverlaptrel_vec:in-zero
!                wakediameterst_vec:out wakecentersyt_vec:out
SUBROUTINE FLORIS_OVERLAP_BV(nturbines, turbinexw, turbineyw, turbineywb&
& , rotordiameter, rotordiameterb, wakediameterst_vec, &
& wakediameterst_vecb, wakecentersyt_vec, wakecentersyt_vecb, cos_spread&
& , wakeoverlaptrel_vecb, cosfac_vecb, &
& nbdirs)
!  USE DIFFSIZES
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), INTENT(IN) :: cos_spread
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: turbinexw, turbineyw, &
& rotordiameter
  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: turbineywb, &
& rotordiameterb
  REAL(dp), DIMENSION(3*nturbines*nturbines), INTENT(IN) :: &
& wakediameterst_vec
  REAL(dp), DIMENSION(nbdirs, 3*nturbines*nturbines), intent(out) :: &
& wakediameterst_vecb
  REAL(dp), DIMENSION(nturbines*nturbines), INTENT(IN) :: &
& wakecentersyt_vec
  REAL(dp), DIMENSION(nbdirs, nturbines*nturbines), intent(out) :: &
& wakecentersyt_vecb
! out
  REAL(dp), DIMENSION(3*nturbines*nturbines) :: wakeoverlaptrel_vec
  REAL(dp), DIMENSION(nbdirs, 3*nturbines*nturbines) :: &
& wakeoverlaptrel_vecb
  REAL(dp), DIMENSION(3*nturbines*nturbines) :: cosfac_vec
  REAL(dp), DIMENSION(nbdirs, 3*nturbines*nturbines) :: cosfac_vecb
! local
  REAL(dp) :: rmax
  REAL(dp), DIMENSION(nbdirs) :: rmaxb
  INTEGER :: turbi, turb, zone
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlaptrel_mat, &
& wakediameterst_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlaptrel_matb, wakediameterst_matb
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: cosfac_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: cosfac_matb
  REAL(dp), DIMENSION(nturbines, nturbines) :: wakecentersyt_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines) :: &
& wakecentersyt_matb
  INTRINSIC KIND
  INTRINSIC DABS
  INTRINSIC COS
  INTEGER :: nd
  INTEGER :: branch
  INTEGER :: nbdirs
  REAL :: tempb0(nbdirs)
  DOUBLE PRECISION :: dabs0b(nbdirs)
  REAL(dp) :: tempb(nbdirs)
  DOUBLE PRECISION :: dabs0
! execute    
! pack wakeDiametersT_v vector into a matrix
  DO turbi=1,nturbines
    wakediameterst_mat(turbi, :, 1) = wakediameterst_vec(3*nturbines*(&
&     turbi-1)+1:3*nturbines*(turbi-1)+nturbines)
    wakediameterst_mat(turbi, :, 2) = wakediameterst_vec(3*nturbines*(&
&     turbi-1)+nturbines+1:3*nturbines*(turbi-1)+2*nturbines)
    wakediameterst_mat(turbi, :, 3) = wakediameterst_vec(3*nturbines*(&
&     turbi-1)+2*nturbines+1:3*nturbines*(turbi-1)+3*nturbines)
  END DO
! pack wakeCentersYT_v vector into a matrix
  DO turbi=1,nturbines
    wakecentersyt_mat(turbi, :) = wakecentersyt_vec(nturbines*(turbi-1)+&
&     1:nturbines*(turbi-1)+nturbines)
  END DO
! calculate cosine factor
  DO turbi=1,nturbines
    DO turb=1,nturbines
      DO zone=1,3
        IF (wakecentersyt_mat(turbi, turb) - turbineyw(turbi) .GE. 0.) &
&       THEN
          CALL PUSHREAL8(dabs0)
          dabs0 = wakecentersyt_mat(turbi, turb) - turbineyw(turbi)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHREAL8(dabs0)
          dabs0 = -(wakecentersyt_mat(turbi, turb)-turbineyw(turbi))
          CALL PUSHCONTROL1B(1)
        END IF
      END DO
    END DO
  END DO
  DO nd=1,nbdirs
    cosfac_matb(nd, :, :, :) = 0.0
    wakeoverlaptrel_matb(nd, :, :, :) = 0.0
  END DO
  DO turbi=nturbines,1,-1
    DO nd=1,nbdirs
      cosfac_matb(nd, turbi, :, 3) = cosfac_matb(nd, turbi, :, 3) + &
&       cosfac_vecb(nd, 3*nturbines*(turbi-1)+2*nturbines+1:3*nturbines*&
&       (turbi-1)+3*nturbines)
      cosfac_vecb(nd, 3*nturbines*(turbi-1)+2*nturbines+1:3*nturbines*(&
&     turbi-1)+3*nturbines) = 0.0
      cosfac_matb(nd, turbi, :, 2) = cosfac_matb(nd, turbi, :, 2) + &
&       cosfac_vecb(nd, 3*nturbines*(turbi-1)+nturbines+1:3*nturbines*(&
&       turbi-1)+2*nturbines)
      cosfac_vecb(nd, 3*nturbines*(turbi-1)+nturbines+1:3*nturbines*(&
&     turbi-1)+2*nturbines) = 0.0
      cosfac_matb(nd, turbi, :, 1) = cosfac_matb(nd, turbi, :, 1) + &
&       cosfac_vecb(nd, 3*nturbines*(turbi-1)+1:3*nturbines*(turbi-1)+&
&       nturbines)
      cosfac_vecb(nd, 3*nturbines*(turbi-1)+1:3*nturbines*(turbi-1)+&
&     nturbines) = 0.0
      wakeoverlaptrel_matb(nd, turbi, :, 3) = wakeoverlaptrel_matb(nd, &
&       turbi, :, 3) + wakeoverlaptrel_vecb(nd, 3*nturbines*(turbi-1)+2*&
&       nturbines+1:3*nturbines*(turbi-1)+3*nturbines)
      wakeoverlaptrel_vecb(nd, 3*nturbines*(turbi-1)+2*nturbines+1:3*&
&     nturbines*(turbi-1)+3*nturbines) = 0.0
      wakeoverlaptrel_matb(nd, turbi, :, 2) = wakeoverlaptrel_matb(nd, &
&       turbi, :, 2) + wakeoverlaptrel_vecb(nd, 3*nturbines*(turbi-1)+&
&       nturbines+1:3*nturbines*(turbi-1)+2*nturbines)
      wakeoverlaptrel_vecb(nd, 3*nturbines*(turbi-1)+nturbines+1:3*&
&     nturbines*(turbi-1)+2*nturbines) = 0.0
      wakeoverlaptrel_matb(nd, turbi, :, 1) = wakeoverlaptrel_matb(nd, &
&       turbi, :, 1) + wakeoverlaptrel_vecb(nd, 3*nturbines*(turbi-1)+1:&
&       3*nturbines*(turbi-1)+nturbines)
      wakeoverlaptrel_vecb(nd, 3*nturbines*(turbi-1)+1:3*nturbines*(&
&     turbi-1)+nturbines) = 0.0
    END DO
  END DO
  DO nd=1,nbdirs
    rotordiameterb(nd, :) = 0.0
    turbineywb(nd, :) = 0.0
    wakediameterst_matb(nd, :, :, :) = 0.0
    wakecentersyt_matb(nd, :, :) = 0.0
  END DO
  DO turbi=nturbines,1,-1
    DO turb=nturbines,1,-1
      DO zone=3,1,-1
        rmax = cos_spread*0.5_dp*(wakediameterst_mat(turbi, turb, 3)+&
&         rotordiameter(turbi))
        DO nd=1,nbdirs
          tempb0(nd) = -(pi*SIN(pi*(dabs0/rmax))*0.5_dp*cosfac_matb(nd, &
&           turbi, turb, zone)/rmax)
          dabs0b(nd) = tempb0(nd)
          rmaxb(nd) = -(dabs0*tempb0(nd)/rmax)
          cosfac_matb(nd, turbi, turb, zone) = 0.0
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(dabs0)
          DO nd=1,nbdirs
            wakecentersyt_matb(nd, turbi, turb) = wakecentersyt_matb(nd&
&             , turbi, turb) + dabs0b(nd)
            turbineywb(nd, turbi) = turbineywb(nd, turbi) - dabs0b(nd)
          END DO
        ELSE
          CALL POPREAL8(dabs0)
          DO nd=1,nbdirs
            turbineywb(nd, turbi) = turbineywb(nd, turbi) + dabs0b(nd)
            wakecentersyt_matb(nd, turbi, turb) = wakecentersyt_matb(nd&
&             , turbi, turb) - dabs0b(nd)
          END DO
        END IF
        DO nd=1,nbdirs
          tempb(nd) = cos_spread*0.5_dp*rmaxb(nd)
          wakediameterst_matb(nd, turbi, turb, 3) = wakediameterst_matb(&
&           nd, turbi, turb, 3) + tempb(nd)
          rotordiameterb(nd, turbi) = rotordiameterb(nd, turbi) + tempb(&
&           nd)
        END DO
      END DO
    END DO
  END DO
  CALL CALCOVERLAPAREAS_BV(nturbines, turbinexw, turbineyw, turbineywb, &
&                    rotordiameter, rotordiameterb, wakediameterst_mat, &
&                    wakediameterst_matb, wakecentersyt_mat, &
&                    wakecentersyt_matb, wakeoverlaptrel_mat, &
&                    wakeoverlaptrel_matb, nbdirs)
  DO nd=1,nbdirs
    wakecentersyt_vecb(nd, :) = 0.0
  END DO
  DO turbi=nturbines,1,-1
    DO nd=1,nbdirs
      wakecentersyt_vecb(nd, nturbines*(turbi-1)+1:nturbines*(turbi-1)+&
&     nturbines) = wakecentersyt_vecb(nd, nturbines*(turbi-1)+1:&
&       nturbines*(turbi-1)+nturbines) + wakecentersyt_matb(nd, turbi, :&
&       )
      wakecentersyt_matb(nd, turbi, :) = 0.0
    END DO
  END DO
  DO nd=1,nbdirs
    wakediameterst_vecb(nd, :) = 0.0
  END DO
  DO turbi=nturbines,1,-1
    DO nd=1,nbdirs
      wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+2*nturbines+1:3*&
&     nturbines*(turbi-1)+3*nturbines) = wakediameterst_vecb(nd, 3*&
&       nturbines*(turbi-1)+2*nturbines+1:3*nturbines*(turbi-1)+3*&
&       nturbines) + wakediameterst_matb(nd, turbi, :, 3)
      wakediameterst_matb(nd, turbi, :, 3) = 0.0
      wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+nturbines+1:3*&
&     nturbines*(turbi-1)+2*nturbines) = wakediameterst_vecb(nd, 3*&
&       nturbines*(turbi-1)+nturbines+1:3*nturbines*(turbi-1)+2*&
&       nturbines) + wakediameterst_matb(nd, turbi, :, 2)
      wakediameterst_matb(nd, turbi, :, 2) = 0.0
      wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+1:3*nturbines*(turbi&
&     -1)+nturbines) = wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+1:3&
&       *nturbines*(turbi-1)+nturbines) + wakediameterst_matb(nd, turbi&
&       , :, 1)
      wakediameterst_matb(nd, turbi, :, 1) = 0.0
    END DO
  END DO
  DO nd=1,nbdirs
    wakeoverlaptrel_vecb(nd, :) = 0.0
  END DO
END SUBROUTINE FLORIS_OVERLAP_BV


!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5902M) - 15 Dec 2015 09:00
!
!  Differentiation of floris_velocity in reverse (adjoint) mode:
!   gradient     of useful results: velocitiesturbines
!   with respect to varying inputs: rotordiameter turbinexw yaw_deg
!                wakeoverlaptrel_v cosfac_v velocitiesturbines
!                ct a_in
!   RW status of diff variables: rotordiameter:out turbinexw:out
!                yaw_deg:out wakeoverlaptrel_v:out cosfac_v:out
!                velocitiesturbines:in-zero ct:out a_in:out
SUBROUTINE FLORIS_VELOCITY_BV(nturbines, wakeoverlaptrel_v, &
& wakeoverlaptrel_vb, cosfac_v, cosfac_vb, ct, ctb, a_in, a_inb, &
& axialindprovided, useaubu, kecorrct, region2ct, ke_in, vinf, &
& kecorrarray, turbinexw, turbinexwb, yaw_deg, yaw_degb, rotordiameter, &
& rotordiameterb, mu, au, bu, velocitiesturbinesb, nbdirs)
!  USE DIFFSIZES
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
! Number of turbines
  INTEGER, INTENT(IN) :: nturbines
! relative overlap and cosine factor (vector of length 3*nTurbines**2)
  REAL(dp), DIMENSION(3*nturbines*nturbines), INTENT(IN) :: &
& wakeoverlaptrel_v, cosfac_v
  REAL(dp), DIMENSION(nbdirs, 3*nturbines*nturbines), intent(out) :: &
& wakeoverlaptrel_vb, cosfac_vb
! thrust coeff, axial induction, turbine yaw
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: ct, a_in
  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: ctb, a_inb
! option logicals
  LOGICAL, INTENT(IN) :: axialindprovided, useaubu
  REAL(dp), INTENT(IN) :: kecorrct, region2ct, ke_in, vinf, kecorrarray
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: turbinexw, yaw_deg
  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: turbinexwb, yaw_degb
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: rotordiameter
  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: rotordiameterb
  REAL(dp), DIMENSION(3), INTENT(IN) :: mu
  REAL(dp), INTENT(IN) :: au, bu
! local
! relative overlap and Cosine factor (NxNx3)
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlaptrel, &
& cosfac
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlaptrelb, cosfacb
  REAL(dp), DIMENSION(nturbines) :: a, ke, kearray, yaw
  REAL(dp), DIMENSION(nbdirs, nturbines) :: ab, keb, kearrayb, yawb
  REAL(dp), DIMENSION(3) :: mmu
  REAL(dp), DIMENSION(nbdirs, 3) :: mmub
  REAL(dp) :: s, wakeeffcoeff, wakeeffcoeffperzone, deltax
  REAL(dp), DIMENSION(nbdirs) :: sb, wakeeffcoeffb, &
& wakeeffcoeffperzoneb, deltaxb
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp
  INTEGER :: turb, turbi, zone
! out
  REAL(dp), DIMENSION(nturbines) :: velocitiesturbines
  REAL(dp), DIMENSION(nbdirs, nturbines) :: velocitiesturbinesb
  INTRINSIC COS
  INTRINSIC KIND
  INTRINSIC SUM
  INTRINSIC SQRT
  INTEGER :: nd
  INTEGER :: branch
  INTEGER :: nbdirs
  REAL(dp) :: temp3
  REAL(dp) :: temp2
  REAL(dp) :: temp1
  REAL(dp) :: temp0
  REAL(dp) :: tempb4(nbdirs)
  REAL(dp) :: tempb3(nbdirs)
  REAL(dp) :: tempb2(nbdirs)
  REAL(dp) :: tempb1(nbdirs)
  REAL(dp) :: tempb0(nbdirs)
  REAL(dp) :: tempb(nbdirs)
  REAL(dp) :: temp
  REAL(dp) :: temp7
  REAL(dp) :: temp6
  REAL(dp) :: temp5
  REAL(dp) :: temp4
! convert yaw from degrees to radians
  yaw = yaw_deg*pi/180.0_dp
! pack overlap and cosine factor vectors into matrices
  DO turbi=1,nturbines
    wakeoverlaptrel(turbi, :, 1) = wakeoverlaptrel_v(3*nturbines*(turbi-&
&     1)+1:3*nturbines*(turbi-1)+nturbines)
    wakeoverlaptrel(turbi, :, 2) = wakeoverlaptrel_v(3*nturbines*(turbi-&
&     1)+nturbines+1:3*nturbines*(turbi-1)+2*nturbines)
    wakeoverlaptrel(turbi, :, 3) = wakeoverlaptrel_v(3*nturbines*(turbi-&
&     1)+2*nturbines+1:3*nturbines*(turbi-1)+3*nturbines)
    cosfac(turbi, :, 1) = cosfac_v(3*nturbines*(turbi-1)+1:3*nturbines*(&
&     turbi-1)+nturbines)
    cosfac(turbi, :, 2) = cosfac_v(3*nturbines*(turbi-1)+nturbines+1:3*&
&     nturbines*(turbi-1)+2*nturbines)
    cosfac(turbi, :, 3) = cosfac_v(3*nturbines*(turbi-1)+2*nturbines+1:3&
&     *nturbines*(turbi-1)+3*nturbines)
  END DO
  IF (axialindprovided) THEN
    a = a_in
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL CTTOAXIALIND(ct, nturbines, a)
    CALL PUSHCONTROL1B(1)
  END IF
! adjust ke to Ct as adjusted to yaw
  ke = ke_in + kecorrct*(ct-region2ct)
  DO turb=1,nturbines
    CALL PUSHREAL4ARRAY(s, dp/4)
    s = SUM(wakeoverlaptrel(turb, :, 1) + wakeoverlaptrel(turb, :, 2))
    kearray(turb) = ke(turb)*(1+s*kecorrarray)
  END DO
! find effective wind speeds at downstream turbines, then predict the power of 
! downstream turbine    
  velocitiesturbines = vinf
  DO turbi=1,nturbines
    CALL PUSHREAL4ARRAY(wakeeffcoeff, dp/4)
    wakeeffcoeff = 0.0_dp
! find overlap-area weighted effect of each wake zone
    DO turb=1,nturbines
      CALL PUSHREAL4ARRAY(wakeeffcoeffperzone, dp/4)
      wakeeffcoeffperzone = 0.0_dp
      deltax = turbinexw(turbi) - turbinexw(turb)
      IF (useaubu) THEN
        CALL PUSHREAL4ARRAY(mmu, dp*3/4)
        mmu = mu/COS(au*pi/180.0_dp+bu*yaw(turb))
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (deltax .GT. 0.0_dp) THEN
        DO zone=1,3
          IF (useaubu) THEN
            wakeeffcoeffperzone = wakeeffcoeffperzone + (rotordiameter(&
&             turb)/(rotordiameter(turb)+2.0_dp*kearray(turb)*mmu(zone)*&
&             deltax)*cosfac(turbi, turb, zone))**2*wakeoverlaptrel(&
&             turbi, turb, zone)
            CALL PUSHCONTROL1B(1)
          ELSE
            wakeeffcoeffperzone = wakeeffcoeffperzone + (rotordiameter(&
&             turb)/(rotordiameter(turb)+2.0_dp*kearray(turb)*mu(zone)*&
&             deltax)*cosfac(turbi, turb, zone))**2*wakeoverlaptrel(&
&             turbi, turb, zone)
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        wakeeffcoeff = wakeeffcoeff + (a(turb)*wakeeffcoeffperzone)**2
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    CALL PUSHREAL4ARRAY(wakeeffcoeff, dp/4)
    wakeeffcoeff = 1.0_dp - 2.0_dp*SQRT(wakeeffcoeff)
! multiply the inflow speed with the wake coefficients to find effective wind 
! speed at turbine
  END DO
  DO nd=1,nbdirs
    rotordiameterb(nd, :) = 0.0
    turbinexwb(nd, :) = 0.0
    kearrayb(nd, :) = 0.0
    wakeoverlaptrelb(nd, :, :, :) = 0.0
    yawb(nd, :) = 0.0
    mmub(nd, :) = 0.0
    cosfacb(nd, :, :, :) = 0.0
    ab(nd, :) = 0.0
  END DO
  DO turbi=nturbines,1,-1
    DO nd=1,nbdirs
      wakeeffcoeffb(nd) = velocitiesturbines(turbi)*velocitiesturbinesb(&
&       nd, turbi)
      velocitiesturbinesb(nd, turbi) = wakeeffcoeff*velocitiesturbinesb(&
&       nd, turbi)
    END DO
    CALL POPREAL4ARRAY(wakeeffcoeff, dp/4)
    DO nd=1,nbdirs
      IF (wakeeffcoeff .EQ. 0.0) THEN
        wakeeffcoeffb(nd) = 0.0
      ELSE
        wakeeffcoeffb(nd) = -(2.0_dp*wakeeffcoeffb(nd)/(2.0*SQRT(&
&         wakeeffcoeff)))
      END IF
    END DO
    DO turb=nturbines,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO nd=1,nbdirs
          deltaxb(nd) = 0.0
        END DO
      ELSE
        DO nd=1,nbdirs
          tempb4(nd) = 2*a(turb)*wakeeffcoeffperzone*wakeeffcoeffb(nd)
          ab(nd, turb) = ab(nd, turb) + wakeeffcoeffperzone*tempb4(nd)
          wakeeffcoeffperzoneb(nd) = a(turb)*tempb4(nd)
        END DO
        deltax = turbinexw(turbi) - turbinexw(turb)
        DO nd=1,nbdirs
          deltaxb(nd) = 0.0
        END DO
        DO zone=3,1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            temp7 = 2.0_dp*mu(zone)
            temp4 = rotordiameter(turb) + temp7*kearray(turb)*deltax
            temp6 = rotordiameter(turb)**2*wakeoverlaptrel(turbi, turb, &
&             zone)
            temp5 = cosfac(turbi, turb, zone)**2
            DO nd=1,nbdirs
              tempb2(nd) = wakeeffcoeffperzoneb(nd)/temp4**2
              tempb3(nd) = -(temp5*temp6*2*tempb2(nd)/temp4)
              cosfacb(nd, turbi, turb, zone) = cosfacb(nd, turbi, turb, &
&               zone) + temp6*2*cosfac(turbi, turb, zone)*tempb2(nd)
              rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&               tempb3(nd) + wakeoverlaptrel(turbi, turb, zone)*temp5*2*&
&               rotordiameter(turb)*tempb2(nd)
              wakeoverlaptrelb(nd, turbi, turb, zone) = wakeoverlaptrelb&
&               (nd, turbi, turb, zone) + temp5*rotordiameter(turb)**2*&
&               tempb2(nd)
              kearrayb(nd, turb) = kearrayb(nd, turb) + temp7*deltax*&
&               tempb3(nd)
              deltaxb(nd) = deltaxb(nd) + temp7*kearray(turb)*tempb3(nd)
            END DO
          ELSE
            temp1 = rotordiameter(turb) + 2.0_dp*kearray(turb)*deltax*&
&             mmu(zone)
            temp3 = rotordiameter(turb)**2*wakeoverlaptrel(turbi, turb, &
&             zone)
            temp2 = cosfac(turbi, turb, zone)**2
            DO nd=1,nbdirs
              tempb(nd) = wakeeffcoeffperzoneb(nd)/temp1**2
              tempb0(nd) = -(temp2*temp3*2*tempb(nd)/temp1)
              tempb1(nd) = 2.0_dp*mmu(zone)*tempb0(nd)
              cosfacb(nd, turbi, turb, zone) = cosfacb(nd, turbi, turb, &
&               zone) + temp3*2*cosfac(turbi, turb, zone)*tempb(nd)
              rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&               tempb0(nd) + wakeoverlaptrel(turbi, turb, zone)*temp2*2*&
&               rotordiameter(turb)*tempb(nd)
              wakeoverlaptrelb(nd, turbi, turb, zone) = wakeoverlaptrelb&
&               (nd, turbi, turb, zone) + temp2*rotordiameter(turb)**2*&
&               tempb(nd)
              kearrayb(nd, turb) = kearrayb(nd, turb) + deltax*tempb1(nd&
&               )
              deltaxb(nd) = deltaxb(nd) + kearray(turb)*tempb1(nd)
              mmub(nd, zone) = mmub(nd, zone) + 2.0_dp*kearray(turb)*&
&               deltax*tempb0(nd)
            END DO
          END IF
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL4ARRAY(mmu, dp*3/4)
        temp0 = au*pi/180.0_dp + bu*yaw(turb)
        temp = COS(temp0)
        DO nd=1,nbdirs
          yawb(nd, turb) = yawb(nd, turb) - SIN(temp0)*bu*SUM(-(mu*mmub(&
&           nd, :)/temp))/temp
        END DO
        DO nd=1,nbdirs
          mmub(nd, :) = 0.0
        END DO
      END IF
      DO nd=1,nbdirs
        turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + deltaxb(nd)
        turbinexwb(nd, turb) = turbinexwb(nd, turb) - deltaxb(nd)
      END DO
      CALL POPREAL4ARRAY(wakeeffcoeffperzone, dp/4)
    END DO
    CALL POPREAL4ARRAY(wakeeffcoeff, dp/4)
  END DO
  DO nd=1,nbdirs
    keb(nd, :) = 0.0
  END DO
  DO turb=nturbines,1,-1
    DO nd=1,nbdirs
      keb(nd, turb) = keb(nd, turb) + (kecorrarray*s+1)*kearrayb(nd, &
&       turb)
      sb(nd) = ke(turb)*kecorrarray*kearrayb(nd, turb)
      kearrayb(nd, turb) = 0.0
      wakeoverlaptrelb(nd, turb, :, 1) = wakeoverlaptrelb(nd, turb, :, 1&
&       ) + sb(nd)
      wakeoverlaptrelb(nd, turb, :, 2) = wakeoverlaptrelb(nd, turb, :, 2&
&       ) + sb(nd)
    END DO
    CALL POPREAL4ARRAY(s, dp/4)
  END DO
  DO nd=1,nbdirs
    ctb(nd, :) = 0.0
    ctb(nd, :) = kecorrct*keb(nd, :)
  END DO
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    DO nd=1,nbdirs
      a_inb(nd, :) = 0.0
      a_inb(nd, :) = ab(nd, :)
    END DO
  ELSE
    CALL CTTOAXIALIND_BV(ct, ctb, nturbines, a, ab, nbdirs)
    DO nd=1,nbdirs
      a_inb(nd, :) = 0.0
    END DO
  END IF
  DO nd=1,nbdirs
    wakeoverlaptrel_vb(nd, :) = 0.0
    cosfac_vb(nd, :) = 0.0
  END DO
  DO turbi=nturbines,1,-1
    DO nd=1,nbdirs
      cosfac_vb(nd, 3*nturbines*(turbi-1)+2*nturbines+1:3*nturbines*(&
&     turbi-1)+3*nturbines) = cosfac_vb(nd, 3*nturbines*(turbi-1)+2*&
&       nturbines+1:3*nturbines*(turbi-1)+3*nturbines) + cosfacb(nd, &
&       turbi, :, 3)
      cosfacb(nd, turbi, :, 3) = 0.0
      cosfac_vb(nd, 3*nturbines*(turbi-1)+nturbines+1:3*nturbines*(turbi&
&     -1)+2*nturbines) = cosfac_vb(nd, 3*nturbines*(turbi-1)+nturbines+1&
&       :3*nturbines*(turbi-1)+2*nturbines) + cosfacb(nd, turbi, :, 2)
      cosfacb(nd, turbi, :, 2) = 0.0
      cosfac_vb(nd, 3*nturbines*(turbi-1)+1:3*nturbines*(turbi-1)+&
&     nturbines) = cosfac_vb(nd, 3*nturbines*(turbi-1)+1:3*nturbines*(&
&       turbi-1)+nturbines) + cosfacb(nd, turbi, :, 1)
      cosfacb(nd, turbi, :, 1) = 0.0
      wakeoverlaptrel_vb(nd, 3*nturbines*(turbi-1)+2*nturbines+1:3*&
&     nturbines*(turbi-1)+3*nturbines) = wakeoverlaptrel_vb(nd, 3*&
&       nturbines*(turbi-1)+2*nturbines+1:3*nturbines*(turbi-1)+3*&
&       nturbines) + wakeoverlaptrelb(nd, turbi, :, 3)
      wakeoverlaptrelb(nd, turbi, :, 3) = 0.0
      wakeoverlaptrel_vb(nd, 3*nturbines*(turbi-1)+nturbines+1:3*&
&     nturbines*(turbi-1)+2*nturbines) = wakeoverlaptrel_vb(nd, 3*&
&       nturbines*(turbi-1)+nturbines+1:3*nturbines*(turbi-1)+2*&
&       nturbines) + wakeoverlaptrelb(nd, turbi, :, 2)
      wakeoverlaptrelb(nd, turbi, :, 2) = 0.0
      wakeoverlaptrel_vb(nd, 3*nturbines*(turbi-1)+1:3*nturbines*(turbi-&
&     1)+nturbines) = wakeoverlaptrel_vb(nd, 3*nturbines*(turbi-1)+1:3*&
&       nturbines*(turbi-1)+nturbines) + wakeoverlaptrelb(nd, turbi, :, &
&       1)
      wakeoverlaptrelb(nd, turbi, :, 1) = 0.0
    END DO
  END DO
  DO nd=1,nbdirs
    yaw_degb(nd, :) = 0.0
    yaw_degb(nd, :) = pi*yawb(nd, :)/180.0_dp
  END DO
  DO nd=1,nbdirs
    velocitiesturbinesb(nd, :) = 0.0
  END DO
END SUBROUTINE FLORIS_VELOCITY_BV
