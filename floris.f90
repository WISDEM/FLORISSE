! Flow field calculations have been intentionally left out to save development time.
! The flow field can be calculated using the pure python version of floris 

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
                            wakeCenters, p_near0, wakeOverlapTRel_m)
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
    real(dp), intent(in) :: p_near0
    
    ! out    
    real(dp), dimension(nTurbines, nTurbines, 3), intent(out) :: wakeOverlapTRel_m
    
    ! local
    integer :: turb, turbI, zone
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp) :: OVdYd, OVr, OVRR, OVL, OVz
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeOverlap
        
    wakeOverlapTRel_m = 0.0_dp
    wakeOverlap = 0.0_dp
    
    do turb = 1, nTurbines
        do turbI = 1, nTurbines
            if (turbineX(turbI) > turbineX(turb) - p_near0*rotorDiameter(turb)) then
                OVdYd = wakeCenters(turbI, turb)-turbineY(turbI)    ! distance between wake center and rotor center
                OVr = rotorDiameter(turbI)/2                        ! rotor diameter
                do zone = 1, 3
                    OVRR = wakeDiameters(turbI, turb, zone)/2.0_dp        ! wake diameter
                    OVdYd = abs(OVdYd)
                    if (OVdYd /= 0.0_dp) then
                        ! calculate the distance from the wake center to the vertical line between
                        ! the two circle intersection points
                        OVL = (-OVr*OVr+OVRR*OVRR+OVdYd*OVdYd)/(2.0_dp*OVdYd)
                    else
                        OVL = 0.0_dp
                    end if

                    OVz = OVRR*OVRR-OVL*OVL

                    ! Finish calculating the distance from the intersection line to the outer edge of the wake zone
                    if (OVz > 0.0_dp) then
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
    


subroutine floris_wcent_wdiam(nTurbines, kd, initialWakeDisplacement, &
                              initialWakeAngle, ke_in, keCorrCT, Region2CT, yaw_deg, Ct, &
                              turbineXw, turbineYw, rotorDiameter, me, wakeCentersYT_vec, &
                              wakeDiametersT_vec)
    
    implicit none
    
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)
    
    ! in
    integer, intent(in) :: nTurbines
    real(dp), intent(in) :: kd, initialWakeDisplacement, initialWakeAngle, ke_in
    real(dp), intent(in) :: keCorrCT, Region2CT
    real(dp), dimension(nTurbines), intent(in) :: yaw_deg, Ct, turbineXw, turbineYw
    real(dp), dimension(nTurbines), intent(in) :: rotorDiameter
    real(dp), dimension(3), intent(in) :: me
    
    ! out
    real(dp), dimension(nTurbines*nTurbines), intent(out) :: wakeCentersYT_vec
    real(dp), dimension(3*nTurbines*nTurbines), intent(out) :: wakeDiametersT_vec
    
    ! local
    real(dp), parameter :: p_unity = 0.25
    real(dp), parameter :: p_center0 = 0.25, p_center1 = 0.25
    real(dp), parameter :: p_near0 = 1, p_near1 = p_unity
    real(dp), parameter :: p_far0 = p_unity
    real(dp), parameter :: p_mix0 = p_unity, p_mix1 = 0.25
    real(dp), dimension(nTurbines) :: ke, yaw
    Integer :: turb, turbI, zone
    real(dp) :: wakeAngleInit
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp) :: deltax, factor, displacement, x, x0, x1, y0, dy0, dx_1, factor_1, y1
    real(dp) :: b, d, dy1_yaw, dy1, wakeDiameter0, x2, y2, dy2
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeDiametersT_mat
    real(dp), dimension(nTurbines, nTurbines) :: wakeCentersYT_mat
    
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
        wakeAngleInit = 0.5_dp*sin(yaw(turb))*Ct(turb)+initialWakeAngle*pi/180.0_dp
        
        do turbI = 1, nTurbines
            deltax = turbineXw(turbI) - turbineXw(turb)
            factor = (2.0_dp*kd*deltax/rotorDiameter(turb)) + 1.0_dp
            
            if (turbineXw(turb)+p_center1*rotorDiameter(turb) < turbineXw(turbI)) then
                wakeCentersYT_mat(turbI, turb) = turbineYw(turb) - initialWakeDisplacement
                ! yaw-induced wake center displacement   
                displacement = (wakeAngleInit*(15.0_dp*(factor*factor*factor*factor) &
                               +(wakeAngleInit*wakeAngleInit))/((30.0_dp*kd* &
                               (factor*factor*factor*factor*factor))/ &
                               rotorDiameter(turb))) - (wakeAngleInit* &
                               rotorDiameter(turb)*(15.0_dp + &
                               (wakeAngleInit*wakeAngleInit))/(30.0_dp*kd)) 
                wakeCentersYT_mat(turbI, turb) = wakeCentersYT_mat(turbI, turb) + displacement
                
            else if (turbineXw(turb)+p_center1*rotorDiameter(turb) >= turbineXw(turbI) &
                    .and. turbineXw(turbI) >= &
                    turbineXw(turb)-p_center0*rotorDiameter(turb)) then
                
                ! set up spline values
                x = turbineXw(turbI)                                    ! point of interest

                x0 = turbineXw(turb) - p_center0*rotorDiameter(turb)    ! upwind point
                x1 = turbineXw(turb) + p_center1*rotorDiameter(turb)    ! downwind point

                y0 = turbineYw(turb)-initialWakeDisplacement            ! upwind point
                dy0 = 0.0_dp                                             ! upwind slope

                dx_1 = x1 - turbineXw(turb)
                factor_1 = (2.0_dp*kd*dx_1/rotorDiameter(turb)) + 1.0_dp
                !print *, 'dx_1, factor_1 = ', dx_1, factor_1
                y1 = turbineYw(turb) - initialWakeDisplacement 
     
                displacement = (wakeAngleInit*(15.0_dp*(factor_1*factor_1*factor_1*factor_1) &
                               +(wakeAngleInit*wakeAngleInit))/((30.0_dp*kd* &
                               (factor_1*factor_1*factor_1*factor_1*factor_1))/ &
                               rotorDiameter(turb))) - (wakeAngleInit* &
                               rotorDiameter(turb)*(15.0_dp + &
                               (wakeAngleInit*wakeAngleInit))/(30.0_dp*kd))
    
                y1 = y1 + displacement                                    ! downwind point
                
                b = 2.0_dp*kd/rotorDiameter(turb)
                d = b*dx_1 + 1.0_dp
                dy1_yaw = -(5.0_dp/(d*d) + (wakeAngleInit*wakeAngleInit)/(3.0_dp* &
                          (d*d*d*d*d*d*d*d))) + 4.0_dp/(d*d)
                dy1 = wakeAngleInit*dy1_yaw                                ! downwind slope
                

                call Hermite_Spline(x, x0, x1, y0, dy0, y1, dy1, wakeCentersYT_mat(turbI, turb))
                
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
        wakeDiameter0 = rotorDiameter(turb)*cos(yaw(turb))
        
        do turbI = 1, nTurbines
        
            ! turbine separation
            deltax = turbineXw(turbI) - turbineXw(turb)            
            
            ! x position of interest
            x = turbineXw(turbI)                         
            
            ! point where all zones have equal diameter
            x1 = turbineXw(turb)-p_unity*rotorDiameter(turb)  
            
            ! diameter at point of unity
            y1 = wakeDiameter0-2.0_dp*ke(turb)*me(2)*p_unity*rotorDiameter(turb) 
            
            ! derivative of diameter at upwind point of spline w.r.t downwind position
            dy1 = 2*ke(turb)*me(2)
            
            zone = 1
            if (turbineXw(turb)+p_near1*rotorDiameter(turb) < turbineXw(turbI)) then
                wakeDiametersT_mat(turbI, turb, zone) = wakeDiameter0+2.0_dp*ke(turb)*me(zone)*deltax            
                    
            else if (turbineXw(turb)+p_near1*rotorDiameter(turb) >= turbineXw(turbI) &
                 .and. turbineXw(turbI) > turbineXw(turb)-p_unity*rotorDiameter(turb)) then
                
                ! calculate spline values
                x2 = turbineXw(turb)+p_near1*rotorDiameter(turb)  ! downwind point              

                ! diameter at downwind point
                y2 = wakeDiameter0+2*ke(turb)*me(zone)*(x2-turbineXw(turb))
                
                ! slope at downwind point
                dy2 = 2.0_dp*ke(turb)*me(zone)
                
                ! solve for the wake zone diameter and its derivative w.r.t. the downwind
                ! location at the point of interest
                call Hermite_Spline(x, x1, x2, y1, dy1, y2, dy2, wakeDiametersT_mat(turbI, turb, zone))
                
            else if (turbineXw(turb)-p_near0*rotorDiameter(turb) <= turbineXw(turbI) &
                    .and. turbineXw(turbI) <= turbineXw(turb)) then
                    
                x0 = turbineXw(turb)-p_near0*rotorDiameter(turb)  ! downwind end point of spline
                
                ! diameter at upwind point of spline
                y0 = 0
                
                ! derivative of diameter at upwind point of spline
                dy0 = 0
                                
                ! solve for the wake zone diameter and its derivative w.r.t. the downwind
                ! location at the point of interest
                call Hermite_Spline(x, x0, x1, y0, dy0, y1, dy1, &
                                    wakeDiametersT_mat(turbI, turb, zone))
            end if
            
            zone = 2            
            if (turbineXw(turb)-p_far0*rotorDiameter(turb) < turbineXw(turbI)) then
                wakeDiametersT_mat(turbI, turb, zone) = wakeDiameter0 + 2.0_dp*ke(turb)*me(zone)*deltax
            else    
                wakeDiametersT_mat(turbI, turb, zone) = wakeDiametersT_mat(turbI, turb, 1)                    
            end if
            
            zone = 3
            if (turbineXw(turb)+p_mix1*rotorDiameter(turb) < turbineXw(turbI)) then
                wakeDiametersT_mat(turbI, turb, zone) = wakeDiameter0+2.0_dp*ke(turb)*me(zone)*deltax
            else if (turbineXw(turb)+p_mix1*rotorDiameter(turb) >= turbineXw(turbI) &
                     .and. turbineXw(turbI) > turbineXw(turb)-p_mix0*rotorDiameter(turb)) then
                
                x2 = turbineXw(turb)+p_mix1*rotorDiameter(turb)  ! downwind end point of spline

                ! diameter at upwind point of spline
                y2 = wakeDiameter0+2.0_dp*ke(turb)*me(zone)*p_mix1*rotorDiameter(turb)
                ! derivative of diameter at upwind point of spline w.r.t downwind position
                dy2 = 2.0_dp*ke(turb)*me(zone)

                ! solve for the wake zone diameter and its derivative w.r.t. the downwind
                ! location at the point of interest
                call Hermite_Spline(x, x1, x2, y1, dy1, y2, dy2, wakeDiametersT_mat(turbI, turb, zone))
                        
            else
                wakeDiametersT_mat(turbI, turb, zone) = wakeDiametersT_mat(turbI, turb, 1)
            
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
                          wakeDiametersT_vec, wakeCentersYT_vec, wakeOverlapTRel_vec)
    implicit none
        
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)
    
    ! in
    integer, intent(in) :: nTurbines
    real(dp), dimension(nTurbines), intent(in) :: turbineXw, turbineYw, rotorDiameter
    real(dp), dimension(3*nTurbines*nTurbines), intent(in) :: wakeDiametersT_vec
    real(dp), dimension(nTurbines*nTurbines), intent(in) :: wakeCentersYT_vec
    
    ! out
    real(dp), dimension(3*nTurbines*nTurbines), intent(out) :: wakeOverlapTRel_vec
    
    ! local
    real(dp), parameter :: p_near0 = 1 
    integer :: turbI
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeOverlapTRel_mat, wakeDiametersT_mat
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
    
    call calcOverlapAreas(nTurbines, turbineXw, turbineYw, rotorDiameter, wakeDiametersT_mat, &
                            wakeCentersYT_mat, p_near0, wakeOverlapTRel_mat)
                            
    wakeOverlapTRel_vec = 0.0_dp
    
    ! pack relative wake overlap into a vector
    do turbI = 1, nTurbines
            wakeOverlapTRel_vec(3*nTurbines*(turbI-1)+1:3*nTurbines*(turbI-1)+nTurbines) &
                                 = wakeOverlapTRel_mat(turbI, :, 1)
            wakeOverlapTRel_vec(3*nTurbines*(turbI-1)+nTurbines+1:3*nTurbines*(turbI-1) &
                                   +2*nTurbines) = wakeOverlapTRel_mat(turbI, :, 2)
            wakeOverlapTRel_vec(3*nTurbines*(turbI-1)+2*nTurbines+1:3*nTurbines*(turbI-1) &
                                   +3*nTurbines) = wakeOverlapTRel_mat(turbI, :, 3)
    end do
    
end subroutine floris_overlap



subroutine floris_power(nTurbines, wakeOverlapTRel_v, Ct, a_in, &
                 axialIndProvided, keCorrCT, Region2CT, ke_in, Vinf, keCorrArray, &
                 turbineXw, p_near0, rotorDiameter, MU, rho, Cp, generator_efficiency, &
                 velocitiesTurbines, wt_power, power)

    implicit none
    
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)


    ! in
    integer, intent(in) :: nTurbines                                            ! Number of turbines
    real(dp), dimension(3*nTurbines*nTurbines), intent(in) :: wakeOverlapTRel_v ! relative overlap (vector of length 3*nTurbines**2)
    real(dp), dimension(nTurbines), intent(in) :: Ct, a_in             ! thrust coeff, axial induction, turbine yaw
    logical, intent(in) :: axialIndProvided                        ! option logicals
    real(dp), intent(in) :: keCorrCT, Region2CT, ke_in, Vinf, keCorrArray
    real(dp), dimension(nTurbines), intent(in) :: turbineXw
    real(dp), intent(in) :: p_near0
    real(dp), dimension(nTurbines), intent(in) :: rotorDiameter
    real(dp), dimension(3), intent(in) :: MU
    real(dp), intent(in) :: rho
    real(dp), dimension(nTurbines), intent(in) :: Cp, generator_efficiency
    
    ! local
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeOverlapTRel    ! relative overlap (NxNx3)
    real(dp), dimension(nTurbines) :: a, ke, keArray
    real(dp) :: s, wakeEffCoeff, wakeEffCoeffPerZone, deltax
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), dimension(nTurbines) :: rotorArea
    integer :: turb, turbI, zone
    
    ! out
    real(dp), dimension(nTurbines), intent(out) :: velocitiesTurbines, wt_power
    real(dp), intent(out) :: power
    

    ! pack overlap vector into a matrix
    do turbI = 1, nTurbines
        wakeOverlapTRel(turbI, :, 1) = wakeOverlapTRel_v(3*nTurbines*(turbI-1)+1: &
                                   3*nTurbines*(turbI-1)+nTurbines)
        wakeOverlapTRel(turbI, :, 2) = wakeOverlapTRel_v(3*nTurbines*(turbI-1)+nTurbines+1:&
                                   3*nTurbines*(turbI-1)+2*nTurbines)
        wakeOverlapTRel(turbI, :, 3) = wakeOverlapTRel_v(3*nTurbines*(turbI-1)+2*nTurbines+1:&
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
            
            if (deltax > -1.0_dp*p_near0*rotorDiameter(turb) .and. turbI /= turb) then
                do zone = 1, 3
                    wakeEffCoeffPerZone = wakeEffCoeffPerZone + &
                    (((rotorDiameter(turb))/(rotorDiameter(turb)+2.0_dp*keArray(turb) &
                    *MU(zone)*deltax))**2)*wakeOverlapTRel(turbI, turb, zone)            
                end do
                wakeEffCoeff = wakeEffCoeff + (a(turb)*wakeEffCoeffPerZone)**2
            end if
        end do
        wakeEffCoeff = 1.0_dp - 2.0_dp*sqrt(wakeEffCoeff)
        
        ! multiply the inflow speed with the wake coefficients to find effective wind 
        ! speed at turbine
        velocitiesTurbines(turbI) = velocitiesTurbines(turbI)*wakeEffCoeff
    end do
    
    rotorArea = pi*rotorDiameter*rotorDiameter/4.0_dp
    
    ! find turbine powers
    wt_power = (velocitiesTurbines*velocitiesTurbines*velocitiesTurbines)* &
               (0.5_dp*rho*rotorArea*Cp)*generator_efficiency
               
    wt_power = wt_power/1000.0_dp
    
    ! calculate total wind farm power
    power = sum(wt_power)

end subroutine floris_power



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                BEGIN TAPENADE GENERATED CODE FOR OBTAINING GRADIENTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!!        Generated by TAPENADE     (INRIA, Tropics team)
!!  Tapenade 3.10 (r5498) - 20 Jan 2015 09:48
!!
!!  Differentiation of cttoaxialind in forward (tangent) mode:
!!   variations   of useful results: axial_induction
!!   with respect to varying inputs: ct
!! Flow field calculations have been intentionally left out to save development time.
!! The flow field can be calculated using the pure python version of floris 
!SUBROUTINE CTTOAXIALIND_DV(ct, ctd, nturbines, axial_induction, &
!& axial_inductiond, nbdirs)
!!  USE DIFFSIZES
!!  Hint: nbdirs should be the maximum number of differentiation directions
!  IMPLICIT NONE
!! define precision to be the standard for a double precision ! on local system
!  INTEGER, PARAMETER :: dp=KIND(0.d0)
!! in
!  INTEGER, INTENT(IN) :: nturbines
!  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: ct
!  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(IN) :: ctd
!! local
!  INTEGER :: i
!! out
!  REAL(dp), DIMENSION(nturbines), INTENT(OUT) :: axial_induction
!  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(OUT) :: &
!& axial_inductiond
!  INTRINSIC KIND
!  INTRINSIC SQRT
!  REAL(dp) :: arg1
!  REAL(dp), DIMENSION(nbdirs) :: arg1d
!  REAL(dp) :: result1
!  REAL(dp), DIMENSION(nbdirs) :: result1d
!  INTEGER :: nd
!  INTEGER :: nbdirs
!  axial_induction = 0.0_dp
!  DO nd=1,nbdirs
!    axial_inductiond(nd, :) = 0.0
!  END DO
!! execute
!  DO i=1,nturbines
!    IF (ct(i) .GT. 0.96) THEN
!      arg1 = 0.0203_dp - 0.6427_dp*(0.889_dp-ct(i))
!      DO nd=1,nbdirs
!! Glauert condition
!        arg1d(nd) = 0.6427_dp*ctd(nd, i)
!        IF (arg1 .EQ. 0.0) THEN
!          result1d(nd) = 0.0
!        ELSE
!          result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
!        END IF
!        axial_inductiond(nd, i) = result1d(nd)
!      END DO
!      result1 = SQRT(arg1)
!      axial_induction(i) = 0.143_dp + result1
!    ELSE
!      arg1 = 1.0_dp - ct(i)
!      DO nd=1,nbdirs
!        arg1d(nd) = -ctd(nd, i)
!        IF (arg1 .EQ. 0.0) THEN
!          result1d(nd) = 0.0
!        ELSE
!          result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
!        END IF
!        axial_inductiond(nd, i) = -(0.5_dp*result1d(nd))
!      END DO
!      result1 = SQRT(arg1)
!      axial_induction(i) = 0.5_dp*(1.0_dp-result1)
!    END IF
!  END DO
!END SUBROUTINE CTTOAXIALIND_DV
!
!
!
!!        Generated by TAPENADE     (INRIA, Tropics team)
!!  Tapenade 3.10 (r5498) - 20 Jan 2015 09:48
!!
!!  Differentiation of floris_power in forward (tangent) mode:
!!   variations   of useful results: wt_power velocitiesturbines
!!                power
!!   with respect to varying inputs: rotordiameter turbinexw wakeoverlaptrel_v
!!                cp ct a_in
!!   RW status of diff variables: rotordiameter:in turbinexw:in
!!                wakeoverlaptrel_v:in wt_power:out velocitiesturbines:out
!!                power:out cp:in ct:in a_in:in
!SUBROUTINE FLORIS_POWER_DV(nturbines, wakeoverlaptrel_v, &
!& wakeoverlaptrel_vd, ct, ctd, a_in, a_ind, axialindprovided, kecorrct, &
!& region2ct, ke_in, vinf, kecorrarray, turbinexw, turbinexwd, p_near0, &
!& rotordiameter, rotordiameterd, mu, rho, cp, cpd, generator_efficiency&
!& , velocitiesturbines, velocitiesturbinesd, wt_power, wt_powerd, power&
!& , powerd, nbdirs)
!!  USE DIFFSIZES
!!  Hint: nbdirs should be the maximum number of differentiation directions
!  IMPLICIT NONE
!! define precision to be the standard for a double precision ! on local system
!  INTEGER, PARAMETER :: dp=KIND(0.d0)
!! in
!! Number of turbines
!  INTEGER, INTENT(IN) :: nturbines, nbdirs
!! relative overlap (vector of length 3*nTurbines**2)
!  REAL(dp), DIMENSION(3*nturbines*nturbines), INTENT(IN) :: &
!& wakeoverlaptrel_v
!  REAL(dp), DIMENSION(nbdirs, 3*nturbines*nturbines), INTENT(IN) :: &
!& wakeoverlaptrel_vd
!! thrust coeff, axial induction, turbine yaw
!  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: ct, a_in
!  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(IN) :: ctd, a_ind
!! option logicals
!  LOGICAL, INTENT(IN) :: axialindprovided
!  REAL(dp), INTENT(IN) :: kecorrct, region2ct, ke_in, vinf, kecorrarray
!  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: turbinexw
!  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(IN) :: turbinexwd
!  REAL(dp), INTENT(IN) :: p_near0
!  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: rotordiameter
!  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(IN) :: &
!& rotordiameterd
!  REAL(dp), DIMENSION(3), INTENT(IN) :: mu
!  REAL(dp), INTENT(IN) :: rho
!  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: cp, generator_efficiency
!  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(IN) :: cpd
!! local
!! relative overlap (NxNx3)
!  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlaptrel
!  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
!& wakeoverlaptreld
!  REAL(dp), DIMENSION(nturbines) :: a, ke, kearray
!  REAL(dp), DIMENSION(nbdirs, nturbines) :: ad, ked, kearrayd
!  REAL(dp) :: s, wakeeffcoeff, wakeeffcoeffperzone, deltax
!  REAL(dp), DIMENSION(nbdirs) :: sd, wakeeffcoeffd, &
!& wakeeffcoeffperzoned, deltaxd
!  REAL(dp), PARAMETER :: pi=3.141592653589793_dp
!  REAL(dp), DIMENSION(nturbines) :: rotorarea
!  REAL(dp), DIMENSION(nbdirs, nturbines) :: rotoraread
!  INTEGER :: turb, turbi, zone
!! out
!  REAL(dp), DIMENSION(nturbines), INTENT(OUT) :: velocitiesturbines, &
!& wt_power
!  REAL(dp), DIMENSION(nbdirs, 3*nturbines*nturbines), INTENT(OUT) :: &
!& velocitiesturbinesd, wt_powerd
!  REAL(dp), INTENT(OUT) :: power
!  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: powerd
!  INTRINSIC KIND
!  INTRINSIC SUM
!  INTRINSIC SQRT
!  REAL(dp), DIMENSION(nturbines) :: arg1
!  REAL(dp), DIMENSION(nbdirs, nturbines) :: arg1d
!  REAL(dp) :: result1
!  REAL(dp), DIMENSION(nbdirs) :: result1d
!  INTEGER :: nd
!!  INTEGER :: nbdirs
!  DO nd=1,nbdirs
!    wakeoverlaptreld(nd, :, :, :) = 0.0
!  END DO
!! pack overlap vector into a matrix
!  DO turbi=1,nturbines
!    DO nd=1,nbdirs
!      wakeoverlaptreld(nd, turbi, :, 1) = wakeoverlaptrel_vd(nd, 3*&
!&       nturbines*(turbi-1)+1:3*nturbines*(turbi-1)+nturbines)
!      wakeoverlaptreld(nd, turbi, :, 2) = wakeoverlaptrel_vd(nd, 3*&
!&       nturbines*(turbi-1)+nturbines+1:3*nturbines*(turbi-1)+2*&
!&       nturbines)
!      wakeoverlaptreld(nd, turbi, :, 3) = wakeoverlaptrel_vd(nd, 3*&
!&       nturbines*(turbi-1)+2*nturbines+1:3*nturbines*(turbi-1)+3*&
!&       nturbines)
!    END DO
!    wakeoverlaptrel(turbi, :, 1) = wakeoverlaptrel_v(3*nturbines*(turbi-&
!&     1)+1:3*nturbines*(turbi-1)+nturbines)
!    wakeoverlaptrel(turbi, :, 2) = wakeoverlaptrel_v(3*nturbines*(turbi-&
!&     1)+nturbines+1:3*nturbines*(turbi-1)+2*nturbines)
!    wakeoverlaptrel(turbi, :, 3) = wakeoverlaptrel_v(3*nturbines*(turbi-&
!&     1)+2*nturbines+1:3*nturbines*(turbi-1)+3*nturbines)
!  END DO
!  IF (axialindprovided) THEN
!    DO nd=1,nbdirs
!      ad(nd, :) = a_ind(nd, :)
!    END DO
!    a = a_in
!  ELSE
!    CALL CTTOAXIALIND_DV(ct, ctd, nturbines, a, ad, nbdirs)
!  END IF
!  DO nd=1,nbdirs
!! adjust ke to Ct as adjusted to yaw
!    ked(nd, :) = kecorrct*ctd(nd, :)
!  END DO
!  ke = ke_in + kecorrct*(ct-region2ct)
!  DO nd=1,nbdirs
!    kearrayd(nd, :) = 0.0
!  END DO
!  DO turb=1,nturbines
!    arg1(:) = wakeoverlaptrel(turb, :, 1) + wakeoverlaptrel(turb, :, 2)
!    s = SUM(arg1(:))
!    DO nd=1,nbdirs
!      arg1d(nd, :) = wakeoverlaptreld(nd, turb, :, 1) + wakeoverlaptreld&
!&       (nd, turb, :, 2)
!      sd(nd) = SUM(arg1d(nd, :))
!      kearrayd(nd, turb) = ked(nd, turb)*(1+s*kecorrarray) + ke(turb)*&
!&       kecorrarray*sd(nd)
!    END DO
!    kearray(turb) = ke(turb)*(1+s*kecorrarray)
!  END DO
!! find effective wind speeds at downstream turbines, then predict the power of 
!! downstream turbine    
!  velocitiesturbines = vinf
!  DO nd=1,nbdirs
!    velocitiesturbinesd(nd, :) = 0.0
!  END DO
!  DO turbi=1,nturbines
!    wakeeffcoeff = 0.0_dp
!    DO nd=1,nbdirs
!      wakeeffcoeffd(nd) = 0.0
!    END DO
!! find overlap-area weighted effect of each wake zone
!    DO turb=1,nturbines
!      wakeeffcoeffperzone = 0.0_dp
!      DO nd=1,nbdirs
!        deltaxd(nd) = turbinexwd(nd, turbi) - turbinexwd(nd, turb)
!      END DO
!      deltax = turbinexw(turbi) - turbinexw(turb)
!      IF (deltax .GT. -(1.0_dp*p_near0*rotordiameter(turb)) .AND. turbi &
!&         .NE. turb) THEN
!        DO nd=1,nbdirs
!          wakeeffcoeffperzoned(nd) = 0.0
!        END DO
!        DO zone=1,3
!          DO nd=1,nbdirs
!            wakeeffcoeffperzoned(nd) = wakeeffcoeffperzoned(nd) + 2*&
!&             rotordiameter(turb)*(rotordiameterd(nd, turb)*(&
!&             rotordiameter(turb)+2.0_dp*kearray(turb)*mu(zone)*deltax)-&
!&             rotordiameter(turb)*(rotordiameterd(nd, turb)+2.0_dp*mu(&
!&             zone)*(kearrayd(nd, turb)*deltax+kearray(turb)*deltaxd(nd)&
!&             )))*wakeoverlaptrel(turbi, turb, zone)/(rotordiameter(turb&
!&             )+2.0_dp*kearray(turb)*mu(zone)*deltax)**3 + rotordiameter&
!&             (turb)**2*wakeoverlaptreld(nd, turbi, turb, zone)/(&
!&             rotordiameter(turb)+2.0_dp*kearray(turb)*mu(zone)*deltax)&
!&             **2
!          END DO
!          wakeeffcoeffperzone = wakeeffcoeffperzone + (rotordiameter(&
!&           turb)/(rotordiameter(turb)+2.0_dp*kearray(turb)*mu(zone)*&
!&           deltax))**2*wakeoverlaptrel(turbi, turb, zone)
!        END DO
!        DO nd=1,nbdirs
!          wakeeffcoeffd(nd) = wakeeffcoeffd(nd) + 2*a(turb)*&
!&           wakeeffcoeffperzone*(ad(nd, turb)*wakeeffcoeffperzone+a(turb&
!&           )*wakeeffcoeffperzoned(nd))
!        END DO
!        wakeeffcoeff = wakeeffcoeff + (a(turb)*wakeeffcoeffperzone)**2
!      END IF
!    END DO
!    DO nd=1,nbdirs
!      IF (wakeeffcoeff .EQ. 0.0) THEN
!        result1d(nd) = 0.0
!      ELSE
!        result1d(nd) = wakeeffcoeffd(nd)/(2.0*SQRT(wakeeffcoeff))
!      END IF
!      wakeeffcoeffd(nd) = -(2.0_dp*result1d(nd))
!    END DO
!    result1 = SQRT(wakeeffcoeff)
!    wakeeffcoeff = 1.0_dp - 2.0_dp*result1
!    DO nd=1,nbdirs
!! multiply the inflow speed with the wake coefficients to find effective wind 
!! speed at turbine
!      velocitiesturbinesd(nd, turbi) = velocitiesturbinesd(nd, turbi)*&
!&       wakeeffcoeff + velocitiesturbines(turbi)*wakeeffcoeffd(nd)
!    END DO
!    velocitiesturbines(turbi) = velocitiesturbines(turbi)*wakeeffcoeff
!  END DO
!  rotorarea = pi*rotordiameter*rotordiameter/4.0_dp
!  DO nd=1,nbdirs
!    rotoraread(nd, :) = pi*(rotordiameterd(nd, :)*rotordiameter+&
!&     rotordiameter*rotordiameterd(nd, :))/4.0_dp
!! find turbine powers
!    wt_powerd(nd, :) = 0.5_dp*rho*generator_efficiency*(((&
!&     velocitiesturbinesd(nd, :)*velocitiesturbines+velocitiesturbines*&
!&     velocitiesturbinesd(nd, :))*cp+velocitiesturbines**2*cpd(nd, :))*&
!&     velocitiesturbines*rotorarea+velocitiesturbines**2*cp*(&
!&     velocitiesturbinesd(nd, :)*rotorarea+velocitiesturbines*rotoraread&
!&     (nd, :)))
!    wt_powerd(nd, :) = wt_powerd(nd, :)/1000.0_dp
!! calculate total wind farm power
!    powerd(nd) = SUM(wt_powerd(nd, :))
!  END DO
!  wt_power = velocitiesturbines*velocitiesturbines*velocitiesturbines*(&
!&   0.5_dp*rho*rotorarea*cp)*generator_efficiency
!  wt_power = wt_power/1000.0_dp
!  power = SUM(wt_power)
!  print *, 'powerd in fortran: ', powerd
!END SUBROUTINE FLORIS_POWER_DV


 !All code generated at approximately 19 Jun 2015 15:30


!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.10 (r5498) - 20 Jan 2015 09:48 
!
!  Differentiation of cttoaxialind in reverse (adjoint) mode:
!   gradient     of useful results: axial_induction ct
!   with respect to varying inputs: ct
! Flow field calculations have been intentionally left out to save development time.
! The flow field can be calculated using the pure python version of floris 
SUBROUTINE CTTOAXIALIND_BV(ct, ctb, nturbines, axial_induction, &
& axial_inductionb, nbdirs)

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



!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.10 (r5498) - 20 Jan 2015 09:48
!
!  Differentiation of calcoverlapareas in reverse (adjoint) mode:
!   gradient     of useful results: wakeoverlaptrel_m
!   with respect to varying inputs: rotordiameter turbiney wakediameters
!                wakecenters
SUBROUTINE CALCOVERLAPAREAS_BV(nturbines, turbinex, turbiney, turbineyb&
& , rotordiameter, rotordiameterb, wakediameters, wakediametersb, &
& wakecenters, wakecentersb, p_near0, wakeoverlaptrel_m, &
& wakeoverlaptrel_mb, nbdirs)
 
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
  REAL(dp), INTENT(IN) :: p_near0
! out    
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlaptrel_m
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlaptrel_mb
! local
  INTEGER :: turb, turbi, zone
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp
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
      IF (turbinex(turbi) .GT. turbinex(turb) - p_near0*rotordiameter(&
&         turb)) THEN
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
          IF (ovdyd .NE. 0.0_dp) THEN
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
          IF (ovz .GT. 0.0_dp) THEN
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
  DO nd=1,nbdirs
    rotordiameterb(nd, :) = 0.0
  END DO
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
  DO nd=1,nbdirs
    turbineyb(nd, :) = 0.0
    wakediametersb(nd, :, :, :) = 0.0
    wakecentersb(nd, :, :) = 0.0
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



!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.10 (r5498) - 20 Jan 2015 09:48
!
!  Differentiation of hermite_spline in reverse (adjoint) mode:
!   gradient     of useful results: dy1 y1 x y x0 x1 dy0 y0
!   with respect to varying inputs: dy1 y1 x x0 x1 dy0 y0
SUBROUTINE HERMITE_SPLINE_BV(x, xb, x0, x0b, x1, x1b, y0, y0b, dy0, dy0b&
& , y1, y1b, dy1, dy1b, y, yb, nbdirs)
  
  IMPLICIT NONE
!dy_dx = c3*3*x**2 + c2*2*x + c1
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: x, x0, x1, y0, dy0, y1, dy1
  REAL(dp), DIMENSION(nbdirs) :: xb, x0b, x1b, y0b, dy0b, y1b, dy1b
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
  REAL(dp) :: temp18
  REAL(dp) :: temp17
  REAL(dp) :: temp16
  REAL(dp) :: temp15
  REAL(dp) :: temp14
  REAL(dp) :: temp13
  REAL(dp) :: temp12
  REAL(dp) :: temp11
  REAL(dp) :: temp10
  REAL(dp) :: tempb(nbdirs)
  REAL(dp) :: tempb34(nbdirs)
  REAL(dp) :: tempb33(nbdirs)
  REAL(dp) :: tempb32(nbdirs)
  REAL(dp) :: tempb31(nbdirs)
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
  temp14 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp13 = 3.0_dp*x0*x1**2 - x1**3
  temp16 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp15 = 3.0_dp*x1*x0**2 - x0**3
  temp17 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp18 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp8 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp7 = x1**2 + 2.0_dp*x0*x1
  temp10 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp9 = x0**2 + 2.0_dp*x1*x0
  temp11 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp12 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
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
    xb(nd) = xb(nd) + (c1+c2*2*x+c3*3*x**2)*yb(nd)
    c2b(nd) = x**2*yb(nd)
    c1b(nd) = x*yb(nd)
    c0b(nd) = yb(nd)
    tempb(nd) = c0b(nd)/temp14
    tempb0(nd) = y0*tempb(nd)
    tempb1(nd) = -(y0*temp13*tempb(nd)/temp14)
    tempb2(nd) = -(c0b(nd)/temp16)
    tempb3(nd) = y1*tempb2(nd)
    tempb4(nd) = -(y1*temp15*tempb2(nd)/temp16)
    tempb5(nd) = -(c0b(nd)/temp17)
    tempb6(nd) = x1**2*tempb5(nd)
    tempb7(nd) = -(x1**2*x0*dy0*tempb5(nd)/temp17)
    tempb8(nd) = -(c0b(nd)/temp18)
    tempb9(nd) = x0**2*tempb8(nd)
    tempb10(nd) = -(x0**2*x1*dy1*tempb8(nd)/temp18)
    tempb33(nd) = c1b(nd)/temp8
    tempb14(nd) = dy0*tempb33(nd)
    tempb15(nd) = -(dy0*temp7*tempb33(nd)/temp8)
    tempb34(nd) = c1b(nd)/temp10
    tempb16(nd) = dy1*tempb34(nd)
    tempb17(nd) = -(dy1*temp9*tempb34(nd)/temp10)
    tempb18(nd) = 6.0_dp*c1b(nd)/temp11
    tempb19(nd) = -(x0*x1*y1*tempb18(nd)/temp11)
    tempb11(nd) = -(6.0_dp*c1b(nd)/temp12)
    tempb20(nd) = -(x0*x1*y0*tempb11(nd)/temp12)
    tempb13(nd) = 3.0_dp*c2b(nd)/temp3
    tempb31(nd) = -(y0*(x0+x1)*tempb13(nd)/temp3)
    tempb30(nd) = -(c2b(nd)/temp4)
    dy1b(nd) = dy1b(nd) + temp9*tempb34(nd) + c3b(nd)/temp2 + (2.0_dp*x0&
&     +x1)*tempb30(nd) + x1*tempb9(nd)
    tempb29(nd) = -(dy1*(2.0_dp*x0+x1)*tempb30(nd)/temp4)
    tempb28(nd) = -(c2b(nd)/temp5)
    dy0b(nd) = dy0b(nd) + temp7*tempb33(nd) + c3b(nd)/temp1 + (x0+2.0_dp&
&     *x1)*tempb28(nd) + x0*tempb6(nd)
    tempb27(nd) = -(dy0*(x0+2.0_dp*x1)*tempb28(nd)/temp5)
    tempb26(nd) = -(3.0_dp*c2b(nd)/temp6)
    tempb25(nd) = -(y1*(x0+x1)*tempb26(nd)/temp6)
    tempb32(nd) = 2.0_dp*c3b(nd)/temp
    y1b(nd) = y1b(nd) + x0*x1*tempb18(nd) + tempb32(nd) + (x0+x1)*&
&     tempb26(nd) + temp15*tempb2(nd)
    tempb21(nd) = -(y1*tempb32(nd)/temp)
    tempb12(nd) = -(2.0_dp*c3b(nd)/temp0)
    y0b(nd) = y0b(nd) + x0*x1*tempb11(nd) + tempb12(nd) + (x0+x1)*&
&     tempb13(nd) + temp13*tempb(nd)
    tempb22(nd) = -(y0*tempb12(nd)/temp0)
    tempb23(nd) = -(dy0*c3b(nd)/temp1**2)
    tempb24(nd) = -(dy1*c3b(nd)/temp2**2)
    x0b(nd) = x0b(nd) + 2.0_dp*x1*tempb14(nd) + (2*x0-2.0_dp*x1)*tempb15&
&     (nd) + (2.0_dp*x1+2*x0)*tempb16(nd) + (2*x0-2.0_dp*x1)*tempb17(nd)&
&     + y1*x1*tempb18(nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*&
&     tempb19(nd) + y0*x1*tempb11(nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*&
&     x0**2)*tempb20(nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb21&
&     (nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb22(nd) + (2*x0-&
&     2.0_dp*x1)*tempb23(nd) + (2*x0-2.0_dp*x1)*tempb24(nd) + (3.0_dp*x1&
&     **2-x1*3.0_dp*2*x0+3*x0**2)*tempb25(nd) + y1*tempb26(nd) + (2*x0-&
&     2.0_dp*x1)*tempb27(nd) + dy0*tempb28(nd) + (2*x0-2.0_dp*x1)*&
&     tempb29(nd) + dy1*2.0_dp*tempb30(nd) + (3.0_dp*x1**2-x1*3.0_dp*2*&
&     x0+3*x0**2)*tempb31(nd) + y0*tempb13(nd) + (2*x0-2.0_dp*x1)*&
&     tempb10(nd) + x1*dy1*2*x0*tempb8(nd) + (2*x0-2.0_dp*x1)*tempb7(nd)&
&     + dy0*tempb6(nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb4(nd&
&     ) + (x1*3.0_dp*2*x0-3*x0**2)*tempb3(nd) + (3.0_dp*x1**2-x1*3.0_dp*&
&     2*x0+3*x0**2)*tempb1(nd) + 3.0_dp*x1**2*tempb0(nd)
    x1b(nd) = x1b(nd) + (2.0_dp*x0+2*x1)*tempb14(nd) + (2*x1-2.0_dp*x0)*&
&     tempb15(nd) + 2.0_dp*x0*tempb16(nd) + (2*x1-2.0_dp*x0)*tempb17(nd)&
&     + y1*x0*tempb18(nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*&
&     tempb19(nd) + y0*x0*tempb11(nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*&
&     x0**2)*tempb20(nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb21&
&     (nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb22(nd) + (2*x1-&
&     2.0_dp*x0)*tempb23(nd) + (2*x1-2.0_dp*x0)*tempb24(nd) + (x0*3.0_dp&
&     *2*x1-3*x1**2-3.0_dp*x0**2)*tempb25(nd) + y1*tempb26(nd) + (2*x1-&
&     2.0_dp*x0)*tempb27(nd) + dy0*2.0_dp*tempb28(nd) + (2*x1-2.0_dp*x0)&
&     *tempb29(nd) + dy1*tempb30(nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0&
&     **2)*tempb31(nd) + y0*tempb13(nd) + (2*x1-2.0_dp*x0)*tempb10(nd) +&
&     dy1*tempb9(nd) + (2*x1-2.0_dp*x0)*tempb7(nd) + x0*dy0*2*x1*tempb5(&
&     nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb4(nd) + 3.0_dp*x0&
&     **2*tempb3(nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb1(nd) &
&     + (x0*3.0_dp*2*x1-3*x1**2)*tempb0(nd)
  END DO
END SUBROUTINE HERMITE_SPLINE_BV


!
!!        Generated by TAPENADE     (INRIA, Tropics team)
!!  Tapenade 3.10 (r5498) - 20 Jan 2015 09:48
!!
!!  Differentiation of floris_wcent_wdiam in reverse (adjoint) mode:
!!   gradient     of useful results: wakediameterst_vec wakecentersyt_vec
!!   with respect to varying inputs: yaw rotordiameter turbinexw
!!                turbineyw wakediameterst_vec ct wakecentersyt_vec
!!   RW status of diff variables: yaw:out rotordiameter:out turbinexw:out
!!                turbineyw:out wakediameterst_vec:in-out ct:out
!!                wakecentersyt_vec:in-out
!SUBROUTINE FLORIS_WCENT_WDIAM_BV(nturbines, kd, initialwakedisplacement&
!& , initialwakeangle, ke_in, kecorrct, region2ct, yaw, yawb, ct, ctb, &
!& turbinexw, turbinexwb, turbineyw, turbineywb, rotordiameter, &
!& rotordiameterb, me, wakecentersyt_vec, wakecentersyt_vecb, &
!& wakediameterst_vec, wakediameterst_vecb, nbdirs)
! 
!  IMPLICIT NONE
!! define precision to be the standard for a double precision ! on local system
!  INTEGER, PARAMETER :: dp=KIND(0.d0)
!! in
!  INTEGER, INTENT(IN) :: nturbines
!  REAL(dp), INTENT(IN) :: kd, initialwakedisplacement, initialwakeangle&
!& , ke_in
!  REAL(dp), INTENT(IN) :: kecorrct, region2ct
!  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: yaw, ct, turbinexw, &
!& turbineyw
!  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: yawb, ctb, turbinexwb, &
!& turbineywb
!  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: rotordiameter
!  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: rotordiameterb
!  REAL(dp), DIMENSION(3), INTENT(IN) :: me
!! out
!  REAL(dp), DIMENSION(nturbines*nturbines), intent(out) :: wakecentersyt_vec
!  REAL(dp), DIMENSION(nbdirs, nturbines*nturbines) :: &
!& wakecentersyt_vecb
!  REAL(dp), DIMENSION(3*nturbines*nturbines), intent(out) :: wakediameterst_vec
!  REAL(dp), DIMENSION(nbdirs, 3*nturbines*nturbines) :: &
!& wakediameterst_vecb
!! local
!  REAL(dp), PARAMETER :: p_unity=0.25
!  REAL(dp), PARAMETER :: p_center0=0.25, p_center1=0.25
!  REAL(dp), PARAMETER :: p_near0=1, p_near1=p_unity
!  REAL(dp), PARAMETER :: p_far0=p_unity
!  REAL(dp), PARAMETER :: p_mix0=p_unity, p_mix1=0.25
!  REAL(dp), DIMENSION(nturbines) :: ke
!  REAL(dp), DIMENSION(nbdirs, nturbines) :: keb
!  INTEGER :: turb, turbi, zone
!  REAL(dp) :: wakeangleinit
!  REAL(dp), DIMENSION(nbdirs) :: wakeangleinitb
!  REAL(dp), PARAMETER :: pi=3.1415926535898
!  REAL(dp) :: deltax, factor, displacement, x, x0, x1, y0, dy0, dx_1, &
!& factor_1, y1
!  REAL(dp), DIMENSION(nbdirs) :: deltaxb, factorb, displacementb, xb&
!& , x0b, x1b, y0b, dy0b, dx_1b, factor_1b, y1b
!  REAL(dp) :: b, d, dy1_yaw, dy1, wakediameter0, x2, y2, dy2
!  REAL(dp), DIMENSION(nbdirs) :: bb, db, dy1_yawb, dy1b, &
!& wakediameter0b, x2b, y2b, dy2b
!  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakediameterst_mat
!  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
!& wakediameterst_matb
!  REAL(dp), DIMENSION(nturbines, nturbines) :: wakecentersyt_mat
!  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines) :: &
!& wakecentersyt_matb
!  INTRINSIC KIND
!  INTRINSIC SIN
!  INTRINSIC COS
!! execute
!!    if (CTcorrected) then
!!        Ct = Ct_in
!!    else
!!        Ct = Ct_in*cos(yaw)*cos(yaw)
!!    end if
!!    
!! calculate y-locations of wake centers in wind ref. frame
!  INTEGER :: nd
!  REAL(dp) :: tmp
!  REAL(dp) :: tmp0
!  INTEGER :: branch
!  INTEGER :: nbdirs
!  REAL(dp) :: temp3
!  REAL(dp) :: temp2
!  REAL(dp) :: temp1
!  REAL(dp) :: temp0
!  REAL(dp) :: tempb9(nbdirs)
!  REAL(dp) :: tempb8(nbdirs)
!  REAL(dp) :: tempb7(nbdirs)
!  REAL(dp) :: tempb6(nbdirs)
!  REAL(dp) :: tempb5(nbdirs)
!  REAL(dp) :: tempb4(nbdirs)
!  REAL(dp) :: tempb3(nbdirs)
!  REAL(dp) :: tempb2(nbdirs)
!  REAL(dp) :: tempb1(nbdirs)
!  REAL(dp) :: tempb0(nbdirs)
!  REAL(dp) :: tmpb(nbdirs)
!  REAL(dp) :: tempb10(nbdirs)
!  REAL(dp) :: tmpb0(nbdirs)
!  REAL(dp) :: tempb(nbdirs)
!  REAL(dp) :: temp
!  REAL(dp) :: temp5
!  REAL(dp) :: temp4
!  DO turb=1,nturbines
!    CALL PUSHREAL4ARRAY(wakeangleinit, dp/4)
!    wakeangleinit = 0.5_dp*SIN(yaw(turb))*ct(turb) + initialwakeangle*pi&
!&     /180.0_dp
!    DO turbi=1,nturbines
!      CALL PUSHREAL4ARRAY(deltax, dp/4)
!      deltax = turbinexw(turbi) - turbinexw(turb)
!      factor = 2.0_dp*kd*deltax/rotordiameter(turb) + 1.0_dp
!      IF (turbinexw(turb) + p_center1*rotordiameter(turb) .LT. turbinexw&
!&         (turbi)) THEN
!        CALL PUSHCONTROL2B(2)
!      ELSE IF (turbinexw(turb) + p_center1*rotordiameter(turb) .GE. &
!&         turbinexw(turbi) .AND. turbinexw(turbi) .GE. turbinexw(turb) -&
!&         p_center0*rotordiameter(turb)) THEN
!! set up spline values
!! point of interest
!! upwind point
!! downwind point
!        x1 = turbinexw(turb) + p_center1*rotordiameter(turb)
!! upwind point
!! upwind slope
!        dx_1 = x1 - turbinexw(turb)
!        factor_1 = 2.0_dp*kd*dx_1/rotordiameter(turb) + 1.0_dp
!!print *, 'dx_1, factor_1 = ', dx_1, factor_1
!        CALL PUSHREAL4ARRAY(y1, dp/4)
!        y1 = turbineyw(turb) - initialwakedisplacement
!        displacement = wakeangleinit*(15.0_dp*(factor_1*factor_1*&
!&         factor_1*factor_1)+wakeangleinit*wakeangleinit)/(30.0_dp*kd*(&
!&         factor_1*factor_1*factor_1*factor_1*factor_1)/rotordiameter(&
!&         turb)) - wakeangleinit*rotordiameter(turb)*(15.0_dp+&
!&         wakeangleinit*wakeangleinit)/(30.0_dp*kd)
!! downwind point
!        y1 = y1 + displacement
!        b = 2.0_dp*kd/rotordiameter(turb)
!        d = b*dx_1 + 1.0_dp
!        CALL PUSHREAL4ARRAY(dy1_yaw, dp/4)
!        dy1_yaw = -(5.0_dp/(d*d)+wakeangleinit*wakeangleinit/(3.0_dp*(d*&
!&         d*d*d*d*d*d*d))) + 4.0_dp/(d*d)
!! downwind slope
!        CALL PUSHCONTROL2B(1)
!      ELSE
!        CALL PUSHCONTROL2B(0)
!      END IF
!    END DO
!  END DO
!!adjust k_e to C_T, adjusted to yaw
!  ke = ke_in + kecorrct*(ct-region2ct)
!  DO turb=1,nturbines
!    wakediameter0 = rotordiameter(turb)*COS(yaw(turb))
!    DO turbi=1,nturbines
!! turbine separation
!      CALL PUSHREAL4ARRAY(deltax, dp/4)
!! x position of interest
!! point where all zones have equal diameter
!! diameter at point of unity
!      CALL PUSHREAL4ARRAY(y1, dp/4)
!      y1 = wakediameter0 - 2.0_dp*ke(turb)*me(2)*p_unity*rotordiameter(&
!&       turb)
!! derivative of diameter at upwind point of spline w.r.t downwind position
!      zone = 1
!      IF (turbinexw(turb) + p_near1*rotordiameter(turb) .LT. turbinexw(&
!&         turbi)) THEN
!        CALL PUSHCONTROL2B(0)
!      ELSE IF (turbinexw(turb) + p_near1*rotordiameter(turb) .GE. &
!&         turbinexw(turbi) .AND. turbinexw(turbi) .GT. turbinexw(turb) -&
!&         p_unity*rotordiameter(turb)) THEN
!! calculate spline values
!! downwind point              
!        x2 = turbinexw(turb) + p_near1*rotordiameter(turb)
!! diameter at downwind point
!        CALL PUSHREAL4ARRAY(y2, dp/4)
!        y2 = wakediameter0 + 2*ke(turb)*me(zone)*(x2-turbinexw(turb))
!! slope at downwind point
!! solve for the wake zone diameter and its derivative w.r.t. the downwind
!! location at the point of interest
!        CALL PUSHCONTROL2B(1)
!      ELSE IF (turbinexw(turb) - p_near0*rotordiameter(turb) .LE. &
!&         turbinexw(turbi) .AND. turbinexw(turbi) .LE. turbinexw(turb)) &
!&     THEN
!        CALL PUSHCONTROL2B(2)
!      ELSE
!        CALL PUSHCONTROL2B(3)
!      END IF
!      IF (turbinexw(turb) - p_far0*rotordiameter(turb) .LT. turbinexw(&
!&         turbi)) THEN
!        CALL PUSHCONTROL1B(0)
!      ELSE
!        CALL PUSHCONTROL1B(1)
!      END IF
!      zone = 3
!      IF (turbinexw(turb) + p_mix1*rotordiameter(turb) .LT. turbinexw(&
!&         turbi)) THEN
!        CALL PUSHCONTROL2B(2)
!      ELSE IF (turbinexw(turb) + p_mix1*rotordiameter(turb) .GE. &
!&         turbinexw(turbi) .AND. turbinexw(turbi) .GT. turbinexw(turb) -&
!&         p_mix0*rotordiameter(turb)) THEN
!! downwind end point of spline
!! diameter at upwind point of spline
!        CALL PUSHREAL4ARRAY(y2, dp/4)
!        y2 = wakediameter0 + 2.0_dp*ke(turb)*me(zone)*p_mix1*&
!&         rotordiameter(turb)
!! derivative of diameter at upwind point of spline w.r.t downwind position
!! solve for the wake zone diameter and its derivative w.r.t. the downwind
!! location at the point of interest
!        CALL PUSHCONTROL2B(1)
!      ELSE
!        CALL PUSHCONTROL2B(0)
!      END IF
!    END DO
!  END DO
!  DO nd=1,nbdirs
!    wakediameterst_matb(nd, :, :, :) = 0.0
!  END DO
!  DO turbi=nturbines,1,-1
!    DO nd=1,nbdirs
!      wakediameterst_matb(nd, turbi, :, 3) = wakediameterst_matb(nd, &
!&       turbi, :, 3) + wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+2*&
!&       nturbines+1:nturbines*(turbi-1)+3*nturbines)
!      wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+2*nturbines+1:&
!&     nturbines*(turbi-1)+3*nturbines) = 0.0
!      wakediameterst_matb(nd, turbi, :, 2) = wakediameterst_matb(nd, &
!&       turbi, :, 2) + wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+&
!&       nturbines+1:3*nturbines*(turbi-1)+2*nturbines)
!      wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+nturbines+1:3*&
!&     nturbines*(turbi-1)+2*nturbines) = 0.0
!      wakediameterst_matb(nd, turbi, :, 1) = wakediameterst_matb(nd, &
!&       turbi, :, 1) + wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+1:3&
!&       *nturbines*(turbi-1)+nturbines)
!      wakediameterst_vecb(nd, 3*nturbines*(turbi-1)+1:3*nturbines*(turbi&
!&     -1)+nturbines) = 0.0
!    END DO
!  END DO
!  DO nd=1,nbdirs
!    yawb(nd, :) = 0.0
!    rotordiameterb(nd, :) = 0.0
!    turbinexwb(nd, :) = 0.0
!    keb(nd, :) = 0.0
!  END DO
!  DO turb=nturbines,1,-1
!    DO nd=1,nbdirs
!      wakediameter0b(nd) = 0.0
!    END DO
!    DO turbi=nturbines,1,-1
!      CALL POPCONTROL2B(branch)
!      IF (branch .EQ. 0) THEN
!        zone = 3
!        dy1 = 2*ke(turb)*me(2)
!        x = turbinexw(turbi)
!        x1 = turbinexw(turb) - p_unity*rotordiameter(turb)
!        DO nd=1,nbdirs
!          tmpb0(nd) = wakediameterst_matb(nd, turbi, turb, zone)
!          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
!          wakediameterst_matb(nd, turbi, turb, 1) = wakediameterst_matb(&
!&           nd, turbi, turb, 1) + tmpb0(nd)
!        END DO
!        DO nd=1,nbdirs
!          dy1b(nd) = 0.0
!          y1b(nd) = 0.0
!          xb(nd) = 0.0
!          x1b(nd) = 0.0
!        END DO
!      ELSE IF (branch .EQ. 1) THEN
!        DO nd=1,nbdirs
!          xb(nd) = 0.0
!          x1b(nd) = 0.0
!          x2b(nd) = 0.0
!          y1b(nd) = 0.0
!          dy1b(nd) = 0.0
!          y2b(nd) = 0.0
!          dy2b(nd) = 0.0
!        END DO
!        dy1 = 2*ke(turb)*me(2)
!        zone = 3
!        dy2 = 2.0_dp*ke(turb)*me(zone)
!        x = turbinexw(turbi)
!        x1 = turbinexw(turb) - p_unity*rotordiameter(turb)
!        x2 = turbinexw(turb) + p_mix1*rotordiameter(turb)
!        CALL HERMITE_SPLINE_BV(x, xb, x1, x1b, x2, x2b, y1, y1b, dy1, &
!&                        dy1b, y2, y2b, dy2, dy2b, wakediameterst_mat(&
!&                        turbi, turb, zone), wakediameterst_matb(1, &
!&                        turbi, turb, zone), nbdirs)
!        CALL POPREAL4ARRAY(y2, dp/4)
!        DO nd=1,nbdirs
!          tempb10(nd) = me(zone)*p_mix1*2.0_dp*y2b(nd)
!          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
!          keb(nd, turb) = keb(nd, turb) + rotordiameter(turb)*tempb10(nd&
!&           ) + me(zone)*2.0_dp*dy2b(nd)
!          wakediameter0b(nd) = wakediameter0b(nd) + y2b(nd)
!          rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + p_mix1*&
!&           x2b(nd) + ke(turb)*tempb10(nd)
!          turbinexwb(nd, turb) = turbinexwb(nd, turb) + x2b(nd)
!        END DO
!      ELSE
!        deltax = turbinexw(turbi) - turbinexw(turb)
!        zone = 3
!        dy1 = 2*ke(turb)*me(2)
!        x = turbinexw(turbi)
!        x1 = turbinexw(turb) - p_unity*rotordiameter(turb)
!        DO nd=1,nbdirs
!          tempb9(nd) = me(zone)*2.0_dp*wakediameterst_matb(nd, turbi, &
!&           turb, zone)
!          wakediameter0b(nd) = wakediameter0b(nd) + wakediameterst_matb(&
!&           nd, turbi, turb, zone)
!          keb(nd, turb) = keb(nd, turb) + deltax*tempb9(nd)
!          deltaxb(nd) = ke(turb)*tempb9(nd)
!          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
!        END DO
!        DO nd=1,nbdirs
!          dy1b(nd) = 0.0
!          y1b(nd) = 0.0
!          xb(nd) = 0.0
!          x1b(nd) = 0.0
!        END DO
!        GOTO 100
!      END IF
!      deltax = turbinexw(turbi) - turbinexw(turb)
!      DO nd=1,nbdirs
!        deltaxb(nd) = 0.0
!      END DO
! 100  CALL POPCONTROL1B(branch)
!      IF (branch .EQ. 0) THEN
!        zone = 2
!        DO nd=1,nbdirs
!          tempb8(nd) = me(zone)*2.0_dp*wakediameterst_matb(nd, turbi, &
!&           turb, zone)
!          wakediameter0b(nd) = wakediameter0b(nd) + wakediameterst_matb(&
!&           nd, turbi, turb, zone)
!          keb(nd, turb) = keb(nd, turb) + deltax*tempb8(nd)
!          deltaxb(nd) = deltaxb(nd) + ke(turb)*tempb8(nd)
!          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
!        END DO
!      ELSE
!        zone = 2
!        DO nd=1,nbdirs
!          tmpb(nd) = wakediameterst_matb(nd, turbi, turb, zone)
!          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
!          wakediameterst_matb(nd, turbi, turb, 1) = wakediameterst_matb(&
!&           nd, turbi, turb, 1) + tmpb(nd)
!        END DO
!      END IF
!      CALL POPCONTROL2B(branch)
!      IF (branch .LT. 2) THEN
!        IF (branch .EQ. 0) THEN
!          zone = 1
!          DO nd=1,nbdirs
!            tempb6(nd) = me(zone)*2.0_dp*wakediameterst_matb(nd, turbi, &
!&             turb, zone)
!            wakediameter0b(nd) = wakediameter0b(nd) + &
!&             wakediameterst_matb(nd, turbi, turb, zone)
!            keb(nd, turb) = keb(nd, turb) + deltax*tempb6(nd)
!            deltaxb(nd) = deltaxb(nd) + ke(turb)*tempb6(nd)
!            wakediameterst_matb(nd, turbi, turb, zone) = 0.0
!          END DO
!        ELSE
!          DO nd=1,nbdirs
!            x2b(nd) = 0.0
!            y2b(nd) = 0.0
!            dy2b(nd) = 0.0
!          END DO
!          zone = 1
!          dy2 = 2.0_dp*ke(turb)*me(zone)
!          x2 = turbinexw(turb) + p_near1*rotordiameter(turb)
!          CALL HERMITE_SPLINE_BV(x, xb, x1, x1b, x2, x2b, y1, y1b, dy1, &
!&                          dy1b, y2, y2b, dy2, dy2b, wakediameterst_mat(&
!&                          turbi, turb, zone), wakediameterst_matb(1, &
!&                          turbi, turb, zone), nbdirs)
!          CALL POPREAL4ARRAY(y2, dp/4)
!          DO nd=1,nbdirs
!            tempb7(nd) = me(zone)*2*y2b(nd)
!            wakediameterst_matb(nd, turbi, turb, zone) = 0.0
!            keb(nd, turb) = keb(nd, turb) + (x2-turbinexw(turb))*tempb7(&
!&             nd) + me(zone)*2.0_dp*dy2b(nd)
!            wakediameter0b(nd) = wakediameter0b(nd) + y2b(nd)
!            x2b(nd) = x2b(nd) + ke(turb)*tempb7(nd)
!            turbinexwb(nd, turb) = turbinexwb(nd, turb) + x2b(nd) - ke(&
!&             turb)*tempb7(nd)
!            rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
!&             p_near1*x2b(nd)
!          END DO
!        END IF
!      ELSE IF (branch .EQ. 2) THEN
!        DO nd=1,nbdirs
!          x0b(nd) = 0.0
!          y0b(nd) = 0.0
!          dy0b(nd) = 0.0
!        END DO
!        x0 = turbinexw(turb) - p_near0*rotordiameter(turb)
!        zone = 1
!        dy0 = 0
!        y0 = 0
!        CALL HERMITE_SPLINE_BV(x, xb, x0, x0b, x1, x1b, y0, y0b, dy0, &
!&                        dy0b, y1, y1b, dy1, dy1b, wakediameterst_mat(&
!&                        turbi, turb, zone), wakediameterst_matb(1, &
!&                        turbi, turb, zone), nbdirs)
!        DO nd=1,nbdirs
!          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
!          turbinexwb(nd, turb) = turbinexwb(nd, turb) + x0b(nd)
!          rotordiameterb(nd, turb) = rotordiameterb(nd, turb) - p_near0*&
!&           x0b(nd)
!        END DO
!      END IF
!      CALL POPREAL4ARRAY(y1, dp/4)
!      CALL POPREAL4ARRAY(deltax, dp/4)
!      DO nd=1,nbdirs
!        tempb5(nd) = -(me(2)*p_unity*2.0_dp*y1b(nd))
!        keb(nd, turb) = keb(nd, turb) + rotordiameter(turb)*tempb5(nd) +&
!&         me(2)*2*dy1b(nd)
!        wakediameter0b(nd) = wakediameter0b(nd) + y1b(nd)
!        rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + ke(turb)*&
!&         tempb5(nd) - p_unity*x1b(nd)
!        turbinexwb(nd, turb) = turbinexwb(nd, turb) + x1b(nd)
!        turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + deltaxb(nd) + xb&
!&         (nd)
!        turbinexwb(nd, turb) = turbinexwb(nd, turb) - deltaxb(nd)
!      END DO
!    END DO
!    DO nd=1,nbdirs
!      rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + COS(yaw(turb&
!&       ))*wakediameter0b(nd)
!      yawb(nd, turb) = yawb(nd, turb) - rotordiameter(turb)*SIN(yaw(turb&
!&       ))*wakediameter0b(nd)
!    END DO
!  END DO
!  DO nd=1,nbdirs
!    ctb(nd, :) = 0.0
!    ctb(nd, :) = kecorrct*keb(nd, :)
!  END DO
!  DO nd=1,nbdirs
!    wakecentersyt_matb(nd, :, :) = 0.0
!  END DO
!  DO turbi=nturbines,1,-1
!    DO nd=1,nbdirs
!      wakecentersyt_matb(nd, turbi, :) = wakecentersyt_matb(nd, turbi, :&
!&       ) + wakecentersyt_vecb(nd, nturbines*(turbi-1)+1:nturbines*(&
!&       turbi-1)+nturbines)
!      wakecentersyt_vecb(nd, nturbines*(turbi-1)+1:nturbines*(turbi-1)+&
!&     nturbines) = 0.0
!    END DO
!  END DO
!  DO nd=1,nbdirs
!    turbineywb(nd, :) = 0.0
!  END DO
!  DO turb=nturbines,1,-1
!    DO nd=1,nbdirs
!      wakeangleinitb(nd) = 0.0
!    END DO
!    DO turbi=nturbines,1,-1
!      CALL POPCONTROL2B(branch)
!      IF (branch .EQ. 0) THEN
!        DO nd=1,nbdirs
!          turbineywb(nd, turb) = turbineywb(nd, turb) + &
!&           wakecentersyt_matb(nd, turbi, turb)
!          wakecentersyt_matb(nd, turbi, turb) = 0.0
!        END DO
!      ELSE IF (branch .EQ. 1) THEN
!        DO nd=1,nbdirs
!          xb(nd) = 0.0
!          x0b(nd) = 0.0
!          x1b(nd) = 0.0
!          y0b(nd) = 0.0
!          dy0b(nd) = 0.0
!          y1b(nd) = 0.0
!          dy1b(nd) = 0.0
!        END DO
!        dy1 = wakeangleinit*dy1_yaw
!        x = turbinexw(turbi)
!        x0 = turbinexw(turb) - p_center0*rotordiameter(turb)
!        x1 = turbinexw(turb) + p_center1*rotordiameter(turb)
!        dy0 = 0.0_dp
!        y0 = turbineyw(turb) - initialwakedisplacement
!        CALL HERMITE_SPLINE_BV(x, xb, x0, x0b, x1, x1b, y0, y0b, dy0, &
!&                        dy0b, y1, y1b, dy1, dy1b, wakecentersyt_mat(&
!&                        turbi, turb), wakecentersyt_matb(1, turbi, turb&
!&                        ), nbdirs)
!        dx_1 = x1 - turbinexw(turb)
!        b = 2.0_dp*kd/rotordiameter(turb)
!        d = b*dx_1 + 1.0_dp
!        temp5 = 3.0_dp*d**8
!        factor_1 = 2.0_dp*kd*dx_1/rotordiameter(turb) + 1.0_dp
!        temp4 = 30.0_dp*kd*factor_1**5
!        temp3 = wakeangleinit*rotordiameter(turb)
!        temp2 = 15.0_dp*factor_1**4 + wakeangleinit**2
!        DO nd=1,nbdirs
!          displacementb(nd) = y1b(nd)
!          tempb3(nd) = displacementb(nd)/temp4
!          tempb2(nd) = -((wakeangleinit**2+15.0_dp)*displacementb(nd)/(&
!&           30.0_dp*kd))
!          factor_1b(nd) = (15.0_dp*temp3*4*factor_1**3-30.0_dp*kd*temp2*&
!&           temp3*5*factor_1**4/temp4)*tempb3(nd)
!          tempb4(nd) = kd*2.0_dp*factor_1b(nd)/rotordiameter(turb)
!          dy1_yawb(nd) = wakeangleinit*dy1b(nd)
!          wakecentersyt_matb(nd, turbi, turb) = 0.0
!          wakeangleinitb(nd) = wakeangleinitb(nd) + rotordiameter(turb)*&
!&           tempb2(nd) - wakeangleinit**2*rotordiameter(turb)*2*&
!&           displacementb(nd)/(30.0_dp*kd) - 2*wakeangleinit*dy1_yawb(nd&
!&           )/temp5 + (temp2*rotordiameter(turb)+temp3*2*wakeangleinit)*&
!&           tempb3(nd) + dy1_yaw*dy1b(nd)
!          db(nd) = (5.0_dp*2/d**3+3.0_dp*wakeangleinit**2*8*d**7/temp5**&
!&           2-4.0_dp*2/d**3)*dy1_yawb(nd)
!          bb(nd) = dx_1*db(nd)
!          dx_1b(nd) = tempb4(nd) + b*db(nd)
!          turbineywb(nd, turb) = turbineywb(nd, turb) + y0b(nd) + y1b(nd&
!&           )
!          x1b(nd) = x1b(nd) + dx_1b(nd)
!          rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + temp2*&
!&           wakeangleinit*tempb3(nd) + wakeangleinit*tempb2(nd) + &
!&           p_center1*x1b(nd) - p_center0*x0b(nd) - dx_1*tempb4(nd)/&
!&           rotordiameter(turb) - kd*2.0_dp*bb(nd)/rotordiameter(turb)**&
!&           2
!          turbinexwb(nd, turb) = turbinexwb(nd, turb) + x1b(nd) + x0b(nd&
!&           ) - dx_1b(nd)
!          turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + xb(nd)
!        END DO
!        CALL POPREAL4ARRAY(dy1_yaw, dp/4)
!        CALL POPREAL4ARRAY(y1, dp/4)
!      ELSE
!        deltax = turbinexw(turbi) - turbinexw(turb)
!        factor = 2.0_dp*kd*deltax/rotordiameter(turb) + 1.0_dp
!        temp1 = 30.0_dp*kd*factor**5
!        temp0 = wakeangleinit*rotordiameter(turb)
!        temp = 15.0_dp*factor**4 + wakeangleinit**2
!        DO nd=1,nbdirs
!          displacementb(nd) = wakecentersyt_matb(nd, turbi, turb)
!          tempb0(nd) = displacementb(nd)/temp1
!          tempb1(nd) = -((wakeangleinit**2+15.0_dp)*displacementb(nd)/(&
!&           30.0_dp*kd))
!          factorb(nd) = (15.0_dp*temp0*4*factor**3-30.0_dp*kd*temp*temp0&
!&           *5*factor**4/temp1)*tempb0(nd)
!          wakeangleinitb(nd) = wakeangleinitb(nd) + rotordiameter(turb)*&
!&           tempb1(nd) - wakeangleinit**2*rotordiameter(turb)*2*&
!&           displacementb(nd)/(30.0_dp*kd) + (temp*rotordiameter(turb)+&
!&           temp0*2*wakeangleinit)*tempb0(nd)
!          rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
!&           wakeangleinit*tempb1(nd) + temp*wakeangleinit*tempb0(nd)
!          turbineywb(nd, turb) = turbineywb(nd, turb) + &
!&           wakecentersyt_matb(nd, turbi, turb)
!          wakecentersyt_matb(nd, turbi, turb) = 0.0
!        END DO
!        GOTO 110
!      END IF
!      DO nd=1,nbdirs
!        factorb(nd) = 0.0
!      END DO
! 110  DO nd=1,nbdirs
!        tempb(nd) = kd*2.0_dp*factorb(nd)/rotordiameter(turb)
!        deltaxb(nd) = tempb(nd)
!        rotordiameterb(nd, turb) = rotordiameterb(nd, turb) - deltax*&
!&         tempb(nd)/rotordiameter(turb)
!        turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + deltaxb(nd)
!        turbinexwb(nd, turb) = turbinexwb(nd, turb) - deltaxb(nd)
!      END DO
!      CALL POPREAL4ARRAY(deltax, dp/4)
!    END DO
!    CALL POPREAL4ARRAY(wakeangleinit, dp/4)
!    DO nd=1,nbdirs
!      yawb(nd, turb) = yawb(nd, turb) + ct(turb)*0.5_dp*COS(yaw(turb))*&
!&       wakeangleinitb(nd)
!      ctb(nd, turb) = ctb(nd, turb) + 0.5_dp*SIN(yaw(turb))*&
!&       wakeangleinitb(nd)
!    END DO
!  END DO
!END SUBROUTINE FLORIS_WCENT_WDIAM_BV
!

!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.10 (r5498) - 20 Jan 2015 09:48
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
& rotordiameterb, me, wakecentersyt_vec, wakecentersyt_vecb, &
& wakediameterst_vec, wakediameterst_vecb, nbdirs)
!  USE DIFFSIZES
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), INTENT(IN) :: kd, initialwakedisplacement, initialwakeangle&
& , ke_in
  REAL(dp), INTENT(IN) :: kecorrct, region2ct
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: yaw_deg, ct, turbinexw, &
& turbineyw
  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: yaw_degb, ctb, turbinexwb&
& , turbineywb
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: rotordiameter
  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: rotordiameterb
  REAL(dp), DIMENSION(3), INTENT(IN) :: me
! out
  REAL(dp), DIMENSION(nturbines*nturbines), intent(out) :: wakecentersyt_vec
  REAL(dp), DIMENSION(nbdirs, nturbines*nturbines) :: &
& wakecentersyt_vecb
  REAL(dp), DIMENSION(3*nturbines*nturbines), intent(out) :: wakediameterst_vec
  REAL(dp), DIMENSION(nbdirs, 3*nturbines*nturbines) :: &
& wakediameterst_vecb
! local
  REAL(dp), PARAMETER :: p_unity=0.25
  REAL(dp), PARAMETER :: p_center0=0.25, p_center1=0.25
  REAL(dp), PARAMETER :: p_near0=1, p_near1=p_unity
  REAL(dp), PARAMETER :: p_far0=p_unity
  REAL(dp), PARAMETER :: p_mix0=p_unity, p_mix1=0.25
  REAL(dp), DIMENSION(nturbines) :: ke, yaw
  REAL(dp), DIMENSION(nbdirs, nturbines) :: keb, yawb
  INTEGER :: turb, turbi, zone
  REAL(dp) :: wakeangleinit
  REAL(dp), DIMENSION(nbdirs) :: wakeangleinitb
  REAL(dp), PARAMETER :: pi=3.1415926535898_dp
  REAL(dp) :: deltax, factor, displacement, x, x0, x1, y0, dy0, dx_1, &
& factor_1, y1
  REAL(dp), DIMENSION(nbdirs) :: deltaxb, factorb, displacementb, xb&
& , x0b, x1b, y0b, dy0b, dx_1b, factor_1b, y1b
  REAL(dp) :: b, d, dy1_yaw, dy1, wakediameter0, x2, y2, dy2
  REAL(dp), DIMENSION(nbdirs) :: bb, db, dy1_yawb, dy1b, &
& wakediameter0b, x2b, y2b, dy2b
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
  REAL(dp) :: tmp
  REAL(dp) :: tmp0
  INTEGER :: branch
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
  REAL(dp) :: tempb3(nbdirs)
  REAL(dp) :: tempb2(nbdirs)
  REAL(dp) :: tempb1(nbdirs)
  REAL(dp) :: tempb0(nbdirs)
  REAL(dp) :: tmpb(nbdirs)
  REAL(dp) :: tempb10(nbdirs)
  REAL(dp) :: tmpb0(nbdirs)
  REAL(dp) :: tempb(nbdirs)
  REAL(dp) :: temp
  REAL(dp) :: temp5
  REAL(dp) :: temp4
! execute
!    if (CTcorrected) then
!        Ct = Ct_in
!    else
!        Ct = Ct_in*cos(yaw)*cos(yaw)
!    end if
! convert yaw from degrees to radians
  yaw(:) = yaw_deg(:)*pi/180
! calculate y-locations of wake centers in wind ref. frame
  DO turb=1,nturbines
    CALL PUSHREAL4ARRAY(wakeangleinit, dp/4)
    wakeangleinit = 0.5_dp*SIN(yaw(turb))*ct(turb) + initialwakeangle*pi&
&     /180.0_dp
    DO turbi=1,nturbines
      CALL PUSHREAL4ARRAY(deltax, dp/4)
      deltax = turbinexw(turbi) - turbinexw(turb)
      factor = 2.0_dp*kd*deltax/rotordiameter(turb) + 1.0_dp
      IF (turbinexw(turb) + p_center1*rotordiameter(turb) .LT. turbinexw&
&         (turbi)) THEN
        CALL PUSHCONTROL2B(2)
      ELSE IF (turbinexw(turb) + p_center1*rotordiameter(turb) .GE. &
&         turbinexw(turbi) .AND. turbinexw(turbi) .GE. turbinexw(turb) -&
&         p_center0*rotordiameter(turb)) THEN
! set up spline values
! point of interest
! upwind point
! downwind point
        x1 = turbinexw(turb) + p_center1*rotordiameter(turb)
! upwind point
! upwind slope
        dx_1 = x1 - turbinexw(turb)
        factor_1 = 2.0_dp*kd*dx_1/rotordiameter(turb) + 1.0_dp
!print *, 'dx_1, factor_1 = ', dx_1, factor_1
        CALL PUSHREAL4ARRAY(y1, dp/4)
        y1 = turbineyw(turb) - initialwakedisplacement
        displacement = wakeangleinit*(15.0_dp*(factor_1*factor_1*&
&         factor_1*factor_1)+wakeangleinit*wakeangleinit)/(30.0_dp*kd*(&
&         factor_1*factor_1*factor_1*factor_1*factor_1)/rotordiameter(&
&         turb)) - wakeangleinit*rotordiameter(turb)*(15.0_dp+&
&         wakeangleinit*wakeangleinit)/(30.0_dp*kd)
! downwind point
        y1 = y1 + displacement
        b = 2.0_dp*kd/rotordiameter(turb)
        d = b*dx_1 + 1.0_dp
        CALL PUSHREAL4ARRAY(dy1_yaw, dp/4)
        dy1_yaw = -(5.0_dp/(d*d)+wakeangleinit*wakeangleinit/(3.0_dp*(d*&
&         d*d*d*d*d*d*d))) + 4.0_dp/(d*d)
! downwind slope
        CALL PUSHCONTROL2B(1)
      ELSE
        CALL PUSHCONTROL2B(0)
      END IF
    END DO
  END DO
!adjust k_e to C_T, adjusted to yaw
  ke = ke_in + kecorrct*(ct-region2ct)
  DO turb=1,nturbines
    wakediameter0 = rotordiameter(turb)*COS(yaw(turb))
    DO turbi=1,nturbines
! turbine separation
      CALL PUSHREAL4ARRAY(deltax, dp/4)
! x position of interest
! point where all zones have equal diameter
! diameter at point of unity
      CALL PUSHREAL4ARRAY(y1, dp/4)
      y1 = wakediameter0 - 2.0_dp*ke(turb)*me(2)*p_unity*rotordiameter(&
&       turb)
! derivative of diameter at upwind point of spline w.r.t downwind position
      zone = 1
      IF (turbinexw(turb) + p_near1*rotordiameter(turb) .LT. turbinexw(&
&         turbi)) THEN
        CALL PUSHCONTROL2B(0)
      ELSE IF (turbinexw(turb) + p_near1*rotordiameter(turb) .GE. &
&         turbinexw(turbi) .AND. turbinexw(turbi) .GT. turbinexw(turb) -&
&         p_unity*rotordiameter(turb)) THEN
! calculate spline values
! downwind point              
        x2 = turbinexw(turb) + p_near1*rotordiameter(turb)
! diameter at downwind point
        CALL PUSHREAL4ARRAY(y2, dp/4)
        y2 = wakediameter0 + 2*ke(turb)*me(zone)*(x2-turbinexw(turb))
! slope at downwind point
! solve for the wake zone diameter and its derivative w.r.t. the downwind
! location at the point of interest
        CALL PUSHCONTROL2B(1)
      ELSE IF (turbinexw(turb) - p_near0*rotordiameter(turb) .LE. &
&         turbinexw(turbi) .AND. turbinexw(turbi) .LE. turbinexw(turb)) &
&     THEN
        CALL PUSHCONTROL2B(2)
      ELSE
        CALL PUSHCONTROL2B(3)
      END IF
      IF (turbinexw(turb) - p_far0*rotordiameter(turb) .LT. turbinexw(&
&         turbi)) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      zone = 3
      IF (turbinexw(turb) + p_mix1*rotordiameter(turb) .LT. turbinexw(&
&         turbi)) THEN
        CALL PUSHCONTROL2B(2)
      ELSE IF (turbinexw(turb) + p_mix1*rotordiameter(turb) .GE. &
&         turbinexw(turbi) .AND. turbinexw(turbi) .GT. turbinexw(turb) -&
&         p_mix0*rotordiameter(turb)) THEN
! downwind end point of spline
! diameter at upwind point of spline
        CALL PUSHREAL4ARRAY(y2, dp/4)
        y2 = wakediameter0 + 2.0_dp*ke(turb)*me(zone)*p_mix1*&
&         rotordiameter(turb)
! derivative of diameter at upwind point of spline w.r.t downwind position
! solve for the wake zone diameter and its derivative w.r.t. the downwind
! location at the point of interest
        CALL PUSHCONTROL2B(1)
      ELSE
        CALL PUSHCONTROL2B(0)
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
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        zone = 3
        dy1 = 2*ke(turb)*me(2)
        x = turbinexw(turbi)
        x1 = turbinexw(turb) - p_unity*rotordiameter(turb)
        DO nd=1,nbdirs
          tmpb0(nd) = wakediameterst_matb(nd, turbi, turb, zone)
          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
          wakediameterst_matb(nd, turbi, turb, 1) = wakediameterst_matb(&
&           nd, turbi, turb, 1) + tmpb0(nd)
        END DO
        DO nd=1,nbdirs
          dy1b(nd) = 0.0
          y1b(nd) = 0.0
          xb(nd) = 0.0
          x1b(nd) = 0.0
        END DO
      ELSE IF (branch .EQ. 1) THEN
        DO nd=1,nbdirs
          xb(nd) = 0.0
          x1b(nd) = 0.0
          x2b(nd) = 0.0
          y1b(nd) = 0.0
          dy1b(nd) = 0.0
          y2b(nd) = 0.0
          dy2b(nd) = 0.0
        END DO
        dy1 = 2*ke(turb)*me(2)
        zone = 3
        dy2 = 2.0_dp*ke(turb)*me(zone)
        x = turbinexw(turbi)
        x1 = turbinexw(turb) - p_unity*rotordiameter(turb)
        x2 = turbinexw(turb) + p_mix1*rotordiameter(turb)
        CALL HERMITE_SPLINE_BV(x, xb, x1, x1b, x2, x2b, y1, y1b, dy1, &
&                         dy1b, y2, y2b, dy2, dy2b, wakediameterst_mat(&
&                         turbi, turb, zone), wakediameterst_matb(1, &
&                         turbi, turb, zone), nbdirs)
        CALL POPREAL4ARRAY(y2, dp/4)
        DO nd=1,nbdirs
          tempb10(nd) = me(zone)*p_mix1*2.0_dp*y2b(nd)
          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
          keb(nd, turb) = keb(nd, turb) + rotordiameter(turb)*tempb10(nd&
&           ) + me(zone)*2.0_dp*dy2b(nd)
          wakediameter0b(nd) = wakediameter0b(nd) + y2b(nd)
          rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + p_mix1*&
&           x2b(nd) + ke(turb)*tempb10(nd)
          turbinexwb(nd, turb) = turbinexwb(nd, turb) + x2b(nd)
        END DO
      ELSE
        deltax = turbinexw(turbi) - turbinexw(turb)
        zone = 3
        dy1 = 2*ke(turb)*me(2)
        x = turbinexw(turbi)
        x1 = turbinexw(turb) - p_unity*rotordiameter(turb)
        DO nd=1,nbdirs
          tempb9(nd) = me(zone)*2.0_dp*wakediameterst_matb(nd, turbi, &
&           turb, zone)
          wakediameter0b(nd) = wakediameter0b(nd) + wakediameterst_matb(&
&           nd, turbi, turb, zone)
          keb(nd, turb) = keb(nd, turb) + deltax*tempb9(nd)
          deltaxb(nd) = ke(turb)*tempb9(nd)
          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
        END DO
        DO nd=1,nbdirs
          dy1b(nd) = 0.0
          y1b(nd) = 0.0
          xb(nd) = 0.0
          x1b(nd) = 0.0
        END DO
        GOTO 100
      END IF
      deltax = turbinexw(turbi) - turbinexw(turb)
      DO nd=1,nbdirs
        deltaxb(nd) = 0.0
      END DO
 100  CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        zone = 2
        DO nd=1,nbdirs
          tempb8(nd) = me(zone)*2.0_dp*wakediameterst_matb(nd, turbi, &
&           turb, zone)
          wakediameter0b(nd) = wakediameter0b(nd) + wakediameterst_matb(&
&           nd, turbi, turb, zone)
          keb(nd, turb) = keb(nd, turb) + deltax*tempb8(nd)
          deltaxb(nd) = deltaxb(nd) + ke(turb)*tempb8(nd)
          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
        END DO
      ELSE
        zone = 2
        DO nd=1,nbdirs
          tmpb(nd) = wakediameterst_matb(nd, turbi, turb, zone)
          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
          wakediameterst_matb(nd, turbi, turb, 1) = wakediameterst_matb(&
&           nd, turbi, turb, 1) + tmpb(nd)
        END DO
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          zone = 1
          DO nd=1,nbdirs
            tempb6(nd) = me(zone)*2.0_dp*wakediameterst_matb(nd, turbi, &
&             turb, zone)
            wakediameter0b(nd) = wakediameter0b(nd) + &
&             wakediameterst_matb(nd, turbi, turb, zone)
            keb(nd, turb) = keb(nd, turb) + deltax*tempb6(nd)
            deltaxb(nd) = deltaxb(nd) + ke(turb)*tempb6(nd)
            wakediameterst_matb(nd, turbi, turb, zone) = 0.0
          END DO
        ELSE
          DO nd=1,nbdirs
            x2b(nd) = 0.0
            y2b(nd) = 0.0
            dy2b(nd) = 0.0
          END DO
          zone = 1
          dy2 = 2.0_dp*ke(turb)*me(zone)
          x2 = turbinexw(turb) + p_near1*rotordiameter(turb)
          CALL HERMITE_SPLINE_BV(x, xb, x1, x1b, x2, x2b, y1, y1b, dy1&
&                           , dy1b, y2, y2b, dy2, dy2b, &
&                           wakediameterst_mat(turbi, turb, zone), &
&                           wakediameterst_matb(1, turbi, turb, zone), &
&                           nbdirs)
          CALL POPREAL4ARRAY(y2, dp/4)
          DO nd=1,nbdirs
            tempb7(nd) = me(zone)*2*y2b(nd)
            wakediameterst_matb(nd, turbi, turb, zone) = 0.0
            keb(nd, turb) = keb(nd, turb) + (x2-turbinexw(turb))*tempb7(&
&             nd) + me(zone)*2.0_dp*dy2b(nd)
            wakediameter0b(nd) = wakediameter0b(nd) + y2b(nd)
            x2b(nd) = x2b(nd) + ke(turb)*tempb7(nd)
            turbinexwb(nd, turb) = turbinexwb(nd, turb) + x2b(nd) - ke(&
&             turb)*tempb7(nd)
            rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&             p_near1*x2b(nd)
          END DO
        END IF
      ELSE IF (branch .EQ. 2) THEN
        DO nd=1,nbdirs
          x0b(nd) = 0.0
          y0b(nd) = 0.0
          dy0b(nd) = 0.0
        END DO
        x0 = turbinexw(turb) - p_near0*rotordiameter(turb)
        zone = 1
        dy0 = 0
        y0 = 0
        CALL HERMITE_SPLINE_BV(x, xb, x0, x0b, x1, x1b, y0, y0b, dy0, &
&                         dy0b, y1, y1b, dy1, dy1b, wakediameterst_mat(&
&                         turbi, turb, zone), wakediameterst_matb(1, &
&                         turbi, turb, zone), nbdirs)
        DO nd=1,nbdirs
          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
          turbinexwb(nd, turb) = turbinexwb(nd, turb) + x0b(nd)
          rotordiameterb(nd, turb) = rotordiameterb(nd, turb) - p_near0*&
&           x0b(nd)
        END DO
      END IF
      CALL POPREAL4ARRAY(y1, dp/4)
      CALL POPREAL4ARRAY(deltax, dp/4)
      DO nd=1,nbdirs
        tempb5(nd) = -(me(2)*p_unity*2.0_dp*y1b(nd))
        keb(nd, turb) = keb(nd, turb) + rotordiameter(turb)*tempb5(nd) +&
&         me(2)*2*dy1b(nd)
        wakediameter0b(nd) = wakediameter0b(nd) + y1b(nd)
        rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + ke(turb)*&
&         tempb5(nd) - p_unity*x1b(nd)
        turbinexwb(nd, turb) = turbinexwb(nd, turb) + x1b(nd)
        turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + deltaxb(nd) + xb&
&         (nd)
        turbinexwb(nd, turb) = turbinexwb(nd, turb) - deltaxb(nd)
      END DO
    END DO
    DO nd=1,nbdirs
      rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + COS(yaw(turb&
&       ))*wakediameter0b(nd)
      yawb(nd, turb) = yawb(nd, turb) - rotordiameter(turb)*SIN(yaw(turb&
&       ))*wakediameter0b(nd)
    END DO
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
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        DO nd=1,nbdirs
          turbineywb(nd, turb) = turbineywb(nd, turb) + &
&           wakecentersyt_matb(nd, turbi, turb)
          wakecentersyt_matb(nd, turbi, turb) = 0.0
        END DO
      ELSE IF (branch .EQ. 1) THEN
        DO nd=1,nbdirs
          xb(nd) = 0.0
          x0b(nd) = 0.0
          x1b(nd) = 0.0
          y0b(nd) = 0.0
          dy0b(nd) = 0.0
          y1b(nd) = 0.0
          dy1b(nd) = 0.0
        END DO
        dy1 = wakeangleinit*dy1_yaw
        x = turbinexw(turbi)
        x0 = turbinexw(turb) - p_center0*rotordiameter(turb)
        x1 = turbinexw(turb) + p_center1*rotordiameter(turb)
        dy0 = 0.0_dp
        y0 = turbineyw(turb) - initialwakedisplacement
        CALL HERMITE_SPLINE_BV(x, xb, x0, x0b, x1, x1b, y0, y0b, dy0, &
&                         dy0b, y1, y1b, dy1, dy1b, wakecentersyt_mat(&
&                         turbi, turb), wakecentersyt_matb(1, turbi, &
&                         turb), nbdirs)
        dx_1 = x1 - turbinexw(turb)
        b = 2.0_dp*kd/rotordiameter(turb)
        d = b*dx_1 + 1.0_dp
        temp5 = 3.0_dp*d**8
        factor_1 = 2.0_dp*kd*dx_1/rotordiameter(turb) + 1.0_dp
        temp4 = 30.0_dp*kd*factor_1**5
        temp3 = wakeangleinit*rotordiameter(turb)
        temp2 = 15.0_dp*factor_1**4 + wakeangleinit**2
        DO nd=1,nbdirs
          displacementb(nd) = y1b(nd)
          tempb3(nd) = displacementb(nd)/temp4
          tempb2(nd) = -((wakeangleinit**2+15.0_dp)*displacementb(nd)/(&
&           30.0_dp*kd))
          factor_1b(nd) = (15.0_dp*temp3*4*factor_1**3-30.0_dp*kd*temp2*&
&           temp3*5*factor_1**4/temp4)*tempb3(nd)
          tempb4(nd) = kd*2.0_dp*factor_1b(nd)/rotordiameter(turb)
          dy1_yawb(nd) = wakeangleinit*dy1b(nd)
          wakecentersyt_matb(nd, turbi, turb) = 0.0
          wakeangleinitb(nd) = wakeangleinitb(nd) + rotordiameter(turb)*&
&           tempb2(nd) - wakeangleinit**2*rotordiameter(turb)*2*&
&           displacementb(nd)/(30.0_dp*kd) - 2*wakeangleinit*dy1_yawb(nd&
&           )/temp5 + (temp2*rotordiameter(turb)+temp3*2*wakeangleinit)*&
&           tempb3(nd) + dy1_yaw*dy1b(nd)
          db(nd) = (5.0_dp*2/d**3+3.0_dp*wakeangleinit**2*8*d**7/temp5**&
&           2-4.0_dp*2/d**3)*dy1_yawb(nd)
          bb(nd) = dx_1*db(nd)
          dx_1b(nd) = tempb4(nd) + b*db(nd)
          turbineywb(nd, turb) = turbineywb(nd, turb) + y0b(nd) + y1b(nd&
&           )
          x1b(nd) = x1b(nd) + dx_1b(nd)
          rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + temp2*&
&           wakeangleinit*tempb3(nd) + wakeangleinit*tempb2(nd) + &
&           p_center1*x1b(nd) - p_center0*x0b(nd) - dx_1*tempb4(nd)/&
&           rotordiameter(turb) - kd*2.0_dp*bb(nd)/rotordiameter(turb)**&
&           2
          turbinexwb(nd, turb) = turbinexwb(nd, turb) + x1b(nd) + x0b(nd&
&           ) - dx_1b(nd)
          turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + xb(nd)
        END DO
        CALL POPREAL4ARRAY(dy1_yaw, dp/4)
        CALL POPREAL4ARRAY(y1, dp/4)
      ELSE
        deltax = turbinexw(turbi) - turbinexw(turb)
        factor = 2.0_dp*kd*deltax/rotordiameter(turb) + 1.0_dp
        temp1 = 30.0_dp*kd*factor**5
        temp0 = wakeangleinit*rotordiameter(turb)
        temp = 15.0_dp*factor**4 + wakeangleinit**2
        DO nd=1,nbdirs
          displacementb(nd) = wakecentersyt_matb(nd, turbi, turb)
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
        GOTO 110
      END IF
      DO nd=1,nbdirs
        factorb(nd) = 0.0
      END DO
 110  DO nd=1,nbdirs
        tempb(nd) = kd*2.0_dp*factorb(nd)/rotordiameter(turb)
        deltaxb(nd) = tempb(nd)
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
    yaw_degb(nd, :) = pi*yawb(nd, :)/180
  END DO
END SUBROUTINE FLORIS_WCENT_WDIAM_BV

!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.10 (r5498) - 20 Jan 2015 09:48
!
!  Differentiation of floris_overlap in reverse (adjoint) mode:
!   gradient     of useful results: wakeoverlaptrel_vec
!   with respect to varying inputs: rotordiameter turbineyw wakeoverlaptrel_vec
!                wakediameterst_vec wakecentersyt_vec
!   RW status of diff variables: rotordiameter:out turbineyw:out
!                wakeoverlaptrel_vec:in-zero wakediameterst_vec:out
!                wakecentersyt_vec:out
SUBROUTINE FLORIS_OVERLAP_BV(nturbines, turbinexw, turbineyw, turbineywb&
& , rotordiameter, rotordiameterb, wakediameterst_vec, &
& wakediameterst_vecb, wakecentersyt_vec, wakecentersyt_vecb, &
& wakeoverlaptrel_vec, wakeoverlaptrel_vecb, nbdirs)
  
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
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
! local
  REAL(dp), PARAMETER :: p_near0=1
  INTEGER :: turbi
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlaptrel_mat, &
& wakediameterst_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlaptrel_matb, wakediameterst_matb
  REAL(dp), DIMENSION(nturbines, nturbines) :: wakecentersyt_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines) :: &
& wakecentersyt_matb
  INTRINSIC KIND
  INTEGER :: nd
  INTEGER :: nbdirs
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
  DO nd=1,nbdirs
    wakeoverlaptrel_matb(nd, :, :, :) = 0.0
  END DO
  DO turbi=nturbines,1,-1
    DO nd=1,nbdirs
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
  CALL CALCOVERLAPAREAS_BV(nturbines, turbinexw, turbineyw, turbineywb, &
&                    rotordiameter, rotordiameterb, wakediameterst_mat, &
&                    wakediameterst_matb, wakecentersyt_mat, &
&                    wakecentersyt_matb, p_near0, wakeoverlaptrel_mat, &
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



!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.10 (r5498) - 20 Jan 2015 09:48
!
!  Differentiation of floris_power in reverse (adjoint) mode:
!   gradient     of useful results: wt_power velocitiesturbines
!                power
!   with respect to varying inputs: rotordiameter turbinexw wakeoverlaptrel_v
!                wt_power velocitiesturbines power cp ct a_in
!   RW status of diff variables: rotordiameter:out turbinexw:out
!                wakeoverlaptrel_v:out wt_power:in-zero velocitiesturbines:in-zero
!                power:in-zero cp:out ct:out a_in:out
SUBROUTINE FLORIS_POWER_BV(nturbines, wakeoverlaptrel_v, &
& wakeoverlaptrel_vb, ct, ctb, a_in, a_inb, axialindprovided, kecorrct, &
& region2ct, ke_in, vinf, kecorrarray, turbinexw, turbinexwb, p_near0, &
& rotordiameter, rotordiameterb, mu, rho, cp, cpb, generator_efficiency&
& , velocitiesturbines, velocitiesturbinesb, wt_power, wt_powerb, power&
& , powerb, nbdirs)
  
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
! Number of turbines
  INTEGER, INTENT(IN) :: nturbines
! relative overlap (vector of length 3*nTurbines**2)
  REAL(dp), DIMENSION(3*nturbines*nturbines), INTENT(IN) :: &
& wakeoverlaptrel_v
  REAL(dp), DIMENSION(nbdirs, 3*nturbines*nturbines), intent(out) :: &
& wakeoverlaptrel_vb
! thrust coeff, axial induction, turbine yaw
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: ct, a_in
  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: ctb, a_inb
! option logicals
  LOGICAL, INTENT(IN) :: axialindprovided
  REAL(dp), INTENT(IN) :: kecorrct, region2ct, ke_in, vinf, kecorrarray
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: turbinexw
  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: turbinexwb
  REAL(dp), INTENT(IN) :: p_near0
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: rotordiameter
  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: rotordiameterb
  REAL(dp), DIMENSION(3), INTENT(IN) :: mu
  REAL(dp), INTENT(IN) :: rho
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: cp, generator_efficiency
  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: cpb
! local
! relative overlap (NxNx3)
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlaptrel
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlaptrelb
  REAL(dp), DIMENSION(nturbines) :: a, ke, kearray
  REAL(dp), DIMENSION(nbdirs, nturbines) :: ab, keb, kearrayb
  REAL(dp) :: s, wakeeffcoeff, wakeeffcoeffperzone, deltax
  REAL(dp), DIMENSION(nbdirs) :: sb, wakeeffcoeffb, &
& wakeeffcoeffperzoneb, deltaxb
  REAL(dp), PARAMETER :: pi=3.141592653589793
  REAL(dp), DIMENSION(nturbines) :: rotorarea
  REAL(dp), DIMENSION(nbdirs, nturbines) :: rotorareab
  INTEGER :: turb, turbi, zone
! out
  REAL(dp), DIMENSION(nturbines), intent(out) :: velocitiesturbines, wt_power
  REAL(dp), DIMENSION(nbdirs, nturbines) :: velocitiesturbinesb, &
& wt_powerb
  REAL(dp), intent(out) :: power
  REAL(dp), DIMENSION(nbdirs) :: powerb
  INTRINSIC KIND
  INTRINSIC SUM
  INTRINSIC SQRT
  INTEGER :: nd
  INTEGER :: branch
  INTEGER :: nbdirs
  REAL(dp) :: temp0
  REAL(dp) :: tempb3(nbdirs, nturbines)
  REAL(dp) :: tempb2(nbdirs, nturbines)
  REAL(dp) :: tempb1(nbdirs)
  REAL(dp) :: tempb0(nbdirs)
  REAL(dp) :: tempb(nbdirs)
  REAL(dp) :: temp
! pack overlap vector into a matrix
  DO turbi=1,nturbines
    wakeoverlaptrel(turbi, :, 1) = wakeoverlaptrel_v(3*nturbines*(turbi-&
&     1)+1:3*nturbines*(turbi-1)+nturbines)
    wakeoverlaptrel(turbi, :, 2) = wakeoverlaptrel_v(3*nturbines*(turbi-&
&     1)+nturbines+1:3*nturbines*(turbi-1)+2*nturbines)
    wakeoverlaptrel(turbi, :, 3) = wakeoverlaptrel_v(3*nturbines*(turbi-&
&     1)+2*nturbines+1:3*nturbines*(turbi-1)+3*nturbines)
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
      IF (deltax .GT. -(1.0_dp*p_near0*rotordiameter(turb)) .AND. turbi &
&         .NE. turb) THEN
        DO zone=1,3
          wakeeffcoeffperzone = wakeeffcoeffperzone + (rotordiameter(&
&           turb)/(rotordiameter(turb)+2.0_dp*kearray(turb)*mu(zone)*&
&           deltax))**2*wakeoverlaptrel(turbi, turb, zone)
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
    CALL PUSHREAL4ARRAY(velocitiesturbines(turbi), dp/4)
    velocitiesturbines(turbi) = velocitiesturbines(turbi)*wakeeffcoeff
  END DO
  rotorarea = pi*rotordiameter*rotordiameter/4.0_dp
! find turbine powers
  wt_power = velocitiesturbines*velocitiesturbines*velocitiesturbines*(&
&   0.5_dp*rho*rotorarea*cp)*generator_efficiency
  wt_power = wt_power/1000.0_dp
! calculate total wind farm power
  DO nd=1,nbdirs
    wt_powerb(nd, :) = wt_powerb(nd, :) + powerb(nd)
    wt_powerb(nd, :) = wt_powerb(nd, :)/1000.0_dp
    cpb(nd, :) = 0.0
    rotorareab(nd, :) = 0.0
    tempb2(nd, :) = rho*0.5_dp*generator_efficiency*wt_powerb(nd, :)
    tempb3(nd, :) = velocitiesturbines**3*tempb2(nd, :)
    velocitiesturbinesb(nd, :) = velocitiesturbinesb(nd, :) + 3*&
&     velocitiesturbines**2*rotorarea*cp*tempb2(nd, :)
    rotorareab(nd, :) = cp*tempb3(nd, :)
    cpb(nd, :) = rotorarea*tempb3(nd, :)
    rotordiameterb(nd, :) = 0.0
    rotordiameterb(nd, :) = pi*2*rotordiameter*rotorareab(nd, :)/4.0_dp
  END DO
  DO nd=1,nbdirs
    turbinexwb(nd, :) = 0.0
    kearrayb(nd, :) = 0.0
    wakeoverlaptrelb(nd, :, :, :) = 0.0
    ab(nd, :) = 0.0
  END DO
  DO turbi=nturbines,1,-1
    CALL POPREAL4ARRAY(velocitiesturbines(turbi), dp/4)
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
          tempb1(nd) = 2*a(turb)*wakeeffcoeffperzone*wakeeffcoeffb(nd)
          ab(nd, turb) = ab(nd, turb) + wakeeffcoeffperzone*tempb1(nd)
          wakeeffcoeffperzoneb(nd) = a(turb)*tempb1(nd)
        END DO
        deltax = turbinexw(turbi) - turbinexw(turb)
        DO nd=1,nbdirs
          deltaxb(nd) = 0.0
        END DO
        DO zone=3,1,-1
          temp0 = 2.0_dp*mu(zone)
          temp = rotordiameter(turb) + temp0*kearray(turb)*deltax
          DO nd=1,nbdirs
            tempb(nd) = wakeeffcoeffperzoneb(nd)/temp**2
            tempb0(nd) = -(rotordiameter(turb)**2*wakeoverlaptrel(turbi&
&             , turb, zone)*2*tempb(nd)/temp)
            rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + tempb0&
&             (nd) + wakeoverlaptrel(turbi, turb, zone)*2*rotordiameter(&
&             turb)*tempb(nd)
            wakeoverlaptrelb(nd, turbi, turb, zone) = wakeoverlaptrelb(&
&             nd, turbi, turb, zone) + rotordiameter(turb)**2*tempb(nd)
            kearrayb(nd, turb) = kearrayb(nd, turb) + temp0*deltax*&
&             tempb0(nd)
            deltaxb(nd) = deltaxb(nd) + temp0*kearray(turb)*tempb0(nd)
          END DO
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
  END DO
  DO turbi=nturbines,1,-1
    DO nd=1,nbdirs
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
    wt_powerb(nd, :) = 0.0
    velocitiesturbinesb(nd, :) = 0.0
    powerb(nd) = 0.0
  END DO
END SUBROUTINE FLORIS_POWER_BV



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						Standard Tapenade Subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!C Pure tangent mode:

      SUBROUTINE PUSHREAL8_D(V,VD)
      REAL*8 V,VD
      call PUSHREAL8(V)
      call PUSHREAL8(VD)
      END

      SUBROUTINE POPREAL8_D(V,VD)
      REAL*8 V,VD
      call POPREAL8(VD)
      call POPREAL8(V)
      END

      SUBROUTINE PUSHREAL8ARRAY_D(V,VD,size)
      INTEGER size
      REAL*8 V(size),VD(size)
      call PUSHREAL8ARRAY(V,size)
      call PUSHREAL8ARRAY(VD,size)
      END

      SUBROUTINE POPREAL8ARRAY_D(V,VD,size)
      INTEGER size
      REAL*8 V(size),VD(size)
      call POPREAL8ARRAY(VD,size)
      call POPREAL8ARRAY(V,size)
      END

!C Multi-directional mode:

      SUBROUTINE PUSHREAL8_DV(V,VD,nbdirs)
      INTEGER nbdirs
      REAL*8 V,VD(nbdirs)
      call PUSHREAL8(V)
      call PUSHREAL8ARRAY(VD,nbdirs)
      END

      SUBROUTINE POPREAL8_DV(V,VD,nbdirs)
      INTEGER nbdirs
      REAL*8 V,VD(nbdirs)
      call POPREAL8ARRAY(VD,nbdirs)
      call POPREAL8(V)
      END

      SUBROUTINE PUSHREAL8ARRAY_DV(V,VD,size,nbdirs)
      INTEGER size,nbdirs
      REAL*8 V(size),VD(size,nbdirs)
      call PUSHREAL8ARRAY(V,size)
      call PUSHREAL8ARRAY(VD,size*nbdirs)
      END

      SUBROUTINE POPREAL8ARRAY_DV(V,VD,size,nbdirs)
      INTEGER size,nbdirs
      REAL*8 V(size),VD(size,nbdirs)
      call POPREAL8ARRAY(VD,size*nbdirs)
      call POPREAL8ARRAY(V,size)
      END