! Flow field calculations have been intentionally left out to save development time.
! The flow field can be calculated using the pure python version of floris 

! This implementation is fully smooth and differentiable with the exception of a 
! discontinuity at the hub of each turbine. The discontinuity only presents issues if
! turbines are place within 1E-15 * rotor diameter of one another, which is extremely 
! unlikely during optimization if the user does not explicitly place them there.

    
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
                            wakeCenters, wakeOverlapTRel_mat)
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
    real(dp), dimension(nTurbines, nTurbines, 3), intent(out) :: wakeOverlapTRel_mat
    
    ! local
    integer :: turb, turbI, zone
    real(dp), parameter :: pi = 3.141592653589793_dp, tol = 0.000001_dp
    real(dp) :: OVdYd, OVr, OVRR, OVL, OVz
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeOverlap
        
    wakeOverlapTRel_mat = 0.0_dp
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
    
    wakeOverlapTRel_mat = wakeOverlap

    do turbI = 1, nTurbines
            wakeOverlapTRel_mat(turbI, :, :) = wakeOverlapTRel_mat(turbI, :, &
                                                         :)/((pi*rotorDiameter(turbI) &
                                                       *rotorDiameter(turbI))/4.0_dp)
    end do
    
    ! do turbI = 1, nTurbines
!         do turb = 1, nTurbines
!             do zone = 1, 3
!                 print *, "wakeOverlapTRel_mat[", turbI, ", ", turb, ", ", zone, "] = ", wakeOverlapTRel_mat(turbI, turb, zone)
!             end do
!         end do
!     end do
        
   
                                    
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
    

subroutine floris(nTurbines, turbineXw, turbineYw, yawDeg, &
                          & rotorDiameter, hubHeight, Vinf, Ct, a_in, ke_in, kd, me, &
                          & initialWakeDisplacement, bd, MU, aU, bU, initialWakeAngle, &
                          & cos_spread, keCorrCT, Region2CT, keCorrArray, useWakeAngle, &
                          & adjustInitialWakeDiamToYaw, axialIndProvided, useaUbU, &
                          & shearCoefficientAlpha, shearZh, &
                          & wtVelocity)
    
    ! independent variables: yawDeg Ct turbineXw turbineYw rotorDiameter a_in    
    ! dependent variables: wtVelocity
    
    implicit none
    
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)
    
    ! in
    integer, intent(in) :: nTurbines
    real(dp), intent(in) :: kd, initialWakeDisplacement, initialWakeAngle, ke_in
    real(dp), intent(in) :: keCorrCT, Region2CT, bd, cos_spread, Vinf, keCorrArray
    real(dp), dimension(nTurbines), intent(in) :: yawDeg, Ct, a_in, turbineXw, turbineYw
    real(dp), dimension(nTurbines), intent(in) :: rotorDiameter, hubHeight
    real(dp), dimension(3), intent(in) :: me, MU
    real(dp), intent(in) :: aU, bU, shearCoefficientAlpha, shearZh
    logical, intent(in) :: useWakeAngle, adjustInitialWakeDiamToYaw, axialIndProvided, &
                           & useaUbU
                           
    ! local (General)
    real(dp), dimension(nTurbines) :: ke, yaw
    real(dp) :: deltax
    Integer :: turb, turbI, zone
    real(dp), parameter :: pi = 3.141592653589793_dp
    
    ! local (Wake centers and diameters)
    real(dp) :: spline_bound ! in rotor diameters    
    real(dp) :: wakeAngleInit, zeroloc
    real(dp) :: factor, displacement, x, x1, x2, y1, y2, dy1, dy2
    real(dp) :: wakeDiameter0
    real(dp), dimension(3) :: wakeDiameters
    real(dp) :: wakeCenterY
    
    ! local (Wake overlap)
    real(dp) :: rmax
    real(dp), dimension(3) :: wakeOverlap
    real(dp) :: tol = 0.000001_dp
    real(dp) :: OVdYd, OVr, OVRR, OVL, OVz

    ! local (Velocity)
    real(dp), dimension(nTurbines) :: a
    real(dp) :: keArray
    real(dp), dimension(3) :: mmU
    real(dp) :: s, cosFac, wakeEffCoeff, wakeEffCoeffPerZone
    
    ! model out
    real(dp), dimension(nTurbines), intent(out) :: wtVelocity
    
    intrinsic cos, atan, max
    
    yaw = yawDeg*pi/180.0_dp

    wtVelocity = Vinf

    !!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Centers and Diameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    spline_bound = 1.0_dp
       
    ! calculate y-locations of wake centers in wind ref. frame
    wakeCenterY = 0.0_dp

    !adjust k_e to C_T, adjusted to yaw
    ke = ke_in + keCorrCT*(Ct-Region2CT)

    ! initialize axial induction values
    if (axialIndProvided) then
        a = a_in
    else
        call CTtoAxialInd(Ct, nTurbines, a)
    end if
    
    do turbI = 1, nTurbines ! downstream turbines

        wakeEffCoeff = 0.0_dp

        do turb = 1, nTurbines  ! upstream turbines

            ! turbine separation
            deltax = turbineXw(turbI) - turbineXw(turb)

            wakeAngleInit = 0.5_dp*sin(yaw(turb))*Ct(turb)

            if (useWakeAngle) then
                wakeAngleInit = wakeAngleInit + initialWakeAngle*pi/180.0_dp
            end if

            if (adjustInitialWakeDiamToYaw) then
                wakeDiameter0 = rotorDiameter(turb)*cos(yaw(turb))
            else
                wakeDiameter0 = rotorDiameter(turb)
            end if

            ! wake center calculations at each turbine
            
            if (turbineXw(turb) < turbineXw(turbI)) then
                factor = (2.0_dp*kd*deltax/rotorDiameter(turb)) + 1.0_dp
                wakeCenterY = turbineYw(turb)
                
                displacement = wakeAngleInit*(wakeAngleInit* &
                                                 & wakeAngleInit + 15.0_dp*factor*factor* &
                                                 factor*factor)/((30.0_dp*kd/ & 
                                                 rotorDiameter(turb))*(factor*factor* &
                                                 & factor*factor*factor))
                                                 
                displacement = displacement - wakeAngleInit*(wakeAngleInit* &
                                              & wakeAngleInit + 15.0_dp)/(30.0_dp*kd/ &
                                              rotorDiameter(turb))

                wakeCenterY = wakeCenterY + initialWakeDisplacement + displacement
  
                if (useWakeAngle .eqv. .false.) then
                    wakeCenterY = wakeCenterY + bd*(deltax)
                end if
                
            end if

            !!!!!!!!!!!!!!!!!!!!!! calculate the wake diameter of each wake at each turbine !!!!!!!!!!!!!!!!!!!!!!!

            ! x position of interest
            x = turbineXw(turbI)

            wakeDiameters = 0.0_dp

            zone = 1

            ! define centerpoint of spline
            zeroloc = turbineXw(turb) - wakeDiameter0/(2.0_dp*ke(turb)*me(zone))

            if (zeroloc + spline_bound*rotorDiameter(turb) < turbineXw(turbI)) then ! check this
                wakeDiameters(zone) = 0.0_dp

            else if (zeroloc - spline_bound*rotorDiameter(turb) < turbineXw(turbI)) then !check this

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
                call Hermite_Spline(x, x1, x2, y1, dy1, y2, dy2, wakeDiameters(zone))

            else if (turbineXw(turb) < turbineXw(turbI)) then
                wakeDiameters(zone) = wakeDiameter0+2.0_dp*ke(turb)*me(zone)*deltax
            end if


            if (turbineXw(turb) < turbineXw(turbI)) then
                zone = 2
                wakeDiameters(zone) = wakeDiameter0 + 2.0_dp*ke(turb)*me(zone)*deltax
                zone = 3
                wakeDiameters(zone) = wakeDiameter0 + 2.0_dp*ke(turb)*me(zone)*deltax
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Overlap !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (turbineXw(turbI) > turbineXw(turb)) then
                OVdYd = wakeCenterY-turbineYw(turbI)    ! distance between wake center and rotor center
                OVr = rotorDiameter(turbI)/2                        ! rotor diameter
                do zone = 1, 3
                    OVRR = wakeDiameters(zone)/2.0_dp        ! wake diameter
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
                            wakeOverlap(zone) = OVRR*OVRR*dacos(OVL/OVRR) + OVr*OVr*dacos((OVdYd-OVL)/OVr) - OVdYd*OVz
                        else if (OVRR > OVr) then
                            wakeOverlap(zone) = pi*OVr*OVr
                        else
                            wakeOverlap(zone) = pi*OVRR*OVRR
                        end if
                    else
                        wakeOverlap(zone) = 0.0_dp
                    end if

                end do

                wakeOverlap(3) = wakeOverlap(3)-wakeOverlap(2)
                wakeOverlap(2) = wakeOverlap(2)-wakeOverlap(1)

            end if
            
            wakeOverlap = wakeOverlap/((pi*rotorDiameter(turbI) &
                                                       *rotorDiameter(turbI))/4.0_dp)
                                                       
! 			print *, wakeOverlap
			! removed array effect option
            !s = sum(wakeOverlap(1) + wakeOverlap(2))
            keArray = ke(turb) !*(1+s*keCorrArray)

            wakeEffCoeffPerZone = 0.0_dp

            if (useaUbU) then
                if (dabs(aU*pi/180.0_dp + bU*yaw(turb)) < 85.0_dp*pi/180.0_dp) then
                        mmU = MU/cos(aU*pi/180.0_dp + bU*yaw(turb))
                else
                        write(*,*)"Made it here...",MU/cos(aU*pi/180.0_dp +bU*yaw(turb))
                        mmU = MU/cos(85.0_dp*pi/180.0_dp)
                end if
            end if

            ! find overlap-area weighted effect of each wake zone
            if (deltax > 0 .and. turbI /= turb) then
                do zone = 1, 3

                    rmax = cos_spread*0.5_dp*(wakeDiameters(3) + rotorDiameter(turbI))
                    cosFac = 0.5_dp*(1.0_dp + cos(pi*dabs(wakeCenterY &
                                     & - turbineYw(turbI))/rmax))

                    if (useaUbU) then
                        wakeEffCoeffPerZone = wakeEffCoeffPerZone + &
                        (((cosFac*rotorDiameter(turb))/(rotorDiameter(turb)+2.0_dp*keArray &
                        *mmU(zone)*deltax))**2)*wakeOverlap(zone)
                    else
                        wakeEffCoeffPerZone = wakeEffCoeffPerZone + &
                        (((cosFac*rotorDiameter(turb))/(rotorDiameter(turb)+2.0_dp*keArray &
                        *MU(zone)*deltax))**2)*wakeOverlap(zone)
                    end if

                end do
                wakeEffCoeff = wakeEffCoeff + (a(turb)*wakeEffCoeffPerZone)**2
            end if

        end do

        ! find effective wind speeds at downstream turbines

        wakeEffCoeff = 1.0_dp - 2.0_dp*sqrt(wakeEffCoeff)

        ! multiply the inflow speed with the wake coefficients to find effective wind speed at turbine
        wtVelocity(turbI) = wtVelocity(turbI)*wakeEffCoeff
        
    end do

end subroutine floris

subroutine floris_visualize(nTurbines, nSamples, turbineXw, turbineYw, yawDeg, &
                          & rotorDiameter, hubHeight, Vinf, Ct, a_in, ke_in, kd, me, &
                          & initialWakeDisplacement, bd, MU, aU, bU, initialWakeAngle, &
                          & cos_spread, keCorrCT, Region2CT, keCorrArray, useWakeAngle, &
                          & adjustInitialWakeDiamToYaw, axialIndProvided, useaUbU, &
                          & wsPositionXYZw, shearCoefficientAlpha, shearZh, &
                          & wtVelocity, wsArray, &
                          & wakeCentersYT_vec, wakeDiametersT_vec, wakeOverlapTRel_vec)
    
    ! independent variables: yawDeg Ct turbineXw turbineYw rotorDiameter a_in    
    ! dependent variables: wtVelocity
    
    implicit none
    
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)
    
    ! in
    integer, intent(in) :: nTurbines
    integer :: nSamples
    real(dp), intent(in) :: kd, initialWakeDisplacement, initialWakeAngle, ke_in
    real(dp), intent(in) :: keCorrCT, Region2CT, bd, cos_spread, Vinf, keCorrArray
    real(dp), dimension(nTurbines), intent(in) :: yawDeg, Ct, a_in, turbineXw, turbineYw
    real(dp), dimension(nTurbines), intent(in) :: rotorDiameter, hubHeight
    real(dp), dimension(3), intent(in) :: me, MU
    real(dp), intent(in) :: aU, bU, shearCoefficientAlpha, shearZh
    logical, intent(in) :: useWakeAngle, adjustInitialWakeDiamToYaw, axialIndProvided, &
                           & useaUbU
    real(dp), dimension(3, nSamples), intent(in) :: wsPositionXYZw
                           
    ! local (General)
    real(dp), dimension(nTurbines) :: ke, yaw
    real(dp) :: deltax
    Integer :: turb, turbI, zone
    real(dp), parameter :: pi = 3.141592653589793_dp
    ! visualization
    Integer :: loc
    real(dp), dimension(nSamples) :: velX, velY, velZ
    
    
    
    ! local (Wake centers and diameters)
    real(dp) :: spline_bound ! in rotor diameters    
    real(dp) :: wakeAngleInit, zeroloc
    real(dp) :: factor, displacement, x, x1, x2, y1, y2, dy1, dy2
    real(dp) :: wakeDiameter0
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeDiametersT_mat
    real(dp), dimension(nTurbines, nTurbines) :: wakeCentersYT_mat
    ! visualization
    real(dp), dimension(nSamples, nTurbines, 3) :: wakeDiameters
    real(dp), dimension(nSamples, nTurbines) :: wakeCentersY, wakeCentersZ
    
    ! local (Wake overlap)
    real(dp) :: rmax
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeOverlapTRel_mat

    ! local (Velocity)
    real(dp), dimension(nTurbines) :: a, keArray
    real(dp), dimension(3) :: mmU
    real(dp) :: s, cosFac, wakeEffCoeff, wakeEffCoeffPerZone
    ! visualization
    real(dp) :: deltay, deltaz, radiusLoc, axialIndAndNearRotor, reductionFactor
    
    ! model out
    real(dp), dimension(nTurbines), intent(out) :: wtVelocity
    
    ! test out
    real(dp), dimension(nTurbines*nTurbines), intent(out) :: wakeCentersYT_vec
    real(dp), dimension(3*nTurbines*nTurbines), intent(out) :: wakeDiametersT_vec    
    real(dp), dimension(3*nTurbines*nTurbines), intent(out) :: wakeOverlapTRel_vec
    
    ! visualization out
    real(dp), dimension(nSamples), intent(out) :: wsArray
    
    intrinsic cos, atan, max
    
    if (nSamples == 1) then
        nSamples = 0
    end if
    
    yaw = yawDeg*pi/180.0_dp
    
    velX(:) = wsPositionXYZw(1, :)
    velY(:) = wsPositionXYZw(2, :)
    velZ(:) = wsPositionXYZw(3, :)
        
            
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Centers and Diameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    spline_bound = 1.0_dp
       
    ! calculate y-locations of wake centers in wind ref. frame
    wakeCentersYT_mat = 0.0_dp
    wakeCentersY = 0.0_dp
    wakeCentersZ = 0.0_dp
    
    do turb = 1, nTurbines
        wakeAngleInit = 0.5_dp*sin(yaw(turb))*Ct(turb)
        
        if (useWakeAngle) then
            wakeAngleInit = wakeAngleInit + initialWakeAngle*pi/180.0_dp
        end if
        
        ! wake center calculations at each turbine
        do turbI = 1, nTurbines            
            
            if (turbineXw(turb) < turbineXw(turbI)) then
                deltax = turbineXw(turbI) - turbineXw(turb)
                factor = (2.0_dp*kd*deltax/rotorDiameter(turb)) + 1.0_dp
                wakeCentersYT_mat(turbI, turb) = turbineYw(turb)
                
                displacement = wakeAngleInit*(wakeAngleInit* &
                                                 & wakeAngleInit + 15.0_dp*factor*factor* &
                                                 factor*factor)/((30.0_dp*kd/ & 
                                                 rotorDiameter(turb))*(factor*factor* &
                                                 & factor*factor*factor))
                                                 
                displacement = displacement - &
                                                 & wakeAngleInit*(wakeAngleInit* &
                                                 & wakeAngleInit + 15.0_dp)/(30.0_dp*kd/ &
                                                 rotorDiameter(turb))

                wakeCentersYT_mat(turbI, turb) = wakeCentersYT_mat(turbI, turb)+ &
                                                 & initialWakeDisplacement + displacement
  
                if (useWakeAngle .eqv. .false.) then
                    wakeCentersYT_mat(turbI, turb) = wakeCentersYT_mat(turbI, turb) + bd*(deltax)
                end if
                
            end if

        end do
        
        ! wake center calculations at each sample point
        do loc = 1, nSamples            
            
            if (turbineXw(turb) < velX(loc)) then
                deltax = velX(loc) - turbineXw(turb)
                factor = (2.0_dp*kd*deltax/rotorDiameter(turb)) + 1.0_dp
                wakeCentersY(loc, turb) = turbineYw(turb)
                
                displacement = wakeAngleInit*(wakeAngleInit* &
                                 & wakeAngleInit + 15.0_dp*factor*factor* &
                                 factor*factor)/((30.0_dp*kd/ & 
                                 rotorDiameter(turb))*(factor*factor* &
                                 & factor*factor*factor))
                                                 
                displacement = displacement - &
                                 & wakeAngleInit*(wakeAngleInit* &
                                 & wakeAngleInit + 15.0_dp)/(30.0_dp*kd/ &
                                 rotorDiameter(turb))

                wakeCentersY(loc, turb) = wakeCentersY(loc, turb)+ &
                                          & initialWakeDisplacement + displacement
  
                if (useWakeAngle .eqv. .false.) then
                    wakeCentersY(loc, turb) = wakeCentersY(loc, turb) + bd*(deltax)
                end if
                
                wakeCentersZ(loc, turb) = hubHeight(turb)
                
            end if

        end do
        
    end do
    
    !adjust k_e to C_T, adjusted to yaw
    ke = ke_in + keCorrCT*(Ct-Region2CT)
    
    ! calculate wake diameters
    wakeDiametersT_mat = 0.0_dp

    do turb = 1, nTurbines
        
        if (adjustInitialWakeDiamToYaw) then
            wakeDiameter0 = rotorDiameter(turb)*cos(yaw(turb))
        else
            wakeDiameter0 = rotorDiameter(turb)        
        end if
        
        ! calculate the wake diameter of each wake at each turbine
        do turbI = 1, nTurbines
        
            ! turbine separation
            deltax = turbineXw(turbI) - turbineXw(turb)            
            
            ! x position of interest
            x = turbineXw(turbI)                         
              
            zone = 1
            
            ! define centerpoint of spline
            zeroloc = turbineXw(turb) - wakeDiameter0/(2.0_dp*ke(turb)*me(zone))
            
            if (zeroloc + spline_bound*rotorDiameter(turb) < turbineXw(turbI)) then ! check this
                wakeDiametersT_mat(turbI, turb, zone) = 0.0_dp
            
            else if (zeroloc - spline_bound*rotorDiameter(turb) < turbineXw(turbI)) then !check this
                               
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
                wakeDiametersT_mat(turbI, turb, zone) = wakeDiameter0 + 2.0_dp*ke(turb)*me(zone)*deltax
            end if      
            
        end do
        
        ! calculate the wake diameter of each wake at each sample point
        do loc = 1, nSamples
        
            ! turbine separation
            deltax = velX(loc) - turbineXw(turb)            
            
            ! x position of interest
            x = velX(loc)                         
              
            zone = 1
            
            ! define centerpoint of spline
            zeroloc = turbineXw(turb) - wakeDiameter0/(2.0_dp*ke(turb)*me(zone))
            
            if (zeroloc + spline_bound*rotorDiameter(turb) < velX(loc)) then ! check this
                wakeDiameters(loc, turb, zone) = 0.0_dp
            
            else if (zeroloc - spline_bound*rotorDiameter(turb) < velX(loc)) then !check this
                               
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
                call Hermite_Spline(x, x1, x2, y1, dy1, y2, dy2, wakeDiameters(loc, turb, zone))
            
            else if (turbineXw(turb) < velX(loc)) then
                wakeDiameters(loc, turb, zone) = wakeDiameter0 + 2.0_dp*ke(turb)*me(zone)*deltax            
            end if
            
                        
            if (turbineXw(turb) < velX(loc)) then
                zone = 2
                wakeDiameters(loc, turb, zone) = wakeDiameter0 + 2.0_dp*ke(turb)*me(zone)*deltax                   
                zone = 3
                wakeDiameters(loc, turb, zone) = wakeDiameter0 + 2.0_dp*ke(turb)*me(zone)*deltax
            end if      
            
        end do
    end do
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Overlap !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! calculate relative overlap
    call calcOverlapAreas(nTurbines, turbineXw, turbineYw, rotorDiameter, &
                          & wakeDiametersT_mat, wakeCentersYT_mat, wakeOverlapTRel_mat)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Velocity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    ! initialize velocities in full flow field (optional)
    wsArray = Vinf
    
    ! apply shear profile to visualization
!     wsArray = wsArray*(velZ/shearZh)**shearCoefficientAlpha
!     print *, wsArray
    ! initialize axial induction values    
    if (axialIndProvided) then
        a = a_in
    else
        call CTtoAxialInd(Ct, nTurbines, a)
    end if
    
    ! adjust ke to Ct as adjusted to yaw
    ke = ke_in + keCorrCT*(Ct-Region2CT)

    do turb = 1, nTurbines
        s = sum(wakeOverlapTRel_mat(turb, :, 1) + wakeOverlapTRel_mat(turb, :, 2))
        keArray(turb) = ke(turb)*(1+s*keCorrArray) 
    end do
    
    ! find effective wind speeds at downstream turbines
    wtVelocity = Vinf
    do turbI = 1, nTurbines
        wakeEffCoeff = 0.0_dp
        
        ! find overlap-area weighted effect of each wake zone
        do turb = 1, nTurbines
            wakeEffCoeffPerZone = 0.0_dp
            deltax = turbineXw(turbI) - turbineXw(turb)
            
            if (useaUbU) then
                if (dabs(aU*pi/180.0_dp + bU*yaw(turb)) < 85.0_dp*pi/180.0_dp) then
                        mmU = MU/cos(aU*pi/180.0_dp + bU*yaw(turb))
                else
                        mmU = MU/cos(85.0_dp*pi/180.0_dp)
                end if
            end if
            
            if (deltax > 0 .and. turbI /= turb) then
                do zone = 1, 3
                
                    rmax = cos_spread*0.5_dp*(wakeDiametersT_mat(turbI, turb, 3) + rotorDiameter(turbI))
                    cosFac = 0.5_dp*(1.0_dp + cos(pi*dabs(wakeCentersYT_mat(turbI, turb) &
                                     & - turbineYw(turbI))/rmax))
                                
                    if (useaUbU) then
                        wakeEffCoeffPerZone = wakeEffCoeffPerZone + &
                        (((cosFac*rotorDiameter(turb))/(rotorDiameter(turb)+2.0_dp*keArray(turb) &
                        *mmU(zone)*deltax))**2)*wakeOverlapTRel_mat(turbI, turb, zone)   
                    else
                        wakeEffCoeffPerZone = wakeEffCoeffPerZone + &
                        (((cosFac*rotorDiameter(turb))/(rotorDiameter(turb)+2.0_dp*keArray(turb) &
                        *MU(zone)*deltax))**2)*wakeOverlapTRel_mat(turbI, turb, zone)   
                    end if                     
                            
                end do
                wakeEffCoeff = wakeEffCoeff + (a(turb)*wakeEffCoeffPerZone)**2
            end if
        end do
        wakeEffCoeff = 1.0_dp - 2.0_dp*sqrt(wakeEffCoeff)
        
        ! multiply the inflow speed with the wake coefficients to find effective wind 
        ! speed at turbine
        wtVelocity(turbI) = wtVelocity(turbI)*wakeEffCoeff
    end do
    
    ! calculate the velocity at the sample points
    do turb = 1, nTurbines
    
       ! if (useaUbU) then
       !     mmU = MU/cos(aU*pi/180.0_dp + bU*yaw(turb))
       ! else
       !     mmu = MU
       ! end if  
       if (useaUbU) then
           if (dabs(aU*pi/180.0_dp + bU*yaw(turb)) < 85.0_dp*pi/180.0_dp) then
                   mmU = MU/cos(aU*pi/180.0_dp + bU*yaw(turb))
           else
                   mmU = MU/cos(85.0_dp*pi/180.0_dp)
           end if
       else
           mmu = MU
       end if  

        do loc = 1, nSamples
            deltax = velX(loc) - turbineXw(turb)
            deltay = velY(loc) - wakeCentersY(loc, turb)
            deltaz = velZ(loc) - wakeCentersZ(loc, turb)
            radiusLoc = sqrt(deltay*deltay+deltaz*deltaz)            
            axialIndAndNearRotor = 2.0_dp*a(turb)
            
            rmax = cos_spread*0.5_dp*(wakeDiameters(loc, turb, 3) + rotorDiameter(turb))
            cosFac = 0.5_dp*(1.0_dp + cos(pi*radiusLoc/rmax))
            
            if (deltax > 0 .and. radiusLoc < wakeDiameters(loc, turb, 1)/2.0_dp) then   ! check if in zone 1
                reductionFactor = axialIndAndNearRotor*&
                                  & (cosFac*rotorDiameter(turb)/(rotorDiameter(turb)+2.0_dp* &
                                  & keArray(turb)*(mmU(1))*deltax))**2
            else if (deltax > 0 .and. radiusLoc < wakeDiameters(loc, turb, 2)/2.0_dp) then  ! check if in zone 2
                reductionFactor = axialIndAndNearRotor* &
                                  & (cosFac*rotorDiameter(turb)/(rotorDiameter(turb)+2.0_dp* &
                                  & keArray(turb)*(mmU(2))*deltax))**2
            else if (deltax > 0 .and. radiusLoc < wakeDiameters(loc, turb, 3)/2.0_dp) then    ! check if in zone 3
                reductionFactor = axialIndAndNearRotor* &
                                  (cosFac*rotorDiameter(turb)/(rotorDiameter(turb)+2.0_dp* &
                                  & keArray(turb)*(mmU(3))*deltax))**2
            ! use this to add upstream turbine influence to visualization
            ! else if (deltax <= 0 .and. radiusLoc < rotorDiameter(turb)/2.0_dp) then     ! check if axial induction zone in front of rotor
!                 reductionFactor = axialIndAndNearRotor*(0.5_dp+atan(2.0_dp*deltax/ &
!                                   & (rotorDiameter(turb)))/pi)
            else
                reductionFactor = 0.0_dp
            end if    
            
            wsArray(loc) = wsArray(loc)*(1.0_dp-reductionFactor)
            
        end do
        
    end do
    
!     print *, "wsArray: ", wsArray
!     print *, "hubHeight: ", hubHeight
!     print *, "velZ: ", velZ
!     print *, wsArray
        
    ! pack desired matrices into vectors for output
    do turbI = 1, nTurbines
        ! wake centers
        wakeCentersYT_vec(nTurbines*(turbI-1)+1:nTurbines*(turbI-1)+nTurbines) &
                                     = wakeCentersYT_mat(turbI, :)
                                     
        ! wake diameters
        wakeDiametersT_vec(3*nTurbines*(turbI-1)+1:3*nTurbines*(turbI-1)+nTurbines) &
                                 = wakeDiametersT_mat(turbI, :, 1)
        wakeDiametersT_vec(3*nTurbines*(turbI-1)+nTurbines+1:3*nTurbines*(turbI-1) &
                                   +2*nTurbines) = wakeDiametersT_mat(turbI, :, 2)
        wakeDiametersT_vec(3*nTurbines*(turbI-1)+2*nTurbines+1:nTurbines*(turbI-1) &
                                   +3*nTurbines) = wakeDiametersT_mat(turbI, :, 3) 
        
        ! relative wake overlap
        wakeOverlapTRel_vec(3*nTurbines*(turbI-1)+1:3*nTurbines*(turbI-1)+nTurbines) &
                             = wakeOverlapTRel_mat(turbI, :, 1)
        wakeOverlapTRel_vec(3*nTurbines*(turbI-1)+nTurbines+1:3*nTurbines*(turbI-1) &
                               +2*nTurbines) = wakeOverlapTRel_mat(turbI, :, 2)
        wakeOverlapTRel_vec(3*nTurbines*(turbI-1)+2*nTurbines+1:3*nTurbines*(turbI-1) &
                               +3*nTurbines) = wakeOverlapTRel_mat(turbI, :, 3)
        
       
    end do
    
end subroutine floris_visualize




!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5902M) - 15 Dec 2015 09:00
!
!  Differentiation of floris_unified in reverse (adjoint) mode:
!   gradient     of useful results: wtVelocity
!   with respect to varying inputs: rotordiameter turbinexw yawDeg
!                wtVelocity turbineyw ct a_in
!   RW status of diff variables: rotordiameter:out turbinexw:out
!                yawDeg:out wtVelocity:in-zero turbineyw:out
!                ct:out a_in:out

!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
!
!  Differentiation of floris in reverse (adjoint) mode:
!   gradient     of useful results: wtvelocity
!   with respect to varying inputs: rotordiameter turbinexw wtvelocity
!                turbineyw yawdeg ct a_in
!   RW status of diff variables: rotordiameter:out turbinexw:out
!                wtvelocity:in-zero turbineyw:out yawdeg:out ct:out
!                a_in:out
SUBROUTINE FLORIS_BV(nturbines, turbinexw, turbinexwb, turbineyw, &
& turbineywb, yawdeg, yawdegb, rotordiameter, rotordiameterb, &
& vinf, ct, ctb, a_in, a_inb, ke_in, kd, me, initialwakedisplacement, &
& bd, mu, au, bu, initialwakeangle, cos_spread, kecorrct, region2ct, &
& kecorrarray, usewakeangle, adjustinitialwakediamtoyaw, &
& axialindprovided, useaubu, &
& wtvelocityb, nbdirsmax)
!   USE DIFFSIZES
!  Hint: nbdirsmax should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), INTENT(IN) :: kd, initialwakedisplacement, initialwakeangle&
& , ke_in
  REAL(dp), INTENT(IN) :: kecorrct, region2ct, bd, cos_spread, vinf, &
& kecorrarray
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: yawdeg, ct, a_in, &
& turbinexw, turbineyw
  REAL(dp), DIMENSION(nbdirsmax, nturbines), intent(out) :: yawdegb, ctb, a_inb, &
& turbinexwb, turbineywb
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: rotordiameter
  REAL(dp), DIMENSION(nbdirsmax, nturbines), intent(out) :: rotordiameterb
  REAL(dp), DIMENSION(3), INTENT(IN) :: me, mu
  REAL(dp), INTENT(IN) :: au, bu
  LOGICAL, INTENT(IN) :: usewakeangle, adjustinitialwakediamtoyaw, &
& axialindprovided, useaubu
! local (General)
  REAL(dp), DIMENSION(nturbines) :: ke, yaw
  REAL(dp), DIMENSION(nbdirsmax, nturbines) :: keb, yawb
  REAL(dp) :: deltax
  REAL(dp), DIMENSION(nbdirsmax) :: deltaxb
  INTEGER :: turb, turbi, zone
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp
! local (Wake centers and diameters)
! in rotor diameters    
  REAL(dp) :: spline_bound
  REAL(dp) :: wakeangleinit, zeroloc
  REAL(dp), DIMENSION(nbdirsmax) :: wakeangleinitb, zerolocb
  REAL(dp) :: factor, displacement, x, x1, x2, y1, y2, dy1, dy2
  REAL(dp), DIMENSION(nbdirsmax) :: factorb, displacementb, xb, x1b, x2b&
& , y1b, dy1b
  REAL(dp) :: wakediameter0
  REAL(dp), DIMENSION(nbdirsmax) :: wakediameter0b
  REAL(dp), DIMENSION(3) :: wakediameters
  REAL(dp), DIMENSION(nbdirsmax, 3) :: wakediametersb
  REAL(dp) :: wakecentery
  REAL(dp), DIMENSION(nbdirsmax) :: wakecenteryb
! local (Wake overlap)
  REAL(dp) :: rmax
  REAL(dp), DIMENSION(nbdirsmax) :: rmaxb
  REAL(dp), DIMENSION(3) :: wakeoverlap
  REAL(dp), DIMENSION(nbdirsmax, 3) :: wakeoverlapb
  REAL(dp), SAVE :: tol=0.000001_dp
  REAL(dp) :: ovdyd, ovr, ovrr, ovl, ovz
  REAL(dp), DIMENSION(nbdirsmax) :: ovdydb, ovrb, ovrrb, ovlb, ovzb
! local (Velocity)
  REAL(dp), DIMENSION(nturbines) :: a
  REAL(dp), DIMENSION(nbdirsmax, nturbines) :: ab
  REAL(dp) :: kearray
  REAL(dp), DIMENSION(nbdirsmax) :: kearrayb
  REAL(dp), DIMENSION(3) :: mmu
  REAL(dp), DIMENSION(nbdirsmax, 3) :: mmub
  REAL(dp) :: s, cosfac, wakeeffcoeff, wakeeffcoeffperzone
  REAL(dp), DIMENSION(nbdirsmax) :: cosfacb, wakeeffcoeffb, &
& wakeeffcoeffperzoneb
! model out
  REAL(dp), DIMENSION(nturbines) :: wtvelocity
  REAL(dp), DIMENSION(nbdirsmax, nturbines) :: wtvelocityb
  INTRINSIC COS, ATAN, MAX
  INTRINSIC KIND
  INTRINSIC SIN
  INTRINSIC ABS
  INTRINSIC SQRT
  INTRINSIC DACOS
  INTRINSIC DABS
  DOUBLE PRECISION :: dabs0
  DOUBLE PRECISION :: dabs1
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: dabs1b
  INTEGER :: nd
  REAL(dp) :: temp
  REAL(dp) :: temp0
  REAL(dp) :: temp1
  REAL(dp), DIMENSION(nbdirsmax) :: tempb
  REAL(dp), DIMENSION(nbdirsmax) :: tempb0
  REAL(dp), DIMENSION(nbdirsmax) :: tempb1
  REAL(dp) :: temp2
  REAL(dp), DIMENSION(nbdirsmax) :: tempb2
  REAL(dp), DIMENSION(nbdirsmax) :: tempb3
  REAL(dp), DIMENSION(nbdirsmax) :: tempb4
  REAL(dp), DIMENSION(nbdirsmax) :: tempb5
  REAL(dp), DIMENSION(nbdirsmax) :: tempb6
  REAL(dp) :: temp3
  REAL(dp) :: temp4
  REAL(dp), DIMENSION(nbdirsmax) :: tempb7
  REAL(dp), DIMENSION(nbdirsmax) :: tempb8
  REAL(dp) :: temp5
  REAL(dp), DIMENSION(nbdirsmax, 3) :: tempb9
  REAL(dp) :: temp6
  REAL(dp) :: temp7
  REAL(dp), DIMENSION(nbdirsmax) :: tempb10
  REAL, DIMENSION(nbdirsmax) :: tempb11
  REAL(dp) :: temp8
  REAL(dp) :: temp9
  REAL(dp) :: temp10
  REAL(dp), DIMENSION(nbdirsmax) :: tempb12
  REAL(dp), DIMENSION(nbdirsmax) :: tempb13
  REAL(dp), DIMENSION(nbdirsmax) :: tempb14
  REAL(dp) :: temp11
  REAL(dp) :: temp12
  REAL(dp) :: temp13
  REAL(dp) :: temp14
  REAL(dp), DIMENSION(nbdirsmax) :: tempb15
  REAL(dp), DIMENSION(nbdirsmax) :: tempb16
  REAL(dp), DIMENSION(nbdirsmax) :: tempb17
  INTEGER :: branch
  INTEGER :: nbdirsmax
  DOUBLE PRECISION :: x3
  yaw = yawdeg*pi/180.0_dp
  wtvelocity = vinf
!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Centers and Diameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
  spline_bound = 1.0_dp
! calculate y-locations of wake centers in wind ref. frame
  wakecentery = 0.0_dp
!adjust k_e to C_T, adjusted to yaw
  ke = ke_in + kecorrct*(ct-region2ct)
! initialize axial induction values
  IF (axialindprovided) THEN
    a = a_in
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL CTTOAXIALIND(ct, nturbines, a)
    CALL PUSHCONTROL1B(0)
  END IF
! downstream turbines
  DO turbi=1,nturbines
    CALL PUSHREAL4ARRAY(wakeeffcoeff, dp/4)
    wakeeffcoeff = 0.0_dp
! upstream turbines
    DO turb=1,nturbines
! turbine separation
      CALL PUSHREAL4ARRAY(deltax, dp/4)
      deltax = turbinexw(turbi) - turbinexw(turb)
      CALL PUSHREAL4ARRAY(wakeangleinit, dp/4)
      wakeangleinit = 0.5_dp*SIN(yaw(turb))*ct(turb)
      IF (usewakeangle) wakeangleinit = wakeangleinit + initialwakeangle&
&         *pi/180.0_dp
      IF (adjustinitialwakediamtoyaw) THEN
        CALL PUSHREAL4ARRAY(wakediameter0, dp/4)
        wakediameter0 = rotordiameter(turb)*COS(yaw(turb))
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHREAL4ARRAY(wakediameter0, dp/4)
        wakediameter0 = rotordiameter(turb)
        CALL PUSHCONTROL1B(1)
      END IF
! wake center calculations at each turbine
      IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        factor = 2.0_dp*kd*deltax/rotordiameter(turb) + 1.0_dp
        wakecentery = turbineyw(turb)
        displacement = wakeangleinit*(wakeangleinit*wakeangleinit+&
&         15.0_dp*factor*factor*factor*factor)/(30.0_dp*kd/rotordiameter&
&         (turb)*(factor*factor*factor*factor*factor))
        displacement = displacement - wakeangleinit*(wakeangleinit*&
&         wakeangleinit+15.0_dp)/(30.0_dp*kd/rotordiameter(turb))
        wakecentery = wakecentery + initialwakedisplacement + &
&         displacement
        IF (usewakeangle .EQV. .false.) THEN
          wakecentery = wakecentery + bd*deltax
          CALL PUSHCONTROL2B(0)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(2)
      END IF
!!!!!!!!!!!!!!!!!!!!!! calculate the wake diameter of each wake at each turbine !!!!!!!!!!!!!!!!!!!!!!!
! x position of interest
      x = turbinexw(turbi)
      wakediameters = 0.0_dp
      CALL PUSHINTEGER4(zone)
      zone = 1
! define centerpoint of spline
      zeroloc = turbinexw(turb) - wakediameter0/(2.0_dp*ke(turb)*me(zone&
&       ))
      IF (zeroloc + spline_bound*rotordiameter(turb) .LT. turbinexw(&
&         turbi)) THEN
! check this
        wakediameters(zone) = 0.0_dp
        CALL PUSHCONTROL2B(0)
      ELSE IF (zeroloc - spline_bound*rotordiameter(turb) .LT. turbinexw&
&         (turbi)) THEN
!check this
!!!!!!!!!!!!!!!!!!!!!! calculate spline values !!!!!!!!!!!!!!!!!!!!!!!!!!
! position of upwind point
        x1 = zeroloc - spline_bound*rotordiameter(turb)
! diameter of upwind point
        y1 = wakediameter0 + 2.0_dp*ke(turb)*me(zone)*(x1-turbinexw(turb&
&         ))
! slope at upwind point
        dy1 = 2.0_dp*ke(turb)*me(zone)
! position of downwind point
        x2 = zeroloc + spline_bound*rotordiameter(turb)
! diameter at downwind point
        y2 = 0.0_dp
! slope at downwind point
        dy2 = 0.0_dp
! solve for the wake zone diameter and its derivative w.r.t. the downwind
! location at the point of interest
        CALL HERMITE_SPLINE(x, x1, x2, y1, dy1, y2, dy2, wakediameters(&
&                     zone))
        CALL PUSHCONTROL2B(1)
      ELSE IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        wakediameters(zone) = wakediameter0 + 2.0_dp*ke(turb)*me(zone)*&
&         deltax
        CALL PUSHCONTROL2B(2)
      ELSE
        CALL PUSHCONTROL2B(3)
      END IF
      IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        CALL PUSHINTEGER4(zone)
        zone = 2
        wakediameters(zone) = wakediameter0 + 2.0_dp*ke(turb)*me(zone)*&
&         deltax
        zone = 3
        wakediameters(zone) = wakediameter0 + 2.0_dp*ke(turb)*me(zone)*&
&         deltax
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Overlap !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (turbinexw(turbi) .GT. turbinexw(turb)) THEN
! distance between wake center and rotor center
        CALL PUSHREAL4ARRAY(ovdyd, dp/4)
        ovdyd = wakecentery - turbineyw(turbi)
! rotor diameter
        ovr = rotordiameter(turbi)/2
        CALL PUSHINTEGER4(zone)
        DO zone=1,3
! wake diameter
          CALL PUSHREAL4ARRAY(ovrr, dp/4)
          ovrr = wakediameters(zone)/2.0_dp
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
              CALL PUSHREAL4ARRAY(wakeoverlap(zone), dp/4)
              wakeoverlap(zone) = ovrr*ovrr*DACOS(ovl/ovrr) + ovr*ovr*&
&               DACOS((ovdyd-ovl)/ovr) - ovdyd*ovz
              CALL PUSHCONTROL2B(3)
            ELSE IF (ovrr .GT. ovr) THEN
              CALL PUSHREAL4ARRAY(wakeoverlap(zone), dp/4)
              wakeoverlap(zone) = pi*ovr*ovr
              CALL PUSHCONTROL2B(2)
            ELSE
              CALL PUSHREAL4ARRAY(wakeoverlap(zone), dp/4)
              wakeoverlap(zone) = pi*ovrr*ovrr
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE
            CALL PUSHREAL4ARRAY(wakeoverlap(zone), dp/4)
            wakeoverlap(zone) = 0.0_dp
            CALL PUSHCONTROL2B(0)
          END IF
        END DO
        wakeoverlap(3) = wakeoverlap(3) - wakeoverlap(2)
        wakeoverlap(2) = wakeoverlap(2) - wakeoverlap(1)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      CALL PUSHREAL4ARRAY(wakeoverlap, dp*3/4)
      wakeoverlap = wakeoverlap/(pi*rotordiameter(turbi)*rotordiameter(&
&       turbi)/4.0_dp)
! 			print *, wakeOverlap
! removed array effect option
!s = sum(wakeOverlap(1) + wakeOverlap(2))
!*(1+s*keCorrArray)
      kearray = ke(turb)
      CALL PUSHREAL4ARRAY(wakeeffcoeffperzone, dp/4)
      wakeeffcoeffperzone = 0.0_dp
      IF (useaubu) THEN
        x3 = au*pi/180.0_dp + bu*yaw(turb)
        IF (x3 .GE. 0.) THEN
          dabs0 = x3
        ELSE
          dabs0 = -x3
        END IF
        IF (dabs0 .LT. 85.0_dp*pi/180.0_dp) THEN
          CALL PUSHREAL4ARRAY(mmu, dp*3/4)
          mmu = mu/COS(au*pi/180.0_dp+bu*yaw(turb))
          CALL PUSHCONTROL2B(0)
        ELSE
          CALL PUSHREAL4ARRAY(mmu, dp*3/4)
          mmu = mu/COS(85.0_dp*pi/180.0_dp)
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(2)
      END IF
! find overlap-area weighted effect of each wake zone
      IF (deltax .GT. 0 .AND. turbi .NE. turb) THEN
        CALL PUSHINTEGER4(zone)
        DO zone=1,3
          CALL PUSHREAL4ARRAY(rmax, dp/4)
          rmax = cos_spread*0.5_dp*(wakediameters(3)+rotordiameter(turbi&
&           ))
          IF (wakecentery - turbineyw(turbi) .GE. 0.) THEN
            CALL PUSHREAL8(dabs1)
            dabs1 = wakecentery - turbineyw(turbi)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREAL8(dabs1)
            dabs1 = -(wakecentery-turbineyw(turbi))
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL4ARRAY(cosfac, dp/4)
          cosfac = 0.5_dp*(1.0_dp+COS(pi*dabs1/rmax))
          IF (useaubu) THEN
            wakeeffcoeffperzone = wakeeffcoeffperzone + (cosfac*&
&             rotordiameter(turb)/(rotordiameter(turb)+2.0_dp*kearray*&
&             mmu(zone)*deltax))**2*wakeoverlap(zone)
            CALL PUSHCONTROL1B(1)
          ELSE
            wakeeffcoeffperzone = wakeeffcoeffperzone + (cosfac*&
&             rotordiameter(turb)/(rotordiameter(turb)+2.0_dp*kearray*mu&
&             (zone)*deltax))**2*wakeoverlap(zone)
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        wakeeffcoeff = wakeeffcoeff + (a(turb)*wakeeffcoeffperzone)**2
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
! find effective wind speeds at downstream turbines
    CALL PUSHREAL4ARRAY(wakeeffcoeff, dp/4)
    wakeeffcoeff = 1.0_dp - 2.0_dp*SQRT(wakeeffcoeff)
! multiply the inflow speed with the wake coefficients to find effective wind speed at turbine
  END DO
  DO nd=1,nbdirsmax
    rotordiameterb(nd, :) = 0.0
    turbinexwb(nd, :) = 0.0
    turbineywb(nd, :) = 0.0
    ctb(nd, :) = 0.0
    yawb(nd, :) = 0.0
    keb(nd, :) = 0.0
    mmub(nd, :) = 0.0
    wakeoverlapb(nd, :) = 0.0
    wakecenteryb(nd) = 0.0
    ab(nd, :) = 0.0
  END DO
  DO turbi=nturbines,1,-1
    DO nd=1,nbdirsmax
      wakeeffcoeffb(nd) = wtvelocity(turbi)*wtvelocityb(nd, turbi)
      wtvelocityb(nd, turbi) = wakeeffcoeff*wtvelocityb(nd, turbi)
    END DO
    CALL POPREAL4ARRAY(wakeeffcoeff, dp/4)
    DO nd=1,nbdirsmax
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
        DO nd=1,nbdirsmax
          kearrayb(nd) = 0.0
          wakediametersb(nd, :) = 0.0
          deltaxb(nd) = 0.0
        END DO
      ELSE
        DO nd=1,nbdirsmax
          tempb17(nd) = 2*a(turb)*wakeeffcoeffperzone*wakeeffcoeffb(nd)
          ab(nd, turb) = ab(nd, turb) + wakeeffcoeffperzone*tempb17(nd)
          wakeeffcoeffperzoneb(nd) = a(turb)*tempb17(nd)
        END DO
        kearray = ke(turb)
        deltax = turbinexw(turbi) - turbinexw(turb)
        DO nd=1,nbdirsmax
          kearrayb(nd) = 0.0
          wakediametersb(nd, :) = 0.0
          deltaxb(nd) = 0.0
        END DO
        DO zone=3,1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            temp14 = 2.0_dp*mu(zone)
            temp11 = rotordiameter(turb) + temp14*kearray*deltax
            temp13 = rotordiameter(turb)**2
            temp12 = cosfac**2*wakeoverlap(zone)
            DO nd=1,nbdirsmax
              tempb15(nd) = wakeeffcoeffperzoneb(nd)/temp11**2
              tempb16(nd) = -(temp12*temp13*2*tempb15(nd)/temp11)
              cosfacb(nd) = wakeoverlap(zone)*temp13*2*cosfac*tempb15(nd&
&               )
              wakeoverlapb(nd, zone) = wakeoverlapb(nd, zone) + temp13*&
&               cosfac**2*tempb15(nd)
              rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&               tempb16(nd) + temp12*2*rotordiameter(turb)*tempb15(nd)
              kearrayb(nd) = kearrayb(nd) + temp14*deltax*tempb16(nd)
              deltaxb(nd) = deltaxb(nd) + temp14*kearray*tempb16(nd)
            END DO
          ELSE
            temp8 = rotordiameter(turb) + 2.0_dp*mmu(zone)*kearray*&
&             deltax
            temp10 = rotordiameter(turb)**2
            temp9 = cosfac**2*wakeoverlap(zone)
            DO nd=1,nbdirsmax
              tempb12(nd) = wakeeffcoeffperzoneb(nd)/temp8**2
              tempb13(nd) = -(temp9*temp10*2*tempb12(nd)/temp8)
              tempb14(nd) = 2.0_dp*mmu(zone)*tempb13(nd)
              cosfacb(nd) = wakeoverlap(zone)*temp10*2*cosfac*tempb12(nd&
&               )
              wakeoverlapb(nd, zone) = wakeoverlapb(nd, zone) + temp10*&
&               cosfac**2*tempb12(nd)
              rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&               tempb13(nd) + temp9*2*rotordiameter(turb)*tempb12(nd)
              mmub(nd, zone) = mmub(nd, zone) + 2.0_dp*kearray*deltax*&
&               tempb13(nd)
              kearrayb(nd) = kearrayb(nd) + deltax*tempb14(nd)
              deltaxb(nd) = deltaxb(nd) + kearray*tempb14(nd)
            END DO
          END IF
          CALL POPREAL4ARRAY(cosfac, dp/4)
          DO nd=1,nbdirsmax
            tempb11(nd) = -(pi*SIN(pi*(dabs1/rmax))*0.5_dp*cosfacb(nd)/&
&             rmax)
            dabs1b(nd) = tempb11(nd)
            rmaxb(nd) = -(dabs1*tempb11(nd)/rmax)
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(dabs1)
            DO nd=1,nbdirsmax
              wakecenteryb(nd) = wakecenteryb(nd) + dabs1b(nd)
              turbineywb(nd, turbi) = turbineywb(nd, turbi) - dabs1b(nd)
            END DO
          ELSE
            CALL POPREAL8(dabs1)
            DO nd=1,nbdirsmax
              turbineywb(nd, turbi) = turbineywb(nd, turbi) + dabs1b(nd)
              wakecenteryb(nd) = wakecenteryb(nd) - dabs1b(nd)
            END DO
          END IF
          CALL POPREAL4ARRAY(rmax, dp/4)
          DO nd=1,nbdirsmax
            tempb10(nd) = cos_spread*0.5_dp*rmaxb(nd)
            wakediametersb(nd, 3) = wakediametersb(nd, 3) + tempb10(nd)
            rotordiameterb(nd, turbi) = rotordiameterb(nd, turbi) + &
&             tempb10(nd)
          END DO
        END DO
        CALL POPINTEGER4(zone)
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL4ARRAY(mmu, dp*3/4)
        temp7 = au*pi/180.0_dp + bu*yaw(turb)
        temp6 = COS(temp7)
        DO nd=1,nbdirsmax
          yawb(nd, turb) = yawb(nd, turb) - SIN(temp7)*bu*SUM(-(mu*mmub(&
&           nd, :)/temp6))/temp6
        END DO
      ELSE IF (branch .EQ. 1) THEN
        CALL POPREAL4ARRAY(mmu, dp*3/4)
      ELSE
        GOTO 100
      END IF
      DO nd=1,nbdirsmax
        mmub(nd, :) = 0.0
      END DO
 100  CALL POPREAL4ARRAY(wakeeffcoeffperzone, dp/4)
      CALL POPREAL4ARRAY(wakeoverlap, dp*3/4)
      temp5 = pi*rotordiameter(turbi)**2
      DO nd=1,nbdirsmax
        keb(nd, turb) = keb(nd, turb) + kearrayb(nd)
        tempb9(nd, :) = 4.0_dp*wakeoverlapb(nd, :)/temp5
        rotordiameterb(nd, turbi) = rotordiameterb(nd, turbi) + pi*2*&
&         rotordiameter(turbi)*SUM(-(wakeoverlap*tempb9(nd, :)/temp5))
        wakeoverlapb(nd, :) = tempb9(nd, :)
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO nd=1,nbdirsmax
          wakeoverlapb(nd, 1) = wakeoverlapb(nd, 1) - wakeoverlapb(nd, 2&
&           )
          wakeoverlapb(nd, 2) = wakeoverlapb(nd, 2) - wakeoverlapb(nd, 3&
&           )
        END DO
        ovr = rotordiameter(turbi)/2
        DO nd=1,nbdirsmax
          ovdydb(nd) = 0.0
          ovrb(nd) = 0.0
        END DO
        DO zone=3,1,-1
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              CALL POPREAL4ARRAY(wakeoverlap(zone), dp/4)
              DO nd=1,nbdirsmax
                wakeoverlapb(nd, zone) = 0.0
              END DO
              DO nd=1,nbdirsmax
                ovlb(nd) = 0.0
                ovrrb(nd) = 0.0
                ovzb(nd) = 0.0
              END DO
              GOTO 110
            ELSE
              CALL POPREAL4ARRAY(wakeoverlap(zone), dp/4)
              DO nd=1,nbdirsmax
                ovrrb(nd) = pi*2*ovrr*wakeoverlapb(nd, zone)
                wakeoverlapb(nd, zone) = 0.0
              END DO
            END IF
          ELSE IF (branch .EQ. 2) THEN
            CALL POPREAL4ARRAY(wakeoverlap(zone), dp/4)
            DO nd=1,nbdirsmax
              ovrb(nd) = ovrb(nd) + pi*2*ovr*wakeoverlapb(nd, zone)
              wakeoverlapb(nd, zone) = 0.0
            END DO
            DO nd=1,nbdirsmax
              ovrrb(nd) = 0.0
            END DO
          ELSE
            CALL POPREAL4ARRAY(wakeoverlap(zone), dp/4)
            temp3 = ovl/ovrr
            temp4 = (ovdyd-ovl)/ovr
            DO nd=1,nbdirsmax
              IF (temp3 .EQ. 1.0 .OR. temp3 .EQ. (-1.0)) THEN
                tempb7(nd) = 0.0
              ELSE
                tempb7(nd) = -(ovrr*wakeoverlapb(nd, zone)/SQRT(1.D0-&
&                 temp3**2))
              END IF
              IF (temp4 .EQ. 1.0 .OR. temp4 .EQ. (-1.0)) THEN
                tempb8(nd) = 0.0
              ELSE
                tempb8(nd) = -(ovr*wakeoverlapb(nd, zone)/SQRT(1.D0-&
&                 temp4**2))
              END IF
              ovrrb(nd) = DACOS(temp3)*2*ovrr*wakeoverlapb(nd, zone) - &
&               temp3*tempb7(nd)
              ovlb(nd) = tempb7(nd) - tempb8(nd)
              ovrb(nd) = ovrb(nd) + DACOS(temp4)*2*ovr*wakeoverlapb(nd, &
&               zone) - temp4*tempb8(nd)
              ovdydb(nd) = ovdydb(nd) + tempb8(nd) - ovz*wakeoverlapb(nd&
&               , zone)
              ovzb(nd) = -(ovdyd*wakeoverlapb(nd, zone))
              wakeoverlapb(nd, zone) = 0.0
            END DO
            GOTO 110
          END IF
          DO nd=1,nbdirsmax
            ovlb(nd) = 0.0
            ovzb(nd) = 0.0
          END DO
 110      CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4ARRAY(ovz, dp/4)
            DO nd=1,nbdirsmax
              IF (ovz .EQ. 0.0) THEN
                ovzb(nd) = 0.0
              ELSE
                ovzb(nd) = ovzb(nd)/(2.0*SQRT(ovz))
              END IF
            END DO
          ELSE
            DO nd=1,nbdirsmax
              ovzb(nd) = 0.0
            END DO
          END IF
          CALL POPREAL4ARRAY(ovz, dp/4)
          DO nd=1,nbdirsmax
            ovrrb(nd) = ovrrb(nd) + 2*ovrr*ovzb(nd)
            ovlb(nd) = ovlb(nd) - 2*ovl*ovzb(nd)
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4ARRAY(ovl, dp/4)
            DO nd=1,nbdirsmax
              tempb6(nd) = ovlb(nd)/(2.0_dp*ovdyd)
              ovrrb(nd) = ovrrb(nd) + 2*ovrr*tempb6(nd)
              ovrb(nd) = ovrb(nd) - 2*ovr*tempb6(nd)
              ovdydb(nd) = ovdydb(nd) + (2*ovdyd-(ovrr**2-ovr**2+ovdyd**&
&               2)/ovdyd)*tempb6(nd)
            END DO
          ELSE
            CALL POPREAL4ARRAY(ovl, dp/4)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4ARRAY(ovdyd, dp/4)
          ELSE
            CALL POPREAL4ARRAY(ovdyd, dp/4)
            DO nd=1,nbdirsmax
              ovdydb(nd) = -ovdydb(nd)
            END DO
          END IF
          CALL POPREAL4ARRAY(ovrr, dp/4)
          DO nd=1,nbdirsmax
            wakediametersb(nd, zone) = wakediametersb(nd, zone) + ovrrb(&
&             nd)/2.0_dp
          END DO
        END DO
        CALL POPINTEGER4(zone)
        DO nd=1,nbdirsmax
          rotordiameterb(nd, turbi) = rotordiameterb(nd, turbi) + ovrb(&
&           nd)/2
          wakecenteryb(nd) = wakecenteryb(nd) + ovdydb(nd)
          turbineywb(nd, turbi) = turbineywb(nd, turbi) - ovdydb(nd)
        END DO
        CALL POPREAL4ARRAY(ovdyd, dp/4)
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        zone = 3
        DO nd=1,nbdirsmax
          tempb4(nd) = me(zone)*2.0_dp*wakediametersb(nd, zone)
          wakediameter0b(nd) = wakediametersb(nd, zone)
          wakediametersb(nd, zone) = 0.0
        END DO
        zone = 2
        DO nd=1,nbdirsmax
          tempb5(nd) = me(zone)*2.0_dp*wakediametersb(nd, zone)
          keb(nd, turb) = keb(nd, turb) + deltax*tempb5(nd) + deltax*&
&           tempb4(nd)
          deltaxb(nd) = deltaxb(nd) + ke(turb)*tempb5(nd) + ke(turb)*&
&           tempb4(nd)
          wakediameter0b(nd) = wakediameter0b(nd) + wakediametersb(nd, &
&           zone)
          wakediametersb(nd, zone) = 0.0
        END DO
        CALL POPINTEGER4(zone)
      ELSE
        DO nd=1,nbdirsmax
          wakediameter0b(nd) = 0.0
        END DO
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          DO nd=1,nbdirsmax
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
          CALL HERMITE_SPLINE_BV(x, xb, x1, x1b, x2, x2b, y1, y1b, dy1&
&                           , dy1b, y2, dy2, wakediameters(zone), &
&                           wakediametersb(1, zone), nbdirsmax)
          DO nd=1,nbdirsmax
            tempb2(nd) = me(zone)*2.0_dp*y1b(nd)
            x1b(nd) = x1b(nd) + ke(turb)*tempb2(nd)
            wakediametersb(nd, zone) = 0.0
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
          DO nd=1,nbdirsmax
            tempb3(nd) = me(zone)*2.0_dp*wakediametersb(nd, zone)
            wakediameter0b(nd) = wakediameter0b(nd) + wakediametersb(nd&
&             , zone)
            keb(nd, turb) = keb(nd, turb) + deltax*tempb3(nd)
            deltaxb(nd) = deltaxb(nd) + ke(turb)*tempb3(nd)
          END DO
        END IF
        DO nd=1,nbdirsmax
          xb(nd) = 0.0
          zerolocb(nd) = 0.0
        END DO
      END IF
      temp2 = 2.0_dp*me(zone)*ke(turb)
      DO nd=1,nbdirsmax
        turbinexwb(nd, turb) = turbinexwb(nd, turb) + zerolocb(nd)
        wakediameter0b(nd) = wakediameter0b(nd) - zerolocb(nd)/temp2
        keb(nd, turb) = keb(nd, turb) + wakediameter0*2.0_dp*me(zone)*&
&         zerolocb(nd)/temp2**2
        turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + xb(nd)
      END DO
      CALL POPINTEGER4(zone)
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        DO nd=1,nbdirsmax
          deltaxb(nd) = deltaxb(nd) + bd*wakecenteryb(nd)
        END DO
      ELSE IF (branch .NE. 1) THEN
        DO nd=1,nbdirsmax
          wakeangleinitb(nd) = 0.0
        END DO
        GOTO 120
      END IF
      factor = 2.0_dp*kd*deltax/rotordiameter(turb) + 1.0_dp
      temp1 = 30.0_dp*kd*factor**5
      temp0 = wakeangleinit*rotordiameter(turb)
      temp = wakeangleinit**2 + 15.0_dp*factor**4
      DO nd=1,nbdirsmax
        displacementb(nd) = wakecenteryb(nd)
        tempb(nd) = -((wakeangleinit**2+15.0_dp)*displacementb(nd)/(&
&         30.0_dp*kd))
        tempb0(nd) = displacementb(nd)/temp1
        wakeangleinitb(nd) = (temp*rotordiameter(turb)+temp0*2*&
&         wakeangleinit)*tempb0(nd) + rotordiameter(turb)*tempb(nd) - &
&         wakeangleinit**2*rotordiameter(turb)*2*displacementb(nd)/(&
&         30.0_dp*kd)
        factorb(nd) = (15.0_dp*temp0*4*factor**3-30.0_dp*kd*temp*temp0*5&
&         *factor**4/temp1)*tempb0(nd)
        turbineywb(nd, turb) = turbineywb(nd, turb) + wakecenteryb(nd)
        tempb1(nd) = kd*2.0_dp*factorb(nd)/rotordiameter(turb)
        rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + temp*&
&         wakeangleinit*tempb0(nd) - deltax*tempb1(nd)/rotordiameter(&
&         turb) + wakeangleinit*tempb(nd)
        deltaxb(nd) = deltaxb(nd) + tempb1(nd)
      END DO
      DO nd=1,nbdirsmax
        wakecenteryb(nd) = 0.0
      END DO
 120  CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL4ARRAY(wakediameter0, dp/4)
        DO nd=1,nbdirsmax
          rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + COS(yaw(&
&           turb))*wakediameter0b(nd)
          yawb(nd, turb) = yawb(nd, turb) - rotordiameter(turb)*SIN(yaw(&
&           turb))*wakediameter0b(nd)
        END DO
      ELSE
        CALL POPREAL4ARRAY(wakediameter0, dp/4)
        DO nd=1,nbdirsmax
          rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&           wakediameter0b(nd)
        END DO
      END IF
      CALL POPREAL4ARRAY(wakeangleinit, dp/4)
      CALL POPREAL4ARRAY(deltax, dp/4)
      DO nd=1,nbdirsmax
        yawb(nd, turb) = yawb(nd, turb) + ct(turb)*0.5_dp*COS(yaw(turb))&
&         *wakeangleinitb(nd)
        ctb(nd, turb) = ctb(nd, turb) + 0.5_dp*SIN(yaw(turb))*&
&         wakeangleinitb(nd)
        turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + deltaxb(nd)
        turbinexwb(nd, turb) = turbinexwb(nd, turb) - deltaxb(nd)
      END DO
    END DO
    CALL POPREAL4ARRAY(wakeeffcoeff, dp/4)
  END DO
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    CALL CTTOAXIALIND_BV(ct, ctb, nturbines, a, ab, nbdirsmax)
    DO nd=1,nbdirsmax
      a_inb(nd, :) = 0.0
    END DO
  ELSE
    DO nd=1,nbdirsmax
      a_inb(nd, :) = 0.0
      a_inb(nd, :) = ab(nd, :)
    END DO
  END IF
  DO nd=1,nbdirsmax
    ctb(nd, :) = ctb(nd, :) + kecorrct*keb(nd, :)
    yawdegb(nd, :) = 0.0
    yawdegb(nd, :) = pi*yawb(nd, :)/180.0_dp
  END DO
  DO nd=1,nbdirsmax
    wtvelocityb(nd, :) = 0.0
  END DO
END SUBROUTINE FLORIS_BV

!  Differentiation of hermite_spline in reverse (adjoint) mode:
!   gradient     of useful results: y
!   with respect to varying inputs: x x0 x1 dy0 y0
! Flow field calculations have been intentionally left out to save development time.
! The flow field can be calculated using the pure python version of floris 
! This implementation is fully smooth and differentiable with the exception of a 
! discontinuity at the hub of each turbine. The discontinuity only presents issues if
! turbines are place within 1E-15 * rotor diameter of one another, which is extremely 
! unlikely during optimization if the user does not explicitly place them there.
SUBROUTINE HERMITE_SPLINE_BV(x, xb, x0, x0b, x1, x1b, y0, y0b, dy0, &
& dy0b, y1, dy1, y, yb, nbdirsmax)
  ! USE DIFFSIZES
!  Hint: nbdirsmax should be the maximum number of differentiation directions
  IMPLICIT NONE
!dy_dx = c3*3*x**2 + c2*2*x + c1
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: x, x0, x1, y0, dy0, y1, dy1
  REAL(dp), DIMENSION(nbdirsmax) :: xb, x0b, x1b, y0b, dy0b
! out
!, dy_dx
  REAL(dp) :: y
  REAL(dp), DIMENSION(nbdirsmax) :: yb
! local
  REAL(dp) :: c3, c2, c1, c0
  REAL(dp), DIMENSION(nbdirsmax) :: c3b, c2b, c1b, c0b
  INTRINSIC KIND
  REAL(dp) :: temp
  REAL(dp) :: temp0
  REAL(dp) :: temp1
  REAL(dp) :: temp2
  REAL(dp) :: temp3
  REAL(dp) :: temp4
  REAL(dp) :: temp5
  REAL(dp) :: temp6
  REAL(dp) :: temp7
  REAL(dp) :: temp8
  REAL(dp) :: temp9
  REAL(dp) :: temp10
  REAL(dp) :: temp11
  REAL(dp) :: temp12
  REAL(dp) :: temp13
  REAL(dp) :: temp14
  REAL(dp) :: temp15
  REAL(dp) :: temp16
  INTEGER :: nd
  REAL(dp), DIMENSION(nbdirsmax) :: tempb
  REAL(dp), DIMENSION(nbdirsmax) :: tempb0
  REAL(dp), DIMENSION(nbdirsmax) :: tempb1
  REAL(dp), DIMENSION(nbdirsmax) :: tempb2
  REAL(dp), DIMENSION(nbdirsmax) :: tempb3
  REAL(dp), DIMENSION(nbdirsmax) :: tempb4
  REAL(dp), DIMENSION(nbdirsmax) :: tempb5
  REAL(dp), DIMENSION(nbdirsmax) :: tempb6
  REAL(dp), DIMENSION(nbdirsmax) :: tempb7
  REAL(dp), DIMENSION(nbdirsmax) :: tempb8
  REAL(dp), DIMENSION(nbdirsmax) :: tempb9
  REAL(dp), DIMENSION(nbdirsmax) :: tempb10
  REAL(dp), DIMENSION(nbdirsmax) :: tempb11
  REAL(dp), DIMENSION(nbdirsmax) :: tempb12
  REAL(dp), DIMENSION(nbdirsmax) :: tempb13
  REAL(dp), DIMENSION(nbdirsmax) :: tempb14
  REAL(dp), DIMENSION(nbdirsmax) :: tempb15
  REAL(dp), DIMENSION(nbdirsmax) :: tempb16
  REAL(dp), DIMENSION(nbdirsmax) :: tempb17
  REAL(dp), DIMENSION(nbdirsmax) :: tempb18
  REAL(dp), DIMENSION(nbdirsmax) :: tempb19
  REAL(dp), DIMENSION(nbdirsmax) :: tempb20
  REAL(dp), DIMENSION(nbdirsmax) :: tempb21
  REAL(dp), DIMENSION(nbdirsmax) :: tempb22
  REAL(dp), DIMENSION(nbdirsmax) :: tempb23
  REAL(dp), DIMENSION(nbdirsmax) :: tempb24
  REAL(dp), DIMENSION(nbdirsmax) :: tempb25
  REAL(dp), DIMENSION(nbdirsmax) :: tempb26
  REAL(dp), DIMENSION(nbdirsmax) :: tempb27
  REAL(dp), DIMENSION(nbdirsmax) :: tempb28
  REAL(dp), DIMENSION(nbdirsmax) :: tempb29
  REAL(dp), DIMENSION(nbdirsmax) :: tempb30
  INTEGER :: nbdirsmax
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
  DO nd=1,nbdirsmax
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

!  Differentiation of cttoaxialind in reverse (adjoint) mode:
!   gradient     of useful results: axial_induction ct
!   with respect to varying inputs: ct
SUBROUTINE CTTOAXIALIND_BV(ct, ctb, nturbines, axial_induction, &
& axial_inductionb, nbdirsmax)
  !USE DIFFSIZES
!  Hint: nbdirsmax should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: ct
  REAL(dp), DIMENSION(nbdirsmax, nturbines) :: ctb
! local
  INTEGER :: i
! out
  REAL(dp), DIMENSION(nturbines) :: axial_induction
  REAL(dp), DIMENSION(nbdirsmax, nturbines) :: axial_inductionb
  INTRINSIC KIND
  INTRINSIC SQRT
  INTEGER :: nd
  INTEGER :: branch
  INTEGER :: nbdirsmax
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
      DO nd=1,nbdirsmax
        IF (.NOT.1.0_dp - ct(i) .EQ. 0.0) ctb(nd, i) = ctb(nd, i) + &
&           0.5_dp*axial_inductionb(nd, i)/(2.0*SQRT(1.0_dp-ct(i)))
        axial_inductionb(nd, i) = 0.0
      END DO
    ELSE
      DO nd=1,nbdirsmax
        IF (.NOT.0.0203_dp - 0.6427_dp*(0.889_dp-ct(i)) .EQ. 0.0) ctb(nd&
&         , i) = ctb(nd, i) + 0.6427_dp*axial_inductionb(nd, i)/(2.0*&
&           SQRT(0.0203_dp-0.6427_dp*(0.889_dp-ct(i))))
        axial_inductionb(nd, i) = 0.0
      END DO
    END IF
  END DO
END SUBROUTINE CTTOAXIALIND_BV


