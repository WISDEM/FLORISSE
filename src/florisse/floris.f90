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

! define precision to be the standard for a double precision! on local system
    integer, parameter :: dp = kind(0.d0)

! in
    real(dp), intent(in) :: x, x0, x1, y0, dy0, y1, dy1

! out
    real(dp), intent(out) :: y
!, dy_dx

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


subroutine calcOverlapAreas(nTurbines, turbineX, turbineY, turbineZ, rotorDiameter, wakeDiameters, &
                            wakeCentersYT, wakeCentersZT, wakeOverlapTRel_mat)
!    calculate overlap of rotors and wake zones (wake zone location defined by wake
!    center and wake diameter)
!   turbineX,turbineY is x,y-location of center of rotor
!
!    wakeOverlap(TURBI,TURB,ZONEI) = overlap area of zone ZONEI of wake of turbine
!     TURB with rotor of downstream turbine
!    TURBI

    implicit none

! define precision to be the standard for a double precision! on local system
    integer, parameter :: dp = kind(0.d0)

! in
    integer, intent(in) :: nTurbines
    real(dp), dimension(nTurbines), intent(in) :: turbineX, turbineY, turbineZ, rotorDiameter
    real(dp), dimension(nTurbines, nTurbines, 3), intent(in) :: wakeDiameters
    real(dp), dimension(nTurbines, nTurbines), intent(in) :: wakeCentersYT, wakeCentersZT

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
                OVdYd = sqrt((wakeCentersYT(turbI, turb)-turbineY(turbI))**2+(wakeCentersZT(turbI, turb)-turbineZ(turbI))**2)
! distance between wake center and rotor center
                OVr = rotorDiameter(turbI)/2
! rotor diameter
                do zone = 1, 3
                    OVRR = wakeDiameters(turbI, turb, zone)/2.0_dp
! wake diameter
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

                    if (OVdYd < (OVr+OVRR)) then
! if the rotor overlaps the wake zone

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

! define precision to be the standard for a double precision! on local system
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
        if (CT(i) > 0.96) then
! Glauert condition
            axial_induction(i) = 0.143_dp + sqrt(0.0203_dp-0.6427_dp*(0.889_dp - CT(i)))
        else
            axial_induction(i) = 0.5_dp*(1.0_dp-sqrt(1.0_dp-CT(i)))
        end if
    end do

end subroutine CTtoAxialInd


subroutine floris(nTurbines, nSamples, turbineXw, turbineYw, turbineZ, yawDeg, &
                          & rotorDiameter, Vinf, shearExp, zref, z0, Ct, a_in, ke_in, kd, me, &
                          & initialWakeDisplacement, bd, MU, aU, bU, initialWakeAngle, &
                          & cos_spread, keCorrCT, Region2CT, keCorrArray, useWakeAngle, &
                          & adjustInitialWakeDiamToYaw, axialIndProvided, useaUbU, &
                          & wsPositionXYZw, shearCoefficientAlpha, shearZh, &
                          & wtVelocity, wsArray, wakeCentersYT_vec, wakeCentersZT_vec, &
                          & wakeDiametersT_vec, wakeOverlapTRel_vec)

! independent variables: yawDeg Ct turbineXw turbineYw rotorDiameter a_in
! dependent variables: wtVelocity

    implicit none

! define precision to be the standard for a double precision! on local system
    integer, parameter :: dp = kind(0.d0)

! in
    integer, intent(in) :: nTurbines
    integer :: nSamples
    real(dp), intent(in) :: kd, initialWakeDisplacement, initialWakeAngle, ke_in
    real(dp), intent(in) :: keCorrCT, Region2CT, bd, cos_spread, Vinf, keCorrArray, shearExp, zref, z0
    real(dp), dimension(nTurbines), intent(in) :: yawDeg, Ct, a_in, turbineXw, turbineYw, turbineZ
    real(dp), dimension(nTurbines), intent(in) :: rotorDiameter
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
    real(dp) :: spline_bound
! in rotor diameters
    real(dp) :: wakeAngleInit, zeroloc
    real(dp) :: factor, displacement, x, x1, x2, y1, y2, dy1, dy2
    real(dp) :: wakeDiameter0
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeDiametersT_mat
    real(dp), dimension(nTurbines, nTurbines) :: wakeCentersYT_mat, wakeCentersZT_mat
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
    real(dp), dimension(nTurbines) :: adjusted_wind_speed
!wind speed at each turbine


! visualization
    real(dp) :: deltay, deltaz, radiusLoc, axialIndAndNearRotor, reductionFactor

! model out
    real(dp), dimension(nTurbines), intent(out) :: wtVelocity

! test out
    real(dp), dimension(nTurbines*nTurbines), intent(out) :: wakeCentersYT_vec, wakecentersZT_vec
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


    do turb = 1, nTurbines
! check that the point of interest is above ground level
      if (turbineZ(turb) >= 0.0) then
! adjusted wind speed for wind shear if point is above ground
          adjusted_wind_speed(turb) = Vinf*((turbineZ(turb)-z0)/(zref-z0))**shearExp
      else
! if the point of interest is below ground, set the wind speed to 0.0
          adjusted_wind_speed(turb) = 0.0_dp
      end if

    end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Centers and Diameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!THESE ARE THE Z CALCULATIONS
                wakeCentersZT_mat(turbI, turb) = turbineZ(turb)

!THESE ARE THE Y CALCULATIONS
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

                wakeCentersZ(loc, turb) = turbineZ(turb)

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

            if (zeroloc + spline_bound*rotorDiameter(turb) < turbineXw(turbI)) then
! check this
                wakeDiametersT_mat(turbI, turb, zone) = 0.0_dp

            else if (zeroloc - spline_bound*rotorDiameter(turb) < turbineXw(turbI)) then
!check this

!!!!!!!!!!!!!!!!!!!!!! calculate spline values!!!!!!!!!!!!!!!!!!!!!!!!!!

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

            if (zeroloc + spline_bound*rotorDiameter(turb) < velX(loc)) then
! check this
                wakeDiameters(loc, turb, zone) = 0.0_dp

            else if (zeroloc - spline_bound*rotorDiameter(turb) < velX(loc)) then
!check this

!!!!!!!!!!!!!!!!!!!!!! calculate spline values!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Overlap!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! calculate relative overlap
    call calcOverlapAreas(nTurbines, turbineXw, turbineYw, turbineZ, rotorDiameter, &
                          & wakeDiametersT_mat, wakeCentersYT_mat, wakeCentersZT_mat, wakeOverlapTRel_mat)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Velocity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialize velocities in full flow field (optional)
    wsArray = minval(adjusted_wind_speed)

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
    wtVelocity = adjusted_wind_speed
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

            if (deltax > 0 .and. radiusLoc < wakeDiameters(loc, turb, 1)/2.0_dp) then
! check if in zone 1
                reductionFactor = axialIndAndNearRotor*&
                                  & (cosFac*rotorDiameter(turb)/(rotorDiameter(turb)+2.0_dp* &
                                  & keArray(turb)*(mmU(1))*deltax))**2
            else if (deltax > 0 .and. radiusLoc < wakeDiameters(loc, turb, 2)/2.0_dp) then
! check if in zone 2
                reductionFactor = axialIndAndNearRotor* &
                                  & (cosFac*rotorDiameter(turb)/(rotorDiameter(turb)+2.0_dp* &
                                  & keArray(turb)*(mmU(2))*deltax))**2
            else if (deltax > 0 .and. radiusLoc < wakeDiameters(loc, turb, 3)/2.0_dp) then
! check if in zone 3
                reductionFactor = axialIndAndNearRotor* &
                                  (cosFac*rotorDiameter(turb)/(rotorDiameter(turb)+2.0_dp* &
                                  & keArray(turb)*(mmU(3))*deltax))**2
! use this to add upstream turbine influence to visualization
! else if (deltax <= 0 .and. radiusLoc < rotorDiameter(turb)/2.0_dp) then! check if axial induction zone in front of rotor
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

        wakeCentersZT_vec(nTurbines*(turbI-1)+1:nTurbines*(turbI-1)+nTurbines) &
                                = wakeCentersZT_mat(turbI, :)

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

end subroutine floris





!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.13 (r6666M) -  1 Mar 2018 15:30
!
!  Differentiation of floris in reverse (adjoint) mode:
!   gradient     of useful results: wtvelocity
!   with respect to varying inputs: rotordiameter turbinez turbinexw
!                wtvelocity vinf turbineyw yawdeg ct a_in
!   RW status of diff variables: rotordiameter:out turbinez:out
!                turbinexw:out wtvelocity:in-zero vinf:out turbineyw:out
!                yawdeg:out ct:out a_in:out
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.13 (r6666M) -  1 Mar 2018 15:30
!
!  Differentiation of floris in reverse (adjoint) mode:
!   gradient     of useful results: wtvelocity
!   with respect to varying inputs: rotordiameter turbinez turbinexw
!                wtvelocity vinf turbineyw yawdeg ct a_in
!   RW status of diff variables: rotordiameter:out turbinez:out
!                turbinexw:out wtvelocity:in-zero vinf:out turbineyw:out
!                yawdeg:out ct:out a_in:out
SUBROUTINE FLORIS_BV(nturbines, turbinexw, turbinexwb, &
& turbineyw, turbineywb, turbinez, turbinezb, yawdeg, yawdegb, &
& rotordiameter, rotordiameterb, vinf, vinfb, shearexp, zref, z0, ct, &
& ctb, a_in, a_inb, ke_in, kd, me, initialwakedisplacement, bd, mu, au, &
& bu, initialwakeangle, cos_spread, kecorrct, region2ct, kecorrarray, &
& usewakeangle, adjustinitialwakediamtoyaw, axialindprovided, useaubu, &
& wtvelocityb, nbdirs)
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), INTENT(IN) :: kd, initialwakedisplacement, initialwakeangle&
& , ke_in
  REAL(dp), INTENT(IN) :: kecorrct, region2ct, bd, cos_spread, vinf, &
& kecorrarray, shearexp, zref, z0
  REAL(dp), DIMENSION(nbdirs), INTENT(out) :: vinfb
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: yawdeg, ct, a_in, &
& turbinexw, turbineyw, turbinez
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(out) :: yawdegb, ctb, a_inb, &
& turbinexwb, turbineywb, turbinezb
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: rotordiameter
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(out) :: rotordiameterb
  REAL(dp), DIMENSION(3), INTENT(IN) :: me, mu
  REAL(dp), INTENT(IN) :: au, bu
  LOGICAL, INTENT(IN) :: usewakeangle, adjustinitialwakediamtoyaw, &
& axialindprovided, useaubu
! local (General)
  REAL(dp), DIMENSION(nturbines) :: ke, yaw
  REAL(dp), DIMENSION(nbdirs, nturbines) :: keb, yawb
  REAL(dp) :: deltax
  REAL(dp), DIMENSION(nbdirs) :: deltaxb
  INTEGER :: turb, turbi, zone
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp
! visualization
  INTEGER :: loc
! local (Wake centers and diameters)
  REAL(dp) :: spline_bound
! in rotor diameters
  REAL(dp) :: wakeangleinit, zeroloc
  REAL(dp), DIMENSION(nbdirs) :: wakeangleinitb, zerolocb
  REAL(dp) :: factor, displacement, x, x1, x2, y1, y2, dy1, dy2
  REAL(dp), DIMENSION(nbdirs) :: factorb, displacementb, xb, x1b, x2b&
& , y1b, dy1b
  REAL(dp) :: wakediameter0
  REAL(dp), DIMENSION(nbdirs) :: wakediameter0b
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakediameterst_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakediameterst_matb
  REAL(dp), DIMENSION(nturbines, nturbines) :: wakecentersyt_mat, &
& wakecenterszt_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines) :: &
& wakecentersyt_matb, wakecenterszt_matb
! local (Wake overlap)
  REAL(dp) :: rmax
  REAL(dp), DIMENSION(nbdirs) :: rmaxb
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlaptrel_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlaptrel_matb
! local (Velocity)
  REAL(dp), DIMENSION(nturbines) :: a, kearray
  REAL(dp), DIMENSION(nbdirs, nturbines) :: ab, kearrayb
  REAL(dp), DIMENSION(3) :: mmu
  REAL(dp), DIMENSION(nbdirs, 3) :: mmub
  REAL(dp) :: s, cosfac, wakeeffcoeff, wakeeffcoeffperzone
  REAL(dp), DIMENSION(nbdirs) :: sb, cosfacb, wakeeffcoeffb, &
& wakeeffcoeffperzoneb
  REAL(dp), DIMENSION(nturbines) :: adjusted_wind_speed
  REAL(dp), DIMENSION(nbdirs, nturbines) :: adjusted_wind_speedb
!wind speed at each turbine
! visualization
  REAL(dp) :: deltay, deltaz, radiusloc, axialindandnearrotor, &
& reductionfactor
! model out
  REAL(dp), DIMENSION(nturbines) :: wtvelocity
  REAL(dp), DIMENSION(nbdirs, nturbines) :: wtvelocityb
! visualization out
  INTRINSIC COS, ATAN, MAX
  INTRINSIC SIN
  INTRINSIC MINVAL
  INTRINSIC SUM
  INTRINSIC DABS
  INTRINSIC SQRT
  DOUBLE PRECISION :: dabs0
  DOUBLE PRECISION :: dabs1
  DOUBLE PRECISION, DIMENSION(nbdirs) :: dabs1b
  DOUBLE PRECISION :: dabs2
  INTEGER :: nd
  REAL(dp) :: temp
  REAL(dp) :: temp0
  REAL(dp) :: temp1
  REAL(dp) :: temp2
  REAL(dp), DIMENSION(nbdirs) :: tempb
  REAL(dp), DIMENSION(nbdirs) :: tempb0
  REAL(dp), DIMENSION(nbdirs) :: tempb1
  REAL(dp) :: temp3
  REAL(dp), DIMENSION(nbdirs) :: tempb2
  REAL(dp), DIMENSION(nbdirs) :: tempb3
  REAL(dp), DIMENSION(nbdirs) :: tempb4
  REAL(dp), DIMENSION(nbdirs) :: tempb5
  REAL(dp) :: temp4
  REAL(dp) :: temp5
  REAL(dp), DIMENSION(nbdirs) :: tempb6
  REAL*8, DIMENSION(nbdirs) :: tempb7
  REAL(dp) :: temp6
  REAL(dp) :: temp7
  REAL(dp) :: temp8
  REAL(dp), DIMENSION(nbdirs) :: tempb8
  REAL(dp), DIMENSION(nbdirs) :: tempb9
  REAL(dp), DIMENSION(nbdirs) :: tempb10
  REAL(dp) :: temp9
  REAL(dp) :: temp10
  REAL(dp) :: temp11
  REAL(dp) :: temp12
  REAL(dp), DIMENSION(nbdirs) :: tempb11
  REAL(dp), DIMENSION(nbdirs) :: tempb12
  REAL(dp), DIMENSION(nbdirs) :: tempb13
  INTEGER :: branch
  INTEGER :: nbdirs
  DOUBLE PRECISION :: x3
  yaw = yawdeg*pi/180.0_dp
  DO turb=1,nturbines
! check that the point of interest is above ground level
    IF (turbinez(turb) .GE. 0.0) THEN
! adjusted wind speed for wind shear if point is above ground
      adjusted_wind_speed(turb) = vinf*((turbinez(turb)-z0)/(zref-z0))**&
&       shearexp
      CALL PUSHCONTROL1B(1)
    ELSE
! if the point of interest is below ground, set the wind speed to 0.0
      adjusted_wind_speed(turb) = 0.0_dp
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Centers and Diameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  spline_bound = 1.0_dp
! calculate y-locations of wake centers in wind ref. frame
  wakecentersyt_mat = 0.0_dp
  DO turb=1,nturbines
    CALL PUSHREAL8(wakeangleinit)
    wakeangleinit = 0.5_dp*SIN(yaw(turb))*ct(turb)
    IF (usewakeangle) wakeangleinit = wakeangleinit + initialwakeangle*&
&       pi/180.0_dp
! wake center calculations at each turbine
    DO turbi=1,nturbines
      IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        CALL PUSHREAL8(deltax)
        deltax = turbinexw(turbi) - turbinexw(turb)
        factor = 2.0_dp*kd*deltax/rotordiameter(turb) + 1.0_dp
!THESE ARE THE Z CALCULATIONS
        wakecenterszt_mat(turbi, turb) = turbinez(turb)
!THESE ARE THE Y CALCULATIONS
        wakecentersyt_mat(turbi, turb) = turbineyw(turb)
        displacement = wakeangleinit*(wakeangleinit*wakeangleinit+&
&         15.0_dp*factor*factor*factor*factor)/(30.0_dp*kd/rotordiameter&
&         (turb)*(factor*factor*factor*factor*factor))
        displacement = displacement - wakeangleinit*(wakeangleinit*&
&         wakeangleinit+15.0_dp)/(30.0_dp*kd/rotordiameter(turb))
        wakecentersyt_mat(turbi, turb) = wakecentersyt_mat(turbi, turb) &
&         + initialwakedisplacement + displacement
        IF (usewakeangle .EQV. .false.) THEN
          wakecentersyt_mat(turbi, turb) = wakecentersyt_mat(turbi, turb&
&           ) + bd*deltax
          CALL PUSHCONTROL2B(2)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(0)
      END IF
    END DO
  END DO
!adjust k_e to C_T, adjusted to yaw
  ke = ke_in + kecorrct*(ct-region2ct)
! calculate wake diameters
  wakediameterst_mat = 0.0_dp
  DO turb=1,nturbines
    IF (adjustinitialwakediamtoyaw) THEN
      CALL PUSHREAL8(wakediameter0)
      wakediameter0 = rotordiameter(turb)*COS(yaw(turb))
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHREAL8(wakediameter0)
      wakediameter0 = rotordiameter(turb)
      CALL PUSHCONTROL1B(0)
    END IF
! calculate the wake diameter of each wake at each turbine
    DO turbi=1,nturbines
! turbine separation
      CALL PUSHREAL8(deltax)
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
! check this
        wakediameterst_mat(turbi, turb, zone) = 0.0_dp
        CALL PUSHCONTROL2B(0)
      ELSE IF (zeroloc - spline_bound*rotordiameter(turb) .LT. turbinexw&
&         (turbi)) THEN
!check this
!!!!!!!!!!!!!!!!!!!!!! calculate spline values!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        CALL HERMITE_SPLINE(x, x1, x2, y1, dy1, y2, dy2, &
&                     wakediameterst_mat(turbi, turb, zone))
        CALL PUSHCONTROL2B(1)
      ELSE IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        wakediameterst_mat(turbi, turb, zone) = wakediameter0 + 2.0_dp*&
&         ke(turb)*me(zone)*deltax
        CALL PUSHCONTROL2B(2)
      ELSE
        CALL PUSHCONTROL2B(3)
      END IF
      IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        CALL PUSHINTEGER4(zone)
        zone = 2
        wakediameterst_mat(turbi, turb, zone) = wakediameter0 + 2.0_dp*&
&         ke(turb)*me(zone)*deltax
        zone = 3
        wakediameterst_mat(turbi, turb, zone) = wakediameter0 + 2.0_dp*&
&         ke(turb)*me(zone)*deltax
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Overlap!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate relative overlap
  CALL CALCOVERLAPAREAS(nturbines, turbinexw, turbineyw, turbinez, &
&                 rotordiameter, wakediameterst_mat, wakecentersyt_mat, &
&                 wakecenterszt_mat, wakeoverlaptrel_mat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Velocity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize velocities in full flow field (optional)
! apply shear profile to visualization
!     wsArray = wsArray*(velZ/shearZh)**shearCoefficientAlpha
!     print *, wsArray
! initialize axial induction values
  IF (axialindprovided) THEN
    a = a_in
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL CTTOAXIALIND(ct, nturbines, a)
    CALL PUSHCONTROL1B(1)
  END IF
! adjust ke to Ct as adjusted to yaw
  CALL PUSHREAL8ARRAY(ke, nturbines)
  ke = ke_in + kecorrct*(ct-region2ct)
  DO turb=1,nturbines
    CALL PUSHREAL8(s)
    s = SUM(wakeoverlaptrel_mat(turb, :, 1) + wakeoverlaptrel_mat(turb, &
&     :, 2))
    kearray(turb) = ke(turb)*(1+s*kecorrarray)
  END DO
! find effective wind speeds at downstream turbines
  wtvelocity = adjusted_wind_speed
  DO turbi=1,nturbines
    CALL PUSHREAL8(wakeeffcoeff)
    wakeeffcoeff = 0.0_dp
! find overlap-area weighted effect of each wake zone
    DO turb=1,nturbines
      CALL PUSHREAL8(wakeeffcoeffperzone)
      wakeeffcoeffperzone = 0.0_dp
      CALL PUSHREAL8(deltax)
      deltax = turbinexw(turbi) - turbinexw(turb)
      IF (useaubu) THEN
        x3 = au*pi/180.0_dp + bu*yaw(turb)
        IF (x3 .GE. 0.) THEN
          dabs0 = x3
        ELSE
          dabs0 = -x3
        END IF
        IF (dabs0 .LT. 85.0_dp*pi/180.0_dp) THEN
          CALL PUSHREAL8ARRAY(mmu, 3)
          mmu = mu/COS(au*pi/180.0_dp+bu*yaw(turb))
          CALL PUSHCONTROL2B(0)
        ELSE
          CALL PUSHREAL8ARRAY(mmu, 3)
          mmu = mu/COS(85.0_dp*pi/180.0_dp)
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(2)
      END IF
      IF (deltax .GT. 0 .AND. turbi .NE. turb) THEN
        CALL PUSHINTEGER4(zone)
        DO zone=1,3
          rmax = cos_spread*0.5_dp*(wakediameterst_mat(turbi, turb, 3)+&
&           rotordiameter(turbi))
          IF (wakecentersyt_mat(turbi, turb) - turbineyw(turbi) .GE. 0.&
&         ) THEN
            CALL PUSHREAL8(dabs1)
            dabs1 = wakecentersyt_mat(turbi, turb) - turbineyw(turbi)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREAL8(dabs1)
            dabs1 = -(wakecentersyt_mat(turbi, turb)-turbineyw(turbi))
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8(cosfac)
          cosfac = 0.5_dp*(1.0_dp+COS(pi*dabs1/rmax))
          IF (useaubu) THEN
            wakeeffcoeffperzone = wakeeffcoeffperzone + (cosfac*&
&             rotordiameter(turb)/(rotordiameter(turb)+2.0_dp*kearray(&
&             turb)*mmu(zone)*deltax))**2*wakeoverlaptrel_mat(turbi, &
&             turb, zone)
            CALL PUSHCONTROL1B(1)
          ELSE
            wakeeffcoeffperzone = wakeeffcoeffperzone + (cosfac*&
&             rotordiameter(turb)/(rotordiameter(turb)+2.0_dp*kearray(&
&             turb)*mu(zone)*deltax))**2*wakeoverlaptrel_mat(turbi, turb&
&             , zone)
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        wakeeffcoeff = wakeeffcoeff + (a(turb)*wakeeffcoeffperzone)**2
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    CALL PUSHREAL8(wakeeffcoeff)
    wakeeffcoeff = 1.0_dp - 2.0_dp*SQRT(wakeeffcoeff)
! multiply the inflow speed with the wake coefficients to find effective wind
! speed at turbine
  END DO
  DO nd=1,nbdirs
    rotordiameterb(nd, :) = 0.0_8
    turbinexwb(nd, :) = 0.0_8
    turbineywb(nd, :) = 0.0_8
    kearrayb(nd, :) = 0.0_8
    yawb(nd, :) = 0.0_8
    wakeoverlaptrel_matb(nd, :, :, :) = 0.0_8
    wakediameterst_matb(nd, :, :, :) = 0.0_8
    wakecentersyt_matb(nd, :, :) = 0.0_8
    mmub(nd, :) = 0.0_8
    ab(nd, :) = 0.0_8
  END DO
  DO turbi=nturbines,1,-1
    DO nd=1,nbdirs
      wakeeffcoeffb(nd) = wtvelocity(turbi)*wtvelocityb(nd, turbi)
      wtvelocityb(nd, turbi) = wakeeffcoeff*wtvelocityb(nd, turbi)
    END DO
    CALL POPREAL8(wakeeffcoeff)
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
          deltaxb(nd) = 0.0_8
        END DO
      ELSE
        DO nd=1,nbdirs
          tempb13(nd) = 2*a(turb)*wakeeffcoeffperzone*wakeeffcoeffb(nd)
          ab(nd, turb) = ab(nd, turb) + wakeeffcoeffperzone*tempb13(nd)
          wakeeffcoeffperzoneb(nd) = a(turb)*tempb13(nd)
        END DO
        deltax = turbinexw(turbi) - turbinexw(turb)
        DO nd=1,nbdirs
          deltaxb(nd) = 0.0_8
        END DO
        DO zone=3,1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            temp12 = 2.0_dp*mu(zone)
            temp9 = rotordiameter(turb) + temp12*kearray(turb)*deltax
            temp11 = wakeoverlaptrel_mat(turbi, turb, zone)
            temp10 = cosfac**2*rotordiameter(turb)**2
            DO nd=1,nbdirs
              tempb11(nd) = wakeeffcoeffperzoneb(nd)/temp9**2
              tempb12(nd) = -(temp10*temp11*2*tempb11(nd)/temp9)
              cosfacb(nd) = rotordiameter(turb)**2*temp11*2*cosfac*&
&               tempb11(nd)
              rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&               tempb12(nd) + cosfac**2*temp11*2*rotordiameter(turb)*&
&               tempb11(nd)
              wakeoverlaptrel_matb(nd, turbi, turb, zone) = &
&               wakeoverlaptrel_matb(nd, turbi, turb, zone) + temp10*&
&               tempb11(nd)
              kearrayb(nd, turb) = kearrayb(nd, turb) + temp12*deltax*&
&               tempb12(nd)
              deltaxb(nd) = deltaxb(nd) + temp12*kearray(turb)*tempb12(&
&               nd)
            END DO
          ELSE
            temp6 = rotordiameter(turb) + 2.0_dp*kearray(turb)*deltax*&
&             mmu(zone)
            temp8 = wakeoverlaptrel_mat(turbi, turb, zone)
            temp7 = cosfac**2*rotordiameter(turb)**2
            DO nd=1,nbdirs
              tempb8(nd) = wakeeffcoeffperzoneb(nd)/temp6**2
              tempb9(nd) = -(temp7*temp8*2*tempb8(nd)/temp6)
              tempb10(nd) = 2.0_dp*mmu(zone)*tempb9(nd)
              cosfacb(nd) = rotordiameter(turb)**2*temp8*2*cosfac*tempb8&
&               (nd)
              rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&               tempb9(nd) + cosfac**2*temp8*2*rotordiameter(turb)*&
&               tempb8(nd)
              wakeoverlaptrel_matb(nd, turbi, turb, zone) = &
&               wakeoverlaptrel_matb(nd, turbi, turb, zone) + temp7*&
&               tempb8(nd)
              kearrayb(nd, turb) = kearrayb(nd, turb) + deltax*tempb10(&
&               nd)
              deltaxb(nd) = deltaxb(nd) + kearray(turb)*tempb10(nd)
              mmub(nd, zone) = mmub(nd, zone) + 2.0_dp*kearray(turb)*&
&               deltax*tempb9(nd)
            END DO
          END IF
          rmax = cos_spread*0.5_dp*(wakediameterst_mat(turbi, turb, 3)+&
&           rotordiameter(turbi))
          CALL POPREAL8(cosfac)
          DO nd=1,nbdirs
            tempb7(nd) = -(pi*SIN(pi*(dabs1/rmax))*0.5_dp*cosfacb(nd)/&
&             rmax)
            dabs1b(nd) = tempb7(nd)
            rmaxb(nd) = -(dabs1*tempb7(nd)/rmax)
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(dabs1)
            DO nd=1,nbdirs
              wakecentersyt_matb(nd, turbi, turb) = wakecentersyt_matb(&
&               nd, turbi, turb) + dabs1b(nd)
              turbineywb(nd, turbi) = turbineywb(nd, turbi) - dabs1b(nd)
            END DO
          ELSE
            CALL POPREAL8(dabs1)
            DO nd=1,nbdirs
              turbineywb(nd, turbi) = turbineywb(nd, turbi) + dabs1b(nd)
              wakecentersyt_matb(nd, turbi, turb) = wakecentersyt_matb(&
&               nd, turbi, turb) - dabs1b(nd)
            END DO
          END IF
          DO nd=1,nbdirs
            tempb6(nd) = cos_spread*0.5_dp*rmaxb(nd)
            wakediameterst_matb(nd, turbi, turb, 3) = &
&             wakediameterst_matb(nd, turbi, turb, 3) + tempb6(nd)
            rotordiameterb(nd, turbi) = rotordiameterb(nd, turbi) + &
&             tempb6(nd)
          END DO
        END DO
        CALL POPINTEGER4(zone)
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8ARRAY(mmu, 3)
        temp5 = au*pi/180.0_dp + bu*yaw(turb)
        temp4 = COS(temp5)
        DO nd=1,nbdirs
          yawb(nd, turb) = yawb(nd, turb) - SIN(temp5)*bu*SUM(-(mu*mmub(&
&           nd, :)/temp4))/temp4
        END DO
      ELSE IF (branch .EQ. 1) THEN
        CALL POPREAL8ARRAY(mmu, 3)
      ELSE
        GOTO 100
      END IF
      DO nd=1,nbdirs
        mmub(nd, :) = 0.0_8
      END DO
 100  CALL POPREAL8(deltax)
      CALL POPREAL8(wakeeffcoeffperzone)
      DO nd=1,nbdirs
        turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + deltaxb(nd)
        turbinexwb(nd, turb) = turbinexwb(nd, turb) - deltaxb(nd)
      END DO
    END DO
    CALL POPREAL8(wakeeffcoeff)
  END DO
  DO nd=1,nbdirs
    adjusted_wind_speedb(nd, :) = 0.0_8
    adjusted_wind_speedb(nd, :) = wtvelocityb(nd, :)
  END DO
  DO nd=1,nbdirs
    keb(nd, :) = 0.0_8
  END DO
  DO turb=nturbines,1,-1
    DO nd=1,nbdirs
      keb(nd, turb) = keb(nd, turb) + (kecorrarray*s+1)*kearrayb(nd, &
&       turb)
      sb(nd) = ke(turb)*kecorrarray*kearrayb(nd, turb)
      kearrayb(nd, turb) = 0.0_8
      wakeoverlaptrel_matb(nd, turb, :, 1) = wakeoverlaptrel_matb(nd, &
&       turb, :, 1) + sb(nd)
      wakeoverlaptrel_matb(nd, turb, :, 2) = wakeoverlaptrel_matb(nd, &
&       turb, :, 2) + sb(nd)
    END DO
    CALL POPREAL8(s)
  END DO
  DO nd=1,nbdirs
    ctb(nd, :) = 0.0_8
    ctb(nd, :) = kecorrct*keb(nd, :)
  END DO
  CALL POPREAL8ARRAY(ke, nturbines)
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    DO nd=1,nbdirs
      a_inb(nd, :) = 0.0_8
      a_inb(nd, :) = ab(nd, :)
    END DO
  ELSE
    CALL CTTOAXIALIND_BV(ct, ctb, nturbines, a, ab, nbdirs)
    DO nd=1,nbdirs
      a_inb(nd, :) = 0.0_8
    END DO
  END IF
  CALL CALCOVERLAPAREAS_BV(nturbines, turbinexw, turbineyw, turbineywb, &
&                    turbinez, turbinezb, rotordiameter, rotordiameterb&
&                    , wakediameterst_mat, wakediameterst_matb, &
&                    wakecentersyt_mat, wakecentersyt_matb, &
&                    wakecenterszt_mat, wakecenterszt_matb, &
&                    wakeoverlaptrel_mat, wakeoverlaptrel_matb, nbdirs)
  DO nd=1,nbdirs
    keb(nd, :) = 0.0_8
  END DO
  DO turb=nturbines,1,-1
    DO nd=1,nbdirs
      wakediameter0b(nd) = 0.0_8
    END DO
    DO turbi=nturbines,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO nd=1,nbdirs
          deltaxb(nd) = 0.0_8
        END DO
      ELSE
        deltax = turbinexw(turbi) - turbinexw(turb)
        zone = 3
        DO nd=1,nbdirs
          tempb4(nd) = me(zone)*2.0_dp*wakediameterst_matb(nd, turbi, &
&           turb, zone)
          wakediameter0b(nd) = wakediameter0b(nd) + wakediameterst_matb(&
&           nd, turbi, turb, zone)
          wakediameterst_matb(nd, turbi, turb, zone) = 0.0_8
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
          wakediameterst_matb(nd, turbi, turb, zone) = 0.0_8
        END DO
        CALL POPINTEGER4(zone)
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          DO nd=1,nbdirs
            wakediameterst_matb(nd, turbi, turb, zone) = 0.0_8
          END DO
          DO nd=1,nbdirs
            xb(nd) = 0.0_8
            zerolocb(nd) = 0.0_8
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
            wakediameterst_matb(nd, turbi, turb, zone) = 0.0_8
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
            wakediameterst_matb(nd, turbi, turb, zone) = 0.0_8
          END DO
        END IF
        DO nd=1,nbdirs
          xb(nd) = 0.0_8
          zerolocb(nd) = 0.0_8
        END DO
      END IF
      temp3 = 2.0_dp*me(zone)*ke(turb)
      DO nd=1,nbdirs
        turbinexwb(nd, turb) = turbinexwb(nd, turb) + zerolocb(nd)
        wakediameter0b(nd) = wakediameter0b(nd) - zerolocb(nd)/temp3
        keb(nd, turb) = keb(nd, turb) + wakediameter0*2.0_dp*me(zone)*&
&         zerolocb(nd)/temp3**2
        turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + deltaxb(nd) + xb&
&         (nd)
        turbinexwb(nd, turb) = turbinexwb(nd, turb) - deltaxb(nd)
      END DO
      CALL POPINTEGER4(zone)
      CALL POPREAL8(deltax)
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(wakediameter0)
      DO nd=1,nbdirs
        rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&         wakediameter0b(nd)
      END DO
    ELSE
      CALL POPREAL8(wakediameter0)
      DO nd=1,nbdirs
        rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + COS(yaw(&
&         turb))*wakediameter0b(nd)
        yawb(nd, turb) = yawb(nd, turb) - rotordiameter(turb)*SIN(yaw(&
&         turb))*wakediameter0b(nd)
      END DO
    END IF
  END DO
  DO nd=1,nbdirs
    ctb(nd, :) = ctb(nd, :) + kecorrct*keb(nd, :)
  END DO
  DO turb=nturbines,1,-1
    DO nd=1,nbdirs
      wakeangleinitb(nd) = 0.0_8
    END DO
    DO turbi=nturbines,1,-1
      CALL POPCONTROL2B(branch)
      IF (branch .NE. 0) THEN
        IF (branch .EQ. 1) THEN
          DO nd=1,nbdirs
            deltaxb(nd) = 0.0_8
          END DO
        ELSE
          deltax = turbinexw(turbi) - turbinexw(turb)
          DO nd=1,nbdirs
            deltaxb(nd) = bd*wakecentersyt_matb(nd, turbi, turb)
          END DO
        END IF
        factor = 2.0_dp*kd*deltax/rotordiameter(turb) + 1.0_dp
        temp2 = 30.0_dp*kd*factor**5
        temp1 = wakeangleinit*rotordiameter(turb)
        temp0 = wakeangleinit**2 + 15.0_dp*factor**4
        DO nd=1,nbdirs
          displacementb(nd) = wakecentersyt_matb(nd, turbi, turb)
          tempb(nd) = -((wakeangleinit**2+15.0_dp)*displacementb(nd)/(&
&           30.0_dp*kd))
          tempb0(nd) = displacementb(nd)/temp2
          wakeangleinitb(nd) = wakeangleinitb(nd) + (temp0*rotordiameter&
&           (turb)+temp1*2*wakeangleinit)*tempb0(nd) + rotordiameter(&
&           turb)*tempb(nd) - wakeangleinit**2*rotordiameter(turb)*2*&
&           displacementb(nd)/(30.0_dp*kd)
          factorb(nd) = (15.0_dp*temp1*4*factor**3-30.0_dp*kd*temp0*&
&           temp1*5*factor**4/temp2)*tempb0(nd)
          turbineywb(nd, turb) = turbineywb(nd, turb) + &
&           wakecentersyt_matb(nd, turbi, turb)
          wakecentersyt_matb(nd, turbi, turb) = 0.0_8
          turbinezb(nd, turb) = turbinezb(nd, turb) + wakecenterszt_matb&
&           (nd, turbi, turb)
          wakecenterszt_matb(nd, turbi, turb) = 0.0_8
          tempb1(nd) = kd*2.0_dp*factorb(nd)/rotordiameter(turb)
          rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + temp0*&
&           wakeangleinit*tempb0(nd) - deltax*tempb1(nd)/rotordiameter(&
&           turb) + wakeangleinit*tempb(nd)
          deltaxb(nd) = deltaxb(nd) + tempb1(nd)
          turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + deltaxb(nd)
          turbinexwb(nd, turb) = turbinexwb(nd, turb) - deltaxb(nd)
        END DO
        CALL POPREAL8(deltax)
      END IF
    END DO
    CALL POPREAL8(wakeangleinit)
    DO nd=1,nbdirs
      yawb(nd, turb) = yawb(nd, turb) + ct(turb)*0.5_dp*COS(yaw(turb))*&
&       wakeangleinitb(nd)
      ctb(nd, turb) = ctb(nd, turb) + 0.5_dp*SIN(yaw(turb))*&
&       wakeangleinitb(nd)
    END DO
  END DO
  DO nd=1,nbdirs
    vinfb(nd) = 0.0_8
  END DO
  DO turb=nturbines,1,-1
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO nd=1,nbdirs
        adjusted_wind_speedb(nd, turb) = 0.0_8
      END DO
    ELSE
      temp = (turbinez(turb)-z0)/(zref-z0)
      DO nd=1,nbdirs
        vinfb(nd) = vinfb(nd) + temp**shearexp*adjusted_wind_speedb(nd, &
&         turb)
        IF (.NOT.(temp .LE. 0.0 .AND. (shearexp .EQ. 0.0 .OR. shearexp &
&           .NE. INT(shearexp)))) turbinezb(nd, turb) = turbinezb(nd, &
&           turb) + shearexp*temp**(shearexp-1)*vinf*&
&           adjusted_wind_speedb(nd, turb)/(zref-z0)
        adjusted_wind_speedb(nd, turb) = 0.0_8
      END DO
    END IF
  END DO
  DO nd=1,nbdirs
    yawdegb(nd, :) = 0.0_8
    yawdegb(nd, :) = pi*yawb(nd, :)/180.0_dp
  END DO
  DO nd=1,nbdirs
    wtvelocityb(nd, :) = 0.0_8
  END DO
END SUBROUTINE FLORIS_BV


!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.13 (r6666M) -  1 Mar 2018 15:30
!
!  Differentiation of cttoaxialind in reverse (adjoint) mode:
!   gradient     of useful results: axial_induction ct
!   with respect to varying inputs: ct
SUBROUTINE CTTOAXIALIND_BV(ct, ctb, nturbines, axial_induction, &
& axial_inductionb, nbdirs)
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  ! INTRINSIC KIND
! define precision to be the standard for a double precision! on local system
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
        axial_inductionb(nd, i) = 0.0_8
      END DO
    ELSE
      DO nd=1,nbdirs
        IF (.NOT.0.0203_dp - 0.6427_dp*(0.889_dp-ct(i)) .EQ. 0.0) ctb(nd&
&         , i) = ctb(nd, i) + 0.6427_dp*axial_inductionb(nd, i)/(2.0*&
&           SQRT(0.0203_dp-0.6427_dp*(0.889_dp-ct(i))))
        axial_inductionb(nd, i) = 0.0_8
      END DO
    END IF
  END DO
END SUBROUTINE CTTOAXIALIND_BV



!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.13 (r6666M) -  1 Mar 2018 15:30
!
!  Differentiation of calcoverlapareas in reverse (adjoint) mode:
!   gradient     of useful results: rotordiameter turbiney wakeoverlaptrel_mat
!                wakediameters wakecentersyt
!   with respect to varying inputs: rotordiameter turbiney turbinez
!                wakediameters wakecentersyt wakecenterszt
SUBROUTINE CALCOVERLAPAREAS_BV(nturbines, turbinex, turbiney, turbineyb&
& , turbinez, turbinezb, rotordiameter, rotordiameterb, wakediameters, &
& wakediametersb, wakecentersyt, wakecentersytb, wakecenterszt, &
& wakecentersztb, wakeoverlaptrel_mat, wakeoverlaptrel_matb, nbdirs)
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! do turbI = 1, nTurbines
!         do turb = 1, nTurbines
!             do zone = 1, 3
!                 print *, "wakeOverlapTRel_mat[", turbI, ", ", turb, ", ", zone, "] = ", wakeOverlapTRel_mat(turbI, turb, zone)
!             end do
!         end do
!     end do
  ! INTRINSIC KIND
! define precision to be the standard for a double precision! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: turbinex, turbiney, &
& turbinez, rotordiameter
  REAL(dp), DIMENSION(nbdirs, nturbines) :: turbineyb, turbinezb, &
& rotordiameterb
  REAL(dp), DIMENSION(nturbines, nturbines, 3), INTENT(IN) :: &
& wakediameters
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakediametersb
  REAL(dp), DIMENSION(nturbines, nturbines), INTENT(IN) :: wakecentersyt&
& , wakecenterszt
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines) :: wakecentersytb&
& , wakecentersztb
! out
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlaptrel_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlaptrel_matb
! local
  INTEGER :: turb, turbi, zone
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp, tol=0.000001_dp
  REAL(dp) :: ovdyd, ovr, ovrr, ovl, ovz
  REAL(dp), DIMENSION(nbdirs) :: ovdydb, ovrb, ovrrb, ovlb, ovzb
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlap
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlapb
  INTRINSIC SQRT
  INTRINSIC ABS
  INTRINSIC DACOS
  REAL(dp) :: temp
  REAL(dp) :: temp0
  INTEGER :: nd
  REAL(dp), DIMENSION(nbdirs) :: tempb
  REAL(dp), DIMENSION(nbdirs) :: tempb0
  REAL(dp), DIMENSION(nbdirs) :: tempb1
  REAL(dp), DIMENSION(nbdirs) :: tempb2
  REAL(dp) :: temp1
  REAL(dp) :: temp2
  REAL(dp), DIMENSION(nbdirs) :: tempb3
  REAL(dp), DIMENSION(nbdirs) :: tempb4
  REAL(dp) :: temp3
  REAL(dp), DIMENSION(nbdirs, nturbines, 3) :: tempb5
  INTEGER :: branch
  INTEGER :: nbdirs
  wakeoverlap = 0.0_dp
  DO turb=1,nturbines
    DO turbi=1,nturbines
      IF (turbinex(turbi) .GT. turbinex(turb)) THEN
        CALL PUSHREAL8(ovdyd)
        ovdyd = SQRT((wakecentersyt(turbi, turb)-turbiney(turbi))**2 + (&
&         wakecenterszt(turbi, turb)-turbinez(turbi))**2)
! distance between wake center and rotor center
        CALL PUSHREAL8(ovr)
        ovr = rotordiameter(turbi)/2
! rotor diameter
        DO zone=1,3
          ovrr = wakediameters(turbi, turb, zone)/2.0_dp
          IF (ovdyd .GE. 0.) THEN
            CALL PUSHREAL8(ovdyd)
            ovdyd = ovdyd
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREAL8(ovdyd)
            ovdyd = -ovdyd
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ovdyd .GE. 0.0_dp + tol) THEN
! calculate the distance from the wake center to the vertical line between
! the two circle intersection points
            CALL PUSHREAL8(ovl)
            ovl = (-(ovr*ovr)+ovrr*ovrr+ovdyd*ovdyd)/(2.0_dp*ovdyd)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREAL8(ovl)
            ovl = 0.0_dp
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8(ovz)
          ovz = ovrr*ovrr - ovl*ovl
! Finish calculating the distance from the intersection line to the outer edge of the wake zone
          IF (ovz .GT. 0.0_dp + tol) THEN
            CALL PUSHREAL8(ovz)
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
  wakeoverlaptrel_mat = wakeoverlap
  DO turbi=nturbines,1,-1
    temp3 = pi*rotordiameter(turbi)**2
    DO nd=1,nbdirs
      tempb5(nd, :, :) = 4.0_dp*wakeoverlaptrel_matb(nd, turbi, :, :)/&
&       temp3
      rotordiameterb(nd, turbi) = rotordiameterb(nd, turbi) + pi*2*&
&       rotordiameter(turbi)*SUM(-(wakeoverlaptrel_mat(turbi, :, :)*&
&       tempb5(nd, :, :)/temp3))
      wakeoverlaptrel_matb(nd, turbi, :, :) = tempb5(nd, :, :)
    END DO
  END DO
  DO nd=1,nbdirs
    wakeoverlapb(nd, :, :, :) = 0.0_8
    wakeoverlapb(nd, :, :, :) = wakeoverlaptrel_matb(nd, :, :, :)
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
    turbinezb(nd, :) = 0.0_8
    wakecentersztb(nd, :, :) = 0.0_8
  END DO
  DO turb=nturbines,1,-1
    DO turbi=nturbines,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        DO nd=1,nbdirs
          ovdydb(nd) = 0.0_8
          ovrb(nd) = 0.0_8
        END DO
        DO zone=3,1,-1
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              DO nd=1,nbdirs
                wakeoverlapb(nd, turbi, turb, zone) = 0.0_8
              END DO
              ovrr = wakediameters(turbi, turb, zone)/2.0_dp
              DO nd=1,nbdirs
                ovlb(nd) = 0.0_8
                ovrrb(nd) = 0.0_8
                ovzb(nd) = 0.0_8
              END DO
              GOTO 100
            ELSE
              ovrr = wakediameters(turbi, turb, zone)/2.0_dp
              DO nd=1,nbdirs
                ovrrb(nd) = pi*2*ovrr*wakeoverlapb(nd, turbi, turb, zone&
&                 )
                wakeoverlapb(nd, turbi, turb, zone) = 0.0_8
              END DO
            END IF
          ELSE IF (branch .EQ. 2) THEN
            DO nd=1,nbdirs
              ovrb(nd) = ovrb(nd) + pi*2*ovr*wakeoverlapb(nd, turbi, &
&               turb, zone)
              wakeoverlapb(nd, turbi, turb, zone) = 0.0_8
            END DO
            ovrr = wakediameters(turbi, turb, zone)/2.0_dp
            DO nd=1,nbdirs
              ovrrb(nd) = 0.0_8
            END DO
          ELSE
            ovrr = wakediameters(turbi, turb, zone)/2.0_dp
            temp1 = ovl/ovrr
            temp2 = (ovdyd-ovl)/ovr
            DO nd=1,nbdirs
              IF (temp1 .EQ. 1.0 .OR. temp1 .EQ. (-1.0)) THEN
                tempb3(nd) = 0.0
              ELSE
                tempb3(nd) = -(ovrr*wakeoverlapb(nd, turbi, turb, zone)/&
&                 SQRT(1.D0-temp1**2))
              END IF
              IF (temp2 .EQ. 1.0 .OR. temp2 .EQ. (-1.0)) THEN
                tempb4(nd) = 0.0
              ELSE
                tempb4(nd) = -(ovr*wakeoverlapb(nd, turbi, turb, zone)/&
&                 SQRT(1.D0-temp2**2))
              END IF
              ovrrb(nd) = DACOS(temp1)*2*ovrr*wakeoverlapb(nd, turbi, &
&               turb, zone) - temp1*tempb3(nd)
              ovlb(nd) = tempb3(nd) - tempb4(nd)
              ovrb(nd) = ovrb(nd) + DACOS(temp2)*2*ovr*wakeoverlapb(nd, &
&               turbi, turb, zone) - temp2*tempb4(nd)
              ovdydb(nd) = ovdydb(nd) + tempb4(nd) - ovz*wakeoverlapb(nd&
&               , turbi, turb, zone)
              ovzb(nd) = -(ovdyd*wakeoverlapb(nd, turbi, turb, zone))
              wakeoverlapb(nd, turbi, turb, zone) = 0.0_8
            END DO
            GOTO 100
          END IF
          DO nd=1,nbdirs
            ovlb(nd) = 0.0_8
            ovzb(nd) = 0.0_8
          END DO
 100      CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(ovz)
            DO nd=1,nbdirs
              IF (ovz .EQ. 0.0) THEN
                ovzb(nd) = 0.0
              ELSE
                ovzb(nd) = ovzb(nd)/(2.0*SQRT(ovz))
              END IF
            END DO
          ELSE
            DO nd=1,nbdirs
              ovzb(nd) = 0.0_8
            END DO
          END IF
          CALL POPREAL8(ovz)
          DO nd=1,nbdirs
            ovrrb(nd) = ovrrb(nd) + 2*ovrr*ovzb(nd)
            ovlb(nd) = ovlb(nd) - 2*ovl*ovzb(nd)
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(ovl)
            DO nd=1,nbdirs
              tempb2(nd) = ovlb(nd)/(2.0_dp*ovdyd)
              ovrrb(nd) = ovrrb(nd) + 2*ovrr*tempb2(nd)
              ovrb(nd) = ovrb(nd) - 2*ovr*tempb2(nd)
              ovdydb(nd) = ovdydb(nd) + (2*ovdyd-(ovrr**2-ovr**2+ovdyd**&
&               2)/ovdyd)*tempb2(nd)
            END DO
          ELSE
            CALL POPREAL8(ovl)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(ovdyd)
          ELSE
            CALL POPREAL8(ovdyd)
            DO nd=1,nbdirs
              ovdydb(nd) = -ovdydb(nd)
            END DO
          END IF
          DO nd=1,nbdirs
            wakediametersb(nd, turbi, turb, zone) = wakediametersb(nd, &
&             turbi, turb, zone) + ovrrb(nd)/2.0_dp
          END DO
        END DO
        CALL POPREAL8(ovr)
        CALL POPREAL8(ovdyd)
        temp0 = wakecenterszt(turbi, turb) - turbinez(turbi)
        temp = wakecentersyt(turbi, turb) - turbiney(turbi)
        DO nd=1,nbdirs
          rotordiameterb(nd, turbi) = rotordiameterb(nd, turbi) + ovrb(&
&           nd)/2
          IF (temp**2 + temp0**2 .EQ. 0.0) THEN
            tempb(nd) = 0.0
          ELSE
            tempb(nd) = ovdydb(nd)/(2.0*SQRT(temp**2+temp0**2))
          END IF
          tempb0(nd) = 2*temp*tempb(nd)
          tempb1(nd) = 2*temp0*tempb(nd)
          wakecentersytb(nd, turbi, turb) = wakecentersytb(nd, turbi, &
&           turb) + tempb0(nd)
          turbineyb(nd, turbi) = turbineyb(nd, turbi) - tempb0(nd)
          wakecentersztb(nd, turbi, turb) = wakecentersztb(nd, turbi, &
&           turb) + tempb1(nd)
          turbinezb(nd, turbi) = turbinezb(nd, turbi) - tempb1(nd)
        END DO
      END IF
    END DO
  END DO
END SUBROUTINE CALCOVERLAPAREAS_BV



!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.13 (r6666M) -  1 Mar 2018 15:30
!
!  Differentiation of hermite_spline in reverse (adjoint) mode:
!   gradient     of useful results: y
!   with respect to varying inputs: x x0 x1 dy0 y0
! Flow field calculations have been intentionally left out to save development time.
! The flow field can be calculated using the pure python version of floris
! This implementation is fully smooth and differentiable with the exception of a
! discontinuity at the hub of each turbine. The discontinuity only presents issues if
! turbines are place within 1E-15 * rotor diameter of one another, which is extremely
! unlikely during optimization if the user does not explicitly place them there.
SUBROUTINE HERMITE_SPLINE_BV(x, xb, x0, x0b, x1, x1b, y0, y0b, dy0, dy0b&
& , y1, dy1, y, yb, nbdirs)
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
!dy_dx = c3*3*x**2 + c2*2*x + c1
  ! INTRINSIC KIND
! define precision to be the standard for a double precision! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: x, x0, x1, y0, dy0, y1, dy1
  REAL(dp), DIMENSION(nbdirs) :: xb, x0b, x1b, y0b, dy0b
! out
  REAL(dp) :: y
  REAL(dp), DIMENSION(nbdirs) :: yb
!, dy_dx
! local
  REAL(dp) :: c3, c2, c1, c0
  REAL(dp), DIMENSION(nbdirs) :: c3b, c2b, c1b, c0b
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
  REAL(dp), DIMENSION(nbdirs) :: tempb
  REAL(dp), DIMENSION(nbdirs) :: tempb0
  REAL(dp), DIMENSION(nbdirs) :: tempb1
  REAL(dp), DIMENSION(nbdirs) :: tempb2
  REAL(dp), DIMENSION(nbdirs) :: tempb3
  REAL(dp), DIMENSION(nbdirs) :: tempb4
  REAL(dp), DIMENSION(nbdirs) :: tempb5
  REAL(dp), DIMENSION(nbdirs) :: tempb6
  REAL(dp), DIMENSION(nbdirs) :: tempb7
  REAL(dp), DIMENSION(nbdirs) :: tempb8
  REAL(dp), DIMENSION(nbdirs) :: tempb9
  REAL(dp), DIMENSION(nbdirs) :: tempb10
  REAL(dp), DIMENSION(nbdirs) :: tempb11
  REAL(dp), DIMENSION(nbdirs) :: tempb12
  REAL(dp), DIMENSION(nbdirs) :: tempb13
  REAL(dp), DIMENSION(nbdirs) :: tempb14
  REAL(dp), DIMENSION(nbdirs) :: tempb15
  REAL(dp), DIMENSION(nbdirs) :: tempb16
  REAL(dp), DIMENSION(nbdirs) :: tempb17
  REAL(dp), DIMENSION(nbdirs) :: tempb18
  REAL(dp), DIMENSION(nbdirs) :: tempb19
  REAL(dp), DIMENSION(nbdirs) :: tempb20
  REAL(dp), DIMENSION(nbdirs) :: tempb21
  REAL(dp), DIMENSION(nbdirs) :: tempb22
  REAL(dp), DIMENSION(nbdirs) :: tempb23
  REAL(dp), DIMENSION(nbdirs) :: tempb24
  REAL(dp), DIMENSION(nbdirs) :: tempb25
  REAL(dp), DIMENSION(nbdirs) :: tempb26
  REAL(dp), DIMENSION(nbdirs) :: tempb27
  REAL(dp), DIMENSION(nbdirs) :: tempb28
  REAL(dp), DIMENSION(nbdirs) :: tempb29
  REAL(dp), DIMENSION(nbdirs) :: tempb30
  INTEGER :: nbdirs
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
