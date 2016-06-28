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
        
    ! define precision to be the standard for a double precision ! on local system
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
                OVdYd = sqrt((wakeCentersYT(turbI, turb)-turbineY(turbI))**2+(wakeCentersZT(turbI, turb)-turbineZ(turbI))**2)     ! distance between wake center and rotor center
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
    

subroutine floris(nTurbines, nSamples, turbineXw, turbineYw, turbineZ, yawDeg, &
                          & rotorDiameter, Vinf, Ct, a_in, ke_in, kd, me, &
                          & initialWakeDisplacement, bd, MU, aU, bU, initialWakeAngle, &
                          & cos_spread, keCorrCT, Region2CT, keCorrArray, useWakeAngle, &
                          & adjustInitialWakeDiamToYaw, axialIndProvided, useaUbU, &
                          & wsPositionXYZw, shearCoefficientAlpha, shearZh, &
                          & wtVelocity, wsArray, &
                          & wakeCentersYT_vec, wakeCentersZT_vec, wakeDiametersT_vec, wakeOverlapTRel_vec) 
    
    ! independent variables: yawDeg Ct turbineXw turbineYw turbineZ  rotorDiameter a_in    
    ! dependent variables: wtVelocity
    
    implicit none
    
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)
    
    ! in
    integer, intent(in) :: nTurbines
    integer :: nSamples
    real(dp), intent(in) :: kd, initialWakeDisplacement, initialWakeAngle, ke_in
    real(dp), intent(in) :: keCorrCT, Region2CT, bd, cos_spread, keCorrArray !TODO move Vinf to len(nTurbines)
    real(dp), dimension(nTurbines), intent(in) :: yawDeg, Ct, a_in, turbineXw, turbineYw, turbineZ, Vinf
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
    real(dp) :: spline_bound ! in rotor diameters    
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
    ! visualization
    real(dp) :: deltay, deltaz, radiusLoc, axialIndAndNearRotor, reductionFactor
    
    ! model out
    real(dp), dimension(nTurbines), intent(out) :: wtVelocity
    
    ! test out
    real(dp), dimension(nTurbines*nTurbines), intent(out) :: wakeCentersYT_vec, wakeCentersZT_vec 
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
       
    ! calculate locations of wake centers in wind ref. frame
    wakeCentersYT_mat = 0.0_dp
    wakeCentersZT_mat = 0.0_dp
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
    call calcOverlapAreas(nTurbines, turbineXw, turbineYw, turbineZ, rotorDiameter, &
                          & wakeDiametersT_mat, wakeCentersYT_mat, wakeCentersZT_mat, wakeOverlapTRel_mat) 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Velocity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    ! initialize velocities in full flow field (optional)
    ! wsArray = Vinf TODO MAY NEED TO FIX THIS?
    
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
                mmU = MU/cos(aU*pi/180.0_dp + bU*yaw(turb))
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
    
        if (useaUbU) then
            mmU = MU/cos(aU*pi/180.0_dp + bU*yaw(turb))
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


!TODO TAPENADE: Run TAPENADE (floris.f90 as source, adStack.c and adBuffer.f as an include) input language "from the files extensions), Name of the top routine: floris. Differentiate in Multiobjective Adjoint Mode.

!TODO After you run, 1.) Change all nbdirs max (all one word) to nbdirs 2.) delete  3.) in florisbv, chnage all the input variables b to intent out 4.) delete wtvelocity from florisbv




