
MODULE epstoluene_module
    IMPLICIT NONE
    CONTAINS

    FUNCTION calculate_epstoluene(wp, i) RESULT(epstoluene)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: wp
        INTEGER, INTENT(IN) :: i
        DOUBLE PRECISION :: epstoluene
        DOUBLE PRECISION, DIMENSION(8) :: lc, lw

        lc = (/ 5.61D-3, 6.97D-2, 8.07D-3, 5.15D-1, 5.74D-1, 9.91D-3, 1.18D-1, 1.01D-2 /)
        lw = (/ 3.40D-2, 9.97D-2, 1.15D0, 7.23D0, 1.50D1, 2.08D1, 2.59D1, 7.65D1 /) * (1.602177D-19/1.0545D-34)

        IF (i > 0) THEN
            epstoluene = 1.D0 + SUM(lc / (1.D0 + (wp / lw)**2.D0))
        ELSE
            epstoluene = 2.395D0
        END IF
    END FUNCTION calculate_epstoluene

    END MODULE epstoluene_module


MODULE epsmagnetite_module
    IMPLICIT NONE
    CONTAINS

    FUNCTION calculate_epsmagnetite(wp, i) RESULT(epsmagnetite)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: wp
        INTEGER, INTENT(IN) :: i
        DOUBLE PRECISION :: epsmagnetite
        DOUBLE PRECISION, DIMENSION(9) :: mc, mw

        mc = (/ 49.338D0, 7.7867D0, 4.7422D0, 3.6331D0, 1.2529D0, 2.2625D0, 1.1828D0, 0.3136D0, 0.0672D0 /)
        mw = (/ 0.0209D0, 0.0772D0, 0.3127D0, 0.6912D0, 2.4508D0, 5.3771D0, 13.0017D0, 30.7646D0, 72.2563D0 /) &
        * (1.602177D-19 / 1.0545D-34)
        IF (i > 0) THEN
            epsmagnetite = 1.D0 + SUM(mc / (1.D0 + (wp / mw)**2.D0))
        ELSE
            epsmagnetite = 54.8D0
        END IF
        END FUNCTION calculate_epsmagnetite

END MODULE epsmagnetite_module



MODULE epsteflon_module
    IMPLICIT NONE
    CONTAINS

    FUNCTION calculate_epsteflon(wp, i) RESULT(epsteflon)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: wp
        INTEGER, INTENT(IN) :: i
        DOUBLE PRECISION :: epsteflon
        DOUBLE PRECISION, DIMENSION(8) :: tc, tw

        tc = (/ 9.30D-3, 1.83D-2, 1.39D-1, 1.12D-1, 1.95D-1, 4.38D-1, 1.06D-1, 3.86D-2 /)
        tw = (/ 3.00D-4, 7.60D-3, 5.57D-2, 1.26D-1, 6.71D0, 1.86D1, 4.21D1, 7.76D1 /) * (1.602177D-19 / 1.0545D-34)

        IF (i > 0) THEN
            epsteflon = 1.D0 + SUM(tc / (1.D0 + (wp / tw)**2.D0))
        ELSE
            epsteflon = 2.1D0
        END IF
    END FUNCTION calculate_epsteflon

END MODULE epsteflon_module


MODULE epspolystyrene_module
    IMPLICIT NONE
    CONTAINS

    FUNCTION calculate_epspolystyrene(wp, i) RESULT(epspolystyrene)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: wp
        INTEGER, INTENT(IN) :: i
        DOUBLE PRECISION :: epspolystyrene
        DOUBLE PRECISION, DIMENSION(8) :: pc, pw

        pc = (/ 1.21D-2, 2.19D-2, 1.79D-2, 3.06D-2, 3.03D-1, 6.23D-1, 3.25D-1, 3.31D-2 /)
        pw = (/ 1.00D-3, 1.32D-2, 3.88D0, 1.31D-1, 5.99D0, 1.02D1, 1.88D1, 5.15D1 /) * (1.602177D-19 / 1.0545D-34)

        IF (i > 0) THEN
            epspolystyrene = 1.D0 + SUM(pc / (1.D0 + (wp / pw)**2.D0))
        ELSE
            epspolystyrene = 2.4D0
        END IF
    END FUNCTION calculate_epspolystyrene

END MODULE epspolystyrene_module





PROGRAM main
  USE epstoluene_module
  USE epsmagnetite_module
  USE epsteflon_module
  USE epspolystyrene_module
  IMPLICIT NONE
  DOUBLE PRECISION :: wp, epstoluene, epsmagnetite, epsteflon, epspolystyrene, eps_mix, G3  !< dilectric functions 
  DOUBLE PRECISION :: dist   !< distance between the particles
  DOUBLE PRECISION :: bsum   !< energy sum
  DOUBLE PRECISION :: qsum, qmax, qmin, dq, q, qbweight, qweight  !< q intrgration variables
  INTEGER :: kmax  !< max k points
  DOUBLE PRECISION :: phi  !< mixing parameter
  DOUBLE PRECISION :: mu2, mu3, mu4 !< magnetic susceptibility parameters
  DOUBLE PRECISION :: diam !< diameter of the nano-particle
  DOUBLE PRECISION :: gamma2, gamma3, gamma4, D234, D232, D134, D132, film1, film2, Gdisp1, Gdisp2, s  !< variables for calculation




INTEGER :: i,dloop,k
DOUBLE PRECISION :: Kb, T, h, pi, c     !< constants

Kb = 1.380649D-23  !< Boltzmann constant in J/K
T = 300.D0         !< Temperature in Kelvin
h = 1.054571800D-34 !< Planck constant in J*s
pi = 3.141592653589793D0 !< pi
c=2.99792458D8     !<velocity of light 


diam=20D0
mu2 = 1.0D0
mu4 = 1.0D0
kmax=8000

OPEN(UNIT=10, FILE='results_20nm_phi0.01.csv', STATUS='REPLACE', ACTION='WRITE')   !< open file for writing
OPEN(UNIT=20, FILE='free_energy_20nm_phi0.01.csv', STATUS='REPLACE', ACTION='WRITE')   !< open file for writing
! WRITE(10, '(A)') 'i,wp,epstoluene,epsmagnetite'
dist=0D-10
DO dloop=1,1500
    dist=dist+10D-10
    bsum=0D0
  DO i = 0, 2000
    wp = (2.D0 * pi * i * Kb * T / h)

    epstoluene = calculate_epstoluene(wp,i)
    epsmagnetite = calculate_epsmagnetite(wp,i)
    epsteflon = calculate_epsteflon(wp,i)
    epspolystyrene = calculate_epspolystyrene(wp,i)
    ! PRINT *, 'For i =', i, 'wp =', wp, 'epstoluene =', epstoluene, 'epsmagnetite =', epsmagnetite, 'epsteflon =', epsteflon, &
    ! 'epspolystyrene =', epspolystyrene
    
     phi = 0.01D0
     mu3=1.D0+(phi/0.05D0)*(1.24D0-1.D0)*((diam/10.D0)**3.D0)  
     
     G3=(phi*(epsmagnetite-1.D0)/(epsmagnetite+2))+(1-phi)*(epstoluene-1.D0)/(epstoluene+2)
     
     eps_mix=(1+2.D0*G3)/(1-G3)
    WRITE(10, '(ES15.8E3,1X,ES15.8E3,1X,ES15.8E3,1X,ES15.8E3)') wp, epsteflon, eps_mix, epspolystyrene
     qsum=0D0
     qmax=(25D0)/(dist)
     qmin=0     
     dq= (qmax-qmin)/DFLOAT(kmax)
    DO k = 1, kmax
        q = qmin + k * dq

        ! weight in Simpson's formula
        IF (k == 0 .OR. k == kmax) THEN
            qweight = 1.D0
            qbweight = 1.D0
        ELSE IF (MOD(k, 2) /= 0) THEN
            qweight = 4.D0
            qbweight = 0.D0
        ELSE
            qweight = 2.D0
            IF (MOD(k, 4) /= 0) THEN
                qbweight = 4.D0
            ELSE
                qbweight = 2.D0
            END IF
        END IF
        
        gamma2=DSQRT(q*q+ epspolystyrene*wp*wp/(c*c))
        gamma3=DSQRT(q*q+ eps_mix*wp*wp/(c*c))
        gamma4=DSQRT(q*q+ epsteflon*wp*wp/(c*c))

        D234=-(epsteflon*gamma3-eps_mix*gamma4)/(epsteflon*gamma3+eps_mix*gamma4) !the 2 indicated that it is TM refl coeff
        D232=-(epspolystyrene*gamma3-eps_mix*gamma2)/(epspolystyrene*gamma3+eps_mix*gamma2)
        
        D134=-(gamma3-gamma4)/(gamma3+gamma4)			!!the 1 indicated that it is TE refl coeff

        D132=-(gamma3-gamma2)/(gamma3+gamma2)
        IF (i==0) THEN
            D134=-(mu4*gamma3-mu3*gamma4)/(mu4*gamma3+mu3*gamma4)			
            D132=-(mu2*gamma3-mu3*gamma2)/(mu2*gamma3+mu3*gamma2)
        END IF

        
        film1=D234*D232         !product of TM reflection coefficient
        film2=D134*D132	 !product of TE reflection coefficient
        
        Gdisp1=1.D0-DEXP(-2.D0*dist*gamma3)*film1	!retardation TM thin films 
        Gdisp2=1.D0-DEXP(-2.D0*dist*gamma3)*film2	!retardation TE thin films

        s=(Kb*T/(2.D0*pi))*(LOG(Gdisp1)+LOG(Gdisp2))  
        IF (i==0) THEN
            s=s/2.D0
        END IF

      qsum=qsum+s*dq*q*qweight/3.D0 

    END DO
     
    bsum=bsum+qsum
    
  END DO
  print*, dist*1D9,' , ', bsum 
  WRITE(20,*) dist*1D9,' ,  ',2.D0*pi*bsum*1D6
CLOSE(UNIT=10)


END DO

CLOSE(UNIT=20)

END PROGRAM main

