PROGRAM Casimir

!skip definitions and go down to the program!

!some of these parameters are no longer in use		
DOUBLE PRECISION :: x1,x2,x,G 	
DOUBLE PRECISION :: vikt,viktb,sum,sumb,sum3W,sum3IN
DOUBLE PRECISION :: dq,q,D1,D2,D3,D2b,D3b
INTEGER :: L,k,dloop,ii,jj,jjj

INTEGER, PARAMETER :: kmax=8000 !steps in q integration vary to test accuracy 
	
DOUBLE PRECISION :: bsum,bsum2,bsum3,qmax
DOUBLE PRECISION :: qsum,qsum2,qsum3,D,qweight,qbweight
DOUBLE PRECISION :: weight,bweight
DOUBLE PRECISION :: qe, pi 
DOUBLE PRECISION :: eps1,eps2,eps3,eps4
DOUBLE PRECISION :: gamma1,gamma2,gamma3,gamma4 
DOUBLE PRECISION :: D234,D241,D232,D134,D141,D132
DOUBLE PRECISION :: film1,film2,bdist 
DOUBLE PRECISION :: Kb,T,qmin,c,s,h,Ng
DOUBLE PRECISION :: dist,D21,D43,D12,D20,D10,Gdisp,Gdisp1,Gdisp2
DOUBLE PRECISION :: wp1,wp2,wp3,DTM12,DTM23,DTE12,DTE23	
DOUBLE PRECISION :: B1,B2, wIR,wUV,CIR,CUV

!parameters for dielectric function of silica
DOUBLE PRECISION :: eps0,epsINF,Wlo,wto,epsII,GGG1, GGG2, GGG3,SiO2epsp3

!ignore these parameters or check which are not used and remove those
!DOUBLE PRECISION :: e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11
DOUBLE PRECISION :: f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11
DOUBLE PRECISION :: g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11
DOUBLE PRECISION :: ee1,ee2,ee3,ee4,ee5,ee6,ee7,ee8,ee9
DOUBLE PRECISION :: ff1,ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9
DOUBLE PRECISION :: gg1,gg2,gg3,gg4,gg5,gg6,gg7,gg8,gg9
DOUBLE PRECISION :: W1,W2,W3,W4,W5,A1,A2,A3,A4,A5
!DOUBLE PRECISION :: cc1,cc2,tau1,tau2,temp,omega
!DOUBLE PRECISION :: ome1,ome2,ome3,ome4,ome5,ome6,ome7,ome8
!DOUBLE PRECISION :: ome9,ome10,ome11,ome12,ome13,ome14,ome15,ome16
!DOUBLE PRECISION :: c1,c2,c3,c4,c5,c6,c7,c8,c9
!DOUBLE PRECISION :: c10,c11,c12,c13,c14,c15,c16
!DOUBLE PRECISION :: amma1,amma2,amma3,amma4,amma5 
!DOUBLE PRECISION :: amma6,amma7,amma8
!DOUBLE PRECISION :: amma9,amma10,amma11,amma12 
!DOUBLE PRECISION :: amma13,amma14,amma15,amma16
DOUBLE PRECISION :: epswater,epsCelluloseIW
DOUBLE PRECISION ::  alphaModal_1, freqModal_1, alphaModal_2, freqModal_2, alphaModal_3, freqModal_3
DOUBLE PRECISION ::  alphaModal_4, freqModal_4, alphaModal_5, freqModal_5, alphaModal_6, freqModal_6
DOUBLE PRECISION ::  alphaModal_7, freqModal_7, alphaModal_8, freqModal_8, alphaModal_9, freqModal_9
DOUBLE PRECISION ::  alphaModal_10, freqModal_10, alphaModal_11, freqModal_11, alphaModal_12, freqModal_12 
DOUBLE PRECISION ::  alphaModal_13, freqModal_13, alphaModal_14, freqModal_14, alphaModal_15, freqModal_15
DOUBLE PRECISION ::  alphaModal_16, freqModal_16, alphaModal_17, freqModal_17, alphaModal_18, freqModal_18
DOUBLE PRECISION ::  alphaModal_19, freqModal_19

DOUBLE PRECISION :: lc1,lc2,lc3,lc4,lc5,lc6,lc7,lc8,lw1,lw2,lw3,lw4,lw5,lw6,lw7,lw8

integer :: i



!Mathias Edits
DOUBLE PRECISION :: cc1,cc2,tau1,tau2,temp,omega
DOUBLE PRECISION :: ome1,ome2,ome3,ome4,ome5,ome6,ome7,ome8
DOUBLE PRECISION :: ome9,ome10,ome11,ome12,ome13,ome14,ome15,ome16
DOUBLE PRECISION :: ome17,ome18,ome19
DOUBLE PRECISION :: c1,c2,c3,c4,c5,c6,c7,c8,c9
DOUBLE PRECISION :: c10,c11,c12,c13,c14,c15,c16
DOUBLE PRECISION :: c17,c18,c19
DOUBLE PRECISION :: amma1,amma2,amma3,amma4,amma5 
DOUBLE PRECISION :: amma6,amma7,amma8
DOUBLE PRECISION :: amma9,amma10,amma11,amma12 
DOUBLE PRECISION :: amma13,amma14,amma15,amma16
DOUBLE PRECISION :: amma17,amma18,amma19

DOUBLE PRECISION :: pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8
DOUBLE PRECISION :: pw1,pw2,pw3,pw4,pw5,pw6,pw7,pw8

DOUBLE PRECISION :: tc1,tc2,tc3,tc4,tc5,tc6,tc7,tc8
DOUBLE PRECISION :: tw1,tw2,tw3,tw4,tw5,tw6,tw7,tw8

DOUBLE PRECISION :: mc1,mc2,mc3,mc4,mc5,mc6,mc7,mc8,mc9
DOUBLE PRECISION :: mw1,mw2,mw3,mw4,mw5,mw6,mw7,mw8,mw9
DOUBLE PRECISION :: epsmagnetite,epstoluene,mu2,mu3,mu4


 
pi=3.14159265359D0
h=1.054571800D-34  !as noted by the numbers this is what is usually called "hbar" 
Kb=1.38064852D-23  !Boltzmann constant
c=2.99792458D8     !velocity of light 
T=300D0  !293.15D0  !273.55D0          !temperature in K
qe=1.60217662D-19 !unit charge in SI units



!phi=0.01455D0

!phi=0.01455 !0%  
!phi=0.02072 !10% 
!phi=0.02702 !20%  
!phi=0.03320 !30%  
!phi=0.03630 !35% 


dist=0D-10    !thickness water layer


DO dloop=1,1500

dist=dist+10D-10  !*sqrt(10.D0)  !thickness condensed water layer

bsum=0.D0
bsum2=0.D0
bsum3=0.D0

DO i=0,2000  !2000  !3998      !Matsubara freq 	

wp=(2.D0*pi*(i)*Kb*T/h)	






 
!PHYSICAL REVIEW A 81, 062502 (2010)
!Repulsive Casimir forces between solid materials with high-refractive-index intervening liquids
!P. J. van Zwol and G. Palasantzas
!Toluene

!5.61[−3] 6.97[−2] 8.07[−3] 5.15[−1] 5.74[−1] 9.91[−3] 1.18[−1] 1.01[−2]

IF (i>0) THEN
lc1=5.61D-3
lc2=6.97D-2
lc3=8.07D-3
lc4=5.15D-1
lc5=5.74D-1
lc6=9.91D-3
lc7=1.18D-1
lc8=1.01D-2


!3.40[−2] 9.97[−2] 1.15[+0] 7.23[+0] 1.50[+1] 2.08[+1] 2.59[+1] 7.65[+1]
lw1=(1.602177D-19/1.0545D-34)*3.40D-2
lw2=(1.602177D-19/1.0545D-34)*9.97D-2
lw3=(1.602177D-19/1.0545D-34)*1.15D0
lw4=(1.602177D-19/1.0545D-34)*7.23D0
lw5=(1.602177D-19/1.0545D-34)*1.50D1
lw6=(1.602177D-19/1.0545D-34)*2.08D1
lw7=(1.602177D-19/1.0545D-34)*2.59D1
lw8=(1.602177D-19/1.0545D-34)*7.65D1



epstoluene=1.D0+(lc1/(1.D0+(wp/lw1)**2.D0))+(lc2/(1.D0+(wp/lw2)**2.D0))+(lc3/(1.D0+(wp/lw3)**2.D0))+&
&(lc4/(1.D0+(wp/lw4)**2.D0))+(lc5/(1.D0+(wp/lw5)**2.D0))+(lc6/(1.D0+(wp/lw6)**2.D0))+&
&(lc7/(1.D0+(wp/lw7)**2.D0))+(lc8/(1.D0+(wp/lw8)**2.D0))

ELSE
epstoluene=2.395D0
END IF
 

IF (i>0) THEN
mc1=49.338D0
mc2=7.7867D0
mc3=4.7422D0
mc4=3.6331D0
mc5=1.2529D0
mc6=2.2625D0
mc7=1.1828D0
mc8=0.3136D0
mc9=0.0672D0

mw1=(1.602177D-19/1.0545D-34)*0.0209D0
mw2=(1.602177D-19/1.0545D-34)*0.0772D0
mw3=(1.602177D-19/1.0545D-34)*0.3127D0
mw4=(1.602177D-19/1.0545D-34)*0.6912D0
mw5=(1.602177D-19/1.0545D-34)*2.4508D0
mw6=(1.602177D-19/1.0545D-34)*5.3771D0
mw7=(1.602177D-19/1.0545D-34)*13.0017D0
mw8=(1.602177D-19/1.0545D-34)*30.7646D0
mw9=(1.602177D-19/1.0545D-34)*72.2563D0




epsmagnetite=1.D0+(mc1/(1.D0+(wp/mw1)**2.D0))+(mc2/(1.D0+(wp/mw2)**2.D0))+(mc3/(1.D0+(wp/mw3)**2.D0))+&
&(mc4/(1.D0+(wp/mw4)**2.D0))+(mc5/(1.D0+(wp/mw5)**2.D0))+(mc6/(1.D0+(wp/mw6)**2.D0))+&
&(mc7/(1.D0+(wp/mw7)**2.D0))+(mc8/(1.D0+(wp/mw8)**2.D0))+(mc9/(1.D0+(wp/mw9)**2.D0))

ELSE
epsmagnetite=54.8D0
END IF
 

phi=0.01D0
mu2=1.D0
mu4=1.D0
diam=5D0
!diam=20D0

mu3=1.D0+(phi/0.05D0)*(1.24D0-1.D0)*((diam/10.D0)**3.D0)  
!mu3=1.D0


G3=(phi*(epsmagnetite-1.D0)/(epsmagnetite+2))+(1-phi)*(epstoluene-1.D0)/(epstoluene+2)

eps3=(1+2.D0*G3)/(1-G3)
!eps3=epstoluene


tc1=9.30D-3	! teflon PTFE van Zwol PRA2010
tc2=1.83D-2
tc3=1.39D-1
tc4=1.12D-1
tc5=1.95D-1
tc6=4.38D-1
tc7=1.06D-1
tc8=3.86D-2

tw1=(1.602177D-19/1.0545D-34)*3.00D-4
tw2=(1.602177D-19/1.0545D-34)*7.60D-3
tw3=(1.602177D-19/1.0545D-34)*5.57D-2
tw4=(1.602177D-19/1.0545D-34)*1.26D-1
tw5=(1.602177D-19/1.0545D-34)*6.71D0
tw6=(1.602177D-19/1.0545D-34)*1.86D1
tw7=(1.602177D-19/1.0545D-34)*4.21D1
tw8=(1.602177D-19/1.0545D-34)*7.76D1


eps4=1.D0+(tc1/(1.D0+(wp/tw1)**2.D0))+(tc2/(1.D0+(wp/tw2)**2.D0))+&
&(tc3/(1.D0+(wp/tw3)**2.D0))+&
&(tc4/(1.D0+(wp/tw4)**2.D0))+(tc5/(1.D0+(wp/tw5)**2.D0))+&
&(tc6/(1.D0+(wp/tw6)**2.D0))+&
&(tc7/(1.D0+(wp/tw7)**2.D0)) +&
&(tc8/(1.D0+(wp/tw8)**2.D0))


IF (i==0) THEN
eps4=2.1D0  ! 2.1 teflon  PTFE van Zwol PRA 2010
END IF

pc1=1.21D-2		! Polystyrene (set 1, 2008) van Zwol PRA 2010
pc2=2.19D-2
pc3=1.79D-2
pc4=3.06D-2
pc5=3.03D-1
pc6=6.23D-1
pc7=3.25D-1
pc8=3.31D-2


pw1=(1.602177D-19/1.0545D-34)*1.00D-3
pw2=(1.602177D-19/1.0545D-34)*1.32D-2
pw3=(1.602177D-19/1.0545D-34)*3.88D0
pw4=(1.602177D-19/1.0545D-34)*1.31D-1
pw5=(1.602177D-19/1.0545D-34)*5.99D0
pw6=(1.602177D-19/1.0545D-34)*1.02D1
pw7=(1.602177D-19/1.0545D-34)*1.88D1
pw8=(1.602177D-19/1.0545D-34)*5.15D1





!pc1=3.12D-2		! Polystyrene (set 2, 1977) van Zwol PRA 2010
!pc2=1.17D-2
!pc3=2.17D-2
!pc4=9.20D-3
!pc5=2.93D-1
!pc6=6.54D-1
!pc7=4.17D-1
!pc8=2.13D-2
!pw1=(1.602177D-19/1.0545D-34)*1.18D-1
!pw2=(1.602177D-19/1.0545D-34)*9.00D-4
!pw3=(1.602177D-19/1.0545D-34)*1.19D-2
!pw4=(1.602177D-19/1.0545D-34)*1.56D0
!pw5=(1.602177D-19/1.0545D-34)*6.12D0
!pw6=(1.602177D-19/1.0545D-34)*1.01D1
!pw7=(1.602177D-19/1.0545D-34)*2.02D1
!pw8=(1.602177D-19/1.0545D-34)*6.86D1


eps2=1.D0+(pc1/(1.D0+(wp/pw1)**2.D0))+(pc2/(1.D0+(wp/pw2)**2.D0))+&
&(pc3/(1.D0+(wp/pw3)**2.D0))+&
&(pc4/(1.D0+(wp/pw4)**2.D0))+(pc5/(1.D0+(wp/pw5)**2.D0))+&
&(pc6/(1.D0+(wp/pw6)**2.D0))+&
&(pc7/(1.D0+(wp/pw7)**2.D0)) +&
&(pc8/(1.D0+(wp/pw8)**2.D0))



IF (i==0) THEN
eps2=2.4D0  ! 2.4-2.5 Polystyrene (set 1, 2008) van Zwol PRA 2010
END IF




!eps2=epsMoistCellulouseIM
!eps3=epswater
!eps4=epsvapor

qsum=0.D0
qsum2=0.D0
qsum3=0.D0

qmax=(25D0)/(dist)
qmin=0     

dq= (qmax-qmin)/DFLOAT(kmax)

 
DO k=1,kmax  !q integration from qmin to qmax using simpson's rule

q=qmin+k*dq

 
! weight in simpsons formula
	
  	
IF (k/2*2/=k) THEN 
  		
qweight=4.D0 
		
qbweight=0.D0
  	
END IF
  
  	
IF (k/2*2==k) THEN 
  		 
qweight=2.D0
		
IF ((k)/4*4/=k) THEN
		
qbweight=4.D0
		
END IF
		
IF ((k)/4*4==k) THEN
		 
qbweight=2.D0
		
END IF
	
END IF


	
IF ((k==0) .OR. (k==kmax) ) THEN 
  		
qweight=1.D0 
		
qbweight=1.D0
  	
END IF





!gamma1=DSQRT(q*q+eps1*wp*wp/(c*c))	!retarded wave vectors in media 1
gamma2=DSQRT(q*q+eps2*wp*wp/(c*c))
gamma3=DSQRT(q*q+eps3*wp*wp/(c*c))
gamma4=DSQRT(q*q+eps4*wp*wp/(c*c))



!print*, 'gamma1=', gamma1

D234=-(eps4*gamma3-eps3*gamma4)/(eps4*gamma3+eps3*gamma4) !the 2 indicated that it is TM refl coeff
!D241=-(eps1*gamma4-eps4*gamma1)/(eps1*gamma4+eps4*gamma1)
D232=-(eps2*gamma3-eps3*gamma2)/(eps2*gamma3+eps3*gamma2)




!print*, 'D234=', D234

D134=-(gamma3-gamma4)/(gamma3+gamma4)			!!the 1 indicated that it is TE refl coeff
!D141=-(gamma4-gamma1)/(gamma4+gamma1)
D132=-(gamma3-gamma2)/(gamma3+gamma2)



IF (i==0) THEN
D134=-(mu4*gamma3-mu3*gamma4)/(mu4*gamma3+mu3*gamma4)			!!the 1 indicated that it is TE refl coeff
!D141=-(gamma4-gamma1)/(gamma4+gamma1)
D132=-(mu2*gamma3-mu3*gamma2)/(mu2*gamma3+mu3*gamma2)
END IF



film1=D234*D232         !product of TM reflection coefficient
film2=D134*D132	 !product of TE reflection coefficient

!print*, 'film1=', film1

!(see the expression for Lifshitz energy contributions from TM and TE gives below)

Gdisp1=1.D0-DEXP(-2.D0*dist*gamma3)*film1	!retardation TM thin films 
Gdisp2=1.D0-DEXP(-2.D0*dist*gamma3)*film2	!retardation TE thin films

s=(Kb*T/(2.D0*pi))*(LOG(Gdisp1)+LOG(Gdisp2) )  
!Expression for free energy/unit area except for constants and sum+integration 
 
!use the derivative of free energy, below, if you want to calculate pressure rather than free energy
!B1=DEXP(-2.D0*dist*gamma3)*film1   !not used for free energy. If you want to calculate force use it
!B2= DEXP(-2.D0*dist*gamma3)*film2 !not used for free energy. If you want to calculate force use it
!s=-(Kb*T/(pi))*gamma3*((B1/Gdisp1)+(B2/Gdisp2)) !pressure



IF (i==0) THEN
	s=s/2.D0	!original sum is from minus to plus infinity so when going over to sum
			!from zero to infinity the zero frequency term should be multiplied by 0.5
END IF
 

qsum=qsum+s*dq*q*qweight/3.D0           !q integration
!qsum2=qsum2+s*2.D0*dq*q*qbweight/3.D0 	!q integration with different steplength to test accuracy
					!check Simpsons rule


END DO		!end of q loop


!qsum3=(qsum+(qsum-qsum2)/15)   !a way to improve the accuracy often just using qsum is ok
qsum3=qsum


bsum=bsum+qsum3 		!summand in frequency summation
 

!WRITE(40,*) wpeV,' ',(eps2-eps3)*(eps4-eps3)/((eps2+eps3)*(eps4+eps3)) 

!WRITE(44,*) wpeV,' ',eps2,eps3,eps4 


END DO		!End of i loop (Matsubara frequency summation)

 

print*, dist*1D9,' , ', bsum 
!WRITE(*,*) dist*1D9,'  , ',bsum   !print on screen dist in nm and free energy in J/m^2
WRITE(120,*) dist*1D9,' ,  ',bsum !print to file dist in nm and free energy in J/m^2 
  
 

 !Always test the dielectric function data to make sure they are generated correctly.
!save data for dieelectric function to data file and then plot often

  
END DO !dloop



END PROGRAM Casimir


