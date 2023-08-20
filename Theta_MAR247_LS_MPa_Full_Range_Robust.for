C Subroutine for Creep analysis using Theta Projection Method
C Updated 2nd July 2023 by VaradhaYamunan KK
C Material constants derived for MAR-M-247 using Least Square regression
C with Robust weights calculated through variances of Theta values
C Creep curve for 800degC&1200 MPa & 800degC&10 MPa added artificially 
C in material constant calculation to prevent very high initial strains 
C Units : Stresses in MPa, Forces in N, lengths in mm, time in s
C      
      SUBROUTINE CREEP(DECRA,DESWA,STATEV,SERD,EC,ESW,P,QTILD,TEMP,
     1DTEMP,PREDEF,DPRED,TIME,DTIME,CMNAME,LEXIMP,LEND,COORDS,NSTATV,
     2NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
	CHARACTER*80 CMNAME
C
	DIMENSION DECRA(5),DESWA(5),PREDEF(*),DPRED(*),TIME(3),STATEV(*),
     1 COORDS(*),EC(2),ESW(2)
C     
      DOUBLE PRECISION A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,A4,B4,C4,D4,
     1 TH1,TH2,TH3,TH4,EB,HC,RC,WC,H,R,W,T1,T2,T3,T4,HCxEB,Z,STRESS,F,
     2 dQTILD, STRESSJUMP,AF,BF,CF,DF
C         
      PARAMETER one=1.0d0
      PARAMETER asmall=1.0d-27
C     
	DO I=1,5
	DECRA(I)=0.0d0
	DESWA(I)=0.0d0
      ENDDO
C
      LEND=0
C
C
	A1= -1.97064763d+01
	B1= 8.13114243d-03
	C1= 1.13154440d-02
	D1= -3.06894657d-06

C      
	A2= -5.51920534d+01
	B2= 2.48375981d-03
	C2= 3.25455041d-02
	D2= 1.43266551d-05
C     
	A3= -1.15705943d+02
	B3= 9.39829012d-02
	C3= 9.16212439d-02
	D3= -7.86600386d-05
C    
	A4= -3.35468858d+01
	B4= -2.60404124d-02
	C4= 1.30591578d-02
	D4= 3.94639082d-05
C
C MATERIAL CONSTANTS FOR CALCULATING TIME TO FRACTURE
C
      AF= -2.88170858d+01
      BF= 3.57331943d-02
      CF= 2.22331983d-02
      DF= -3.11971782d-05
C      
C CORRECTION FOR STRESS SINGULARITY DUE TO ELEMENT DAMAGE      
C      
C      dQTILD = QTILD-STATEV(5)
C      
C      IF(STATEV(5).EQ.0.0d0) GO TO 1
C      IF(dQTILD.EQ.0.0d0) GO TO 1
C      
C      STRESSJUMP = ABS(dQTILD/STATEV(5))
C      
C      IF(STRESSJUMP.GT.1.0d0) THEN
C      STRESS = ABS(QTILD-dQTILD)
C      GO TO 2
C      ENDIF
C                 
    1 STRESS = QTILD
C     
C DEFINE LOWER & UPPER LIMITS OF STRESS
C      
      IF((STRESS.LE.10.0d0).OR.(STRESS.GE.1.5d3)) THEN
          DECRA(1)=asmall
          STATEV(4)=Z+DECRA(1)
      RETURN
      ENDIF
C            
C THETA CALCULATIONS
C      
    2 TH1=EXP(A1+(B1*STRESS)+(C1*TEMP)+(D1*TEMP*STRESS))
      TH2=EXP(A2+(B2*STRESS)+(C2*TEMP)+(D2*TEMP*STRESS))
      TH3=EXP(A3+(B3*STRESS)+(C3*TEMP)+(D3*TEMP*STRESS))
      TH4=EXP(A4+(B4*STRESS)+(C4*TEMP)+(D4*TEMP*STRESS))
C
C FRACTURE STRAIN TO DEFINE ELEMENT DELETION
C      
      F=EXP(AF+(BF*STRESS)+(CF*TEMP)+(DF*STRESS*TEMP))
C
      IF(F.LT.1.0d-2)THEN
          F=1.0d-2
      ENDIF
C 
C CONSTITUTIVE EQUATIONS
C (WITH APPROXIMATIONS TO PREVENT MATHEMATICAL SINGULARITY)   
C
      EB=(TH1*TH2)+(TH3*TH4)
C                 
      IF(EB.eq.0.0d0) THEN
         EB=asmall
      ENDIF
C 
   	HC=(TH2/EB)
C               
	RC=(HC)*(TH3*TH4)
C
      IF(TH3.eq.0.0d0) THEN
         TH3=asmall
      ENDIF
C      
      IF(TH2.eq.0.0d0) THEN
         TH2=asmall
      ENDIF
C           
      WC=(one/TH3)
C
C STATE VARIABLES DECLARATION
C
      H=STATEV(1)
      R=STATEV(2)
      W=STATEV(3)
      Z=STATEV(4)
C    
C CREEP STRAIN INCREMENT CALCULATION
C      
      DECRA(1) = (EB*(one+H+R+(W*RC/(TH2)))*DTIME)
C     
      IF(DECRA(1).LT.0.0d0) THEN
          DECRA(1)=asmall
      ENDIF
C
      Z=EC(1)
C      
C CALCULATION OF DECRA(2) & DECRA(5) FOR IMPLICIT INTEGRATION
C     
      IF(LEXIMP.EQ.1) THEN
C
      DECRA(2) = ((((-HC*EB)*(one+H+R+W))+RC+((RC*WC/HC)*
     1(one+H+R+((RC*W)/(HC*EB)))))*DTIME)/
     2(one+H+R+(RC*W/(HC*EB)))
C      
      T1 = B1+(D1*TEMP)
      T2 = B2+(D2*TEMP)
      T3 = B3+(D3*TEMP)
      T4 = B4+(D4*TEMP)
C      
      DECRA(5) = (((T1+T2*(one+LOG(one+H+R)))*(EB-(RC/HC))*(one+H+R))+
     1 ((T3+T4*(one+LOG(one+H+R+W)))*(RC/HC)*(one+H+R+W)))*DTIME
C
      ENDIF
C
C UPDATING STATE VARIABLES
C
    3	STATEV(1) = H-((one+H+R)*HC*EB*DTIME)
	STATEV(2) = R+(RC*DTIME)
      STATEV(3) = W+((one+W)*(WC*RC/HC)*DTIME)
      STATEV(4) = Z+DECRA(1)
      STATEV(5) = STRESS
C      
      IF (STATEV(4).GE.F) THEN
          STATEV(4)=0.0d0
      ENDIF
C
    4 RETURN
      END