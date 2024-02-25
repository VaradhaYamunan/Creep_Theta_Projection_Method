C Subroutine for Creep analysis using Theta Projection Method
C Updated 14th Dec 2023 by VaradhaYamunan KK
C Material constants derived for MAR-M-247 using Unique solution (Inversion)
C Experimental curves taken : 800degC_500 MPa, 900_250, 950_300, 1000_200
C Units : Stresses in MPa, Forces in N, lengths in mm, time in s
C
*USER SUBROUTINE
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
     1 TH1,TH2,TH3,TH4,EB,HC,RC,WC,H,R,W,T1,T2,T3,T4,HCxEB,STRESS,TF,
     2 dQTILD,STRESSJUMP,AF,BF,CF,DF,CREEP_RATE,dEB_dS,dH_dS,dR_dS,dW_dS
     3 WF,CLF
C         
      PARAMETER one=1.0d0
      PARAMETER asmall=1.0d-27
C     
	DO I=1,5
	DECRA(I)=0.0d0
	DESWA(I)=0.0d0
      ENDDO
C
	A1= -2.9949480806d+01
	B1= -1.2575371687d-03
	C1= 1.7711238000d-02
	D1= 1.1085235857d-05
C      
	A2= -6.3920693872d+01
	B2= -1.3215089272d-02
	C2= 3.9407616148d-02
	D2= 3.2218635290d-05
C     
	A3= 6.9504924093d+00
	B3= -1.3207327259d-01
	C3= -1.4542686231d-02
	D3= 1.2525706687d-04
C    
	A4= -5.1855226072d+01
	B4= 4.6840626209d-03
	C4= 2.8571788751d-02
	D4= 1.2513241184d-05
C
C MATERIAL CONSTANTS FOR CALCULATING TIME TO FRACTURE
C
      AF= 5.17133890d+01
      BF= 2.14571905d-02
      CF= -2.65877477d-02
      DF= -3.84887680d-05
C            
    1 STRESS = QTILD
C     
C DEFINE LOWER & UPPER LIMITS OF STRESS
C      
      IF(STRESS.LE.1d1)THEN
          STRESS=1d1
      ELSEIF(STRESS.GE.1d3)THEN
          STRESS=1d3
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
      TF=EXP(AF+(BF*STRESS)+(CF*TEMP)+(DF*STRESS*TEMP))
C
      IF(TF.LT.1.0d-2)THEN
          TF=1.0d-2
      ENDIF
      IF(TF.GE.1d30)THEN
          TF=1d30
      ENDIF
C      
      WF = EXP(TF*TH4)-1d0
      IF(WF.GE.1d30)THEN
          WF=1d30
      ENDIF
C      
      IF(WF.LE.1d-30)THEN
          WF=1d-30
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
C    
C CREEP STRAIN INCREMENT CALCULATION
C 
      CREEP_RATE=EB*(1d0+H+R+(W*RC/TH2))
C
      DECRA(1) = CREEP_RATE*DTIME
C     
      IF(DECRA(1).LT.0.0d0) THEN
          DECRA(1)=asmall
      ENDIF
C      
C CALCULATION OF DECRA(2) & DECRA(5) FOR IMPLICIT INTEGRATION
C     
      IF(LEXIMP.EQ.1) THEN
C
      DECRA(2)=((-TH2*CREEP_RATE)+(TH2*TH3*TH4)+((1d0+W)*TH3*(TH4**2)))
     1*(DTIME/CREEP_RATE)
C            
      dEB_dS=((B1+B2+(TEMP*(D1+D2)))*TH1*TH2)+
     1((B3+B4+(TEMP*(D3+D4)))*TH3*TH4)
C
      dR_dS=(RC*TIME(2)/EB)*((EB*(B2+B3+B4+(TEMP*(D2+D3+D4))))-dEB_dS)
C
      dH_dS=(-(1d0+R)*(B2+(D2*TEMP))*TH2*TIME(2)*(exp(-TH2*TIME(2))))
     1+(dR_dS*((exp(-TH2*TIME(2)))-1d0))
C
      dW_dS=(B4+(D4*TEMP))*TH4*TIME(2)*(exp(TH4*TIME(2)))
C
      DECRA(5)=((dEB_dS*(1d0+H+R))+(EB*(dH_dS+dR_dS))+(TH3*TH4*dW_dS)
     1+(W*TH4*TH3*(B3+(D3*TEMP)))+(W*TH3*TH4*(B4+(D4*TEMP))))*DTIME
C
      ENDIF
C      
       IF(KSTEP.LE.1)THEN
          IF(KINC.LE.2)THEN
                  DECRA(2)=0d0
                  DECRA(5)=0d0
          ENDIF
      ENDIF
C
C UPDATING STATE VARIABLES
C
    3	STATEV(1) = H-((one+H+R)*HC*EB*DTIME)
	STATEV(2) = R+(RC*DTIME)
      STATEV(3) = W+((one+W)*(WC*RC/HC)*DTIME)
C
C CREEP LIFE FRACTION CALCULATION (AS PER MINER'S RULE)
C      
      CLF = STATEV(4)+(((1d0+W)*(WC*RC/HC)*DTIME)/WF)
      STATEV(4)=CLF
C      
      IF(CLF.LT.1d-10) THEN
          CLF=1d-10
          STATEV(4)=1d-10
          STATEV(5)=1d10
      ENDIF     
C
C LIFE PREDICTION (IN HOURS)
C      
      IF(ABS(CLF).GE.1d-10) THEN
          STATEV(5) = (TIME(2))/((3.6d3)*ABS(STATEV(4)))
      ENDIF
C
C ELEMENT DELETION BASED ON CREEP LIFE FRACTION
C      
      IF (STATEV(3).GE.WF) THEN
          STATEV(5)=0d0
      ENDIF
C
      
    4 RETURN
      END