c      subroutine isz
      PARAMETER(PI=3.1415926535897932384626433832795,DNO=180./PI)
      ALLOCATABLE TOBS(:),X(:),V(:),XIN(:),VIN(:),ALF(:),DEL(:)
      ALLOCATABLE ALFC(:),DELC(:),RO(:),DISO(:,:,:)
      DIMENSION DNORM(6,6),RTL(6),OBR(6,6),POPR(6),ED(6,6),COV(6,6)
      DIMENSION XOP(6),EL(6),DER(6,6),PROM(6,6),COVEL(6,6),SIGEL(6)
      COMMON /PERT/IMON,ISUN,IGARM /UPBORDER/NM,MN /TJD0/TEPJD,TEPSEC
      COMMON /NOBS/NOBS /ATM/IATM,SM,CD,HBURN/NEXCL/NEXCL
      COMMON /NUMS/NUMS(1000) /INTERV/INTERV /NT/NT /TIDE/ITIDE
      COMMON /RELAT/IREL /RELCUST/IGPR,ISCHW,ILTER /LIGHT/ILGT,QOTR
      COMMON /IREGIM/IREGIM /TOPC/DXY,DZ,DLAM,TOPC(3) /HXV/HX,HV 
      COMMON /DA/DA
      DATA SUT/86400./ ZERO,ONE/0.,1./ NI/2/ NCLASS/2/ NEXCL/-1/
      CHARACTER*12 OBSFILE
      CHARACTER*4 ELNM(6)
      DATA ELNM/'a','e','i','Node','Peri','M'/

      OPEN(3,FILE='ISZM_PUC.IN')
      OPEN(4,FILE='EPH.OUT') 

      READ(3,*) IREGIM
      IF(IREGIM==1) THEN; NT=1; ELSE
      IF(IREGIM==2) THEN; NT=7; ELSE
      STOP 'Неверно задан режим, допустимы значения 1 и 2'
      ENDIF; ENDIF
      
      
      READ(3,*) IYEAR,IMONTH,DAY,IHR,IMIN,SEC
      CALL EPOCH(IYEAR,IMONTH,DAY,IHR,IMIN,SEC,TEPJD,TEPSEC)
C      TEPSEC=TEPSEC-UTCMTT(TEPJD+TEPSEC/SUT)
      READ(3,*) NSP
      IF(IREGIM==1) THEN; NT=NSP; ENDIF 
      DO I=1,NT; NUMS(I)=I; ENDDO 
      CALL INITCONST(NT)
      NV=3*NT
      ALLOCATE(X(NV),V(NV),XIN(NV),VIN(NV))
      DO ISP=1,NSP
      II=3*ISP-3
      READ(3,*) (X(II+K),K=1,3)
      READ(3,*) (V(II+K),K=1,3)
      ENDDO
      XIN=X; VIN=V
      READ(3,*); READ(3,*)
      READ(3,*) IYEAR,IMONTH,DAY,IHR,IMIN,SEC
      CALL EPOCH(IYEAR,IMONTH,DAY,IHR,IMIN,SEC,T0JD,T0SEC)
      READ(3,*) IYEAR,IMONTH,DAY,IHR,IMIN,SEC
      CALL EPOCH(IYEAR,IMONTH,DAY,IHR,IMIN,SEC,TFJD,TFSEC)
      READ(3,*) STEP
      READ(3,*) DA
      READ(3,*); READ(3,*)
      READ(3,*) OBSFILE
      READ(3,*) NOBSIMP
      READ(3,*) TOPC
      DLAM=ATAN2(TOPC(2),TOPC(1))
      DXY=SQRT(TOPC(1)**2+TOPC(2)**2)
      DZ=TOPC(3)
      READ(3,*) EPSIMP
      READ(3,*) HX,HV
      READ(3,*) HIMP
      READ(3,*) INAC
      READ(3,*) ROIN
      READ(3,*); READ(3,*)
      READ(3,*) XL
      READ(3,*) NOR
      READ(3,*) LL
      READ(3,*) INTERV
      READ(3,*); READ(3,*)
      READ(3,*) NM,MN
      IF(NM<2) THEN; IGARM=0; ELSE; IGARM=1; ENDIF
      IF(MN>NM) 
     * STOP 'Неверные входные данные: для гармоник должно быть N>=M'
      NM=NM+2
      MN=MN+1
      READ(3,*) IMON
      READ(3,*) ISUN
      READ(3,*) ILGT
      READ(3,*) IGPR,ISCHW,ILTER
      READ(3,*) ITIDE
      READ(3,*) IATM
      READ(3,*) HBURN
      READ(3,*); READ(3,*); READ(3,*) DMASS
      READ(3,*) SM; SM=SM/DMASS
      READ(3,*) CD
      READ(3,*) QOTR
      CLOSE(3)
      IREL=0
      IF(IGPR/=0.OR.ISCHW/=0.OR.ILTER/=0) IREL=1

      IF(IREGIM==2)  THEN
        OPEN(7,FILE=OBSFILE)
        NOBS=NOBSIMP-1
        ALLOCATE(TOBS(0:NOBS),ALF(0:NOBS),DEL(0:NOBS))
        DO I=0,NOBS
        READ(7,*) IYEAR,IMONTH,DAY,IHR,IMIN,SEC,
     p    ALH,ALM,ALS,DELG,DELM,DELS
        ZDEL=SIGNUM(DELG)
        ALF50=(ALH+(ALM+ALS/60.)/60.)*15./DNO
        DEL50=ZDEL*(ABS(DELG)+(DELM+DELS/60.)/60.)/DNO
        CALL P1950_2000(ALF50,DEL50,ALF(I),DEL(I))
C        WRITE(*,'(4F15.12)') ALF50,DEL50,ALF(I),DEL(I)
C        ALF(I)=ALF50; DEL(I)=DEL50
        CALL EPOCH(IYEAR,IMONTH,DAY,IHR,IMIN,SEC,TCURJD,TCURSEC)
        dt=UTCMTT(TCURJD+TCURSEC/SUT)
        TCURSEC=TCURSEC-dt
        IF(I==0) THEN
          T0JD=TCURJD; T0SEC=TCURSEC-0.1
          IF(INAC==2) THEN
            TEPJD=TCURJD; TEPSEC=TCURSEC
            ENDIF
          ENDIF
        IF(I==NOBS) THEN
          TFJD=TCURJD; TFSEC=TCURSEC+0.1
          ENDIF
        TOBS(I)=DIFJD(TCURJD,TEPJD)*SUT+(TCURSEC-TEPSEC)
        ENDDO
        ENDIF
                                
      ENEP=TEPJD+TEPSEC/SUT                       
      ENT0=T0JD+T0SEC/SUT                         
      ENTF=TFJD+TFSEC/SUT

      IF(IREGIM==1) THEN
      NOBS=CEILING((ENTF-ENT0)*SUT/STEP)
      ALLOCATE(TOBS(0:NOBS))
      ENDIF

      ALLOCATE(ALFC(0:NOBS),DELC(0:NOBS),RO(0:NOBS),DISO(0:NOBS,3,6))       

      IF(IREGIM==1) THEN
      TOBS(0)=DIFJD(T0JD,TEPJD)*SUT+(T0SEC-TEPSEC)
      DO I=1,NOBS                                 
       TOBS(I)=TOBS(0)+FLOAT(I)*STEP             
       ENDDO
      ENDIF                          
      
      IF(IREGIM==1) THEN
      CALL PROGNOS(ENEP,ENT0,ENTF,TEPJD,T0JD,TFJD,TEPSEC,T0SEC,TFSEC,X,
     :V,XL,LL,NV,NI,NF,NS,NCLASS,NOR,NOBS,TOBS,ALFC,DELC,RO,NT,NEXCL,
     :DISO)
      ELSE                      
      OPEN(8,FILE='IMPROVE.OUT')
      IF(INAC==2) CALL NACPRIBL(NOBS,TOBS,ALF,DEL,ROIN,NV,XIN,VIN) 
      RAZ=1.; IT=0              
      DO WHILE(RAZ>EPSIMP)
      IT=IT+1
      WRITE(8,*) 'Итерация:',it
      WRITE(*,*) 'Итерация:',it
      X=XIN; V=VIN 
      CALL VARINIT(NV,X,V,HX,HV)   
      CALL PROGNOS(ENEP,ENT0,ENTF,TEPJD,T0JD,TFJD,TEPSEC,T0SEC,TFSEC,X,
     :V,XL,LL,NV,NI,NF,NS,NCLASS,NOR,NOBS,TOBS,ALFC,DELC,RO,NT,NEXCL,
     :DISO) 
      CALL NORMMATR(NOBS,ALF,DEL,ALFC,DELC,RO,DISO,DNORM,RTL,SIGMA,
     :TOBS) 
      
      CALL OBRSING(6,DNORM,OBR)

      COV=OBR*SIGMA

      POPR=0.
       DO I=1,6
       DO K=1,6; POPR(I)=POPR(I)+OBR(I,K)*RTL(K); ENDDO
       ENDDO
      RAZ=SQRT(POPR(1)**2+POPR(2)**2+POPR(3)**2)
      XIN(1:3)=XIN(1:3)-HIMP*POPR(1:3)
      VIN(1:3)=VIN(1:3)-HIMP*POPR(4:6)
      WRITE(8,*) 'Улучшенные координаты и скорости:'
      WRITE(8,'(3G27.19)') XIN(1:3),VIN(1:3)

      XOP(1:3)=XIN(1:3); XOP(4:6)=VIN(1:3)
      CALL COOREL(0,0.,0.,XOP,EL)
      CALL DERIVELC(XOP,1E-4,DER)
      ONX=1.-EL(2)**2
      DER(1,1:6)=(DER(1,1:6)*ONX+2.*EL(2)*EL(1)*DER(2,1:6))/ONX**2
      CALL UMATRN(6,0,0,DER,COV,PROM)
      CALL UMATRN(6,0,1,PROM,DER,COVEL)
       DO I=1,6
       IF(COVEL(I,I)>=0) THEN; 
       SIGEL(I)=SQRT(COVEL(I,I))*3.
       IF(I>=3) SIGEL(I)=SIGEL(I)*DNO
       ELSE; SIGEL(I)=1E300
       ENDIF
       ENDDO

      EL(1)=EL(1)/ONX
      EL(3:6)=EL(3:6)*DNO
      WRITE(8,*) 'Оскулирующие элементы:'
      WRITE(8,1) (ELNM(I),EL(I),SIGEL(I),I=1,6)
   1  FORMAT(A4,G27.19,' +- ',G14.7)
      SIGMA=SQRT(SIGMA)*DNO*3600.
      WRITE(8,'(A7,F10.2,A1)') 'sigma= ',SIGMA,'"'
      WRITE(*,*) 'sigma= ',SIGMA,' POPR=',RAZ
      ENDDO
      ENDIF
      CLOSE(4)
      CLOSE(3)
      CLOSE(7)
      END
************************************************************************
      SUBROUTINE NACPRIBL(NOBS,TOBS,ALF,DEL,ROIN,NV,X,V)
      ! Приближенное определение скорости при улучшении
      DIMENSION TOBS(0:NOBS),ALF(0:NOBS),DEL(0:NOBS),X(NV),V(NV)
      DIMENSION X2(3)
      COMMON /TOPC/DXY,DZ,DLAM /TJD0/TEPJD,TEPSEC
      S=SID2000(TEPJD+(TEPSEC+TOBS(0))/86400.)+DLAM  
      H=ROIN-6378.
      X(1)=H*COS(DEL(0))*COS(ALF(0))+DXY*COS(S) 
      X(2)=H*COS(DEL(0))*SIN(ALF(0))+DXY*SIN(S) 
      X(3)=H*SIN(DEL(0))+DZ              
      X2(1)=H*COS(DEL(1))*COS(ALF(1))+DXY*COS(S) 
      X2(2)=H*COS(DEL(1))*SIN(ALF(1))+DXY*SIN(S) 
      X2(3)=H*SIN(DEL(1))+DZ
      TIME=TOBS(1)-TOBS(0)
      V(1:3)=(X2(1:3)-X(1:3))/TIME 
      END
************************************************************************
      SUBROUTINE PROGNOS(ENEP,ENT0,ENTF,TEPJD,T0JD,TFJD,TEPSEC,T0SEC,
      ! Численное интегрирование
     :TFSEC,X,V,XL,LL,NV,NI,NF,NS,NCLASS,NOR,NOBS,TOBS,ALFC,DELC,RO,NT,
     :NEXCL,DISO)
      PARAMETER(SUT=86400.)
      DIMENSION X(NV),V(NV),XIN(NV),VIN(NV),TOBS(0:NOBS),ALFC(0:NOBS)
      DIMENSION DELC(0:NOBS),RO(0:NOBS),DISO(0:NOBS,3,6)
      DATA ZERO/0./
      XIN=X; VIN=V; NTIN=NT
      TI=ZERO
      IF(ENEP<=ENT0) THEN
        TF=DIFJD(TFJD,TEPJD)*SUT+TFSEC-TEPSEC+0.1 
   11   CALL RADA39(X,V,TI,TF,XL,LL,NV,NI,NF,NS,NCLASS,NOR,NOBS,TOBS,
     :  ALFC,DELC,RO,DISO)
        IF(NT==0) STOP 'Не осталось объектов'
        IF(NEXCL/=-1.AND.NT>=1) THEN; NEXCL=-1; GOTO 11; ENDIF
        GOTO 1
        ENDIF
      IF(ENEP>=ENTF) THEN
        TF=DIFJD(T0JD,TEPJD)*SUT+T0SEC-TEPSEC-0.1
   12   CALL RADA39(X,V,TI,TF,XL,LL,NV,NI,NF,NS,NCLASS,NOR,NOBS,TOBS,
     :  ALFC,DELC,RO,DISO)
        IF(NT==0) STOP 'Не осталось объектов'
        IF(NEXCL/=-1.AND.NT>=1) THEN; NEXCL=-1; GOTO 12; ENDIF
        GOTO 1
        ENDIF
      IF(ENEP>ENT0.AND.ENEP<ENTF) THEN
        TF=DIFJD(T0JD,TEPJD)*SUT+T0SEC-TEPSEC-0.1
   13   CALL RADA39(X,V,TI,TF,XL,LL,NV,NI,NF,NS,NCLASS,NOR,NOBS,TOBS,
     :  ALFC,DELC,RO,DISO)
        IF(NT==0) STOP 'Не осталось объектов'
        IF(NEXCL/=-1.AND.NT>=1) THEN; NEXCL=-1; GOTO 13; ENDIF
        X=XIN; V=VIN; NT=NTIN; NV=3*NT
        TF=DIFJD(TFJD,TEPJD)*SUT+TFSEC-TEPSEC+0.1
   14   CALL RADA39(X,V,TI,TF,XL,LL,NV,NI,NF,NS,NCLASS,NOR,NOBS,TOBS,
     :  ALFC,DELC,RO,DISO)
        IF(NT==0) STOP 'Не осталось объектов'
        IF(NEXCL/=-1.AND.NT>=1) THEN; NEXCL=-1; GOTO 14; ENDIF
        ENDIF                                              
    1 CONTINUE
      END
************************************************************************
      SUBROUTINE NORMMATR(NOBS,ALF,DEL,ALFC,DELC,RO,DISO,DNORM,RTL,
     :SIGMA,TOBS)
      ! Вычисление нормальной матрицы и ср. кв. ошибки
      PARAMETER(PI=3.1415926535897932384626433832795,DNO=180./PI)
      DIMENSION ALF(0:NOBS),DEL(0:NOBS),ALFC(0:NOBS),DELC(0:NOBS)
      DIMENSION RO(0:NOBS),DISO(0:NOBS,3,6),DNORM(6,6),RTL(6),DIS(2,6)
      DIMENSION TOBS(0:NOBS)
      
      DNORM=0.; RTL=0.; SIGMA=0.;
      DO ICN=0,NOBS
      SIALF=SIN(ALFC(ICN)); SIDEL=SIN(DELC(ICN))
      COALF=COS(ALFC(ICN)); CODEL=COS(DELC(ICN)) 
       DO I=1,6
       DIS(1,I)=(-SIALF*DISO(ICN,1,I)+COALF*DISO(ICN,2,I))/RO(ICN)
       DIS(2,I)=(-SIDEL*COALF*DISO(ICN,1,I)-SIDEL*SIALF*DISO(ICN,2,I)
     :           +CODEL*DISO(ICN,3,I))/RO(ICN)
       ENDDO
      DEALF=ANGDIF(ALFC(ICN),ALF(ICN))*CODEL
      DEDEL=DELC(ICN)-DEL(ICN)
      WRITE(*,'(I3,2F10.2)') ICN+1,DEALF*DNO*3600.,DEDEL*DNO*3600.
      WRITE(8,'(I3,2F10.2)') ICN+1,DEALF*DNO*3600.,DEDEL*DNO*3600.
      SIGMA=SIGMA+DEALF**2+DEDEL**2
       DO I=1,6
       RTL(I)=RTL(I)+DIS(1,I)*DEALF+DIS(2,I)*DEDEL
       DO K=1,6
       DNORM(I,K)=DNORM(I,K)+DIS(1,I)*DIS(1,K)+DIS(2,I)*DIS(2,K)
       ENDDO
       ENDDO
      ENDDO
      SIGMA=SIGMA/(2*NOBS-4)
      END
************************************************************************
      FUNCTION ANGDIF(ANG1,ANG2)
      ! Разность углов
      PARAMETER(PI=3.1415926535897932384626433832795,PI2=PI*2.)
      DIF=ANG1-ANG2
      IF(ABS(DIF)>PI) THEN
      IF(DIF>0) DIF=DIF-PI2
      IF(DIF<0) DIF=DIF+PI2
      ENDIF
      ANGDIF=DIF
      END
************************************************************************
      SUBROUTINE VARINIT(NV,X,V,HX,HV)
      ! Вариация начальных параметров 
      DIMENSION X(NV),V(NV)
      DO I=1,6
      II=3*I
      X(II+1:II+3)=X(1:3)
      V(II+1:II+3)=V(1:3)
      IF(I<=3) THEN 
      X(II+I)=X(I)+HX
      ELSE 
      V(II+I-3)=V(I)+HV 
      ENDIF
      ENDDO
      END   
************************************************************************
      FUNCTION SIGNUM(A)
      ! Знак числа
      IF(A>=0.) THEN; SIGNUM=1.; ELSE; SIGNUM=-1.; ENDIF
      END
************************************************************************
      FUNCTION DIFJD(A,B)
      PARAMETER(FAC=10.**7)
      CHARACTER*50 S
      DIF=(A-B)*FAC
      WRITE(S,'(F50.0)') DIF
      READ(S,*) DIF
      DIFJD=DIF/FAC
      END
************************************************************************
      SUBROUTINE EPOCH(IYEAR,IMONTH,DAY,IHR,IMIN,SEC,DJD,TSEC)
      ! YMDHMS -> JD Sec
      DJD=DJDATE(IYEAR,IMONTH,DAY)
      TSEC=FLOAT(IHR)*3600.+FLOAT(IMIN)*60.+SEC
      END
************************************************************************
      FUNCTION DJDATE(IYEAR,IMONTH,DAY)
      ! YMD -> JD
      DATA HAN/100./ C4/1E4/
      REAL M1,JD
      DATE=FLOAT(IYEAR)+FLOAT(IMONTH)/HAN+DAY/C4
      I=DATE
      M1=(DATE-FLOAT(I))*HAN
      ME=M1
      D=(M1-FLOAT(ME))*HAN
      IF(ME>2) GOTO 1
      I=I-1                                                    
      ME=ME+12                                                 
    1 JD=AINT(365.25*I)+AINT(30.6001*(ME+1))+D+1720994.5 
      IF(DATE<1582.1015) THEN                                     
        DJDATE=JD                                             
        RETURN                                                    
        ENDIF                                                     
      JA=I/100
      JB=2-JA+JA/4                                      
      JD=JD+FLOAT(JB)
      DJDATE=JD
      END
************************************************************************
      SUBROUTINE WRITEDATE(SEC)
      ! Запись даты
      COMMON /TJD0/TEPJD,TEPSEC                     
      DATA SUT/86400./
      DJD=TEPJD+(TEPSEC+SEC)/SUT
      CALL CALDAT(DJD,IYEAR,IMONTH,IDAY,IHR,IMN,DSEC)
      WRITE(*,1) IYEAR,IMONTH,IDAY,IHR,IMN,DSEC
    1 FORMAT(I5,2I3,I5,'h',I3,'m',F10.6,'s')
      END
************************************************************************
      SUBROUTINE CALDAT(DJD,IYEAR,IMONTH,IDAY,IHR,IMN,DSEC)
      ! JD -> YMDHMS
      DJD0=AINT(DJD+0.5)
      IF(DJD0<2299161.) THEN
      IB=0
      C=DJD0+1524.
      ELSE
      IB=(DJD0-1867216.25)/36524.25
      C=DJD0+FLOAT(IB-IB/4)+1525.
      ENDIF
      ID=(C-122.1)/365.25
      E=365.*FLOAT(ID)+ID/4
      KF=(C-E)/30.6001
      IDAY=AINT(C-E+0.5)-AINT(30.6001*KF)
      IMONTH=KF-1-12*(KF/14)
      IYEAR=ID-4715-(7+IMONTH)/10
      HOUR=24.*(DJD+0.5-DJD0)
      IHR=HOUR
      DMN=(HOUR-FLOAT(IHR))*60.
      IMN=DMN
      DSEC=(DMN-FLOAT(IMN))*60.
      END
************************************************************************
      SUBROUTINE FORCE(NV,X,V,T,F)
      ! Пр. части диф. ур-ний
      DIMENSION X(NV), V(NV), F(NV), XM(6), XS(6), CONX(3), DFATM(3)
      DIMENSION DFREL(3),SP(3)
      COMMON /NT/NT /GR/GE,GM,GS /PERT/IMON,ISUN,IGARM/TJD0/TJD0,TEPSEC
      COMMON /ATM/IATM,SM,CD /NOBJ/I /TIDE/ITIDE /RELAT/IREL
      COMMON /LIGHT/ILGT,QOTR
      DATA SUT/86400./ ZERO/0./
      IF(IMON/=0.OR.ISUN/=0.OR.ITIDE/=0.OR.IREL/=0.OR.ILGT/=0) 
     :CALL READ405(T,XM,XS) ! Координаты Луны (и Солнца)
      IF(IGARM/=0) THEN
        TJD=TJD0+(TEPSEC+T)/SUT
        CALL CELTOTER(TJD,ZERO)        
        IF(ITIDE/=0) CALL TIDES(XS,XM)
        ENDIF
      DO I=1,NT
      ! Кеплеровский член
      II=3*I-3
      R2=X(II+1)**2+X(II+2)**2+X(II+3)**2
      R3=R2*SQRT(R2)
      F(II+1:II+3)=-GE*X(II+1:II+3)/R3
      ! Возмущения от Луны
      IF(IMON/=0) THEN 
       R2=XM(1)**2+XM(2)**2+XM(3)**2
       R3=R2*SQRT(R2)
       D2=ZERO
       DO K=1,3; D2=D2+(XM(K)-X(II+K))**2; ENDDO
       D3=D2*SQRT(D2)
       F(II+1:II+3)=F(II+1:II+3)+GM*((XM(1:3)-X(II+1:II+3))/D3
     : -XM(1:3)/R3)
       ENDIF
      ! Возмущения от Солнца
      IF(ISUN/=0) THEN 
       R2=XS(1)**2+XS(2)**2+XS(3)**2
       R3=R2*SQRT(R2)
       D2=ZERO
       DO K=1,3; D2=D2+(XS(K)-X(II+K))**2; ENDDO
       D3=D2*SQRT(D2)
       F(II+1:II+3)=F(II+1:II+3)+GS*((XS(1:3)-X(II+1:II+3))/D3
     : -XS(1:3)/R3)
       ENDIF
      ! Возмущения от несферичности Земли
      IF(IGARM/=0) THEN
        CALL GEOPOT(X(II+1),CONX)
        F(II+1:II+3)=F(II+1:II+3)+CONX(1:3)
        ENDIF
      ! Возмущения от сопротивления атмосферы
      IF(IATM/=0) THEN
        CALL ATMOS(X(II+1),V(II+1),SM,CD,DFATM)
        F(II+1:II+3)=F(II+1:II+3)+DFATM(1:3)
        ENDIF
      IF(IREL/=0) THEN
        CALL RLTVT(X(II+1),V(II+1),XS,DFREL)
        F(II+1:II+3)=F(II+1:II+3)+DFREL(1:3)
        ENDIF
      IF(ILGT/=0) THEN
        CALL SUN_LITE(X(II+1),V(II+1),XS,QOTR,SM,SP)
        F(II+1:II+3)=F(II+1:II+3)+SP(1:3)
        ENDIF
      ENDDO
      END
************************************************************************
      SUBROUTINE RADA39(X,V,TI,TF,XL,LL,NV,NI,NF,NS,NCLASS,NOR,
     :NOBS,TOBS,ALFC,DELC,RO,DISO)
      ! Интегратор Эверхарта до 39-го порядка
      DIMENSION TOBS(0:NOBS),ALFC(0:NOBS),DELC(0:NOBS),RO(0:NOBS)
      DIMENSION DISO(0:NOBS,3,6)
      DIMENSION X(NV),V(NV),C(171),D(171),R(171),XI(171),HH(99),H(20)
      DIMENSION W(19),U(19),NW(20),MC(18),NXI(171)
      DIMENSION F1(NV),FJ(NV),Y(NV),Z(NV)
      DIMENSION BE(19,NV),BT(19,NV),B(19,NV)
      LOGICAL J2,NPQ,NSF,NPER,NCL,NES
      COMMON /TM/TM,T  /NEXCL/NEXCL /INTERV/INTERV
      DATA NW
     *  /0,0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153,171/
      DATA MC
     * /1,19,36,52,67,81,94,106,117,127,136,144,151,157,162,166,169,171/
      DATA EPS,C001,C04,C05,C08,ZERO,ONE,SR12,SR15,TEN,ELEVEN,C10000
     */1E-10,1E-2,4E-1,5E-1,8E-1,0.,1.,1.2,1.5,10.,11.,1E4/
      DATA NXI/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, 3,6,10,15,
     1 21,28,36,45,55,66,78,91,105,120,136,153,171,4,10,20,35,56,84,
     2 120,165,220,286,364,455,560,680,816,969, 5,15,35,70,126,210,330,
     3 495,715,1001,1365,1820,2380,3060,3876, 6,21,56,126,252,462,792,
     4 1287,2002,3003,4368,6188,8568,11628, 7,28,84,210,462,924,1716,
     5 3003,5005,8008,12376,18564,27132, 8,36,120,330,792,1716,3432,
     6 6435,11440,19448,31824,50388, 9,45,165,495,1287,3003,6435,12870,
     7 24310,43758,75582, 10,55,220,715,2002,5005,11440,24310,48620,
     8 92378, 11,66,286,1001,3003,8008,19448,43758,92378, 12,78,364,
     9 1365,4368,12376,31824,75582, 13,91,455,1820,6188,18564,50388,
     A 14,105,560,2380,8568,27132, 15,120,680,3060,11628, 16,136,816,
     B 3876, 17,153,969, 18,171, 19/
      DATA  HH/
     1 0.21234053823915294397475811012400038,  ! 7 порядок
     2 0.59053313555926528913507374793117011,
     3 0.91141204048729605260445385623054380,
     1 0.098535085798826426123498897887752690, ! 11 порядок
     2 0.30453572664636390548538517627883284,
     3 0.56202518975261385599498747999477115,
     4 0.80198658212639182746420786320470376,
     5 0.96019014294853125765919330990666683,
     1 0.056262560536922146465652191032311176, ! 15 порядок
     2 0.18024069173689236498757994280918178,
     3 0.35262471711316963737390777017124120,
     4 0.54715362633055538300144855765234885,
     5 0.73421017721541053152321060830661000,
     6 0.88532094683909576809035976293248537,
     7 0.97752061356128750189117450042915494,
     1 0.036257812883209460941164300768081452, ! 19 порядок
     2 0.11807897878999870019228511199473522,
     3 0.23717698481496038531730669285327447,
     4 0.38188276530470597536077024839649552,
     5 0.53802959891898906511685689131943569,
     6 0.69033242007236218294037953277051595,
     7 0.82388334383700471813682425392743141,
     8 0.92561261029080395536408181404400141,
     9 0.98558759035112345136717325918918676,
     1 0.025273620397520349753331186461629102, ! 23 порядок
     2 0.083041613447405146706865372981972582,
     3 0.16917510037718142596943345609434330,
     4 0.27779671510903207443667869219538909,
     5 0.40150272023286081677227928632695641,
     6 0.53186238691041595791688961924224641,
     7 0.65999184208533481176639476610298247,
     8 0.77715939295616214449216854654263724,
     9 0.87538077485555692626470041273609077,
     A 0.94796454887281944741645730422703568,
     B 0.98998171953831959415697527013219522,
     1 0.018610365010987851439719377840287688, ! 27 порядок
     2 0.061475540899268987602366613234700648,
     3 0.12630517869331058063228543286825611,
     4 0.20984297172656251444713666750062407,
     5 0.30789899828039834310295804831233027,
     6 0.41555603597865954449577915218908929,
     7 0.52741561399588227482490535732140273,
     8 0.63786860271776119959131870177267844,
     9 0.74137645929423748341020926717730900,
     A 0.83274898860844226850447752124064288,
     B 0.90740477530099736471710862456138364,
     C 0.96160186126032164962316747513585831,
     D 0.99263534897391067834930850158617766,
     1 0.014269454736825774734099366940870758, ! 31 порядок
     2 0.047299590094166685661955792475737906,
     3 0.097713299320621973368761495337990920,
     4 0.16356903939438987602444091434581680,
     5 0.24233526096865728800292572225971468,
     6 0.33098480497004012346130436094686043,
     7 0.42611083909331411932854614476247356,
     8 0.52405769153676513942741100798415032,
     9 0.62106131135302196189347099085722615,
     A 0.71339391374247294001597395451560359,
     B 0.79750724494989595243178001167976956,
     C 0.87016897444640894402874546190571009,
     D 0.92858704688484115994521609825326508,
     E 0.97051770135205751336835901528199629,
     F 0.99435931102748829024249353342055580,
     1 0.011285959367756248582878597747039905, ! 35 порядок
     2 0.037498802176649958749539224611494439,
     3 0.077756965407701591512892671769618231,
     4 0.13083394043258638079555350198884246,
     5 0.19511642147503350800082229860297450,
     6 0.26865104726594861786697850443702412,
     7 0.34920344154661773526851259267685342,
     8 0.43432603311336566809669302399035436,
     9 0.52143240314692041121628953714063214,
     A 0.60787586527790434196638053666809691,
     B 0.69102988224780409668247573871373687,
     C 0.76836787422282781705320277426164125,
     D 0.83753999657186345001837045608589682,
     E 0.89644456706480293658650877303868760,
     F 0.94329203123092176046728666844347936,
     G 0.97665990587622572603297847179231688,
     H 0.99554200643221260824549377088845360,
     1 0.0091481947290443146469051647798827689, ! 39 порядок
     2 0.030447362919779114528071213767733995,
     3 0.063304151925634919456094497918054973,
     4 0.10690686501815502924901884543508515,
     5 0.16018138429128892216136370557071174,
     6 0.22181577702323838579070035295792412,
     7 0.29029234834609760115625321688836449,
     8 0.36392495512072910530890894290183523,
     9 0.44090050735096801215307202412091851,
     A 0.51932360642144809736677864751729374,
     B 0.59726321383492290487695298889512237,
     C 0.67280019923618828427157003472433978,
     D 0.74407459651789703974795260270471668,
     E 0.80933140509523664891252533511280621,
     F 0.86696381144190122457591257602478368,
     G 0.91555277721579017940015209954883045,
     H 0.95390206695156891642607522916454091,
     I 0.98106816296841178274973547568305076,
     J 0.99638835718144310696554678987374802/
      FFLOAT(IIIIII)=FLOAT(IIIIII)
      FSQRT(XXXXXX)=SQRT(XXXXXX)
      FABS(XXXXXX)=ABS(XXXXXX)
      IF(TF==TI) RETURN
      EPS=10.**(-15)
      KD=(NOR-3)/2
      KD2=KD/2
      KE=KD+1
      KF=KD+2
      PW=ONE/FFLOAT(KD+3)
      NPER=.FALSE.
      NSF=.FALSE.
      NCL=NCLASS.EQ.1
      NPQ=NCLASS.LT.2
      SR=SR15
      IF(NV.EQ.1) SR=SR12
      NES=LL.LT.0
      TDIF=TF-TI
      DIR=TDIF/FABS(TDIF)
      IF(NES) XL=FABS(XL)*DIR
      NCLASS=IABS(NCLASS)
      LA=KD2*KD2-1
      DO 14 N=2,KF
      LA=LA+1
      H(N)=HH(LA)
      W(N-1)=ONE/FFLOAT(N+N**2*(NCLASS-1))
   14 U(N-1)=N+1
      DO 22 K=1,NV
      IF(NCL) V(K)=ZERO
      DO 22 L=1,KE
      BT(L,K)=ZERO
   22 B(L,K)=ZERO
      W1=ONE/FFLOAT(NCLASS)
      DO 939 J=1,KD
      M=MC(J)
      JD=J+1
      DO 939 L=JD,KE
      XI(M)=FFLOAT(NXI(M))*W(J)/W(L)
  939 M=M+1
      C(1)=-H(2)*W(1)
      D(1)= H(2)/W(2)
      R(1)=ONE/(H(3)-H(2))
      LA=1
      LC=1
      DO 73 K=3,KE
      LB=LA
      LA=LC+1
      LC=NW(K+1)
      JD=LC-LA
      C(LA)=-H(K)*C(LB)
      C(LC)=(C(LA-1)/W(JD)-H(K))*W(JD+1)
      D(LA)=H(2)*D(LB)*W(K-1)/W(K)
      D(LC)=(D(LA-1)*W(K-1)+H(K))/W(K)
      R(LA)=ONE/(H(K+1)-H(2))
      R(LC)=ONE/(H(K+1)-H(K))
      IF(K.EQ.3) GO TO 73
      DO 72 L=4,K
      LD=LA+L-3
      LE=LB+L-4
      JDM=LD-LA
      C(LD)=W(JDM+1)*C(LE)/W(JDM)-H(K)*C(LE+1)
      D(LD)=(D(LE)+H(L-1)*D(LE+1))*W(K-1)/W(K)
   72 R(LD)=ONE/(H(K+1)-H(L-1))
   73 CONTINUE
      SS=TEN**(-LL)
      NL=NI+30
      TP=((FFLOAT(NOR)/ELEVEN)*C05**(C04*FFLOAT(LL)))*DIR*C05
      IF(NES) TP=XL
      IF(TP/TDIF.GT.C05) TP=C05*TDIF
      NF=0
      NCOUNT=0
 4000 NS=0
      TM=TI
      SM=C10000
      CALL FORCE(NV,X,V,TM,F1)
      NF=NF+1
  722 DO 58 K=1,NV
      BE(KE,K)=B(KE,K)/W(KE)
      DO 58 J=1,KD
      JD=J+1
      BE(J,K)=B(J,K)/W(J)
      DO 58 L=JD,KE
      N=NW(L)+J
   58 BE(J,K)=BE(J,K)+D(N)*B(L,K)
      T=TP
      TVAL=FABS(T)
      T2=T**NCLASS
      DO 175 M=1,NL
      J2=.TRUE.
      DO 174 J=2,KF
      JD=J-1
      LA=NW(JD)
      JDM=J-2
      S=H(J)
      Q=S**(NCLASS-1)
      IF(NPQ) GO TO 5100
      DO 1300 K=1,NV
      RES=B(KE,K)
      TEMP=RES*U(KE)
      DO 7340 L=1,KD
      JR=KE-L
      RES=B(JR,K)+S*RES
 7340 TEMP=B(JR,K)*U(JR)+S*TEMP
      Y(K)=X(K)+Q*(T*V(K)+T2*S*(F1(K)*W1+S*RES))
 1300 Z(K)=V(K)+S*T*(F1(K)+S*TEMP)
      GO TO 5200
 5100 DO 1400 K=1,NV
      RES=B(KE,K)
      DO 2340 L=1,KD
      JR=KE-L
 2340 RES=B(JR,K)+S*RES
 1400 Y(K)=X(K)+Q*(T*V(K)+T2*S*(F1(K)*W1+S*RES))
 5200 CONTINUE
      TTT=TM+S*T
      CALL FORCE(NV,Y,Z,TTT,FJ)
      NF=NF+1
      IF(J2) GO TO 702
      DO 471 K=1,NV
      TEMP=BE(JD,K)
      RES=(FJ(K)-F1(K))/S
      N=LA
      DO 134 L=1,JDM
      N=N+1
  134 RES=(RES-BE(L,K))*R(N)
      BE(JD,K)=RES
      TEMP=RES-TEMP
      B(JD,K)=B(JD,K)+TEMP*W(JD)
      N=LA
      DO 471 L=1,JDM
      N=N+1
  471 B(L,K)=B(L,K)+C(N)*TEMP
      GO TO 174
  702 J2=.FALSE.
      DO 271 K=1,NV
      TEMP=BE(1,K)
      RES=(FJ(K)-F1(K))/S
      BE(1,K)=RES
  271 B(1,K)=B(1,K)+(RES-TEMP)*W(1)
  174 CONTINUE
      IF(M.LT.NI) GO TO 175
      HSUM=ZERO
       VAL=TVAL**(-KE)
      DO 635 K=1,NV
      BDUBL=B(KE,K)
  635 HSUM=HSUM+BDUBL**2
      HSUM=VAL*FSQRT(HSUM)
      IF(NSF) GO TO 175
      IF(FABS((HSUM-SM)/HSUM).LT.C001) GO TO 176
      SM=HSUM
  175 CONTINUE
  176 IF(NSF) GO TO 180
      TP=(SS/HSUM)**PW*DIR
      IF(NES) TP=XL
      IF(NES) GO TO 170
      IF(TP/T.GT.ONE) GO TO 170
    8 FORMAT(2X,I2,2G18.10)
      TP=C08*TP
      NCOUNT=NCOUNT+1
      IF(NCOUNT.GT.10) GOTO 999
      GO TO 4000
  170 NSF=.TRUE.
  180 CALL FINDINTERP(NV,NCLASS,NOR,NCL,T,X,V,F1,B,U,W1,NOBS,TOBS,
     :ALFC,DELC,RO,DISO)
      DO 35 K=1,NV
      RES=B(KE,K)
      DO 34 L=1,KD
   34 RES=RES+B(L,K)
      X(K)=X(K)+V(K)*T+T2*(F1(K)*W1+RES)
      IF(NCL) GO TO35
      RES=B(KE,K)*U(KE)
      DO 33 L=1,KD
   33 RES=RES+B(L,K)*U(L)
      V(K)=V(K)+T*(F1(K)+RES)
   35 CONTINUE
      TM=TM+T
      NS=NS+1
   74 IF(NPER) GOTO 999
      GO TO 2222
      NPER=.FALSE.
 2222 CALL FORCE(NV,X,V,TM,F1)
      NF=NF+1

        IF(NEXCL/=-1) THEN
        CALL EXCLUDE(NV,X,V)
        TI=TM
        RETURN
        ENDIF

        IF (MOD(NS,INTERV).EQ.0) CALL WRITEDATE(TM)

      IF(NES) GO TO 341
      TP=((SS/HSUM)**PW)*DIR
      IF(TP/T.GT.SR) TP=SR*T
  341 IF(NES) TP=XL
      IF(DIR*(TM+TP).LT.DIR*TF-EPS) GO TO 77
      TP=TF-TM
      NPER=.TRUE.
   77 Q=TP/T
      DO 39 K=1,NV
      RES=ONE
      DO 39 J=1,KE
      IF(NS.GT.1) BT(J,K)=B(J,K)-BT(J,K)
      IF(J.EQ.KE) GO TO 740
      M=MC(J)
      JD=J+1
      DO 40 L=JD,KE
      B(J,K)=B(J,K)+XI(M)*B(L,K)
   40 M=M+1
  740 RES=RES*Q
      TEMP=RES*B(J,K)
      B(J,K)=TEMP+BT(J,K)
   39 BT(J,K)=TEMP
      NL=NI
      GO TO 722
  999 CONTINUE
      END
************************************************************************      
      SUBROUTINE FINDINTERP(NV,NCLASS,NOR,NCL,T,X,V,F1,B,U,W1,NOBS,TOBS,
     :ALFC,DELC,RO,DISO)
      ! Результаты интегрирования
      PARAMETER(PI=3.1415926535897932384626433832795)
      DIMENSION TOBS(0:NOBS),ALFC(0:NOBS),DELC(0:NOBS),RO(0:NOBS)
      DIMENSION DISO(0:NOBS,3,6),XOP(6),EL(6),RT2C(3,3)
      DIMENSION X(NV),V(NV),XT(NV),VT(NV),F1(NV),B(19,NV),U(19)
      LOGICAL NCL
      COMMON /TM/TM,T_ /TJD0/TEPJD,TEPSEC /NT/NT /IREGIM/IREGIM 
      COMMON /NUMS/NUMS(1000) /GR/GE,GM_,GS_ /RE/RE /DA/DA
      IF(T>0) THEN; IIN=0; IEND=NOBS; ISTEP=1
              ELSE; IIN=NOBS; IEND=0; ISTEP=-1
      ENDIF
      DO I=IIN,IEND,ISTEP
      IF((TOBS(I)>=TM.AND.TOBS(I)<TM+T).OR.
     :  (TOBS(I)<=TM.AND.TOBS(I)>TM+T)) THEN
        S=(TOBS(I)-TM)/T
        CALL INTERPTOTIME(NV,NCLASS,NOR,NCL,S,T,X,V,XT,VT,F1,B,U,W1)
        TJD=TEPJD+(TEPSEC+TOBS(I))/86400.
        CALL ALFDEL(TJD,XT,ALFC(I),DELC(I),RO(I))
        IF(IREGIM==2) CALL ISOCHRON(I,NV,XT,VT,NOBS,DISO)
        CALL TERTOCEL(TJD,0.,RT2C)
        CALL CALDAT(TJD,IYEAR,IMONTH,IDAY,IHR,IMN,DSEC)
        IF(IREGIM==1) THEN
        WRITE(4,3) TEPJD,TEPSEC+TOBS(I),IYEAR,IMONTH,IDAY,IHR,IMN,DSEC
        ENDIF
    3   FORMAT(F10.1,G30.19,' (',I5,2I3,I5,'h',I3,'m',F10.6,'s)')
        DO J=1,NT
        JJ=3*J-3
        XOP(1:3)=XT(JJ+1:JJ+3); XOP(4:6)=VT(JJ+1:JJ+3)
        CALL COOREL(0,0.,0.,XOP,EL)
        A=EL(1)/(1.-EL(2)**2)
        DN=SQRT(GE/A**3)
        DSIG=1.5*DA/(A-RE)*DN*ABS(TEPSEC+TOBS(I))/PI*180*3600.
        IF(IREGIM==1) THEN
        WRITE(4,1) NUMS(J),(XT(JJ+K),K=1,3) 
        WRITE(4,2) (VT(JJ+K),K=1,3)
        WRITE(4,'(3G27.19)') ((RT2C(K,L),L=1,3),K=1,3)
        WRITE(4,4) DSIG
        ENDIF
    1   FORMAT(I6,3G27.19)
    2   FORMAT(6X,3G27.19)
    4   FORMAT('Прогнозируемая угловая ошибка: ',F13.2,'"'/)
        ENDDO
        ENDIF
      ENDDO
      END
************************************************************************
      SUBROUTINE INTERPTOTIME(NV,NCLASS,NOR,NCL,S,T,X,V,XT,VT,F1,B,U,W1)
      ! Вычисление решения по интерпол. ф-ле метода Эверхарта
      DIMENSION X(NV),V(NV),XT(NV),VT(NV),F1(NV),B(19,NV),U(19)
      LOGICAL NCL
      KD=(NOR-3)/2
      KE=KD+1
      T2=T**NCLASS
      Q=S**(NCLASS-1)
      DO 1300 K=1,NV
      RES=B(KE,K)
      IF(.NOT.NCL) TEMP=RES*U(KE)
      DO 7340 L=1,KD
      JR=KE-L
      RES=B(JR,K)+S*RES
      IF(.NOT.NCL) TEMP=B(JR,K)*U(JR)+S*TEMP
 7340 CONTINUE
      XT(K)=X(K)+Q*(T*V(K)+T2*S*(F1(K)*W1+S*RES))
      IF(.NOT.NCL) VT(K)=V(K)+S*T*(F1(K)+S*TEMP)
 1300 CONTINUE
      END
************************************************************************
      SUBROUTINE ALFDEL(TJD,X,ALFC,DELC,RO)
      ! Вычисление alpha, delta и rho
      PARAMETER(PI2=3.1415926535897932384626433832795*2.)
      DIMENSION X(3),XTOP(3),TOPCN(3)
      COMMON /TOPC/DXY,DZ,DLAM,TOPC(3) /AR/AR(3,3)
C      UTA=TJD+UT1MTT(TJD)/86400.
C      S=SID2000(UTA)+DLAM
C      S=GMST2000(UTA,0.,TJD,0.)+DLAM
C      XTOP(1)=X(1)-DXY*COS(S)
C      XTOP(2)=X(2)-DXY*SIN(S)
C      XTOP(3)=X(3)-DZ
      CALL CELTOTER(TJD,0.)
      TOPCN=0.
      DO I=1,3
      DO K=1,3; TOPCN(I)=TOPCN(I)+AR(K,I)*TOPC(K); ENDDO
      ENDDO
      XTOP(1:3)=X(1:3)-TOPCN(1:3)
      ALFC=ATAN2(XTOP(2),XTOP(1))
      IF(ALFC<0) ALFC=ALFC+PI2
      DELC=ATAN(XTOP(3)/SQRT(XTOP(1)**2+XTOP(2)**2))
      RO=SQRT(XTOP(1)**2+XTOP(2)**2+XTOP(3)**2)
      END
************************************************************************
      FUNCTION SID2000(DJD)
      ! Звездное время
      PARAMETER(PI=3.1415926535897932384626433832795,PI2=PI*2.)
      PARAMETER(DJ2000=2451545.,DJYEAR=36525.,SUT=86400.)
      DM=DJD-AINT(DJD)+0.5
      IF(DM>=1.) DM=DM-0.5
      D=DJD-DM-DJ2000
      T=(D+DM)/DJYEAR
      DMM=DM*SUT
      S=(24110.54841+DMM+236.555367908*(D+DM)+
     : (0.093104*T-6.2E-6*T**2)*T)/SUT*PI2
      B=S-AINT(S/PI2)*PI2
      IF(B<0) B=B+PI2
      SID2000=B
      END
************************************************************************
      SUBROUTINE P1950_2000(AL50,DE50,ALF,DEL)
      ! B1950 -> J2000
      PARAMETER(PI=3.1415926535897932384626433832795,PI2=PI*2.)
      DIMENSION ROTT(3,3),X50(3),X(3)
      DATA ROTT/
     :0.9999256794956877,-0.0111814832204662,-0.0048590038153592,
     :0.0111814832391717, 0.9999374848933135,-0.0000271625947142,
     :0.0048590037723143,-0.0000271702937440, 0.9999881946023742/
      X50(1)=COS(DE50)*COS(AL50)
      X50(2)=COS(DE50)*SIN(AL50)
      X50(3)=SIN(DE50)
      DO I=1,3
      X(I)=0.
      DO J=1,3; X(I)=X(I)+ROTT(J,I)*X50(J); ENDDO
      ENDDO
      ALF=ATAN2(X(2),X(1))
      IF(ALF<0.) ALF=ALF+PI2
      DEL=ATAN(X(3)/SQRT(X(1)**2+X(2)**2))
      END
************************************************************************
      SUBROUTINE ISOCHRON(ICN,NV,XT,VT,NOBS,DISO)
      ! Вычисление изохр. производных
      DIMENSION XT(NV),VT(NV),DISO(0:NOBS,3,6)
      COMMON /HXV/HX,HV
      DO I=1,3
      DO K=1,3
      DISO(ICN,I,K)=(XT(3*K+I)-XT(I))/HX
      DISO(ICN,I,K+3)=(XT(3*(K+3)+I)-XT(I))/HV
      ENDDO
      ENDDO
      END
************************************************************************
      SUBROUTINE EXCLUDE(NV,X,V)
      ! Исключение объекта за счет сгорания в атмосфере
      DIMENSION X(NV),V(NV)
      COMMON/NT/NT /NUMS/NUMS(1)/NEXCL/NEXCL 
      WRITE(*,1) NUMS(NEXCL),NT-1
      WRITE(4,1) NUMS(NEXCL),NT-1
    1 FORMAT('Объект N',I6,' вошел в плотные слои, осталось объектов:',
     :I6)
      DO I=NEXCL,NT-1
      II=3*I-3
      X(II+1:II+3)=X(II+4:II+6)
      V(II+1:II+3)=V(II+4:II+6)
      NUMS(I)=NUMS(I+1)
      ENDDO      
      NT=NT-1
      NV=3*NT
      IF(NT<1) WRITE(*,*) 'Не осталось объектов для интегрирования'
      END
************************************************************************
      SUBROUTINE ATMOS(X,V,SM,CD,P)
      ! Сопротивление атмосферы
      DIMENSION X(3),V(3),P(3),VA(3)
      ! X и V - полож. и скор. в ин. сист. координат (км, км/сек)
      ! VA - скорость во вращ. сист. координат, связ. с Землей
      ! SM - отнош-е площади мид. сеч. к массе спутника (м^2/кг)
      ! CD - безразм. коэф. сопр. атмосферы (2-2.5)
      ! P - ускорение от сопротивления атмосферы (км/сек^2)
      DATA W/7.292115E-5/ ! Частота вращения Земли (сек^-1)
      R2=X(1)**2+X(2)**2+X(3)**2 
      R=SQRT(R2) ! Радиус-вектор
      CALL DENSIT(R,RHO) ! RHO - плотность атмосферы (кг/м^3)
      VA(1)=(V(1)+W*X(2))
      VA(2)=(V(2)-W*X(1))
      VA(3)= V(3)
      VAIR=VA(1)**2+VA(2)**2+VA(3)**2
      CDAM=-0.5E3*CD*SM
      VAIR=SQRT(VAIR)*CDAM*RHO
      P=VAIR*VA
      END
************************************************************************
      SUBROUTINE DENSIT(R,RHO)
      ! Плотность атмосферы
      PARAMETER (RE=6402E0) ! Средний радиус Земли (км)
      COMMON/NEXCL/NEXCL/NOBJ/N /ATM/IATM_,SM_,CD_,HBURN
      DIMENSION HSC(9),RHOSC(9),OH(9),OHH(9)
      DATA HSC/5E3,9E2,6E2,3E2,15E1,1E2,6E1,2E1,0./
      DATA RHOSC/0.,0.58038E-14,0.11495E-12,0.19019E-10,0.20474E-8,
     *0.54733E-6,0.31655E-3,0.91907E-1,1.2522/
      DATA OH/0.,0.39247E-2,0.14474E-1,0.19885E-1,0.45825E-1,
     *0.17527,0.12378,0.16739,0.90764E-1/
      DATA OHH/0.,0.,0.15127E-4,0.97266E-5,0.10167E-3,0.12870E-2,
     *-0.86999E-3,0.62669E-3,-0.20452E-2/
      H=R-RE ! Высота орбиты спутника
      IF(H<HBURN) NEXCL=N
      I=0
    1 I=I+1
      DH=H-HSC(I) 
      IF(DH) 1,2,2
    2 RHO=RHOSC(I)
      IF(I.EQ.1) GO TO 3
      RHO=RHO*EXP((OHH(I)*DH-OH(I))*DH)
    3 RETURN
      END
************************************************************************
      SUBROUTINE RLTVT(X,V,XS,P)
      ! Релятивистские эффекты
      ! DMU  - гравитационный параметр Земли
      ! X    - геоцентр. вектор положения спутника (км)
      ! V    - вектор скорости спутника отн. геоцентра (км/с)
      ! DMUS - гравитационный параметр Солнца
      ! XS   - гелиоцентр. вектор положения и скорости Земли (км) 
      PARAMETER (RC2=1./299792.458**2)
      PARAMETER (W=9.8E2)
      DIMENSION X(3),V(3),XS(6),P(3),PHI_0(3),PHI_1(3),PHI_2(3)
      COMMON /GR/DMU,GM_,DMUS /RELCUST/IGPR,ISCHW,ILTER
      XS=-XS
      PHI_0=0.; PHI_1=0.; PHI_2=0.

      R2=X(1)**2+X(2)**2+X(3)**2; R=SQRT(R2); R3=R2*R
      V2=V(1)**2+V(2)**2+V(3)**2
      RV=X(1)*V(1)+X(2)*V(2)+X(3)*V(3)

      ! Геодезическая прецессия
      IF(IGPR/=0) THEN
      RS2=XS(1)**2+XS(2)**2+XS(3)**2; RS=SQRT(RS2)
      COEF1=3.*DMUS/RS2/RS
      C1=XS(2)*XS(6)-XS(3)*XS(5)
      C2=XS(3)*XS(4)-XS(1)*XS(6)
      C3=XS(1)*XS(5)-XS(2)*XS(4)
      PHI_0(1)=COEF1*(C2*V(3)-C3*V(2))
      PHI_0(2)=COEF1*(C3*V(1)-C1*V(3))
      PHI_0(3)=COEF1*(C1*V(2)-C2*V(1))
      ENDIF

      ! Шварцшильдовский эффект
      IF(ISCHW/=0) THEN
      COEF1=DMU/R3
      PHI_1=COEF1*((4.*DMU/R-V2)*X(1:3)+4.*RV*X(1:3))
      ENDIF

      ! Эффект Лензе-Терринга
      IF(ILTER/=0) THEN
      COEF1=2.*DMU*W/R3; COEF2=3.*X(3)/R2
      PHI_2(1)=COEF1*(COEF2*(X(2)*V(3)-X(3)*V(2))+V(2))
      PHI_2(2)=COEF1*(COEF2*(X(3)*V(1)-X(1)*V(3))-V(1))
      PHI_2(3)=COEF1*(COEF2*(X(1)*V(2)-X(2)*V(1)))
      ENDIF

      P=(PHI_0+PHI_1+PHI_2)*RC2
      RETURN
      END 
************************************************************************
      SUBROUTINE SUN_LITE(X,V,XS,QPR,SM,SP)
      ! Световое давление и эффект Поинтинга-Робертсона
      ! X,XS,V,VS - геоцентр. коор-ты и скорости 
      ! спутника и Солнца (км,км/с)
      ! QPR - отраж. способность спутника;
      ! SM - отнош-е площади мид. сеч. к массе спутника (м^2/кг);
      PARAMETER (OPI=2./6.283185307179586477)
      PARAMETER (ERAD=6402.,SRAD=696000.) ! Радиусы Земли и Солнца (км)
      PARAMETER (AU2=149597870691E-3**2)  ! Квадрат AU (км^2)
      PARAMETER (C=299792.458)            ! Скорость света (км/с)
      REAL, DIMENSION(3) :: X,V,SP,DXS,DVS
      DIMENSION XS(6)
      
      DELTA2=0.; CPHI=0.; XV=0.
      DXS(1:3)=X(1:3)-XS(1:3); DVS(1:3)=V(1:3)-XS(4:6) ! Гелиоцентр. коор-ты и скорости
      DO I=1,3
      DELTA2=DELTA2+DXS(I)**2
      CPHI=CPHI+DXS(I)*X(I)
      XV=XV+DXS(I)*DVS(I)
      END DO
      
      BETA=4.56E-6*QPR*SM/1E3
      FP=BETA*AU2/DELTA2              ! Коэф. F_P
      R=SQRT(X(1)**2+X(2)**2+X(3)**2) ! Расстояние спутник-Земля
      DELTA=SQRT(DELTA2)              ! Расстояние спутник-Солнце
      PHI=ACOS(CPHI/R/DELTA)          ! Угол напр. от спут. на Землю и Солнце
      GS=ASIN(SRAD/DELTA)             ! Угл. радиус Солнца
      GE=ASIN(ERAD/R)                 ! Угл. радиус Земли
      
      IF(PHI.LT.GE+GS) THEN
      IF(PHI.GT.GE-GS) THEN ! Полутень
      COEF=(GE/GS)**2
      CFS=(GS**2-GE**2+PHI**2)/(2.*GS*PHI)
      CFE=(GE**2-GS**2+PHI**2)/(2.*GE*PHI)
      FS=ACOS(CFS); FE=ACOS(CFE)
      SFS=SQRT(1.-CFS**2); SFE=SQRT(1.-CFE**2)
      PSI=1.-OPI*(FS-CFS*SFS+COEF*(FE-CFE*SFE))
      SP=FP*PSI*(DXS*(1./DELTA-XV/DELTA2/C)-DVS/C)
      ELSE
      SP=0. ! Тень
      END IF
      ELSE
      SP=FP*(DXS*(1./DELTA-XV/DELTA2/C)-DVS/C)
      END IF
      RETURN
      END
************************************************************************
      SUBROUTINE READ405(TSEC,XM,XS)
      ! Фонд DE405
      PARAMETER (RATL=1./(81.30056+1.))
      DIMENSION XM(6), XS(6), XB(6), XE(6), BUF(459)
      COMMON/BUF/BUF /PERT/IMON,ISUN,IGARM_ /TJD0/TJD0,TEPSEC 
      COMMON /TIDE/ITIDE /RELAT/IREL /LIGHT/ILGT
      DATA IOP/0/NRC/-1/TMIN/2436112.5/TMAX/2525008.5/SUT/86400./
      T=TJD0+(TEPSEC+TSEC)/SUT
      IF(T<TMIN.OR.T>TMAX)STOP'Дата за пределами фонда (1957-2200 гг.)'
       IF(IOP==0) THEN
       OPEN(30,FILE='MODTOOLS\\57200EMS.405',FORM='FORMATTED',
     : ACCESS='DIRECT',RECL=80)
       IOP=1
       ENDIF
      NR=(T-TMIN)/32
      NR=NR*153+1
       IF(NR/=NRC) THEN
       NRC=NR
       READ(30,REC=NR,FMT='(3E26.18)') BUF
       IF(T<BUF(1).OR.T>BUF(2)) STOP 'Прочитана не та запись в фонде'
       ENDIF    
      CALL COOR(2,T,XM)   
      IF(ISUN/=0.OR.ITIDE/=0.OR.IREL/=0.OR.ILGT/=0) THEN
      CALL COOR(1,T,XB)
      CALL COOR(3,T,XS)
      IF(IREL/=0.OR.ILGT/=0) THEN
       XE=XB-RATL*XM
       XE=XE-XS
       XS=-XE
       ELSE
       XE(1:3)=XB(1:3)-RATL*XM(1:3)
       XE(1:3)=XE(1:3)-XS(1:3)
       XS(1:3)=-XE(1:3)
       ENDIF
      ENDIF
      END
************************************************************************
      SUBROUTINE CHEB(A,B,T,IST,TC,TCP)
      ! Полиномы Чебышева
      DIMENSION TC(0:15),TCP(0:15)
      COMMON /RELAT/IREL /LIGHT/ILGT
      DATA TWO/2./ONE/1./FOUR/4./
      DLIN=B-A
      TAU=TWO*(T-A)/DLIN-ONE
      TAU2=TAU*TWO
      TC(0)=ONE; TC(1)=TAU
      DO I=2,IST
      TC(I)=TAU2*TC(I-1)-TC(I-2)
      ENDDO
      IF(IREL/=0.OR.ILGT/=0) THEN
       RAT=FOUR/DLIN
       TCP(0)=0.; TCP(1)=TWO/DLIN
       DO I=2,IST
       TCP(I)=RAT*TC(I-1)+TAU2*TCP(I-1)-TCP(I-2)
       ENDDO
       ENDIF
      END
************************************************************************
      SUBROUTINE COOR(I,T,XC)
      ! Вычисление координат и скоростей Луны и Солнца
      DIMENSION TC(0:15), TCP(0:15), XC(6)
      DIMENSION IPOW(3), NPER(3), NCF(3), SUT(3)
      COMMON/BUF/BUF(459) /RELAT/IREL /LIGHT/ILGT
      SAVE TC,TCP
      DATA IPOW/12,12,10/NPER/3,81,393/ZERO/0./NCF/39,39,33/
      DATA SUT/16.,4.,16./ SECSUT/86400./
      NINT=(T-BUF(1))/SUT(I)
      A=BUF(1)+FLOAT(NINT)*SUT(I)
      IST=IPOW(I)+1
      B=A+SUT(I)
      IF(I/=3) CALL CHEB(A,B,T,IPOW(I),TC,TCP)
      JN=NPER(I)+NINT*NCF(I)
      DO K=1,3
      XC(K)=ZERO
      JJ=JN+(K-1)*IST
       DO J=0,IPOW(I)
       XC(K)=XC(K)+TC(J)*BUF(JJ+J)
       ENDDO
      ENDDO    
      IF(IREL/=0.OR.ILGT/=0) THEN
       DO K=1,3
       XC(K+3)=ZERO
       JJ=JN+(K-1)*IST
        DO J=1,IPOW(I)
        XC(K+3)=XC(K+3)+TCP(J)*BUF(JJ+J)
        ENDDO
       ENDDO
       XC(4:6)=XC(4:6)/SECSUT
       ENDIF
      END
************************************************************************
      SUBROUTINE GEOPOT(X,CONX)
      ! Влияние от несферичности Земли
      PARAMETER(NGRM=360, NG=NGRM+2)
      DIMENSION CC(2),CP(2),CQ(2)
      DIMENSION X(3),CONX(3), XR(3),VX(3),GF(3),F(6)
      !****** CC=CC(1)+i*CC(2) CP=CP(1)+i*CP(2) CQ=CQ(1)+i*CQ(2) ******!
      EQUIVALENCE(XR(1),XR1),(XR(2),XR2),(XR(3),XR3)      
      COMMON /UPBORDER/NM,MN /GR/GE,GM_,GS_ /RE/RE /AR/AR(3,3)
      COMMON /CVQXZ/CV(2,NG,NG),Q(NG,NG),QX(NG,NG),QZ(NG,NG),GM  
      DATA  ZERO,ONE,TWO,THREE,C05/0.,1.,2.,3.,0.5/ 

      R2=ZERO
      VV=ZERO
      DO 5 I=1,3
      CONX(I)=ZERO
      VX(I)=ZERO
      GF(I)=ZERO
    5 R2=R2+(X(I)/RE)**2
      R=SQRT(R2)
      SINF=ONE/R
      OR2=ONE/R2
      CV(1,1,1)=SINF ! Функция V_00

      DO 1 I=1,3
      VC=ZERO
      DO 121 J=1,3
  121 VC=VC+AR(I,J)*X(J)
    1 XR(I)=VC*OR2/RE     ! Преобразование координат: XR=AR*X/R^2
                     
      DO I=2,NM
      IM1=I-1
      CP(1:2)=CV(1:2,IM1,IM1)
      SQ2=Q(I,I)
      CQ(1)=SQ2*(XR1*CP(1)-XR2*CP(2))
      CQ(2)=SQ2*(XR1*CP(2)+XR2*CP(1))
      CV(1:2,I,I)=CQ(1:2)   ! Вычисление диагональных функций типа V_NN
      ENDDO
      CV(1,2,1)=Q(2,1)*XR3*SINF ! Шаровая функция V_10
      !*********** CC=CNM+i*SNM CP=CPR+i*CPI CQ=CR1+i*CI1 ***********!
      
      DO 3 I=3,NM
      IM1=I-1
      DO 3 J=1,IM1
      CP(1:2)=CV(1:2,IM1,J)
      CC(1:2)=CV(1:2,IM1-1,J)
      SQ1=Q(I,J)*XR3
      SQ2=Q(J,I)*OR2
      CQ(1)=SQ1*CP(1)-SQ2*CC(1)
      CQ(2)=SQ1*CP(2)-SQ2*CC(2)
    3 CV(1:2,I,J)=CQ(1:2)   ! Вычисление остальных функций V_NM (M<N)

      !********************************************************************!
      !**** Замечание. Коэффициенты C_NM, S_NM находятся в верхнем ********!
      !* треугольнике матрицы CV, значение шаровых функций V_NM в нижнем. *!
      !**** Связь значений коэффициентов и шаровых функций задается *******!
      !******************** следующим образом: ****************************!
      !********************************************************************!

                       ! C_NM -> CV(M+1,N+2)
                       ! V_NM -> CV(N+1,M+1)
                       ! CV(I,J) -> C_{J-2,I-1} (J>I)
                       ! CV(I,J) -> V_{I-1,J-1} (J<I или J=I)
      NM1=NM-1
      DO 4 I=2,NM1
      IP1=I+1
      CC(1:2)=CV(1:2,1,IP1)       ! Коэффициент C_{I-1,0}
      SQ1=QX(I,1)*CC(1)
      SQ2=CC(1)
      SQ3=QZ(I,1)*CC(1)
      CP(1:2)=CV(1:2,IP1,2)        ! Шаровая функция V_{I,1}
      CQ(1:2)=CV(1:2,IP1,1)        ! Шаровая функция V_{I,0}
      CC(1:2)=CV(1:2,I,1)          ! Шаровая функция V_{I-1,0}
      VX(1)=VX(1)-SQ1*CP(1)  ! -QX(I,1)*Real(C_{I-1,0})*Real(V_{I,1})
      VX(2)=VX(2)-SQ1*CP(2)  ! -QX(I,1)*Real(C_{I-1,0})*Imag(V_{I,1})
      VX(3)=VX(3)-SQ3*CQ(1) ! -QZ(I,1)*Real(C_{I-1,0})*Real(V_{I,0})
      !*********** CC=CNM+i*SNM CP=CPR+i*CPI CQ=CR1+i*CI1 ***********!
      VV=VV+SQ2*CC(1)
    4 CONTINUE     ! Цикл 4 - вычисление зональной части разложения:
                   ! VX - геопотенциал
      
      DO 6 I=2,NM1
      IP1=I+1
      MM=I
      DO 16 J=2,MM
      IF(J>MN.AND.I==NM1) GOTO 6  
      JP1=J+1
      JM1=J-1
      CC(1:2)=CV(1:2,J,IP1)              ! Коэффициент C_{I-1,J-1}
      CS=CC(1)                    ! Real(C_{I-1,J-1})
      SS=CC(2)                    ! Imag(C_{I-1,J-1})
    7 CP(1:2)=CV(1:2,IP1,JP1)            ! Функция V_{I,J}
      CC(1:2)=CV(1:2,IP1,J)              ! Функция V_{I,J-1}
      CQ(1:2)=CV(1:2,IP1,JM1)            ! Функция V_{I,J-2}
      SQ1=QX(I,J)*C05
      SQ2=QX(JM1,I)*C05
      SQ3=QZ(I,J)
      CP(1:2)=CP(1:2)*SQ1
      CQ(1:2)=CQ(1:2)*SQ2

      GF(1)=GF(1)+(CS*(CQ(1)-CP(1))+SS*(CQ(2)-CP(2)))
      GF(2)=GF(2)+(SS*(CQ(1)+CP(1))-CS*(CQ(2)+CP(2)))
      GF(3)=GF(3)-SQ3*(CS*CC(1)+SS*CC(2))
      CC(1:2)=CV(1:2,J,IP1)
      CP(1:2)=CV(1:2,I,J)
      VV=VV+CC(1)*CP(1)+CC(2)*CP(2)

   16 CONTINUE
    6 CONTINUE ! Вычисление секториальной и тессеральной части разложения
      
      DO 11 I=1,3
      CNM=ZERO
      DO 123 J=1,3
  123 CNM=CNM+(VX(J)+GF(J))*AR(J,I)
   11 CONX(I)=CNM*GM ! Представление потенциальных сил
      END
************************************************************************
      SUBROUTINE TIDES(XS,XM)
      ! Приливные деформации
      ! XS и XM - геоцентрические координаты Солнца и Луны (км)
      COMMON /GLS/GSM(2) ! Отношение масс Солнца и Луны к массе Земли
      PARAMETER(NGRM=360,NG=NGRM+2)
      PARAMETER (RE=6378.1363)
      COMMON /CVQXZ/CV(2,NG,NG) /AR/AR(3,3)
      COMMON /CVNST/CNST(2:4,0:3),SNST(2:4,0:3)

      REAL, DIMENSION(6)   :: XS,XM
      REAL, DIMENSION(3,2) :: X
      DIMENSION DL(2:3,0:3,1:2),DNORM(2:3,0:3),DKLOVE(2:4,0:3)
      DIMENSION DSLAT(2),DLONG(2),RSM(2),DC(2:4,0:3),DS(2:4,0:3)
      SAVE DNORM,DKLOVE
      DATA ZERO/0./KEY/0/

      IF(KEY.EQ.0) THEN; KEY=1
      DNORM(2,0)=SQRT(5.) ! Коэффициенты нормализации
      DNORM(2,1)=SQRT(5./3.)
      DNORM(2,2)=SQRT(5./3.)/2.
      DNORM(3,0)=SQRT(7.)
      DNORM(3,1)=SQRT(7./6.)
      DNORM(3,2)=SQRT(7./15.)/2.
      DNORM(3,3)=SQRT(7./10.)/6.
      DKLOVE(2,0)=0.29525 ! Коэффициенты Лява
      DKLOVE(2,1)=0.29470
      DKLOVE(2,2)=0.29801
      DKLOVE(3,0)=0.093
      DKLOVE(3,1)=0.093
      DKLOVE(3,2)=0.093
      DKLOVE(3,3)=0.094
      DKLOVE(4,0)=-0.00087
      DKLOVE(4,1)=-0.00079
      DKLOVE(4,2)=-0.00057
      END IF

      X=ZERO
      DO I=1,3 ! Получ-е коорд-т Солнца и Луны в сист., связан. с Землей
      DO J=1,3; X(I,1)=X(I,1)+AR(I,J)*XS(J); END DO
      DO J=1,3; X(I,2)=X(I,2)+AR(I,J)*XM(J); END DO
      END DO

      ! Геоцентрические расстояния до Солнца и Луны
      RSM(1)=SQRT(X(1,1)**2+X(2,1)**2+X(3,1)**2)
      RSM(2)=SQRT(X(1,2)**2+X(2,2)**2+X(3,2)**2)

      ! Синусы широт и долготы Солнца и Луны
      DSLAT(1)=X(3,1)/RSM(1); DSLAT(2)=X(3,2)/RSM(2)
      DLONG(1)=ATAN2(X(2,1),X(1,1))
      DLONG(2)=ATAN2(X(2,2),X(1,2))

      DO J=1,2 ! Нормализованные полиномы Лежандра
      Z=DSLAT(J)
      Z2=Z**2; OMZ2=1.-Z2; SOMZ2=SQRT(OMZ2)
      DL(2,0,J)=DNORM(2,0)*(-1.+3.*Z2)/2.   ! 2_0
      DL(2,1,J)=DNORM(2,1)*(-3.*Z*SOMZ2)    ! 2_1
      DL(2,2,J)=DNORM(2,2)*(3.*OMZ2)        ! 2_2
      DL(3,0,J)=DNORM(3,0)*(Z*(-3.+5.*Z2)/2.)            ! 3_0
      DL(3,1,J)=DNORM(3,1)*(-3.*SOMZ2*(-1.+5.*Z2)/2.)    ! 3_1
      DL(3,2,J)=DNORM(3,2)*(15.*Z*OMZ2)                  ! 3_2
      DL(3,3,J)=DNORM(3,3)*(-15.*OMZ2*SOMZ2)             ! 3_3
      END DO

      DC=ZERO; DS=ZERO

      DO J=1,2 ! Поправки 2 и 3 степени
      DO N=2,3
        COEF1=GSM(J)*(RE/RSM(J))**(N+1)/FLOAT(2*N+1)
        DO M=0,N
          COEF2=COEF1*DKLOVE(N,M)*DL(N,M,J)
          ANGL=FLOAT(M)*DLONG(J)
          DC(N,M)=DC(N,M)+COEF2*COS(ANGL)
          DS(N,M)=DS(N,M)+COEF2*SIN(ANGL)
        END DO
      END DO
      END DO

      DO J=1,2 ! Поправки 4 степени
        COEF1=GSM(J)*(RE/RSM(J))**3/5.
        DO M=0,2
          COEF2=COEF1*DKLOVE(4,M)*DL(2,M,J)
          ANGL=FLOAT(M)*DLONG(J)
          DC(4,M)=DC(4,M)+COEF2*COS(ANGL)
          DS(4,M)=DS(4,M)+COEF2*SIN(ANGL)
        END DO
      END DO

      DO N=2,3; DO M=0,N
      CV(1,M+1,N+2)=CNST(N,M)+DC(N,M)
      CV(2,M+1,N+2)=SNST(N,M)+DS(N,M)
      ENDDO; ENDDO
      DO M=0,2
      CV(1,M+1,6)=CNST(4,M)+DC(4,M)
      CV(2,M+1,6)=SNST(4,M)+DS(4,M)
      ENDDO
      END
************************************************************************
      SUBROUTINE COOREL(KOD,T0,T,X,EL)
      ! Переход от координат к элементам
      PARAMETER(PI=3.1415926535897932384626433832795,C2PI=PI*2.)
      DIMENSION X(0:5),EL(0:5)
      COMMON /GR/DMU,GM_,GS_ /V/V                                              
      DATA EPSDIF/1E-8/ 
      
      C1=X(1)*X(5)-X(2)*X(4)                                           
      C2=X(2)*X(3)-X(0)*X(5)                                           
      C3=X(0)*X(4)-X(1)*X(3)
      CC=C1**2+C2**2+C3**2
      C=SQRT(CC)
      EL(0)=CC/DMU  

      R=SQRT(X(0)**2+X(1)**2+X(2)**2)
      VV=SQRT(X(3)**2+X(4)**2+X(5)**2)
      DMU_R=-DMU/R
      DL1=DMU_R*X(0)+C3*X(4)-C2*X(5)
      DL2=DMU_R*X(1)+C1*X(5)-C3*X(3)
      DL3=DMU_R*X(2)+C2*X(3)-C1*X(4)
      DL=SQRT(DL1**2+DL2**2+DL3**2)
      EL(1)=DL/DMU  

      COI=C3/C
      EL(2)=ACOS(COI)                                                       

      SII=SIN(EL(2))
      SIOM= C1/C/SII
      COOM=-C2/C/SII
      EL(3)=ATAN2(SIOM,COOM)                                             
      IF(EL(3)<0.) EL(3)=EL(3)+C2PI

      IF(ABS(COI)>1E-200) THEN
      SIO=(-DL1*SIOM+DL2*COOM)/DL/COI
      COO=( DL1*COOM+DL2*SIOM)/DL  
      ELSE
      SIO=DL3/DL; COO=DL2/DL/SIOM
      ENDIF
      EL(4)=ATAN2(SIO,COO)         
      IF(EL(4)<0.) EL(4)=EL(4)+C2PI

      IF(ABS(COI)>1E-200) THEN
      SIU=(-X(0)*SIOM+X(1)*COOM)/R/COI
      COU=( X(0)*COOM+X(1)*SIOM)/R
      ELSE
      SIU=X(2)/R; COU=X(1)/R/SIOM
      ENDIF
      U=ATAN2(SIU,COU)
      V=U-EL(4)              

      IF(EL(1)<=1.-EPSDIF.AND.EL(1)>=0.) THEN                                     
      CALL ELLIPSE(KOD,T0,T,EL) 
      ELSEIF(EL(1)>=1.+EPSDIF) THEN
      CALL HYPERBOLA(KOD,T0,T,EL) 
      ELSEIF(EL(1)>1.-EPSDIF.AND.EL(1)<1.+EPSDIF) THEN
      CALL PARABOLA(KOD,T0,T,EL)
      ELSE; STOP 'эксцентриситет < 0'
      ENDIF
      IF(EL(5)<0.) EL(5)=EL(5)+C2PI
      END                                                                         
      FUNCTION TAU(KOD,T0,T)
      COMMON /DM/DM /DNN/DNN
      SELECT CASE(KOD)
      CASE(0); TAU=DM+DNN*(T-T0)
      CASE(1); TAU=T-DM/DNN
      CASE DEFAULT
      WRITE(*,*) 'kod должен быть 0 или 1'
      STOP
      END SELECT
      END
      SUBROUTINE ELLIPSE(KOD,T0,T,EL)
      PARAMETER(PI=3.1415926535897932384626433832795,C2PI=PI*2.)
      DIMENSION EL(0:5)
      COMMON /DM/DM /DMU/DMU /DNN/DNN /V/V
      A=EL(0)/(1.-EL(1)**2) 
      DNN=SQRT(DMU/A**3)    
      E=2.*ATAN(SQRT((1.-EL(1))/(1.+EL(1)))*TAN(V/2.)) 
      DM=E-EL(1)*SIN(E)   
      DM=DM-AINT(DM/C2PI)*C2PI 
      EL(5)=TAU(KOD,T0,T)      
      END
      SUBROUTINE HYPERBOLA(KOD,T0,T,EL)
      PARAMETER(PI=3.1415926535897932384626433832795,PI4=PI/4.)
      DIMENSION EL(0:5)
      COMMON /DM/DM /DMU/DMU /DNN/DNN /V/V
      A=EL(0)/(EL(1)**2-1.)
      DNN=SQRT(DMU/A**3)                                                                         
      F=2.*ATAN(SQRT((EL(1)-1.)/(EL(1)+1.))*TAN(V/2.))
      DM=EL(1)*TAN(F)-LOG(TAN(F/2.+PI4))                                          
      EL(5)=TAU(KOD,T0,T)
      END                                                                      
      SUBROUTINE PARABOLA(KOD,T0,T,EL)
      DIMENSION EL(0:5)
      COMMON /DM/DM /DMU/DMU /DNN/DNN /V/V
      Q=EL(0)/2.                                                        
      DNN=SQRT(DMU/2./Q**3)                                                    
      S=TAN(V/2.)
      DM=S*(1.+S**2/3.)                                                              
      EL(5)=TAU(KOD,T0,T)
      END

      SUBROUTINE DERIVELC(X,HX,DER)
      ! Производные от элементов по координатам
      DIMENSION X(6),DER(6,6),XP(6),XM(6),ELP(6),ELM(6)
      HV=HX*1E-3
      DO I=1,6; DO J=1,6
      XP=X; XM=X
      IF(J<=3) THEN; XP(J)=XP(J)+HX; XM(J)=XM(J)-HX
               ELSE; XP(J)=XP(J)+HV; XM(J)=XM(J)-HV
      ENDIF
      CALL COOREL(0,0.,0.,XP,ELP)
      CALL COOREL(0,0.,0.,XM,ELM)
       SELECT CASE(I)
       CASE(1,2); DIF=ELP(I)-ELM(I)
       CASE(3:6); DIF=ANGDIF(ELP(I),ELM(I))
       END SELECT
      IF(J<=3) THEN; DER(I,J)=DIF/(2.*HX)
               ELSE; DER(I,J)=DIF/(2.*HV)
      ENDIF
      ENDDO; ENDDO
      END
************************************************************************
      SUBROUTINE OBRSING(N,A,OBR)
      ! Обращение матрицы
      DIMENSION A(N,N), OBR(N,N), SIGMA(N), WORK(N), V(N,N), U(N,N),
     :B(N,N), C(N,N)
      CALL SVD(N,N,A,SIGMA,U,V,WORK,.TRUE.,.TRUE.,IERR)
      B=0.
      DO I=1,N; B(I,I)=1./SIGMA(I); ENDDO
      C=0.
      CALL UMATRN(N,0,1,V,B,C)
      CALL UMATRN(N,0,1,C,U,OBR)
      END
************************************************************************
      SUBROUTINE UMATRN(N,KOD1,KOD2,A,B,C)
      DIMENSION C(N,N),A(N,N),B(N,N)
      IF(KOD1/=1.AND.KOD2/=1) THEN
        DO I=1,N; DO J=1,N
        S=0.
        DO K=1,N; S=S+A(I,K)*B(K,J); C(I,J)=S; ENDDO
        ENDDO; ENDDO
      RETURN
      ENDIF

      IF(KOD1==1.AND.KOD2==1) THEN
        DO I=1,N; DO J=1,N
        S=0.
        DO K=1,N; S=S+A(K,I)*B(J,K); C(I,J)=S; ENDDO
        ENDDO; ENDDO
      RETURN
      ENDIF

      IF(KOD1==1.AND.KOD2/=1) THEN
        DO I=1,N; DO J=1,N
        S=0.
        DO K=1,N; S=S+A(K,I)*B(K,J); C(I,J)=S; ENDDO
        ENDDO; ENDDO
      RETURN
      ENDIF

      IF(KOD1/=1.AND.KOD2==1) THEN
        DO I=1,N; DO J=1,N
        S=0.
        DO K=1,N; S=S+A(I,K)*B(J,K); C(I,J)=S; ENDDO
        ENDDO; ENDDO
      ENDIF
      END
************************************************************************
      SUBROUTINE SVD(M,N,A,W,U,V,RV1,MATU,MATV,IERR)
      ! Сингулярное разложение матрицы A
      INTEGER, INTENT(IN) :: M,N
      REAL, INTENT(IN) :: A(M,N)
      DIMENSION W(N), U(M,M), V(N,N), RV1(M)
      LOGICAL, INTENT(IN) :: MATU,MATV
      U=0.
      V=0.
      IERR=0.
      U(1:M,1:N)=A(1:M,1:N)
      G=0.
      SCALE=0.
      ANORM=0.
        DO I=1,N 
        L=I+1
        RV1(I)=SCALE*G
        G=0.
        S=0.
        SCALE=0.
        IF(I>M) GOTO 210
        DO K=I,M; SCALE=SCALE+ABS(U(K,I)); ENDDO
        IF(SCALE==0) GOTO 210
          DO K=I,M 
          U(K,I)=U(K,I)/SCALE
          S=S+U(K,I)**2
          ENDDO
        F=U(I,I)
        G=-SIGN(SQRT(S),F)
        H=F*G-S
        U(I,I)=F-G
        IF(I==N) GOTO 190
          DO J=L,N 
          S=0.
          DO K=I,M; S=S+U(K,I)*U(K,J); ENDDO
          F=S/H
          U(I:M,J)=U(I:M,J)+F*U(I:M,I)
          ENDDO
  190   U(I:M,I)=SCALE*U(I:M,I)
  210   W(I)=SCALE*G
        G=0.
        S=0.
        SCALE=0.
        IF(I>M.OR.I==N) GOTO 290
        DO K=L,N; SCALE=SCALE+ABS(U(I,K)); ENDDO
        IF(SCALE==0) GOTO 290
          DO K=L,N
          U(I,K)=U(I,K)/SCALE
          S=S+U(I,K)**2
          ENDDO
        F=U(I,L)
        G=-SIGN(SQRT(S),F)
        H=F*G-S
        U(I,L)=F-G
        RV1(L:N)=U(I,L:N)/H
        IF(I==M) GOTO 270
          DO J=L,M
          S=0.
          DO K=L,N; S=S+U(J,K)*U(I,K); ENDDO
          U(J,L:N)=U(J,L:N)+S*RV1(L:N)
          ENDDO
  270   U(I,L:N)=SCALE*U(I,L:N)
  290   ANORM=MAX(ANORM,ABS(W(I))+ABS(RV1(I)))
        ENDDO
      IF(.NOT.MATV) GOTO 410
        DO II=1,N
        I=N+1-II
        IF(I==N) GOTO 390
        IF(G==0.) GOTO 360
        V(L:N,I)=(U(I,L:N)/U(I,L))/G
          DO J=L,N 
          S=0.
          DO K=L,N; S=S+U(I,K)*V(K,J); ENDDO
          V(L:N,J)=V(L:N,J)+S*V(L:N,I)
          ENDDO
  360   V(I,L:N)=0.; V(L:N,I)=0.
  390   V(I,I)=1.
        G=RV1(I)
        L=I
        ENDDO
  410 IF(.NOT.MATU) GOTO 510
      MN=N
      IF(M<N) MN=M
        DO II=1,MN
        I=MN+1-II
        L=I+1
        G=W(I)
        IF(I==N) GOTO 430
        U(I,L:N)=0.
  430   IF(G==0.) GOTO 475
        IF(I==MN) GOTO 460
          DO J=L,N
          S=0.
          DO K=L,M; S=S+U(K,I)*U(K,J); ENDDO
          F=(S/U(I,I))/G
          U(I:M,J)=U(I:M,J)+F*U(I:M,I)
          ENDDO
  460   U(I:M,I)=U(I:M,I)/G
        GOTO 490
  475   U(I:M,I)=0.
  490   U(I,I)=U(I,I)+1.
        ENDDO
  510   DO KK=1,N
        K1=N-KK
        K=K1+1
        ITS=0
  520     DO LL=1,K
          L1=K-LL
          L=L1+1
          IF((ABS(RV1(L))+ANORM)==ANORM) GOTO 565
          IF((ABS(W(L1))+ANORM)==ANORM)  GOTO 540
          ENDDO
  540   C=0.
        S=1.
          DO I=L,K
          F=S*RV1(I)
          RV1(I)=C*RV1(I)
          IF((ABS(F)+ANORM)==ANORM) GOTO 565
          G=W(I)
          H=SQRT(F**2+G**2)
          W(I)=H
          C=G/H
          S=-F/H
          IF(.NOT.MATU) CYCLE
            DO J=1,M 
            Y=U(J,L1)
            Z=U(J,I)
            U(J,L1)=Y*C+Z*S
            U(J,I)=-Y*S+Z*C
            ENDDO
          ENDDO
  565   Z=W(K)
        IF(L==K) GOTO 650
        IF(ITS==1000) GOTO 1000
        ITS=ITS+1
        X=W(L)
        Y=W(K1)
        G=RV1(K1)
        H=RV1(K)
        F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2*H*Y)
        G=SQRT(F**2+1)
        F=((X-Z)*(X+Z)+H*(Y/(F+SIGN(G,F))-H))/X
        C=1.
        S=1.
          DO I1=L,K1
          I=I1+1
          G=RV1(I)
          Y=W(I)
          H=S*G
          G=C*G
          Z=SQRT(F**2+H**2)
          RV1(I1)=Z
          C=F/Z
          S=H/Z
          F=X*C+G*S
          G=-X*S+G*C
          H=Y*S
          Y=Y*C
          IF(.NOT.MATV) GOTO 575
            DO J=1,N
            X=V(J,I1)
            Z=V(J,I)
            V(J,I1)=X*C+Z*S
            V(J,I)=-X*S+Z*C
            ENDDO
  575     Z=SQRT(F**2+H**2)
          W(I1)=Z
          IF(Z==0) GOTO 580
          C=F/Z
          S=H/Z
  580     F=C*G+S*Y
          X=-S*G+C*Y
          IF(.NOT.MATU) CYCLE
            DO J=1,M
            Y=U(J,I1)
            Z=U(J,I)
            U(J,I1)=Y*C+Z*S
            U(J,I)=-Y*S+Z*C
            ENDDO
          ENDDO
        RV1(L)=0.
        RV1(K)=F
        W(K)=X
        GOTO 520
  650   IF(Z>=0) CYCLE
        W(K)=-Z
        IF(.NOT.MATV) CYCLE
        V(1:N,K)=-V(1:N,K)
        ENDDO
      RETURN
 1000 IERR=K
      END
************************************************************************
      SUBROUTINE CELTOTER (DATE1,DATE2)
      ! Матрица перехода Cel -> Ter
      DIMENSION RT2C(3,3) 
      COMMON /AR/RC2T(3,3)
      CALL TERTOCEL(DATE1,DATE2,RT2C)
      DO I=1,3; DO J=1,3
      RC2T(I,J)=RT2C(J,I)
      ENDDO; ENDDO
      END 
************************************************************************
      SUBROUTINE TERTOCEL (DATE1,DATE2,RT2C)
      ! Матрица перехода Ter -> Cel
      DIMENSION RPOM(3,3), RBPN(3,3), RT2C(3,3)
      DATA SUT/86400./
      DJ=DATE1+DATE2
      UTA=DJ+UT1MTT(DJ)/SUT
      UTB=0.
      CALL NU2000B (DATE1,DATE2,DPSI,DEPS)
      THETA=GST2000 (UTA,UTB,DATE1,DATE2,DPSI)
      CALL CBPN2000 (DATE1,DATE2,DPSI,DEPS,RBPN)
      CALL POLAR (DATE1,DATE2,RPOM)
      CALL T2C2000 (RPOM,THETA,RBPN,RT2C)
      END
************************************************************************
      SUBROUTINE T2C2000 ( RPOM, THETA, RBPN, RT2C )
      DIMENSION RPOM(3,3), RBPN(3,3), RT2C(3,3), R(3,3)
      ! Polar motion.
      CALL iau_CR ( RPOM, R )
      ! Earth rotation.
      CALL iau_RZ ( -THETA, R )
      ! CIP motion.
      CALL iau_RXR ( RBPN, R, RT2C )
      END
************************************************************************
      SUBROUTINE POLAR(DATE1,DATE2,RPOM)
      ! Arcseconds to radians
      PARAMETER ( DAS2R = 4.848136811095359935899141E-6 )
      PARAMETER ( D2PI = 6.283185307179586476925287)
      DIMENSION RPOM(3,3)
      DATA DM/2400000.5/
      DMJD=DATE1+DATE2-DM
      A = D2PI*(DMJD-52340.)/365.25 
      C = D2PI*(DMJD-52340.)/435.           
      sinA=SIN(A)
      cosA=COS(A)
      sinC=SIN(C)
      cosC=COS(C)
      XP=(.0551-.1087*cosA+.0469*sinA-.0577*cosC+.1302*sinC)*DAS2R   
      YP=(.3270+.0422*cosA+.0994*sinA+.1302*cosC+.0577*sinC)*DAS2R      
      SP=-0.047*(DMJD-51544.5)/36525.*DAS2R
      CALL iau_IR ( RPOM )
      CALL iau_RX ( YP, RPOM )
      CALL iau_RY ( XP, RPOM )
      CALL iau_RZ ( -SP, RPOM )
      END
************************************************************************
      FUNCTION GST2000 ( UTA, UTB, TTA, TTB, DPSI )
      REAL iau_ANP

      ! Greenwich Sidereal Time, IAU 2000.
      GST2000 = iau_ANP ( GMST2000 ( UTA, UTB, TTA, TTB ) +
     :                    EE2000 ( TTA, TTB, DPSI ) )

      END
************************************************************************
      FUNCTION GMST2000 ( UTA, UTB, TTA, TTB )
      ! Arcseconds to radians
      PARAMETER ( DAS2R = 4.848136811095359935899141E-6 )
      ! Reference epoch (J2000), JD
      PARAMETER ( DJ0 = 2451545. )
      ! Days per Julian century
      PARAMETER ( DJC = 36525. )
      REAL iau_ANP

      ! TT Julian centuries since J2000.0.
      T = ( ( TTA-DJ0 ) + TTB ) / DJC
      ! Greenwich Mean Sidereal Time, IAU 2000.
      GMST2000 = iau_ANP ( ERA2000 ( UTA, UTB ) +
     :                        (    0.014506   +
     :                        ( 4612.15739966 +
     :                        (  + 1.39667721 +
     :                        (  - 0.00009344 +
     :                        (  + 0.00001882 )
     :                                 * T ) * T ) * T ) * T ) * DAS2R )

      END
************************************************************************
      FUNCTION ERA2000 ( DJ1, DJ2 )
      ! 2Pi
      PARAMETER ( D2PI = 6.283185307179586476925287)
      ! Reference epoch (J2000), JD
      PARAMETER ( DJ0 = 2451545. )
      REAL iau_ANP

      ! Days since fundamental epoch.
      IF ( DJ1 .LT. DJ2 ) THEN
         D1 = DJ1
         D2 = DJ2
      ELSE
         D1 = DJ2
         D2 = DJ1
      END IF
      T = D1 + ( D2-DJ0 )

      ! Fractional part of T (days).
      F = MOD ( D1, 1. ) + MOD ( D2, 1. )

      ! Earth rotation angle at this UT1.
      ERA2000 = iau_ANP ( D2PI * ( F + 0.7790572732640
     :                               + 0.00273781191135448 * T ) )

      END
************************************************************************
      FUNCTION EE2000 ( DATE1, DATE2, DPSI )
      ! Arcseconds to radians
      PARAMETER ( DAS2R = 4.848136811095359935899141E-6 )
      ! Reference epoch (J2000), JD
      PARAMETER ( DJ0 = 2451545. )
      ! Days per Julian century
      PARAMETER ( DJC = 36525. )
      ! J2000 obliquity (Lieske et al. 1977)
      PARAMETER ( EPS0 = 84381.448 * DAS2R )

      ! Interval between fundamental epoch J2000.0 and given date (JC).
      T = ( ( DATE1-DJ0 ) + DATE2 ) / DJC

      ! Mean obliquity from Chapter 5, expression (32).
      EPSA = EPS0 + (  -46.8402 +
     :              (   -0.00059 +
     :              (    0.001813 ) * T ) * T ) * T * DAS2R

      ! Equation of the equinoxes.
      EE2000 = DPSI * COS(EPSA) + EECT2000 ( DATE1, DATE2 )

      END
************************************************************************
      FUNCTION EECT2000 ( DATE1, DATE2 )
      ! 2Pi
      PARAMETER ( D2PI = 6.283185307179586476925287 )
      ! Arcseconds to radians
      PARAMETER ( DAS2R = 4.848136811095359935899141E-6 )
      ! Reference epoch (J2000), JD
      PARAMETER ( DJ0 = 2451545. )
      ! Days per Julian century
      PARAMETER ( DJC = 36525. )
      REAL iau_ANPM
      ! Fundamental arguments
      DIMENSION FA(14)

*  -----------------------------------------
*  The series for the EE complementary terms
*  -----------------------------------------

      ! Number of terms in the series
      INTEGER NE0, NE1
      PARAMETER ( NE0=  33, NE1=  1 )
      ! Coefficients of l,l',F,D,Om,LMe,LVe,LE,LMa,LJu,LSa,LU,LN,pA
      INTEGER KE0 ( 14, NE0 ),
     :        KE1 ( 14, NE1 )
      ! Sine and cosine coefficients
      DIMENSION SE0 ( 2, NE0 ),
     :                 SE1 ( 2, NE1 )
      ! Argument coefficients for t^0
      DATA ( ( KE0(I,J), I=1,14), J =    1,   10 ) /
     :  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA ( ( KE0(I,J), I=1,14), J =   11,   20 ) /
     :  1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  4, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA ( ( KE0(I,J), I=1,14), J =   21,   30 ) /
     :  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1, -2,  2, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1, -2,  2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0, -1,
     :  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  2,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  0, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  4, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA ( ( KE0(I,J), I=1,14), J =   31,  NE0 ) /
     :  0,  0,  2, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0, -2,  0, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      ! Argument coefficients for t^1
      DATA ( ( KE1(I,J), I=1,14), J =    1,  NE1 ) /
     :  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      ! Sine and cosine coefficients for t^0
      DATA ( ( SE0(I,J), I=1,2), J =    1,   10 ) /
     :            +2640.96E-6,          -0.39E-6,
     :              +63.52E-6,          -0.02E-6,
     :              +11.75E-6,          +0.01E-6,
     :              +11.21E-6,          +0.01E-6,
     :               -4.55E-6,          +0.00E-6,
     :               +2.02E-6,          +0.00E-6,
     :               +1.98E-6,          +0.00E-6,
     :               -1.72E-6,          +0.00E-6,
     :               -1.41E-6,          -0.01E-6,
     :               -1.26E-6,          -0.01E-6 /
      DATA ( ( SE0(I,J), I=1,2), J =   11,   20 ) /
     :               -0.63E-6,          +0.00E-6,
     :               -0.63E-6,          +0.00E-6,
     :               +0.46E-6,          +0.00E-6,
     :               +0.45E-6,          +0.00E-6,
     :               +0.36E-6,          +0.00E-6,
     :               -0.24E-6,          -0.12E-6,
     :               +0.32E-6,          +0.00E-6,
     :               +0.28E-6,          +0.00E-6,
     :               +0.27E-6,          +0.00E-6,
     :               +0.26E-6,          +0.00E-6 /
      DATA ( ( SE0(I,J), I=1,2), J =   21,   30 ) /
     :               -0.21E-6,          +0.00E-6,
     :               +0.19E-6,          +0.00E-6,
     :               +0.18E-6,          +0.00E-6,
     :               -0.10E-6,          +0.05E-6,
     :               +0.15E-6,          +0.00E-6,
     :               -0.14E-6,          +0.00E-6,
     :               +0.14E-6,          +0.00E-6,
     :               -0.14E-6,          +0.00E-6,
     :               +0.14E-6,          +0.00E-6,
     :               +0.13E-6,          +0.00E-6 /
      DATA ( ( SE0(I,J), I=1,2), J =   31,  NE0 ) /
     :               -0.11E-6,          +0.00E-6,
     :               +0.11E-6,          +0.00E-6,
     :               +0.11E-6,          +0.00E-6 /
      ! Sine and cosine coefficients for t^1
      DATA ( ( SE1(I,J), I=1,2), J =    1,  NE1 ) /
     :               -0.87E-6,          +0.00E-6 /

      ! Interval between fundamental epoch J2000.0 and current date (JC).
      T = ( ( DATE1-DJ0 ) + DATE2 ) / DJC

      ! Fundamental Arguments (from IERS Conventions 2000)

      ! Mean Anomaly of the Moon.
      FA(1) = iau_ANPM ( ( 485868.249036 +
     :                   ( 715923.2178 +
     :                   (     31.8792 +
     :                   (      0.051635 +
     :                   (     -0.00024470 )
     :                   * T ) * T ) * T ) * T ) * DAS2R
     :                   + MOD ( 1325.*T, 1. ) * D2PI )

      ! Mean Anomaly of the Sun.
      FA(2) = iau_ANPM ( ( 1287104.793048 +
     :                   ( 1292581.0481 +
     :                   (      -0.5532 +
     :                   (      +0.000136 +
     :                   (      -0.00001149 )
     :                   * T ) * T ) * T ) * T ) * DAS2R
     :                   + MOD ( 99.*T, 1. ) * D2PI )

      ! Mean Longitude of the Moon minus Mean Longitude of the Ascending
      ! Node of the Moon.
      FA(3) = iau_ANPM ( (  335779.526232 +
     :                   (  295262.8478 +
     :                   (     -12.7512 +
     :                   (      -0.001037 +
     :                   (       0.00000417 )
     :                   * T ) * T ) * T ) * T ) * DAS2R
     :                   + MOD ( 1342.*T, 1. ) * D2PI )

      ! Mean Elongation of the Moon from the Sun.
      FA(4) = iau_ANPM ( ( 1072260.703692 +
     :                   ( 1105601.2090 +
     :                   (      -6.3706 +
     :                   (       0.006593 +
     :                   (      -0.00003169 )
     :                   * T ) * T ) * T ) * T ) * DAS2R
     :                   + MOD ( 1236.*T, 1. ) * D2PI )

      ! Mean Longitude of the Ascending Node of the Moon.
      FA(5) = iau_ANPM ( (  450160.398036 +
     :                   ( -482890.5431 +
     :                   (       7.4722 +
     :                   (       0.007702 +
     :                   (      -0.00005939 )
     :                   * T ) * T ) * T ) * T ) * DAS2R
     :                   + MOD ( -5.*T, 1. ) * D2PI )

      FA( 6) = iau_ANPM ( 4.402608842 + 2608.7903141574 * T )
      FA( 7) = iau_ANPM ( 3.176146697 + 1021.3285546211 * T )
      FA( 8) = iau_ANPM ( 1.753470314 +  628.3075849991 * T )
      FA( 9) = iau_ANPM ( 6.203480913 +  334.0612426700 * T )
      FA(10) = iau_ANPM ( 0.599546497 +   52.9690962641 * T )
      FA(11) = iau_ANPM ( 0.874016757 +   21.3299104960 * T )
      FA(12) = iau_ANPM ( 5.481293872 +    7.4781598567 * T )
      FA(13) = iau_ANPM ( 5.311886287 +    3.8133035638 * T )
      FA(14) =          ( 0.024381750 +   0.00000538691 * T ) * T

      ! Evaluate the EE complementary terms.
      S0 = 0.
      S1 = 0.

      DO I = NE0,1,-1
         A = 0.
         DO J=1,14
            A = A + FLOAT(KE0(J,I))*FA(J)
         END DO
         S0 = S0 + ( SE0(1,I)*SIN(A) + SE0(2,I)*COS(A) )
      END DO
      DO I = NE1,1,-1
         A = 0.
         DO J=1,14
            A = A + FLOAT(KE1(J,I))*FA(J)
         END DO
         S1 = S1 + ( SE1(1,I)*SIN(A) + SE1(2,I)*COS(A) )
      END DO
      EECT2000 = ( S0 + S1 * T ) * DAS2R
      END
************************************************************************
      FUNCTION UT1MTT(DJD)
      UT1MTT=UT1CTAI(DJD,2)-32.184
      END
************************************************************************
      FUNCTION UTCMTT(DJD)
      UT1_UTC=UT1CTAI(DJD,1)
      UT1_TAI=UT1CTAI(DJD,2)
      UTCMTT=UT1_TAI-UT1_UTC-32.184
      END
************************************************************************
      FUNCTION UT1CTAI(DJD,K12)
      PARAMETER(NTI=30000)
      COMMON /UTTAI/NTAI,DMJD(NTI),DUT(2,NTI)
      DATA DM/2400000.5/
C      WRITE(*,*)'NTAI=',NTAI
C      DO I=1,NTAI
C      WRITE(*,'(F10.0,2F15.7)') DMJD(I),DUT(1,I),DUT(2,I)
C      ENDDO
C      STOP

      T=DJD-DM
      IF(T<=DMJD(1)) THEN
        UT1CTAI=DUT(K12,1)
        RETURN
        ENDIF
      IF(T>=DMJD(NTAI)) THEN
        UT1CTAI=DUT(K12,NTAI)
        RETURN
        ENDIF
      NR=T-DMJD(1)+1
      COEF=DUT(K12,NR+1)-DUT(K12,NR)
      UT1CTAI=DUT(K12,NR)+COEF*(T-DMJD(NR))
      END
************************************************************************
      SUBROUTINE NU2000B ( DATE1, DATE2, DPSI, DEPS )
      ! Arcseconds to radians
      PARAMETER ( DAS2R = 4.848136811095359935899141E-6 )
      ! Milliarcseconds to radians
      PARAMETER ( DMAS2R = DAS2R / 1E3 )
      ! Arc seconds in a full circle
      PARAMETER ( TURNAS = 1296000. )
      ! 2Pi
      PARAMETER ( D2PI = 6.283185307179586476925287 )
      ! Units of 0.1 microarcsecond to radians
      PARAMETER ( U2R = DAS2R/1E7 )
      ! Reference epoch (J2000), JD
      PARAMETER ( DJ0 = 2451545. )
      ! Days per Julian century
      PARAMETER ( DJC = 36525. )

      ! Miscellaneous

*  -------------------------
*  Luni-Solar nutation model
*  -------------------------

      ! Number of terms in the luni-solar nutation model
      PARAMETER ( NLS = 77 )
      ! Coefficients for fundamental arguments
      INTEGER NALS(5,NLS)
      ! Longitude and obliquity coefficients
      DIMENSION CLS(6,NLS)

*  ------------------
*  Planetary nutation (radians)
*  ------------------

      PARAMETER ( DPPLAN = - 0.135 * DMAS2R,
     :            DEPLAN = + 0.388 * DMAS2R )
      ! n.b.  The above fixed terms account for the omission of the
      ! long-period planetary terms in the truncated model.

*  ----------------------------------------
*  Tables of argument and term coefficients
*  ----------------------------------------

      ! Luni-Solar argument multipliers:
      !              L     L'    F     D     Om

      DATA ( ( NALS(I,J), I=1,5 ), J= 1,10 ) /
     :          0,    0,    0,    0,    1,
     :          0,    0,    2,   -2,    2,
     :          0,    0,    2,    0,    2,
     :          0,    0,    0,    0,    2,
     :          0,    1,    0,    0,    0,
     :          0,    1,    2,   -2,    2,
     :          1,    0,    0,    0,    0,
     :          0,    0,    2,    0,    1,
     :          1,    0,    2,    0,    2,
     :          0,   -1,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=11,20 ) /
     :          0,    0,    2,   -2,    1,
     :         -1,    0,    2,    0,    2,
     :         -1,    0,    0,    2,    0,
     :          1,    0,    0,    0,    1,
     :         -1,    0,    0,    0,    1,
     :         -1,    0,    2,    2,    2,
     :          1,    0,    2,    0,    1,
     :         -2,    0,    2,    0,    1,
     :          0,    0,    0,    2,    0,
     :          0,    0,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=21,30 ) /
     :          0,   -2,    2,   -2,    2,
     :         -2,    0,    0,    2,    0,
     :          2,    0,    2,    0,    2,
     :          1,    0,    2,   -2,    2,
     :         -1,    0,    2,    0,    1,
     :          2,    0,    0,    0,    0,
     :          0,    0,    2,    0,    0,
     :          0,    1,    0,    0,    1,
     :         -1,    0,    0,    2,    1,
     :          0,    2,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=31,40 ) /
     :          0,    0,   -2,    2,    0,
     :          1,    0,    0,   -2,    1,
     :          0,   -1,    0,    0,    1,
     :         -1,    0,    2,    2,    1,
     :          0,    2,    0,    0,    0,
     :          1,    0,    2,    2,    2,
     :         -2,    0,    2,    0,    0,
     :          0,    1,    2,    0,    2,
     :          0,    0,    2,    2,    1,
     :          0,   -1,    2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=41,50 ) /
     :          0,    0,    0,    2,    1,
     :          1,    0,    2,   -2,    1,
     :          2,    0,    2,   -2,    2,
     :         -2,    0,    0,    2,    1,
     :          2,    0,    2,    0,    1,
     :          0,   -1,    2,   -2,    1,
     :          0,    0,    0,   -2,    1,
     :         -1,   -1,    0,    2,    0,
     :          2,    0,    0,   -2,    1,
     :          1,    0,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=51,60 ) /
     :          0,    1,    2,   -2,    1,
     :          1,   -1,    0,    0,    0,
     :         -2,    0,    2,    0,    2,
     :          3,    0,    2,    0,    2,
     :          0,   -1,    0,    2,    0,
     :          1,   -1,    2,    0,    2,
     :          0,    0,    0,    1,    0,
     :         -1,   -1,    2,    2,    2,
     :         -1,    0,    2,    0,    0,
     :          0,   -1,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=61,70 ) /
     :         -2,    0,    0,    0,    1,
     :          1,    1,    2,    0,    2,
     :          2,    0,    0,    0,    1,
     :         -1,    1,    0,    1,    0,
     :          1,    1,    0,    0,    0,
     :          1,    0,    2,    0,    0,
     :         -1,    0,    2,   -2,    1,
     :          1,    0,    0,    0,    2,
     :         -1,    0,    0,    1,    0,
     :          0,    0,    2,    1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=71,77 ) /
     :         -1,    0,    2,    4,    2,
     :         -1,    1,    0,    1,    1,
     :          0,   -2,    2,   -2,    1,
     :          1,    0,    2,    2,    1,
     :         -2,    0,    2,    2,    2,
     :         -1,    0,    0,    0,    2,
     :          1,    1,    2,   -2,    2 /

      ! Luni-Solar nutation coefficients, unit 1e-7 arcsec:
      ! longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
      ! Each row of coefficients in CLS belongs with the corresponding row of
      ! fundamental-argument multipliers in NALS.

      DATA ( ( CLS(I,J), I=1,6 ), J= 1,10 ) /
     : -172064161., -174666.,  33386., 92052331.,  9086., 15377.,
     :  -13170906.,   -1675., -13696.,  5730336., -3015., -4587.,
     :   -2276413.,    -234.,   2796.,   978459.,  -485.,  1374.,
     :    2074554.,     207.,   -698.,  -897492.,   470.,  -291.,
     :    1475877.,   -3633.,  11817.,    73871.,  -184., -1924.,
     :    -516821.,    1226.,   -524.,   224386.,  -677.,  -174.,
     :     711159.,      73.,   -872.,    -6750.,     0.,   358.,
     :    -387298.,    -367.,    380.,   200728.,    18.,   318.,
     :    -301461.,     -36.,    816.,   129025.,   -63.,   367.,
     :     215829.,    -494.,    111.,   -95929.,   299.,   132. /
      DATA ( ( CLS(I,J), I=1,6 ), J=11,20 ) /
     :     128227.,     137.,    181.,   -68982.,    -9.,    39.,
     :     123457.,      11.,     19.,   -53311.,    32.,    -4.,
     :     156994.,      10.,   -168.,    -1235.,     0.,    82.,
     :      63110.,      63.,     27.,   -33228.,     0.,    -9.,
     :     -57976.,     -63.,   -189.,    31429.,     0.,   -75.,
     :     -59641.,     -11.,    149.,    25543.,   -11.,    66.,
     :     -51613.,     -42.,    129.,    26366.,     0.,    78.,
     :      45893.,      50.,     31.,   -24236.,   -10.,    20.,
     :      63384.,      11.,   -150.,    -1220.,     0.,    29.,
     :     -38571.,      -1.,    158.,    16452.,   -11.,    68. /
      DATA ( ( CLS(I,J), I=1,6 ), J=21,30 ) /
     :      32481.,       0.,      0.,   -13870.,     0.,     0.,
     :     -47722.,       0.,    -18.,      477.,     0.,   -25.,
     :     -31046.,      -1.,    131.,    13238.,   -11.,    59.,
     :      28593.,       0.,     -1.,   -12338.,    10.,    -3.,
     :      20441.,      21.,     10.,   -10758.,     0.,    -3.,
     :      29243.,       0.,    -74.,     -609.,     0.,    13.,
     :      25887.,       0.,    -66.,     -550.,     0.,    11.,
     :     -14053.,     -25.,     79.,     8551.,    -2.,   -45.,
     :      15164.,      10.,     11.,    -8001.,     0.,    -1.,
     :     -15794.,      72.,    -16.,     6850.,   -42.,    -5. /
      DATA ( ( CLS(I,J), I=1,6 ), J=31,40 ) /
     :      21783.,       0.,     13.,     -167.,     0.,    13.,
     :     -12873.,     -10.,    -37.,     6953.,     0.,   -14.,
     :     -12654.,      11.,     63.,     6415.,     0.,    26.,
     :     -10204.,       0.,     25.,     5222.,     0.,    15.,
     :      16707.,     -85.,    -10.,      168.,    -1.,    10.,
     :      -7691.,       0.,     44.,     3268.,     0.,    19.,
     :     -11024.,       0.,    -14.,      104.,     0.,     2.,
     :       7566.,     -21.,    -11.,    -3250.,     0.,    -5.,
     :      -6637.,     -11.,     25.,     3353.,     0.,    14.,
     :      -7141.,      21.,      8.,     3070.,     0.,     4. /
      DATA ( ( CLS(I,J), I=1,6 ), J=41,50 ) /
     :      -6302.,     -11.,      2.,     3272.,     0.,     4.,
     :       5800.,      10.,      2.,    -3045.,     0.,    -1.,
     :       6443.,       0.,     -7.,    -2768.,     0.,    -4.,
     :      -5774.,     -11.,    -15.,     3041.,     0.,    -5.,
     :      -5350.,       0.,     21.,     2695.,     0.,    12.,
     :      -4752.,     -11.,     -3.,     2719.,     0.,    -3.,
     :      -4940.,     -11.,    -21.,     2720.,     0.,    -9.,
     :       7350.,       0.,     -8.,      -51.,     0.,     4.,
     :       4065.,       0.,      6.,    -2206.,     0.,     1.,
     :       6579.,       0.,    -24.,     -199.,     0.,     2. /
      DATA ( ( CLS(I,J), I=1,6 ), J=51,60 ) /
     :       3579.,       0.,      5.,    -1900.,     0.,     1.,
     :       4725.,       0.,     -6.,      -41.,     0.,     3.,
     :      -3075.,       0.,     -2.,     1313.,     0.,    -1.,
     :      -2904.,       0.,     15.,     1233.,     0.,     7.,
     :       4348.,       0.,    -10.,      -81.,     0.,     2.,
     :      -2878.,       0.,      8.,     1232.,     0.,     4.,
     :      -4230.,       0.,      5.,      -20.,     0.,    -2.,
     :      -2819.,       0.,      7.,     1207.,     0.,     3.,
     :      -4056.,       0.,      5.,       40.,     0.,    -2.,
     :      -2647.,       0.,     11.,     1129.,     0.,     5. /
      DATA ( ( CLS(I,J), I=1,6 ), J=61,70 ) /
     :      -2294.,       0.,    -10.,     1266.,     0.,    -4.,
     :       2481.,       0.,     -7.,    -1062.,     0.,    -3.,
     :       2179.,       0.,     -2.,    -1129.,     0.,    -2.,
     :       3276.,       0.,      1.,       -9.,     0.,     0.,
     :      -3389.,       0.,      5.,       35.,     0.,    -2.,
     :       3339.,       0.,    -13.,     -107.,     0.,     1.,
     :      -1987.,       0.,     -6.,     1073.,     0.,    -2.,
     :      -1981.,       0.,      0.,      854.,     0.,     0.,
     :       4026.,       0.,   -353.,     -553.,     0.,  -139.,
     :       1660.,       0.,     -5.,     -710.,     0.,    -2. /
      DATA ( ( CLS(I,J), I=1,6 ), J=71,77 ) /
     :      -1521.,       0.,      9.,      647.,     0.,     4.,
     :       1314.,       0.,      0.,     -700.,     0.,     0.,
     :      -1283.,       0.,      0.,      672.,     0.,     0.,
     :      -1331.,       0.,      8.,      663.,     0.,     4.,
     :       1383.,       0.,     -2.,     -594.,     0.,    -2.,
     :       1405.,       0.,      4.,     -610.,     0.,     2.,
     :       1290.,       0.,      0.,     -556.,     0.,     0. /
                                                          
      ! Interval between fundamental epoch J2000.0 and given date (JC).
      T = ( ( DATE1-DJ0 ) + DATE2 ) / DJC

*  -------------------
*  Luni-Solar Nutation
*  -------------------

      ! Fundamental (Delaunay) arguments from Simon et al. (1994)

      ! Mean anomaly of the Moon.
      EL  = MOD (         485868.249036 +
     :            T * 1717915923.2178, TURNAS ) * DAS2R

      ! Mean anomaly of the Sun.
      ELP = MOD (        1287104.79305 +
     :            T *  129596581.0481, TURNAS ) * DAS2R

      ! Mean argument of the latitude of the Moon.
      F   = MOD (         335779.526232 +
     :            T * 1739527262.8478, TURNAS ) * DAS2R

      ! Mean elongation of the Moon from the Sun.
      D   = MOD (        1072260.70369 +
     :            T * 1602961601.2090, TURNAS ) * DAS2R

      ! Mean longitude of the ascending node of the Moon.
      OM  = MOD (         450160.398036 -
     :            T *    6962890.5431, TURNAS ) * DAS2R

      ! Initialize the nutation values.
      DP = 0.
      DE = 0.

      ! Summation of luni-solar nutation series (in reverse order).
      DO 100 I = NLS, 1, -1

      !    Argument and functions.
         ARG = MOD ( FLOAT ( NALS(1,I) ) * EL  +
     :               FLOAT ( NALS(2,I) ) * ELP +
     :               FLOAT ( NALS(3,I) ) * F   +
     :               FLOAT ( NALS(4,I) ) * D   +
     :               FLOAT ( NALS(5,I) ) * OM, D2PI )
         SARG = SIN(ARG)
         CARG = COS(ARG)

      !    Term.
         DP = DP + ( CLS(1,I) + CLS(2,I) * T ) * SARG
     :           +   CLS(3,I)                  * CARG
         DE = DE + ( CLS(4,I) + CLS(5,I) * T ) * CARG
     :           +   CLS(6,I)                  * SARG

 100  CONTINUE

      ! Convert from 0.1 microarcsec units to radians.
      DPSILS = DP * U2R
      DEPSLS = DE * U2R

*  ------------------
*  Planetary Nutation
*  ------------------

      ! Fixed terms to allow for long-period nutation.
      DPSIPL = DPPLAN
      DEPSPL = DEPLAN

*  -----
*  Total
*  -----

      ! Add planetary and luni-solar components.
      DPSI = DPSIPL + DPSILS
      DEPS = DEPSPL + DEPSLS

      END
************************************************************************
      SUBROUTINE CBPN2000 ( DATE1, DATE2, DPSI, DEPS, RBPNC )
      DIMENSION RBPNC(3,3)
      ! Arcseconds to radians
      PARAMETER ( DAS2R = 4.848136811095359935899141E-6 )
      ! Reference epoch (J2000), JD
      PARAMETER ( DJ0 = 2451545. )
      ! Days per Julian century
      PARAMETER ( DJC = 36525. )
      ! J2000 obliquity (Lieske et al. 1977)
      PARAMETER ( EPS0 = 84381.448 * DAS2R )
      ! The ICRS RA of the J2000 equinox (Chapront et al., 2002)
      PARAMETER ( DRA0 = -0.0146 * DAS2R )
      ! The precession and obliquity corrections (radians per century)
      PARAMETER ( PRECOR = -0.29965 * DAS2R,
     :            OBLCOR = -0.02524 * DAS2R )
      ! The frame bias corrections in longitude and obliquity
      PARAMETER ( DPSIBI = -0.041775 * DAS2R,
     :            DEPSBI = -0.0068192 * DAS2R )

      ! Interval between fundamental epoch J2000.0 and given date (JC).
      T = ( ( DATE1-DJ0 ) + DATE2 ) / DJC

      ! Precession rate contributions with respect to IAU 1976/80.
      DPSIPR = PRECOR * T
      DEPSPR = OBLCOR * T

      ! IAU 1980 mean obliquity of date.
      EPSA80 = EPS0 + (  -46.8150 +
     :                (   -0.00059 +
     :                (    0.001813 ) * T ) * T ) * T * DAS2R

      ! Precession angles (Lieske et al. 1977)
      PSIA77 =        ( 5038.7784 +
     :                (   -1.07259 +
     :                (   -0.001147 ) * T ) * T ) * T * DAS2R
      OMA77  = EPS0 + (
     :                (    0.05127 +
     :                (   -0.007726 ) * T ) * T ) * T * DAS2R
      CHIA   =        (   10.5526 +
     :                (   -2.38064 +
     :                (   -0.001125 ) * T ) * T ) * T * DAS2R

      ! Apply IAU 2000A precession corrections.
      PSIA = PSIA77 + DPSIPR
      OMA  = OMA77  + DEPSPR
      EPSA = EPSA80 + DEPSPR

      ! Initialize the true-to-celestial matrix.
      CALL iau_IR ( RBPNC )

      ! Remove IAU 2000A nutation (pure: luni-solar and planetary).
      CALL iau_RX ( EPSA+DEPS, RBPNC )
      CALL iau_RZ ( DPSI, RBPNC )
      CALL iau_RX ( -EPSA, RBPNC )

      ! Remove precession (Lieske et al. 1977 plus corrections).
      CALL iau_RZ ( -CHIA, RBPNC )
      CALL iau_RX ( OMA, RBPNC )
      CALL iau_RZ ( PSIA, RBPNC )
      CALL iau_RX ( -EPS0, RBPNC )

      ! Remove frame bias.
      CALL iau_RX ( DEPSBI, RBPNC )
      CALL iau_RY ( -DPSIBI*SIN(EPS0), RBPNC )
      CALL iau_RZ ( -DRA0, RBPNC )

      END
************************************************************************
      SUBROUTINE iau_ZR ( R )
      DIMENSION R(3,3)
      DO 2 J=1,3
         DO 1 I=1,3
            R(I,J) = 0.
 1       CONTINUE
 2    CONTINUE
      END
************************************************************************
      SUBROUTINE iau_IR ( R )
      DIMENSION R(3,3)
      CALL iau_ZR ( R )
      DO 1 I=1,3
         R(I,I) = 1.
 1    CONTINUE
      END
************************************************************************
      SUBROUTINE iau_RXR ( A, B, ATB )
      ! Multiply two r-matrices.
      DIMENSION A(3,3), B(3,3), ATB(3,3), WM(3,3)
      DO 3 I=1,3
         DO 2 J=1,3
            W = 0.
            DO 1 K=1,3
               W = W + A(I,K)*B(K,J)
 1          CONTINUE
            WM(I,J) = W
 2       CONTINUE
 3    CONTINUE
      CALL iau_CR ( WM, ATB )
      END
************************************************************************
      SUBROUTINE iau_CR ( R, C )
      ! Copy an r-matrix.
      DIMENSION R(3,3), C(3,3)
      DO 1 I=1,3
         CALL iau_CP ( R(1,I), C(1,I) )
 1    CONTINUE
      END
************************************************************************
      SUBROUTINE iau_CP ( P, C )
      ! Copy a p-vector.
      DIMENSION P(3), C(3)
      DO 1 I=1,3
         C(I) = P(I)
 1    CONTINUE
      END
************************************************************************
      SUBROUTINE iau_RX ( PHI, R )
      DIMENSION R(3,3), A(3,3), W(3,3)
      ! Matrix representing new rotation.
      S = SIN(PHI)
      C = COS(PHI)
      CALL iau_IR ( A )
      A(2,2) = C
      A(3,2) = -S
      A(2,3) = S
      A(3,3) = C
      ! Rotate.
      CALL iau_RXR ( A, R, W )
      ! Return result.
      CALL iau_CR ( W, R )
      END
************************************************************************
      SUBROUTINE iau_RY ( THETA, R )
      DIMENSION R(3,3), A(3,3), W(3,3)
      ! Matrix representing new rotation.
      S = SIN(THETA)
      C = COS(THETA)
      CALL iau_IR ( A )
      A(1,1) = C
      A(3,1) = S
      A(1,3) = -S
      A(3,3) = C
      ! Rotate.
      CALL iau_RXR ( A, R, W )
      ! Return result.
      CALL iau_CR ( W, R )
      END
************************************************************************
      SUBROUTINE iau_RZ ( PSI, R )
      DIMENSION R(3,3), A(3,3), W(3,3)
      ! Matrix representing new rotation.
      S = SIN(PSI)
      C = COS(PSI)
      CALL iau_IR ( A )
      A(1,1) = C
      A(2,1) = -S
      A(1,2) = S
      A(2,2) = C
      ! Rotate.
      CALL iau_RXR ( A, R, W )
      ! Return result.
      CALL iau_CR ( W, R )
      END
************************************************************************
      REAL FUNCTION iau_ANP ( A )
      ! Normalize angle into the range 0 <= A < 2pi.
      ! 2Pi
      PARAMETER ( D2PI = 6.283185307179586476925287 )
      W = MOD(A,D2PI)
      IF ( W .LT. 0. ) W = W + D2PI
      iau_ANP = W
      END
************************************************************************
      REAL FUNCTION iau_ANPM ( A )
      ! Normalize angle into the range -pi <= A < +pi.
      ! Pi
      PARAMETER ( DPI = 3.141592653589793238462643 )
      ! 2Pi
      PARAMETER ( D2PI = 6.283185307179586476925287 )
      W = MOD(A,D2PI)
      IF ( ABS(W) .GE. DPI ) W = W - SIGN(D2PI,A)
      iau_ANPM = W
      END
************************************************************************
      SUBROUTINE INITCONST(NN)
      PARAMETER(NGRM=360,NG=NGRM+2,NTI=30000)
      COMMON /NT/NT /GR/GE,GM,GS /RE/RE
      COMMON /UTTAI/NTAI,DMJD(NTI),DUT(2,NTI) /UPBORDER/NM,MN
      COMMON /CVQXZ/CV(2,NG,NG),Q(NG,NG),QX(NG,NG),QZ(NG,NG),GNORM
      COMMON /GLS/GLS1,GLS2
      COMMON /CVNST/CNST(2:4,0:3),SNST(2:4,0:3)
      DATA AU/149597870.691/  EMRAT/81.30056/ 
      DATA SUT/86400./ SBRAT/328900.5614/ DE3/1E-3/
      DATA ZERO,ONE,TWO,THREE/0.,1.,2.,3./
      NT=NN
      RE=6378.1363
      GE=398600.4415
      GM=GE/EMRAT
      GS=SBRAT*(1.+1./EMRAT)*GE
      GLS1=GS/GE
      GLS2=GM/GE
      
      OPEN(9,FILE='MODTOOLS\\GARM360.IN')
      CV=0.
      DO N=0,NGRM
      DO M=0,N
      READ(9,*) I,J,CG,SG
      CV(1,M+1,N+2)=CG
      CV(2,M+1,N+2)=SG
      ENDDO
      ENDDO
      CLOSE(9)
      DO N=2,4
      DO M=0,3
      CNST(N,M)=CV(1,M+1,N+2)
      SNST(N,M)=CV(2,M+1,N+2)
      ENDDO
      ENDDO
      
      OPEN(31,FILE='MODTOOLS\\UT-TAI.IN')
      READ(31,*) NTAI
      DO I=1,NTAI
      READ(31,*) DMJD(I),DUT(1,I),DUT(2,I)
      ENDDO
      DUT=DUT*DE3
      CLOSE(31)

      Q=ZERO
      QX=ZERO
      QZ=ZERO
      DO I=1,NG
      N=I-1
      CN=FLOAT(N)
      C2N=CN+CN
      CNM1=CN-ONE  
      C2NP1=C2N+ONE
      C2NM1=C2N-ONE
      C2NP3=C2N+THREE
      C2NM3=C2N-THREE
      DO J=1,I
      M=J-1
      MM1=M-1
      CM=FLOAT(M)
      CNPM=CN+CM
      CNMM=CN-CM
      CN2MM2=CNPM*CNMM
      DM=ONE
      IF(M.EQ.0) DM=TWO
      QX(I,J)=SQRT(C2NP1*(CNPM+TWO)*(CNPM+ONE)/(DM*C2NP3))
      IF(MM1) 109,107,108
  107 DM=TWO
  108 QX(M,I)=SQRT(DM*C2NP1*(CNMM+TWO)*(CNMM+ONE)/C2NP3)
  109 IF(N.EQ.M) GOTO 110
      Q(I,J)=SQRT(C2NP1*C2NM1/CN2MM2)
      Q(J,I)=SQRT(C2NP1*(CNM1+CM)*(CNM1-CM)/(CN2MM2*C2NM3))
      GOTO 111
  110 IF(N.NE.0) Q(I,J)=SQRT(DM*C2NP1/C2N)
  111 QZ(I,J)=SQRT(C2NP1*(CNPM+ONE)*(CNMM+ONE)/C2NP3)
      ENDDO
      ENDDO
      GNORM=GE/RE**2
      END

