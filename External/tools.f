
C Special functions, taken from CERNLIB/MATHLIB
C M.Boonekamp, 2 June 2004

      FUNCTION DBESI0(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LEX
      CHARACTER NAME0*(*),NAME1*(*),NAME0E*(*),NAME1E*(*)
      CHARACTER*80 ERRTXT
      DIMENSION CI(0:24,0:1),CK(0:16,0:1)

      PARAMETER (NAME0 = 'BESK0/DBESK0', NAME0E = 'EBESK0/DBESK0')
      PARAMETER (NAME1 = 'BESK1/DBESK1', NAME1E = 'EBESK1/DBESK1')
      PARAMETER (EPS=1D-15)
      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (CE = 0.57721 56649 01532 86D0)
      PARAMETER (PIH = PI/2, RPIH = 2/PI, RPI2 = 1/(2*PI))

      DATA CI( 0,0) /+1.00827 92054 58740 032D0/
      DATA CI( 1,0) /+0.00844 51226 24920 943D0/
      DATA CI( 2,0) /+0.00017 27006 30777 567D0/
      DATA CI( 3,0) /+0.00000 72475 91099 959D0/
      DATA CI( 4,0) /+0.00000 05135 87726 878D0/
      DATA CI( 5,0) /+0.00000 00568 16965 808D0/
      DATA CI( 6,0) /+0.00000 00085 13091 223D0/
      DATA CI( 7,0) /+0.00000 00012 38425 364D0/
      DATA CI( 8,0) /+0.00000 00000 29801 672D0/
      DATA CI( 9,0) /-0.00000 00000 78956 698D0/
      DATA CI(10,0) /-0.00000 00000 33127 128D0/
      DATA CI(11,0) /-0.00000 00000 04497 339D0/
      DATA CI(12,0) /+0.00000 00000 01799 790D0/
      DATA CI(13,0) /+0.00000 00000 00965 748D0/
      DATA CI(14,0) /+0.00000 00000 00038 604D0/
      DATA CI(15,0) /-0.00000 00000 00104 039D0/
      DATA CI(16,0) /-0.00000 00000 00023 950D0/
      DATA CI(17,0) /+0.00000 00000 00009 554D0/
      DATA CI(18,0) /+0.00000 00000 00004 443D0/
      DATA CI(19,0) /-0.00000 00000 00000 859D0/
      DATA CI(20,0) /-0.00000 00000 00000 709D0/
      DATA CI(21,0) /+0.00000 00000 00000 087D0/
      DATA CI(22,0) /+0.00000 00000 00000 112D0/
      DATA CI(23,0) /-0.00000 00000 00000 012D0/
      DATA CI(24,0) /-0.00000 00000 00000 018D0/

      DATA CI( 0,1) /+0.97580 06023 26285 926D0/
      DATA CI( 1,1) /-0.02446 74429 63276 385D0/
      DATA CI( 2,1) /-0.00027 72053 60763 829D0/
      DATA CI( 3,1) /-0.00000 97321 46728 020D0/
      DATA CI( 4,1) /-0.00000 06297 24238 640D0/
      DATA CI( 5,1) /-0.00000 00659 61142 154D0/
      DATA CI( 6,1) /-0.00000 00096 13872 919D0/
      DATA CI( 7,1) /-0.00000 00014 01140 901D0/
      DATA CI( 8,1) /-0.00000 00000 47563 167D0/
      DATA CI( 9,1) /+0.00000 00000 81530 681D0/
      DATA CI(10,1) /+0.00000 00000 35408 148D0/
      DATA CI(11,1) /+0.00000 00000 05102 564D0/
      DATA CI(12,1) /-0.00000 00000 01804 409D0/
      DATA CI(13,1) /-0.00000 00000 01023 594D0/
      DATA CI(14,1) /-0.00000 00000 00052 678D0/
      DATA CI(15,1) /+0.00000 00000 00107 094D0/
      DATA CI(16,1) /+0.00000 00000 00026 120D0/
      DATA CI(17,1) /-0.00000 00000 00009 561D0/
      DATA CI(18,1) /-0.00000 00000 00004 713D0/
      DATA CI(19,1) /+0.00000 00000 00000 829D0/
      DATA CI(20,1) /+0.00000 00000 00000 743D0/
      DATA CI(21,1) /-0.00000 00000 00000 080D0/
      DATA CI(22,1) /-0.00000 00000 00000 117D0/
      DATA CI(23,1) /+0.00000 00000 00000 011D0/
      DATA CI(24,1) /+0.00000 00000 00000 019D0/

      DATA CK( 0,0) /+0.98840 81742 30825 800D0/
      DATA CK( 1,0) /-0.01131 05046 46928 281D0/
      DATA CK( 2,0) /+0.00026 95326 12762 724D0/
      DATA CK( 3,0) /-0.00001 11066 85196 665D0/
      DATA CK( 4,0) /+0.00000 06325 75108 500D0/
      DATA CK( 5,0) /-0.00000 00450 47337 641D0/
      DATA CK( 6,0) /+0.00000 00037 92996 456D0/
      DATA CK( 7,0) /-0.00000 00003 64547 179D0/
      DATA CK( 8,0) /+0.00000 00000 39043 756D0/
      DATA CK( 9,0) /-0.00000 00000 04579 936D0/
      DATA CK(10,0) /+0.00000 00000 00580 811D0/
      DATA CK(11,0) /-0.00000 00000 00078 832D0/
      DATA CK(12,0) /+0.00000 00000 00011 360D0/
      DATA CK(13,0) /-0.00000 00000 00001 727D0/
      DATA CK(14,0) /+0.00000 00000 00000 275D0/
      DATA CK(15,0) /-0.00000 00000 00000 046D0/
      DATA CK(16,0) /+0.00000 00000 00000 008D0/

      DATA CK( 0,1) /+1.03595 08587 72358 331D0/
      DATA CK( 1,1) /+0.03546 52912 43331 114D0/
      DATA CK( 2,1) /-0.00046 84750 28166 889D0/
      DATA CK( 3,1) /+0.00001 61850 63810 053D0/
      DATA CK( 4,1) /-0.00000 08451 72048 124D0/
      DATA CK( 5,1) /+0.00000 00571 32218 103D0/
      DATA CK( 6,1) /-0.00000 00046 45554 607D0/
      DATA CK( 7,1) /+0.00000 00004 35417 339D0/
      DATA CK( 8,1) /-0.00000 00000 45757 297D0/
      DATA CK( 9,1) /+0.00000 00000 05288 133D0/
      DATA CK(10,1) /-0.00000 00000 00662 613D0/
      DATA CK(11,1) /+0.00000 00000 00089 048D0/
      DATA CK(12,1) /-0.00000 00000 00012 726D0/
      DATA CK(13,1) /+0.00000 00000 00001 921D0/
      DATA CK(14,1) /-0.00000 00000 00000 305D0/
      DATA CK(15,1) /+0.00000 00000 00000 050D0/
      DATA CK(16,1) /-0.00000 00000 00000 009D0/

      NU=0
      LEX=.FALSE.
      GO TO 6

      ENTRY DEBSI0(X)
      NU=0
      LEX=.TRUE.
      GO TO 6

      ENTRY DBESI1(X)
      NU=1
      LEX=.FALSE.
      GO TO 6

      ENTRY DEBSI1(X)
      NU=1
      LEX=.TRUE.

    6 V=ABS(X)
      IF(V .LT. 8) THEN
       Y=(HF*V)**2
       XL=NU+2
       A0=1
       A1=1+2*Y/((XL+1)*(XL-1))
       A2=1+Y*(4+3*Y/((XL+2)*XL))/((XL+3)*(XL-1))
       B0=1
       B1=1-Y/(XL+1)
       B2=1-Y*(1-Y/(2*(XL+2)))/(XL+3)
       W1=3+XL
       V1=3-XL
       V3=XL-1
       V2=V3+V3
       C=0
       DO 3 N = 3,30
       C0=C
       FN=N
       W1=W1+2
       W2=W1-1
       W3=W2-1
       W4=W3-1
       W5=W4-1
       W6=W5-1
       V1=V1+1
       V2=V2+1
       V3=V3+1
       U1=FN*W4
       E=V3/(U1*W3)
       U2=E*Y
       F1=1+Y*V1/(U1*W1)
       F2=(1+Y*V2/(V3*W2*W5))*U2
       F3=-Y*Y*U2/(W4*W5*W5*W6)
       A=F1*A2+F2*A1+F3*A0
       B=F1*B2+F2*B1+F3*B0
       C=A/B
       IF(ABS(C0-C) .LT. EPS*ABS(C)) GO TO 4
       A0=A1
       A1=A2
       A2=A
       B0=B1
       B1=B2
       B2=B
    3  CONTINUE
    4  H=C
       IF(NU .EQ. 1) H=HF*X*H
       IF(LEX) H=EXP(-V)*H
      ELSE
       R=1/V
       H=16*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 24,0,-1
       B0=CI(I,NU)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=SQRT(RPI2*R)*(B0-H*B2)
       IF(NU*X .LT. 0) H=-H
       IF(.NOT.LEX) H=EXP(V)*H
      ENDIF
      GO TO 9

      ENTRY DBESK0(X)
      NU=0
      LEX=.FALSE.
      GO TO 8

      ENTRY DEBSK0(X)
      NU=0
      LEX=.TRUE.
      GO TO 8

      ENTRY DBESK1(X)
      NU=1
      LEX=.FALSE.
      GO TO 8

      ENTRY DEBSK1(X)
      NU=1
      LEX=.TRUE.

    8 IF(X .LE. 0) THEN
       H=0
       WRITE(ERRTXT,101) X
c       IF(NU .EQ. 0 .AND. .NOT.LEX) CALL MTLPRT(NAME0 ,'C313.1',ERRTXT)
c       IF(NU .EQ. 0 .AND.      LEX) CALL MTLPRT(NAME0E,'C313.1',ERRTXT)
c       IF(NU .EQ. 1 .AND. .NOT.LEX) CALL MTLPRT(NAME1 ,'C313.1',ERRTXT)
c       IF(NU .EQ. 1 .AND.      LEX) CALL MTLPRT(NAME1E,'C313.1',ERRTXT)
      ELSEIF(X .LT. 1) THEN
       B=HF*X
       BK=-(LOG(B)+CE)
       F=BK
       P=HF
       Q=HF
       C=1
       D=B**2
       BK1=P
       DO 11 N = 1,15
       FN=N
       RFN=1/FN
       P=P*RFN
       Q=Q*RFN
       F=(F+P+Q)*RFN
       C=C*D*RFN
       G=C*(P-FN*F)
       H=C*F
       BK=BK+H
       BK1=BK1+G
       IF(BK1*H+ABS(G)*BK .LE. EPS*BK*BK1) GO TO 12
   11  CONTINUE
   12  H=BK
       IF(NU .EQ. 1) H=BK1/B
       IF(LEX) H=EXP(X)*H
      ELSEIF(X .LE. 5) THEN
       XN=4*NU**2
       A=9-XN
       B=25-XN
       C=768*X**2
       C0=48*X
       A0=1
       A1=(16*X+7+XN)/A
       A2=(C+C0*(XN+23)+XN*(XN+62)+129)/(A*B)
       B0=1
       B1=(16*X+9-XN)/A
       B2=(C+C0*B)/(A*B)+1
       C=0
       DO 24 N = 3,30
       C0=C
       FN=N
       FN2=FN+FN
       FN1=FN2-1
       FN3=FN1/(FN2-3)
       FN4=12*FN**2-(1-XN)
       FN5=16*FN1*X
       RAN=1/((FN2+1)**2-XN)
       F1=FN3*(FN4-20*FN)+FN5
       F2=28*FN-FN4-8+FN5
       F3=FN3*((FN2-5)**2-XN)
       A=(F1*A2+F2*A1+F3*A0)*RAN
       B=(F1*B2+F2*B1+F3*B0)*RAN
       C=A/B
       IF(ABS(C0-C) .LT. EPS*ABS(C)) GO TO 25
       A0=A1
       A1=A2
       A2=A
       B0=B1
       B1=B2
       B2=B
   24  CONTINUE
   25  H=C/SQRT(RPIH*X)
       IF(.NOT.LEX) H=EXP(-X)*H
      ELSE
       R=1/X
       H=10*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 23 I = 16,0,-1
       B0=CK(I,NU)+ALFA*B1-B2
       B2=B1
   23  B1=B0
       H=SQRT(PIH*R)*(B0-H*B2)
       IF(.NOT.LEX) H=EXP(-X)*H
      ENDIF
    9 CONTINUE
      DBESI0=H

      RETURN
  101 FORMAT(' NON-POSITIVE ARGUMENT X = ',1P,E15.6)
      END

      FUNCTION DGAGNC(A,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C     Calculates the complementary incomplete gamma function G(A,X)
C     as defined in Ref. 1. Based on
C     1. W. Gautschi, ALGORITHM 542 Incomplete Gamma Functions,
C        ACM Trans. Math. Software 5 (1979) 482-489
C     2. W. Gautschi, A computational procedure for incomplete gamma
C        functions, ACM Trans. Math. Software 5 (1979) 466-481

      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'RGAGNC/DGAGNC')

      PARAMETER (EPS = 5D-14)
      PARAMETER (ALH = -0.69314 71805 59945 31D0)
      PARAMETER (Z1 = 1, HALF = Z1/2, QUAR = Z1/4)
      PARAMETER (C1 = 3*Z1/2, KMAX = 600, EPS1 = EPS/100)

      DIMENSION C(25)

      DATA C( 1) / 0.57721 56649 01532 86D0/
      DATA C( 2) /-0.65587 80715 20253 88D0/
      DATA C( 3) /-0.04200 26350 34095 24D0/
      DATA C( 4) / 0.16653 86113 82291 49D0/
      DATA C( 5) /-0.04219 77345 55544 34D0/
      DATA C( 6) /-0.00962 19715 27876 97D0/
      DATA C( 7) / 0.00721 89432 46663 10D0/
      DATA C( 8) /-0.00116 51675 91859 07D0/
      DATA C( 9) /-0.00021 52416 74114 95D0/
      DATA C(10) / 0.00012 80502 82388 12D0/
      DATA C(11) /-0.00002 01348 54780 79D0/
      DATA C(12) /-0.00000 12504 93482 14D0/
      DATA C(13) / 0.00000 11330 27231 98D0/
      DATA C(14) /-0.00000 02056 33841 70D0/
      DATA C(15) / 0.00000 00061 16095 10D0/
      DATA C(16) / 0.00000 00050 02007 64D0/
      DATA C(17) /-0.00000 00011 81274 57D0/
      DATA C(18) / 0.00000 00001 04342 67D0/
      DATA C(19) / 0.00000 00000 07782 26D0/
      DATA C(20) /-0.00000 00000 03696 81D0/
      DATA C(21) / 0.00000 00000 00510 04D0/
      DATA C(22) /-0.00000 00000 00020 58D0/
      DATA C(23) /-0.00000 00000 00005 35D0/
      DATA C(24) / 0.00000 00000 00001 23D0/
      DATA C(25) /-0.00000 00000 00000 12D0/

      GLGAMA(V)=DLGAMA(V)

      H=0
      IF(X .LT. 0) THEN
       WRITE(ERRTXT,101) X
c       CALL MTLPRT(NAME,'C334.1',ERRTXT)
       GO TO 99
      ELSEIF(X .EQ. 0) THEN
       IF(A .LT. 0) THEN
        H=-1/A
       ELSEIF(A .EQ. 0) THEN
c       CALL MTLPRT(NAME,'C334.2','ILLEGAL ARGUMENTS A = X = 0')
       ELSE
        H=1
       ENDIF
       GO TO 99
      ELSE
       ALX=LOG(X)
      ENDIF
      IF(X .LT. QUAR) THEN
       ALFA=ALH/ALX
      ELSE
       ALFA=X+QUAR
      ENDIF
      MA=HALF-A
      AEPS=A+MA

      IF(MA .GT. 0) THEN
       IF(AEPS .NE. 0) THEN
        ALGP1=GLGAMA(1+AEPS)-LOG(ABS(AEPS))
        IF(MA .NE. 1) ALGP1=ALGP1+GLGAMA(1-AEPS)-GLGAMA(MA-AEPS)
       ELSE
        ALGP1=0
       ENDIF
      ELSE
       ALGP1=GLGAMA(1+A)
      ENDIF
      IF(A .GT. ALFA) THEN
       TERM=1
       SUM=1
       DO 1 K = 1,KMAX
       TERM=X*TERM/(A+K)
       SUM=SUM+TERM
       IF(ABS(TERM) .LE. EPS*SUM) GO TO 2
    1  CONTINUE
       GO TO 98
    2  H=1-EXP(A*ALX-X+LOG(SUM)-ALGP1)
      ELSEIF(X .GT. C1) THEN
       P=0
       S=1-A
       Q=(X+S)*(X-1-A)
       R=4*(X+S)
       TERM=1
       SUM=1
       RHO=0
       DO 3 K = 2,KMAX
       P=P+S
       Q=Q+R
       R=R+8
       S=S+2
       T=P*(1+RHO)
       RHO=T/(Q-T)
       TERM=RHO*TERM
       SUM=SUM+TERM
       IF(ABS(TERM) .LE. EPS*SUM) GO TO 4
    3  CONTINUE
       GO TO 98
    4  IF(A .LE. 0) THEN
        H=SUM/(X+1-A)
       ELSE
        H=EXP(A*ALX-X+LOG(A*SUM/(X+1-A))-ALGP1)
       ENDIF
      ELSE
       AE=A
       IF(A .LT. HALF) THEN
        IF(A .LT. -HALF) AE=AEPS
        SUM=C(25)
        DO 12 K = 24,1,-1
   12   SUM=AE*SUM+C(K)
        GA=-SUM/(1+AE*SUM)
        Y=AE*ALX
        IF(ABS(Y) .GE. 1) THEN
         U=GA-(EXP(Y)-1)/AE
        ELSE
         SUM=1
         TERM=1
         DO 7 K = 2,KMAX
         TERM=Y*TERM/K
         SUM=SUM+TERM
         IF(ABS(TERM) .LE. EPS1*SUM) GO TO 8
    7    CONTINUE
         GO TO 98
    8    U=GA-SUM*ALX
        ENDIF
       ELSE
        U=EXP(GLGAMA(A))-X**A/A
       ENDIF
       P=AE*X
       Q=AE+1
       R=AE+3
       TERM=1
       SUM=1
       DO 9 K = 2,KMAX
       P=P+X
       Q=Q+R
       R=R+2
       TERM=-P*TERM/Q
       SUM=SUM+TERM
       IF(ABS(TERM) .LE. EPS1*SUM) GO TO 10
    9  CONTINUE
       GO TO 98
   10  H=U+SUM*X**(AE+1)/(AE+1)
       IF(A .LT. -HALF) THEN
        H=H*EXP(X-AE*ALX)
        DO 13 J = 1,MA
   13   H=(1-X*H)/(J-AE)
       ELSEIF(A .LE. 0) THEN
        H=H*EXP(X-A*ALX)
       ELSE
        H=A*H*EXP(-ALGP1)
       ENDIF
      ENDIF
   99 DGAGNC=H
      RETURN

   98 WRITE(ERRTXT,103) A,X
c      CALL MTLPRT(NAME,'C334.3',ERRTXT)
      GO TO 99
  101 FORMAT('ILLEGAL ARGUMENT  X = ',1P,E15.6,' < 0')
  103 FORMAT('PROBLEMS WITH CONVERGENCE, A = ',1P,E15.8,'  X = ',E15.6)
      END


      FUNCTION DLGAMA(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) NAME
      PARAMETER(NAME='ALGAMA/DLGAMA')
C
      DIMENSION P1(7),Q1(7),P2(7),Q2(7),P3(7),Q3(7),C(5)

      PARAMETER (Z1 = 1, HF = Z1/2, HF1 = 1+HF)
      CHARACTER*80 ERRTXT
      DATA P1
     1/+3.84287 36567 45991D+0, +5.27068 93753 00983D+1,
     2 +5.55840 45723 51531D+1, -2.15135 13573 72570D+2,
     3 -2.45872 61722 29242D+2, -5.75008 93603 04123D+1,
     4 -2.33590 98949 51284D+0/
      DATA Q1
     1/+1.00000 00000 00000D+0, +3.37330 47907 07074D+1,
     2 +1.93877 84034 37713D+2, +3.08829 54973 42428D+2,
     3 +1.50068 39064 89095D+2, +2.01068 51344 33395D+1,
     4 +4.57174 20282 50299D-1/
      DATA P2
     1/+4.87402 01396 83863 6D+0, +2.48845 25168 57407 6D+2,
     2 +2.17973 66058 89591 5D+3, +3.79751 24011 52511 8D+3,
     3 -1.97780 70769 84164 6D+3, -3.69298 34005 59128 2D+3,
     4 -5.60177 73537 80387 7D+2/
      DATA Q2
     1/+1.00000 00000 00000 0D+0, +9.50999 17418 20893 8D+1,
     2 +1.56120 45277 92863 5D+3, +7.23400 87928 94807 1D+3,
     3 +1.04595 76594 05895 9D+4, +4.16994 15153 20023 1D+3,
     4 +2.76785 83623 80410 1D+2/
      DATA P3
     1/-6.88062 40094 59425D+3, -4.30699 69819 57098D+5,
     2 -4.75045 94653 43956D+6, -2.94234 45930 32234D+6,
     3 +3.63218 04931 54257D+7, -3.35677 82814 54576D+6,
     4 -2.48043 69488 28593D+7/
      DATA Q3
     1/+1.00000 00000 00000D+0, -1.42168 29839 65146D+3,
     2 -1.55528 90280 85353D+5, -3.41525 17108 01107D+6,
     3 -2.09696 23255 80444D+7, -3.45441 75093 34395D+7,
     4 -9.16055 82863 71317D+6/
      DATA C
     1/ 1.12249 21356 561D-1,  7.95916 92961 204D-2,
     1 -1.70877 94611 020D-3,  9.18938 53320 467D-1,
     2  1.34699 05627 879D+0/


      ENTRY DLOGAM(X)

      IF(X .LE. 0) THEN
       H=0
       WRITE(ERRTXT,101) X
c       CALL MTLPRT(NAME,'C304.1',ERRTXT)
      ELSE IF(X .EQ. 1 .OR. X .EQ. 2) THEN
       H=0
      ELSE IF(X .LE. HF) THEN
       Y=X+1
       AP=P1(1)
       AQ=Q1(1)
       DO 2 I = 2,7
       AP=P1(I)+Y*AP
    2  AQ=Q1(I)+Y*AQ
       H=-LOG(X)+X*AP/AQ
      ELSE IF(X .LE. HF1) THEN
       AP=P1(1)
       AQ=Q1(1)
       DO 3 I = 2,7
       AP=P1(I)+X*AP
    3  AQ=Q1(I)+X*AQ
       H=(X-1)*AP/AQ
      ELSE IF(X .LE. 4) THEN
       AP=P2(1)
       AQ=Q2(1)
       DO 4 I = 2,7
       AP=P2(I)+X*AP
    4  AQ=Q2(I)+X*AQ
       H=(X-2)*AP/AQ
      ELSE IF(X .LE. 12) THEN
       AP=P3(1)
       AQ=Q3(1)
       DO 5 I = 2,7
       AP=P3(I)+X*AP
    5  AQ=Q3(I)+X*AQ
       H=AP/AQ
      ELSE
       Y=1/X**2
       H=(X-HF)*LOG(X)-X+C(4)+(C(1)+Y*(C(2)+Y*C(3)))/
     1                                        ((C(5)+Y)*X)
      ENDIF
      DLGAMA=H
      RETURN
  101 FORMAT('NON-POSITIVE ARGUMENT  X = ',1P,E15.6)
      END
      SUBROUTINE VEGAS(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z )
      COMMON/BVEG2/NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
     +,D(50,10),DI(50,10),NXI(50,10)
      DIMENSION XIN(50),R(50),DX(10),IA(10),KG(10),DT(10)
      DIMENSION XL(10),XU(10),QRAN(10),X(10)
      COMMON/RESULT/S1,S2,S3,S4
      EXTERNAL FXN
      DATA XL,XU/10*0.D+00,10*1.0D+00/
      DATA NDMX/50/,ALPH/1.5/,ONE/1./,MDS/1/
      IPR=1
      IF(NPRN.GT.0)IPR=0
      NDO=1
      DO 1 J=1,NDIM
1     XI(1,J)=ONE
      ENTRY VEGAS1(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
      NOW=IGRAPH
      IF(IGRAPH.GT.0)CALL INPLOT(NOW,F1,W)
      IT=0
      SI=0.
      SI2=SI
      SWGT=SI
      SCHI=SI
      SCALLS=SI
      ENTRY VEGAS2(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
      ND=NDMX
      NG=1
      IF(MDS.EQ.0) GO TO 2
      NG=(NCALL*0.5D+00)**(1.D+00/NDIM)
      MDS=1
      IF((2*NG-NDMX).LT.0) GO TO 2
      MDS=-1
      NPG=NG/NDMX+1
      ND=NG/NPG
      NG=NPG*ND
2     K=NG**NDIM
      NPG=NCALL/K
      IF(NPG.LT.2)NPG=2
      CALLS=NPG*K
      DXG=ONE/NG
      DV2G=DXG**(2*NDIM)/NPG/NPG/(NPG-ONE)
      XND=ND
      NDM=ND-1
      DXG=DXG*XND
      XJAC=ONE
      DO 3 J=1,NDIM
      DX(J)=XU(J)-XL(J)
3     XJAC=XJAC*DX(J)
      IF(ND.EQ.NDO) GO TO 8
      RC=NDO/XND
      DO 7 J=1,NDIM
      K=0
      XN=0.
      DR=XN
      I=K
4     K=K+1
      DR=DR+ONE
      XO=XN
      XN=XI(K,J)
5     IF(RC.GT.DR) GO TO 4
      I=I+1
      DR=DR-RC
      XIN(I)=XN-(XN-XO)*DR
      IF(I.LT.NDM) GO TO 5
      DO 6  I=1,NDM
6     XI(I,J)=XIN(I)
7     XI(ND,J)=ONE
      NDO=ND
      ACC=BCC
8     IF(NPRN.NE.0.AND.NPRN.NE.10)PRINT 200,NDIM,CALLS,IT,ITMX
     1,ACC,MDS,ND
      IF(NPRN.EQ.10)PRINT 290,NDIM,CALLS,ITMX,ACC,MDS,ND
      ENTRY VEGAS3(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
9     IT=IT+1
      TI=0.
      TSI=TI
      IF(IGRAPH.GT.0)CALL REPLOT(NOW,F1,W)
      DO 10 J=1,NDIM
      KG(J)=1
      DO 10 I=1,ND
      NXI(I,J)=0
      D(I,J)=TI
10    DI(I,J)=TI
11    FB=0.D0
      F2B=FB
      K=0
12    K=K+1
      DO 121 J=1,NDIM
121   QRAN(J)=RANF(0.0d0)
      WGT=XJAC
      DO 15 J=1,NDIM
      XN=(KG(J)-QRAN(J))*DXG+ONE
      IA(J)=XN
      IAJ=IA(J)
      IAJ1=IAJ-1
      IF(IAJ.GT.1) GO TO 13
      XO=XI(IAJ,J)
      RC=(XN-IAJ)*XO
      GO TO 14
13    XO=XI(IAJ,J)-XI(IAJ1,J)
      RC=XI(IAJ1,J)+(XN-IAJ)*XO
14    X(J)=XL(J)+RC*DX(J)
15    WGT=WGT*XO*XND
      F=FXN(X)*WGT
      F1=F/CALLS
      W=WGT/CALLS
      IF(IGRAPH.GT.0)CALL XPLOT(NOW,F1,W)
      F2=F*F
      FB=FB+F
      F2B=F2B+F2
      DO 16 J=1,NDIM
      IAJ=IA(J)
      NXI(IAJ,J)=NXI(IAJ,J)+1
      DI(IAJ,J)=DI(IAJ,J)+F/CALLS
16    IF(MDS.GE.0)  D(IAJ,J)=D(IAJ,J)+F2
      IF(K.LT.NPG) GO TO 12
      F2B=F2B*NPG
      F2B=DSQRT(F2B)
      F2B=(F2B-FB)*(F2B+FB)
      TI=TI+FB
      TSI=TSI+F2B
      IF(MDS.GE.0) GO TO 18
      DO 17 J=1,NDIM
      IAJ=IA(J)
17    D(IAJ,J)=D(IAJ,J)+F2B
18    K=NDIM
19    KG(K)=MOD(KG(K),NG)+1
      IF(KG(K).NE.1) GO TO 11
      K=K-1
      IF(K.GT.0) GO TO 19
      TI=TI/CALLS
      TSI=TSI*DV2G
      TI2=TI*TI
      WGT=TI2/TSI
      SI=SI+TI*WGT
      SI2=SI2+TI2
      SWGT=SWGT+WGT
      SCHI=SCHI+TI2*WGT
      SCALLS=SCALLS+CALLS
      AVGI=SI/SWGT
      SD=SWGT*IT/SI2
      CHI2A=0.
      IF(IT.GT.1)CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(IT-1)
      SD=ONE/SD
      SD=DSQRT(SD)
      IF(NPRN.EQ.0) GO TO 21
      TSI=DSQRT(TSI)
      IF(NPRN.NE.10)PRINT 201,IPR,IT,TI,TSI,AVGI,SD,CHI2A
      IF(NPRN.EQ.10)PRINT 203,IT,TI,TSI,AVGI,SD,CHI2A
      IF(NPRN.GE.0) GO TO 21
      DO 20 J=1,NDIM
      PRINT 202,J
20    PRINT 204,(XI(I,J),DI(I,J),D(I,J),I=1,ND)
21    IF(DABS(SD/AVGI).LE.DABS(ACC).OR.IT.GE.ITMX)NOW=2
      S1=AVGI
      S2=SD
      S3=TI
      S4=TSI
      IF(IGRAPH.GT.0)CALL PLOTIT(NOW,F1,W)
C      DO 23 J=1,NDIM
C      XO=D(1,J)
C      XN=D(2,J)
C      D(1,J)=(XO+XN)*0.5
C      DT(J)=D(1,J)
C      DO 22 I=2,NDM
C      D(I,J)=XO+XN
C      XO=XN
C      XN=D(I+1,J)
C      D(I,J)=(D(I,J)+XN)/3.
C22    DT(J)=DT(J)+D(I,J)
C      D(ND,J)=(XN+XO)*0.5
C23    DT(J)=DT(J)+D(ND,J)
C-----THIS PART OF THE VEGAS-ALGORITHM IS UNSTABLE
C-----IT SHOULD BE REPLACED BY
      DO 23 J=1,NDIM
      DT(J)=0.
      DO 23 I=1,ND
      IF(NXI(I,J).GT.0)D(I,J)=D(I,J)/NXI(I,J)
23    DT(J)=DT(J)+D(I,J)
      DO 28 J=1,NDIM
      RC=0.
      DO 24 I=1,ND
      R(I)=0.
      IF(D(I,J).LE.0.)GO TO 24
      XO=DT(J)/D(I,J)
      R(I)=((XO-ONE)/XO/DLOG(XO))**ALPH
24    RC=RC+R(I)
      RC=RC/XND
      K=0
      XN=0.
      DR=XN
      I=K
25    K=K+1
      DR=DR+R(K)
      XO=XN
      XN=XI(K,J)
26    IF(RC.GT.DR) GO TO 25
      I=I+1
      DR=DR-RC
      XIN(I)=XN-(XN-XO)*DR/R(K)
      IF(I.LT.NDM) GO TO 26
      DO 27 I=1,NDM
27    XI(I,J)=XIN(I)
28    XI(ND,J)=ONE
      IF(IT.LT.ITMX.AND.DABS(ACC).LT.DABS(SD/AVGI))GO TO 9
200   FORMAT(35H0INPUT PARAMETERS FOR VEGAS   NDIM=,I3,
     18H  NCALL=,F8.0/28X,5H  IT=,I5,8H  ITMX =,I5/28X,
     26H  ACC=,G9.3/28X,6H  MDS=,I3,6H   ND=,I4//)
290   FORMAT(13H0VEGAS  NDIM=,I3,8H  NCALL=,F8.0,8H  ITMX =,I5,
     16H  ACC=,G9.3,6H  MDS=,I3,6H   ND=,I4)
201   FORMAT(/I1,20HINTEGRATION BY VEGAS/13H0ITERATION NO,I3,
     114H.   INTEGRAL =,G14.8/20X,10HSTD DEV  =,G10.4/
     234H ACCUMULATED RESULTS.   INTEGRAL =,G14.8/
     324X,10HSTD DEV  =,G10.4 / 24X,18HCHI**2 PER ITN   =,G10.4)
202   FORMAT(14H0DATA FOR AXIS,I2 / 7X,1HX,7X,10H  DELT I  ,
     12X,11H CONVCE    ,11X,1HX,7X,10H  DELT I  ,2X,11H CONVCE
     2,11X,1HX,7X,10H  DELT I  ,2X,11H CONVCE     /)
204   FORMAT(1X,3G12.4,5X,3G12.4,5X,3G12.4)
203   FORMAT(1H ,I3,G20.8,G12.4,G20.8,G12.4,G12.4)
      S1=AVGI
      S2=SD
      S3=CHI2A
      RETURN
      END
 
      SUBROUTINE INPLOT(NOW,FF,PDX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/LPLOT/XL(25),V1(10),V2(10),AV(10)
      DIMENSION ZAV(10),YAV(10),ZSV(10),YSV(10),ZTV(10)
      DIMENSION XLMAX(25),XLMIN(25),NLP(25),LTOP(25),TEXT(8,25)
     1,LL(25)
      DIMENSION NUMB(12)
      DIMENSION XLS(42,25),YLS(42,25),NLSN(42,25),MLSN(42,25)
     1,DLS(25)
     1,XLAV(25),XLSQ(25),XLAVA(25),SXA(25),TLIM(6),TOP(25),XLTQ(25)
      DIMENSION NBIN(41),NLOG(41),SLOG(41),TLOG(41),HV(12)
      DIMENSION V1MAX(10),V1MIN(10),V2MAX(10),V2MIN(10),NV1(10)
     1,NV2(10),VTEXT(6,10)
      DIMENSION VM(12,12,10),NVM(12,12,10),BIN1(10),BIN2(10),VOL(10)
     1,WM(12,12,10),MVM(12,12,10)
      COMMON/RESULT/Y,SI,U,V
      CHARACTER*1 HMIN,HPLUS,HBLANK,HSTAR,CHAR(40)
      SAVE
      DATA TLIM/1.6D+00,2.5D+00,4.0D+00,6.666666667D+00,10.D+00,
     +16.D+00/
      DATA HMIN,HPLUS,HBLANK,HSTAR/'-','+',' ','*'/
      DATA MLS,MAV,NDMAX/25,10,10/
      DATA NGRAPH/0/
C      PRINT '('' IGRAPH = '',I8)',NOW
      IGRAPH=NOW
      NOW=0
      KK=0
      ITT=0
      IF(IGRAPH.LE.0) GO TO 800
      IF(IGRAPH.NE.NGRAPH)READ (12,810)NLS
      IF(IGRAPH.NE.NGRAPH)PRINT 814,NLS
      IF(NLS.LT.0) NLS=0
      IF(NLS.EQ.0) GO TO 802
      IF(NLS.GT.MLS) GO TO 807
      IF(IGRAPH.NE.NGRAPH)PRINT 815
      DO 801 I=1,NLS
      IF(IGRAPH.NE.NGRAPH)
     1READ (12,811)XLMIN(I),XLMAX(I),NLP(I),LTOP(I),LL(I),
     2(TEXT(J,I),J=1,8)
      IF(IGRAPH.NE.NGRAPH) PRINT 816,I,XLMIN(I),XLMAX(I)
     1,NLP(I),LTOP(I),LL(I),(TEXT(J,I),J=1,8)
      IF(NLP(I).LT.1)NLP(I)=1
      IF(NLP(I).GT.40)NLP(I)=40
      DLS(I)=(XLMAX(I)-XLMIN(I))/NLP(I)
      NLPS=NLP(I)+2
      DO 300 J=1,NLPS
      YLS(J,I)=0
300   MLSN(J,I)=0
801   CONTINUE
802   IF(IGRAPH.NE.NGRAPH) READ (12,810)NDD
      IF(IGRAPH.NE.NGRAPH) PRINT 817,NDD
      IF(NDD.LT.0) NDD=0
      IF(NDD.EQ.0) GO TO 804
      IF(NDD.GT.NDMAX) GO TO 807
      IF(IGRAPH.NE.NGRAPH) PRINT 818
      DO 803 I=1,NDD
      IF(IGRAPH.NE.NGRAPH)
     1READ (12,812)V1MIN(I),V1MAX(I),NV1(I),V2MIN(I),V2MAX(I)
     1,NV2(I),(VTEXT(J,I),J=1,6)
      IF(IGRAPH.NE.NGRAPH) PRINT 819,I,V1MIN(I),V1MAX(I),NV1(I)
     1,V2MIN(I),V2MAX(I),NV2(I),(VTEXT(J,I),J=1,6)
      IF(NV1(I).LT.1)NV1(I)=1
      IF(NV2(I).LT.1)NV2(I)=1
      IF(NV1(I).GT.10)NV1(I)=10
      IF(NV2(I).GT.10)NV2(I)=10
      BIN1(I)=(V1MAX(I)-V1MIN(I))/NV1(I)
      BIN2(I)=(V2MAX(I)-V2MIN(I))/NV2(I)
      VOL(I)=BIN1(I)*BIN2(I)
803   CONTINUE
      WTOW=0.
      DO 805 I=1,NDD
      DO 805 J=1,12
      DO 805 K=1,12
      WM(K,J,I)=0.
805   MVM(K,J,I)=0
804   CONTINUE
      IF(IGRAPH.NE.NGRAPH)READ (12,810)NAVE
      IF(IGRAPH.NE.NGRAPH)PRINT 820,NAVE
      IF(NAVE.LT.0)NAVE=0
      IF(NAVE.GT.MAV)GO TO 807
      DO 11 I=1,MAV
      YAV(I)=0.
11    YSV(I)=0.
      KT=0
      GO TO 808
800   NAVE=0
      NLS=0
      NDD=0
      GO TO 808
807   PRINT 813,MAV,MLS,NDMAX
      STOP
808   NGRAPH=IGRAPH
      RETURN
      ENTRY REPLOT(NOW,FF,PDX)
      IF(NAVE.EQ.0) GO TO 49
      DO 62 I=1,NAVE
      ZAV(I)=0.
      ZTV(I)=0.
62    ZSV(I)=0.
49    FSQA=0.
      KT=KT+1
      IF(NLS.EQ.0) GO TO 303
      DO 302 I=1,NLS
      NLPS=NLP(I)+2
      XLAV(I)=0
      XLTQ(I)=0.
      XLSQ(I)=0
      DO 302 J=1,NLPS
      XLS(J,I)=0
302   NLSN(J,I)=0
303   CONTINUE
      IF(NDD.EQ.0) GO TO 403
      DO 402 I=1,NDD
      N1=NV1(I)+2
      N2=NV2(I)+2
      DO 402 I1=1,N1
      DO 402 I2=1,N2
      VM(I1,I2,I)=0
402   NVM(I1,I2,I)=0
403   CONTINUE
      RETURN
      ENTRY XPLOT(NOW,FF,PDX)
      FSQA=FSQA+FF*FF/PDX
      ITT=ITT+1
      IF(NLS.EQ.0) GO TO 305
      DO 304  I=1,NLS
      NLPS=(XL(I)-XLMIN(I))/DLS(I)+1
      IF(NLPS.LT.0)NLPS=0
      IF(NLPS.GT.NLP(I))NLPS=NLP(I)+1
      NLPS=NLPS+1
      XLS(NLPS,I)=XLS(NLPS,I)+FF/DLS(I)
      NLSN(NLPS,I)=NLSN(NLPS,I)+1
      XLAV(I)=XLAV(I)+FF*XL(I)
      XLTQ(I)=XLTQ(I)+FF*FF*XL(I)/PDX
304   XLSQ(I)=XLSQ(I)+(FF*XL(I))**2/PDX
305   CONTINUE
      IF(NDD.EQ.0)GO TO 405
      DO 404 I=1,NDD
      I1=(V1(I)-V1MIN(I))/BIN1(I)+2
      IF(I1.LT.1) I1=1
      IF(I1.GT.NV1(I)+2) I1=NV1(I)+2
      I2=(V2(I)-V2MIN(I))/BIN2(I)+2
      IF(I2.LT.1) I2=1
      IF(I2.GT.NV2(I)+2) I2=NV2(I)+2
      VM(I1,I2,I)=VM(I1,I2,I)+FF/VOL(I)
404   NVM(I1,I2,I)=NVM(I1,I2,I)+1
405   CONTINUE
      IF(NAVE.EQ.0)GO TO 99
      DO 22 I=1,NAVE
      ZAV(I)=ZAV(I)+AV(I)*FF
      ZTV(I)=ZTV(I)+FF*FF*AV(I)/PDX
22    ZSV(I)=ZSV(I)+(AV(I)*FF)**2/PDX
99    RETURN
      ENTRY PLOTIT(NOW,FF,PDX)
      IF(NLS.EQ.0)GO TO 315
      IF(KK.GT.0)GO TO 307
      DO 306 I=1,NLS
      NLPS=NLP(I)+2
      DO 306 J=1,NLPS
      MLSN(J,I)=NLSN(J,I)
306   YLS(J,I)=XLS(J,I)
      GO TO 310
307   VBEF=VTOT
      VU=(V/U)**2
      DO 309 I=1,NLS
      NLPS=NLP(I)+2
      DO 309 J=1,NLPS
      IF(NLSN(J,I).EQ.0)GO TO 309
      IF(MLSN(J,I).EQ.0)GO TO 308
      AL1=VU/NLSN(J,I)
      AL2=VBEF/MLSN(J,I)
      MLSN(J,I)=MLSN(J,I)+NLSN(J,I)
      YLS(J,I)=(AL2*XLS(J,I)+AL1*YLS(J,I))/(AL1+AL2)
      GO TO 309
308   MLSN(J,I)=NLSN(J,I)
      YLS(J,I)=XLS(J,I)
309   CONTINUE
310   CONTINUE
      DO 311 I=1,NLS
      SXF=XLSQ(I)-XLAV(I)*XLAV(I)
      SXT=XLTQ(I)-XLAV(I)*U
      SX2=XLSQ(I)/XLAV(I)**2+FSQA/U**2-2.*XLTQ(I)/(XLAV(I)*U)
      SX2=SX2*(XLAV(I)/U)**2
      IF(KT.NE.1)GO TO 312
      XLAVA(I)=XLAV(I)/U
      SXA(I)=SX2
      GO TO 311
312   XHELP=SX2+SXA(I)
      IF(XHELP.EQ.0)GO TO 311
      XLAVA(I)=(XLAV(I)*SXA(I)/U+XLAVA(I)*SX2)/XHELP
      SXA(I)=SXA(I)*SX2/XHELP
311   CONTINUE
      VTOT=(SI/Y)**2
      IF(NOW.NE.2)GO TO 315
      DO 341 I=1,NLS
      TOP(I)=0.
      NLPS=NLP(I)+1
      DO 341 J=2,NLPS
      XLS(J,I)=YLS(J,I)/Y
      IF(XLS(J,I).GT.TOP(I))TOP(I)=XLS(J,I)
341   CONTINUE
      DO 342 I=1,NLS
      IF(LTOP(I).LE.0)LTOP(I)=I
      LTO=LTOP(I)
      IF(TOP(I).GT.TOP(LTO))TOP(LTO)=TOP(I)
342   CONTINUE
      YLOG=0.5*DLOG10(Y*Y)
      DO 314 I=1,NLS
      PRINT 321,I
      NLPS=NLP(I)+1
      LTO=LTOP(I)
      TOP(I)=TOP(LTO)
      IF(TOP(I).EQ.0)TOP(I)=1.
      AN1=DLOG10(TOP(I))
      N1=AN1
      IF(N1.GT.AN1)N1=N1-1
      Z1=TOP(I)*10.**(-N1)
      DO 343 L=1,4
      IF(Z1.LT.TLIM(L))GO TO 344
343   CONTINUE
      L=5
344   IF(TOP(I).LT.1.6/(XLMAX(I)-XLMIN(I)))L=L+1
      TOPM=TLIM(L)*10.**N1
      DO 345 J=2,NLPS
      NBIN(J)=XLS(J,I)*40./TOPM+1.5
      IF(LL(I).LT.0)NBIN(J)=0
      IF(XLS(J,I).LE.0) GO TO 346
      TLOG(J)=DLOG10(XLS(J,I))
      SLOG(J)=TLOG(J)+YLOG
      NLOG(J)=(TLOG(J)-N1)*8.+33.5
      IF(LL(I).GT.0)NLOG(J)=0
      GO TO 345
346   SLOG(J)=0
      TLOG(J)=0
      NLOG(J)=0
345   CONTINUE
      PRINT 322,(TEXT(J,I),J=1,8)
      N1P1=N1+1
      N1M4=N1-4
      PRINT 323,TLIM(L),N1,N1P1,N1M4
      DO 347 L=1,40
      CHAR(L)=HMIN
      IF(NLOG(L+1).EQ.41)CHAR(L)=HPLUS
      IF(NBIN(L+1).EQ.41)CHAR(L)=HSTAR
347   CONTINUE
      XMIN=XLMIN(I)
      XMAX=XMIN+DLS(I)
      PRINT 324,XMIN,XMAX,YLS(2,I),SLOG(2),XLS(2,I),TLOG(2)
     1,MLSN(2,I),CHAR
      DO 348 J=3,NLPS
      XMIN=XMAX
      XMAX=XMIN+DLS(I)
      DO 349 L=1,40
      CHAR(L)=HBLANK
      IF(NLOG(L+1).EQ.43-J)CHAR(L)=HPLUS
      IF(NBIN(L+1).EQ.43-J)CHAR(L)=HSTAR
349   CONTINUE
      PRINT 324,XMIN,XMAX,YLS(J,I),SLOG(J),XLS(J,I),TLOG(J)
     1,MLSN(J,I),CHAR
348   CONTINUE
      NLPS1=NLPS+1
      IF(NLPS.EQ.41)GO TO 352
      DO 351 J=NLPS1,41
      DO 350 L=1,40
      CHAR(L)=HBLANK
      IF(NLOG(L+1).EQ.43-J)CHAR(L)=HPLUS
      IF(NBIN(L+1).EQ.43-J)CHAR(L)=HSTAR
350   CONTINUE
351   PRINT 325,CHAR
352   DO 353 L=1,40
      CHAR(L)=HMIN
      IF(NLOG(L+1).EQ.1)CHAR(L)=HPLUS
      IF(NBIN(L+1).EQ.1)CHAR(L)=HSTAR
353   CONTINUE
      PRINT 326,CHAR
      EL1=YLS(1,I)*DLS(I)
      EL2=EL1/Y
      PRINT 327,EL1,EL2,MLSN(1,I)
      EL1=YLS(42,I)*DLS(I)
      EL2=EL1/Y
      PRINT 328,EL1,EL2,MLSN(NLPS1,I)
      SXSQ=DSQRT(SXA(I)/ITT)
      PRINT 329,XLAVA(I),SXSQ
314   CONTINUE
315   CONTINUE
      IF(NDD.EQ.0)GO TO 409
      WBEF=WTOT
      DO 500 I=1,NDD
      NX=NV1(I)+2
      NY=NV2(I)+2
      IF(KK.GT.0)GO TO 502
      DO 501 J=1,NX
      DO 501 K=1,NY
      WM(J,K,I)=VM(J,K,I)
501   MVM(J,K,I)=NVM(J,K,I)
      GO TO 500
502   VU=(V/U)**2
      DO 503 J=1,NX
      DO 503 K=1,NY
      IF(NVM(J,K,I).EQ.0)GO TO 503
      IF(MVM(J,K,I).EQ.0)GO TO 504
      AL1=VU/NVM(J,K,I)
      AL2=VBEF/MVM(J,K,I)
      MVM(J,K,I)=MVM(J,K,I)+NVM(J,K,I)
      WM(J,K,I)=(AL2*VM(J,K,I)+AL1*WM(J,K,I))/(AL1+AL2)
      GO TO 503
504   MVM(J,K,I)=NVM(J,K,I)
      WM(J,K,I)=VM(J,K,I)
503   CONTINUE
500   CONTINUE
      WTOT=(SI/Y)**2
      IF(NOW.NE.2) GO TO 409
      DO 408 I=1,NDD
      PRINT 481,I,(VTEXT(J,I),J=1,6)
      VVV=V2MAX(I)
      MVV=NV1(I)+2
      NVV=NV2(I)+1
      SIZE=VOL(I)/Y
      DO 406 I2=1,NVV
      J2=NVV+2-I2
      DO 410 I1=1,MVV
410   NUMB(I1)=1000.*WM(I1,J2,I)*SIZE+.5
      PRINT 486,(NUMB(I1),I1=1,MVV)
      PRINT 483,(WM(I1,J2,I),I1=1,MVV)
      PRINT 486,(MVM(I1,J2,I),I1=1,MVV)
      PRINT 484,VVV
      VVV=VVV-BIN2(I)
      IF(DABS(VVV/BIN2(I)).LT.1.D-10)VVV=0.
406   CONTINUE
      DO 411 I1=1,MVV
411   NUMB(I1)=1000.*WM(I1,1,I)*SIZE+.5
      PRINT 486,(NUMB(I1),I1=1,MVV)
      PRINT 483,(WM(I1,1,I),I1=1,MVV)
      PRINT 486,(MVM(I1,1,I),I1=1,MVV)
      PRINT 482
      MVV=MVV-1
      DO 407 I1=1,MVV
      HV(I1)=V1MIN(I)+(I1-1)*BIN1(I)
      IF(DABS(HV(I1)/BIN1(I)).LT.1.D-10)HV(I1)=0.
407   CONTINUE
      PRINT 485,(HV(I1),I1=1,MVV)
408   CONTINUE
409   CONTINUE
      IF(NAVE.EQ.0) GO TO 23
      IF(NOW.EQ.2) PRINT 26
      DO 24 I=1,NAVE
      SXF=ZSV(I)-ZAV(I)*ZAV(I)
      SXT=ZSV(I)/ZAV(I)**2+FSQA/U**2-2.*ZTV(I)/(ZAV(I)*U)
      SX2=SXT*(ZAV(I)/U)**2
      IF(KT.NE.1) GO TO 21
      YAV(I)=ZAV(I)/U
      YSV(I)=SX2
      GO TO 30
21    XHELP=SX2+YSV(I)
      IF(XHELP.EQ.0)GO TO 30
      YAV(I)=(YSV(I)*ZAV(I)/U+YAV(I)*SX2)/XHELP
      YSV(I)=YSV(I)*SX2/XHELP
30    YSSQ=DSQRT(YSV(I)/ITT)
      IF(NOW.EQ.2)PRINT 27,I,YAV(I),YSSQ
24    CONTINUE
23    NOW=1
      KK=KK+1
      RETURN
27    FORMAT(12X,I2,9X,D15.5,5X,D15.3)
26    FORMAT(1HI,10X,46HTHE FOLLOWING ARE AVERAGES WITH ERROR ESTIMATE/)
321   FORMAT(1HI,40X,
     140HSINGLE DIFFERENTIAL CROSS-SECTION NUMBER,I3///)
322   FORMAT(38H SINGLE DIFFERENTIAL CROSS SECTION OF ,8A4/)
323   FORMAT(11X,6HLIMITS,9X,1HI,16X,
     119HACCUMULATED RESULTS,15X,1HI,24X,
     29HUPPER BIN,6X,9HLOWER BIN/26X,1HI,50X,
     320HI * LINEAR      PLOT,F8.2,5H*10**,I3,8X,1H0/5X,
     45HLOWER,7X,5HUPPER,4X,1HI,5X,5HDS/DX,4X,
     556HALOG10   (DS/DX)/S  ALOG10  POINTS  I + LOGARITHMIC PLOT,
     66X,4H10**,I3,8X,4H10**,I3/2X,24(1H-),1HI,50(1H-),1HI)
324   FORMAT(D12.4,D12.4,3H  I,2(D12.4,F8.2),I8,3H  I,
     14X,1HI,40A1,1HI)
325   FORMAT(26X,1HI,50X,1HI,4X,1HI,40A1,1HI)
326   FORMAT(2X,24(1H-),1HI,50(1H-),1HI,4X,1HI,40A1,1HI)
327   FORMAT(7X,15HTOTAL UNDERFLOW,4X,1HI,D12.4,
     1D20.4,I16,2X,1HI)
328   FORMAT(7X,15HTOTAL  OVERFLOW,4X,1HI,D12.4,D20.4,
     1I16,2X,1HI)
329   FORMAT(//19X,21HACCUMULATED AVERAGE =,D12.5
     1/19X,21HESTIMATED ERROR     =,D12.5)
C---We commented the paragraph between C481 and 482,
C---since it was giving some format errors and replaced
C---it by an empty format statement. 06/24/2003, T.K.
 481  FORMAT()
C481   FORMAT(1HI,45X,40HDOUBLE DIFFERENTIAL CROSS-SECTION 
C     1NUMBER,I3//60X,7HX-AXIS ,3A4/60X,7HY-AXIS ,3A4/)
482   FORMAT(20X,11(1HI,9X))
483   FORMAT(11X,D9.3,11(1HI,D9.3))
484   FORMAT(1X,D10.3,  9H---------,11(10HI---------))
485   FORMAT(1H0,14X,11D10.3)
486   FORMAT(11X,I8,1X,11(1HI,I8,1X))
810   FORMAT(I2)
811   FORMAT(2D12.4,3I2,8A4)
812   FORMAT(2(2D10.3,I4),6A4)
C---We commented the paragraph between C813 and 814,
C---since it was giving some format errors and replaced
C---it by an empty format statement. 06/24/2003, T.K.
 813  FORMAT()
C813   FORMAT(12H1***ERROR***,10X,24HTOO MANY PLOTS REQUESTED//
C     122X,20HTHE UPPER LIMITS ARE //19X,I2,9H 
C     2AVERAGES//19X,I2,22H ONE 
C     3DMENSIONAL PLOTS//19X,I2,22H 
C     4TWO DIMENSIONAL PLOTS////22X,25H***EXE
C     5CUTION IS HALTED***   )
814   FORMAT(37H1NUMBER OF SINGLE DIFFERENTIAL CROSS
     1,20HSECTIONS REQUESTED =,I3/)
815   FORMAT(30H INFORMATION ON THE DATA CARDS//
     13H  I,10X,5HXLMIN,12X,5HXLMAX,7X,
     224HBINS  CORRELLATION  TYPE,19X,4HTEXT/)
816   FORMAT(I3,2E17.4,I8,I9,I10,5X,8A4)
817   FORMAT(37H0NUMBER OF DOUBLE DIFFERENTIAL CROSS
     1,20HSECTIONS REQUESTED =,I3/)
818   FORMAT(30H INFORMATION ON THE DATA CARDS//
     13H  I,8X,5HV1MIN,10X,5HV1MAX,4X,4HBINS,8X,5HV2MIN
     2,10X,5HV2MAX,4X,4HBINS,8X,6HTEXT 1,8X,6HTEXT 2/)
819   FORMAT(I3,2D15.3,I5,1X,2D15.3,I5,6X,3A4,2X,3A4)
820   FORMAT(31H0NUMBER OF AVERAGES REQUESTED =,I3)
      END

       SUBROUTINE IN55(IA,IX)
       PARAMETER (MODULO=1000000000)
       INTEGER IA(55)
       IA(55)=IX
       J=IX
       K=1
       DO 10 I=1,54
       II=MOD(21*I,55)
       IA(II)=K
       K=J-K
       IF(K.LT.0)K=K+MODULO
       J=IA(II)
   10  CONTINUE
       DO 20 I=1,10
       CALL IRN55(IA)
   20  CONTINUE
       END

       SUBROUTINE IRN55(IA)
       PARAMETER (MODULO=1000000000)
       INTEGER IA(55)
       DO 10 I=1,24
       J=IA(I)-IA(I+31)
       IF(J.LT.0)J=J+MODULO
       IA(I)=J
   10  CONTINUE
       DO 20 I=25,55
       J=IA(I)-IA(I-24)
       IF(J.LT.0)J=J+MODULO
       IA(I)=J
   20  CONTINUE
       END

       DOUBLE PRECISION FUNCTION RANF(DUMMY)
*
*      RANDOM NUMBER FUNCTION TAKEN FROM KNUTH
*      (SEMINUMERICAL ALGORITHMS).
*      METHOD IS X(N)=MOD(X(N-55)-X(N-24),1/FMODUL)
*      NO PROVISION YET FOR CONTROL OVER THE SEED NUMBER.
*
*      RANF GIVES ONE RANDOM NUMBER BETWEEN 0 AND 1.
*      IRN55 GENERATES 55 RANDOM NUMBERS BETWEEN 0 AND 1/FMODUL.
*      IN55  INITIALIZES THE 55 NUMBERS AND WARMS UP THE SEQUENCE.
*
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (FMODUL=1.D-09)
       INTEGER IA(55)
       SAVE IA
       DATA NCALL/0/
       DATA MCALL/55/
       IF( NCALL.EQ.0 ) THEN
           CALL IN55 ( IA,234612947 )
           NCALL = 1
       ENDIF
       IF ( MCALL.EQ.0 ) THEN
           CALL IRN55(IA)
           MCALL=55
       ENDIF
       RANF=IA(MCALL)*FMODUL
       MCALL=MCALL-1
       END

      SUBROUTINE SAVE(NDIM,NTAPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/BVEG2/NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
     1 ,D(50,10),DI(50,10),NXI(50,10)
C
C   STORES VEGAS DATA (UNIT NTAPE) FOR LATER RE-INITIALIZATION
      WRITE(NTAPE,200) NDO,IT,SI,SI2,SWGT,SCHI,
     1      ((XI(I,J),I=1,NDO),J=1,NDIM)
     2     ,((DI(I,J),I=1,NDO),J=1,NDIM)
      RETURN
      ENTRY RESTR(NDIM,NTAPE)
C
C   ENTERS INITIALIZATION DATA FOR VEGAS
      READ(NTAPE,200) NDO,IT,SI,SI2,SWGT,SCHI,
     1      ((XI(I,J),I=1,NDO),J=1,NDIM)
     2     ,((DI(I,J),I=1,NDO),J=1,NDIM)
200   FORMAT(2I8,4Z16/(5Z16))
      RETURN
      END


