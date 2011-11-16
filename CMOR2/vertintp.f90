!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program drive_vinth2p
  implicit none
  !
  integer,parameter::imax = 256,nlat = 128,nlevi = 26,nlevip1 = nlevi+1,nlevo = 17
  real,dimension(imax,nlat,nlevi)::dati,revdati
!  real*8,dimension(imax,nlat,nlevo)::dato
  real,dimension(imax,nlat,nlevo)::dato
  real,dimension(nlevi)::plevi,hbcofa,hbcofb
  real,dimension(nlevo)::plevo,rplevo
!  real,dimension(imax,nlat)::psfc,tbot,phis
  real,dimension(imax,nlat)::psfc
  !
!  real*8::p0,spvl
  real::p0,spvl
  integer::intyp,ilev,kxtrp,varflg,i,j,k,kp,r
  !
  plevi = (/ 3.544638,   7.3888135,  13.967214,  23.944625,  37.23029 , &
            53.114605,  70.05915  ,  85.439115, 100.514695, 118.250335, &
           139.115395, 163.66207  , 192.539935, 226.513265, 266.481155, &
           313.501265, 368.81798  , 433.895225, 510.455255, 600.5242  , &
           696.79629 , 787.70206  , 867.160760, 929.648875, 970.55483 , &
           992.5561  /)
  !
  plevo = (/1000.,925.,850.,700.,600.,500.,400.,300.,250.,200.,150.,100.,70.,50.,30.,20.,10./)
  !
  hbcofa = (/0.003544638, 0.0073888135, 0.013967214, 0.023944625 , 0.03723029 , &
             0.053114605, 0.07005915  , 0.07791257 , 0.07660701  , 0.075071085, &
             0.07326415 , 0.071138385 , 0.068637535, 0.065695415 , 0.062234155, &
             0.058162165, 0.05337168  , 0.047735925, 0.041105755 , 0.0333057  , &
             0.02496844 , 0.01709591  , 0.01021471 , 0.004803175 , 0.00126068 , &
             0 /)
  !
  hbcofb = (/0.         , 0.          , 0.          , 0.         , 0.         , &
             0.         , 0.          , 0.007526545 , 0.023907685, 0.04317925 , &
             0.065851245, 0.092523685 , 0.1239024   , 0.16081785 , 0.204247   , &
             0.2553391  , 0.3154463   , 0.3861593   , 0.4693495  , 0.5672185  , &
             0.67182785 , 0.77060615  , 0.85694605  , 0.9248457  , 0.96929415 , &
             0.9925561  /)
  !
!  write(*,*) hbcofa
!  write(*,*) hbcofb
  rplevo = plevo(nlevo:1:-1)
  !
  open (10,file='T.bin',form='unformatted',access='direct',recl=imax*nlat*nlevi*4)
  read (10,rec=1) dati
  close(10)
!  do k = 1,nlevi
!     r = nlevi - k + 1
!     revdati(:,:,r) = dati(:,:,k)
!  enddo
!  do k = 1,nlevi
!     write(*,'('' DATI '',2f12.4)') minval(dati(:,:,k)),maxval(dati(:,:,k))
!     write(*,'(''      '',i12)')    count(mask=dati(:,:,k)/=spvl)
!  enddo
!  tbot = dati(:,:,nlevi)
!  write(*,'('' TBOT '',2f12.4)') minval(tbot),maxval(tbot)
!  write(*,'(''      '',i12)')    count(mask=tbot/=spvl)
!  open (10,file='PHIS.bin',form='unformatted',access='direct',recl=imax*nlat*4)
!  read (10,rec=1) phis
!  close(10)
!  write(*,'('' PHIS '',2f12.4)') minval(phis),maxval(phis)
!  write(*,'(''      '',i12)')    count(mask=phis/=spvl)
  open (10,file='PS.bin',form='unformatted',access='direct',recl=imax*nlat*4)
  read (10,rec=1) psfc
  close(10)
!  write(*,'('' PSFC '',2f12.4)') minval(psfc),maxval(psfc)
  psfc = psfc / 100.
!  write(*,'('' PSFC '',2f12.4)') minval(psfc),maxval(psfc)
!  write(*,'(''      '',i12)')    count(mask=psfc/=spvl)
  !
  P0     = 1000.
  INTYP  = 1
  ILEV   = 0
  SPVL   = -9999.
  KXTRP  = 0
  VARFLG = 1
  DATO   = SPVL
  !
  call VINTH2PECMWF(DATI,DATO,&
       HBCOFA,HBCOFB,&
       P0,&
       PLEVI,RPLEVO,&
       PSFC,&
       SPVL,&
       IMAX,NLAT,&
       NLEVI,NLEVIP1,NLEVO)
!  write(*,'('' DATO '',2f12.4)') minval(dato),maxval(dato)
!  write(*,'('' DATO '',2f12.4)') minval(dato,mask=dato/=spvl),maxval(dato,mask=dato/=spvl)
!  write(*,'(''      '',i12)')    count(mask=dato/=spvl)
!  do k = 1,nlevo
!     do j = 1,nlat
!        write(*,'(256(f8.2))') (real(dato(i,j,k)),i=1,imax)
!     enddo
!  enddo
  open(20,file='TEST_VINT.da',access='direct',recl=imax*nlat*nlevo*4)
  write(20,rec=1) dato(:,:,nlevo:1:-1)
  close(20)
end program drive_vinth2p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCLFORTSTART
      SUBROUTINE VINTH2PECMWF(DATI,DATO,HBCOFA,HBCOFB,P0,PLEVI,PLEVO, &
                              PSFC,SPVL,IMAX,NLAT, &
                              NLEVI,NLEVIP1,NLEVO)
!
!****     THIS ROUTINE INTERPLOATES CCM2/3 HYBRID COORDINATE DATA
!****     TO PRESSURE COORDINATES USING PRESSURE SURFACES AS THE
!****     COORDINATE SURFACE WHERE THE INTERPOLATION IS DONE.  THE
!****     TYPE OF INTERPOLATION IS CURRENTLY A VARIANT OF TRANSFORMED
!****     PRESSURE COORDINATES WITH THE  INTERPOLATION TYPE
!****     SPECIFIED BY INTYP.  ALL HYBRID COORDINATE VALUES ARE
!****     TRANSFORMED TO PRESSURE VALUES. WHERE THE
!****     FORMULA FOR THE PRESSURE OF A HYBRID SURFACE IS;
!****          P(K) = HBCOFA(LEVH,K)*P0 + HBCOFB(LEVH,K)*PSFC
!****     WHERE,
!****          HBCOFA - IS THE "A" OR PRESSURE HYBRID COEF
!****          LEVH   - IS THE LAYER SURFACE (INTERFACE=1 MIDPOINT=2)
!****          P0     - IS THE BASE PRESSURE IN MB
!****          K      - THE LEVEL INDEX (RUNNING FROM TOP TO BOTTOM)
!****          HBCOFB - IS THE "B" OR SIGMA COEFICIENT
!****          P(K)   - IS THE PRESSURE OF A HYBRID SURFACE IN MB.
!****          PSFC   - IS THE SURFACE PRESSURE IN PASCALS
!****                   (MB = .01*PASCALS
!
!****     FOR HYBRID DATA AT LEVEL INTERFACES SINCE THERE IS ONE
!****     MORE VERTICAL LEVEL FOR INTERFACES THAN FOR LEVEL MIDPOINTS
!****     IT IS ASSUNMED THAT THE FIRST INTERFACE LEVEL WITH A DATA
!****     VALUE IS THE SECOND LEVEL FROM THE TOP.
!
!****     ON INPUT-
!****        DATI    - 3 DIMENSIONAL ARRAY (I,J,KI) CONTAINING DATA
!****                  ON HYBRID SURFACES  WHERE I IS LONGTIUDE, J
!****                  IS LATITUDE AND K IS THE VERTICAL HYBRID
!****                  COORDINATE.  THE VERTICAL DATA RUN TOP TO BOTTOM.
!****                  SIGMA DATA WITH THE DATA ORDERED TOP TO BOTTOM.
!****        HBCOFA  - 2 DIMENSIONAL ARRAY CONTAINING "A" OR PRESSURE
!****                  COEFICIENTS FOR COMPUTING PRESSURE AT A LEVEL.
!****                  ARRAY IS 2XNLEVIP1.  THE 1ST INDEX TAKES ON
!****                  THE VALUE OF EITHER
!****                   1 - FOR LEVEL INTERFACES (OR 1/2 LEVELS) OR;
!****                   2 - FOR LEVEL MIDPOINTS  (OR FULL LEVELS WHERE
!****                       VIRTUALLY ALL VARIABLES ARE LOCATED)
!****                  NOTE THAT COEFICIENTS ARE SCALED TO YIELD A
!****                  PRESSURE IN MB.  THEY ARE ORDERED FROM TOP
!****                  OF THE MODEL TO THE BOTTOM.
!****        HBCOFB  - SAME AS HCOFA BUT FOR THE "B" OR SIGMA COEFICIENT
!****        P0      - BASE PRESSURE IN MB FOR COMPUTING PRESSURE
!****                  OF A HYBRID COORDINATE LEVEL
!****        PLEVI -  1 DIMENSIONAL ARRAY TO HOLD PRESSURE VALUES
!****                  OF HYBRID SURFACES FOR A VERTICAL COLUMN
!****                  SLICE
!****        PLEVO   - LIST OF OUTPUT PRESSURE SURFACES IN MB
!****                  LOW TO HIGH PRESSURE
!****        INTYP   - A FLAG INDICATING INTERPOLATION FOR EACH
!****                  FIELD (1 - LINEAR,2 - LOG ,3 - LOG LOG)
!****                  WHERE EACH INTERPOLATION IS DONE IN TRANSFORMED
!****                  PRESSURE COORDINATES.
!****        ILEV    - FLAG TO SHOW WHETHER FIELD IS ON LEVEL INTERFACE
!****                  1/2 LEVEL WHICH HAS A VALUE OF 1 OR A LEVEL
!****                  MIDPOINT A FULL LEVEL (WHERE ALMOST ALL VARIABLES
!****                  ARE LOCATED
!****        PSFC    - MODEL SFC PRESSURE IN PASCALS (WILL BE CONVERTED
!****                  TO MB)
!****        VCOLI   - ARRAY TO STORE A LONGITUDINAL VERTICAL SLICE OF
!****                  INPUT DATA (IMAX BY NLEVI).
!****        VCOLO   - SAME BUT FOR OUTPUT DATA (IMAX BY NLEVO)
!****        IMAX    - LONGITUDINAL DIMENSION OF THE DATA.
!****        NLAT    - LATITUDINAL DIMENSION OF THE DATA.
!****        NLEVI   - NO. OF LEVELS FOR THE HYBRID DATA
!****        NLEVIP1 - NLEVI + 1
!****        NLEVO   - NUMBER OF OUTPUT LEVELS FOR PRESSURE DATA
!****        KXTRP   - FLAG WHICH INDICATES WHETHER OR NOT
!****                  EXTRAPOLATION WILL BE USED WHEN THE OUTPUT
!****                  PRESSURE SURFACE IS BELOW THE LOWEST LEVEL
!****                  OF THE MODEL.
!****                     0 - DON'T EXTRAPOLATE USE SPECIAL VALUE SPVL
!****                     1 - EXTRAPOLATE DATA using ECMWF formulation
!****                         below PSFC
!****        SPVL    - SPECIAL VALUE TO USE WHEN DATA IS NOT
!****                  EXTRAPOLATED
!****        varflg  - flag which indicates the name of the variable
!****                  -1 means geopotential (Z)
!****                  +1 means geopotential (T)
!****                   0 any other variable
!****        tbot    - temperature at level closest to ground
!****        phis    - surface geopotential
!
!****     ON OUTPUT-
!****        DATO  - 3 DIMENSIONAL ARRAY TO HOLD DATA INTERPOLATED
!****                TO PRESSURE SURFACES.
!
      IMPLICIT NONE
!      DOUBLE PRECISION SPVL,DATI,DATO,HBCOFA,HBCOFB,PLEVI,PLEVO,PSFC,P0
!      DOUBLE PRECISION A2LN,A1
      REAL SPVL,DATI,DATO,HBCOFA,HBCOFB,PLEVI,PLEVO,PSFC,P0
      INTEGER IMAX,NLAT,NLEVI,NLEVIP1,NLEVO
      INTEGER I,J,K,KP,KPI
!
      DIMENSION DATI(IMAX,NLAT,NLEVI),DATO(IMAX,NLAT,NLEVO), &
               HBCOFA(NLEVIP1),HBCOFB(NLEVIP1),PLEVI(NLEVIP1), &
               PLEVO(NLEVO),PSFC(IMAX,NLAT)
! NCLEND
!****
!
      DO 70 J = 1,NLAT
          DO 60 I = 1,IMAX

!
!****     GET PRESSURE VALUES FOR HYBRID SURFACES FOR THIS POINT
!****     AND FOR THE TYPE OF MODEL SURFACE THE DATA IS ON.
!****     INTERFACE DATA STARTS AT THE SECOND INTERFACE LEVEL SO
!****     IF THE DATA IS ON THOSE LEVELS START THE
!
              DO K = 1,NLEVI
                  KPI = K
                  PLEVI(K) = (HBCOFA(KPI)*P0) + HBCOFB(KPI)* (PSFC(I,J))
              END DO
!
!****     CALL P2HBD TO PERFORM VERTICAL INTERP. THEN TRANSFER DATA TO
!****     THE OUTPUT ARRAY
!
              DO 50 K = 1,NLEVO
!
!****     CHECK FOR BRACKETING LEVEL KP WILL BE THE INPUT LEVEL THAT
!****     IS THE UPPER PORTION OF 2 INPUT BRACKETING LEVELS.
!
!****     IF BRANCH FOR MODEL TOP
!
                  IF (PLEVO(K).LE.PLEVI(1)) THEN
                      KP = 1
                      GO TO 30
!
!****     IF BRANCH FOR LEVEL BELOW LOWEST HYBRID LEVEL
!
                  ELSE IF (PLEVO(K).GT.PLEVI(NLEVI)) THEN
                     DATO(I,J,K) = SPVL
                     GO TO 40
!
!****     IF BRANCH FOR TO CHECK IF OUTPUT LEVEL IN BETWEEN
!****     2 LOWEST HYBRID LEVELS
!
                  ELSE IF (PLEVO(K).GE.PLEVI(NLEVI-1)) THEN
                      KP = NLEVI - 1
                      GO TO 30
!
!****     IF BRANCH FOR MODEL INTERIOR
!****     LOOP THROUGH INPUT LEVELS TILL YOU ARE BRACKETING
!****     OUTPUT LEVEL
!
                  ELSE
                      KP = 0
   20                 CONTINUE
                      KP = KP + 1
                      IF (PLEVO(K).LE.PLEVI(KP+1)) GO TO 30
                      GO TO 20
                  END IF
   30             CONTINUE
!
!****     LEVEL BRACKETED PICK TYPE OF INTERP.
!
!
!****     LINEAR INTERP.
!
                  DATO(I,J,K) = DATI(I,J,KP) + &
                       (DATI(I,J,KP+1)-DATI(I,J,KP))* &
                       (PLEVO(K)-PLEVI(KP))/ &
                       (PLEVI(KP+1)-PLEVI(KP))
   40             CONTINUE
   50         CONTINUE
   60     CONTINUE
   70 CONTINUE
      RETURN
      END
