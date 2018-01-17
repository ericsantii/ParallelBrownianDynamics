SUBROUTINE ClusterStats
  USE Globals
  USE VectorFunctions
  IMPLICIT NONE
  INTEGER :: i, j
  REAL(8) :: RIJ2
  
  IF ( mod( gOutClusCounter, gOutClust ) .eq. 0 ) THEN 

     CALL ResetValues
     DO i = 1, gNPart - 1
        DO j = i + 1, gNPart
           RIJ2 = sqrdist(i,j)  
           IF ( RIJ2 .LE. RB2 ) CALL ClusterIndex(i,j)
        END DO
     END DO
     
     CALL CalculateNumClusters       
     CALL ClusterMass               
     CALL ClusterDist
     CALL CenterOfMass
     CALL EffRad               
     CALL RadiusGyration             
     CALL CalculateMeanValues
     CALL WriteOutput

     IF ( gT .GE. gAvgTiming ) THEN
        CALL CalcEquilValues 
        CALL RadDist
        gEQLCOUNT = gEQLCOUNT + 1 
     END IF

  END IF
  gOutClusCounter = gOutClusCounter + 1
END SUBROUTINE ClusterStats


SUBROUTINE ClusterIndex(I,J)
  USE Globals
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER :: OLDSMALLERINDEX, OLDSMALLERINDEX2,  NEWOLDSMALLERINDEX
  
  IF( CLUSTER_INDEX_LIST(J) .EQ. J ) THEN
     CLUSTER_INDEX_LIST(J) = CLUSTER_INDEX_LIST(I)
  ELSE
     OLDSMALLERINDEX    = CLUSTER_INDEX_LIST(J)
     OLDSMALLERINDEX2   = CLUSTER_INDEX_LIST(I)
     NEWOLDSMALLERINDEX = min0( OLDSMALLERINDEX, OLDSMALLERINDEX2, CLUSTER_INDEX_LIST(OLDSMALLERINDEX) )
     
     DO K = 1, gNPart
        IF ( (CLUSTER_INDEX_LIST(K) .EQ. OLDSMALLERINDEX) .OR. (CLUSTER_INDEX_LIST(K) .EQ. OLDSMALLERINDEX2) ) THEN
           CLUSTER_INDEX_LIST(K) = NEWOLDSMALLERINDEX
        ENDIF
     ENDDO
  ENDIF
  
END SUBROUTINE ClusterIndex


SUBROUTINE CalculateNumClusters
  USE Globals
  IMPLICIT NONE
  INTEGER(8) :: I, J
  INTEGER(8) ::  COUNTER
  
  COUNTER = 1
  
  DO i = 1, gNPart - 1
     DO j = i + 1, gNPart
        IF( CLUSTER_INDEX_LIST(j) - COUNTER == 1 ) THEN
           COUNTER = COUNTER + 1
        ELSEIF( CLUSTER_INDEX_LIST(j) - COUNTER .GT. 1 ) THEN
           CLUSTER_INDEX_LIST(j) = CLUSTER_INDEX_LIST(j) - 1
        END IF
     END DO
  END DO
  
  NC = MAXVAL(CLUSTER_INDEX_LIST)
  
END SUBROUTINE CalculateNumClusters

SUBROUTINE ClusterMass
  USE GLOBALS
  IMPLICIT NONE
  INTEGER :: i, CLUSTER_LABEL
  
  DO i = 1, gNPart
     CLUSTER_LABEL = CLUSTER_INDEX_LIST(i)
     NpC(CLUSTER_LABEL) = NpC(CLUSTER_LABEL) + 1
  END DO
  
END SUBROUTINE ClusterMass

SUBROUTINE CLUSTERDIST
  USE GLOBALS
  IMPLICIT NONE
  INTEGER :: i, j, MAXCLUST
  REAL(8) :: NOSINGLETS
  
  MAXCLUST = MAXVAL( NpC )
  
  DO i = 1, MAXCLUST
     DO j = 1, NC
        IF ( NpC(j) == i ) CLUSTDIST(i) = CLUSTDIST(i) + 1
     END DO
  END DO
  
  SINGLETS = CLUSTDIST(1)
  !DOUBLETS = CLUSTDIST(2)
  !TRIPLETS = CLUSTDIST(3)
  !QUADRUPLETS = CLUSTDIST(4)
  !PENTUPLETS = CLUSTDIST(5)
! density of populations
 !nSINGLETS = dble(SINGLETS)/gAreaVolumeNoDim 
 !nDOUBLETS = dble(DOUBLETS)/gAreaVolumeNoDim
 !nTRIPLETS = dble(TRIPLETS)/gAreaVolumeNoDim
 !nQUADRUPLETS = dble(QUADRUPLETS)/gAreaVolumeNoDim
 !nPENTUPLETS = dble(PENTUPLETS)/gAreaVolumeNoDim 
 NOSINGLETS = NC - SINGLETS
WRITE(5, 1205) gT, NOSINGLETS
1205 FORMAT( 2(1X, E15.7E3) )
  
END SUBROUTINE CLUSTERDIST

SUBROUTINE CenterOfMass
  USE GLOBALS
  IMPLICIT NONE
  INTEGER :: i, j, CLUSTER_LABEL
  
  DO i = 1, gNPart
     CLUSTER_LABEL = CLUSTER_INDEX_LIST(i)
     RCM(CLUSTER_LABEL,:) = RCM(CLUSTER_LABEL,:) + gR(i,:)
  END DO
  
  DO i = 1, gNPart
     IF ( NpC(i) .NE. 0 ) RCM(i,:) = RCM(i,:) / DFLOAT( NpC(i) )
  END DO
  
END SUBROUTINE CenterOfMass


SUBROUTINE RadiusGyration
  USE GLOBALS
  USE VectorFunctions
  IMPLICIT NONE
  INTEGER :: i, CLUSTER_LABEL
  REAL(8) :: DR2
  REAL(8), DIMENSION(gND) :: DR
  
  DO i = 1, gNPart
     CLUSTER_LABEL = CLUSTER_INDEX_LIST(i)
     DR  = gR(i,:) - RCM(CLUSTER_LABEL,:)
     DR2 = vdot(DR,DR)
     RG2(CLUSTER_LABEL) = RG2(CLUSTER_LABEL) + DR2
  END DO
  
  DO I = 1,NC
     RG2(I) = RG2(I) / DFLOAT( NpC(I) )
  END DO
  !RG2(:) = DSQRT(RG2(:))
  
END SUBROUTINE RadiusGyration


SUBROUTINE CalculateMeanValues
  USE GLOBALS
  IMPLICIT NONE
  INTEGER :: I, S, MAXCLUST
  REAL(8) ::  SUM_NPC2, SUM_RG2NPC, SUM_REFFNPC, SUM_NPC
  REAL(8) :: SUM0, SUM1, SUM2
  
  SUM0 = 0.D0
  SUM1 = 0.D0
  SUM2 = 0.D0

  MAXCLUST    = MAXVAL( NpC )
  SUM_NPC     = DFLOAT(SUM( NpC ))
  SUM_NPC2    = SUM( NpC * NpC )
  SUM_RG2NPC  = SUM( NpC * RG2 )
  SUM_REFFNPC = SUM( NpC * REFF)

  DO S = 1, MAXCLUST
     SUM0 = SUM0 + DFLOAT( CLUSTDIST( S ) )
     SUM2 = SUM2 + DFLOAT( S * S ) * DFLOAT( CLUSTDIST( S ) )
     SUM1 = SUM1 + DFLOAT( S ) * DFLOAT( CLUSTDIST( S ) )
  END DO  
  
  PROB(:)   = DFLOAT(CLUSTDIST(:)) / SUM0   ! Probability Disctribution of clusters: P(s,t)
  NpC_MEAN  = SUM_NPC2 / SUM_NPC    ! Weight averaged mean cluster size  ( Doi-Chen def.)
  RG2_MEAN  = SUM_RG2NPC  / SUM_NPC ! Weight averaged mean Rad. gyr.    ( Doi-Chen def.)
  REFF_MEAN = SUM_REFFNPC / SUM_NPC ! Mean effective( hydrodynamic ) radius
  S_MEAN    = SUM2 / SUM1           ! Mean cluster size
  L_MEAN    = SUM1 / SUM0           ! Mean cluster lenght

  EPS   = 1.D0 - DFLOAT(NC) / DFLOAT(gNPart)
  THETA = 1.D0 - 1.D0 / S_MEAN
  
END SUBROUTINE CalculateMeanValues


SUBROUTINE ResetValues
  USE GLOBALS
  IMPLICIT NONE
  INTEGER :: i

  NC           = 0
  CLUSTDIST(:) = 0
  NpC(:)       = 0
  RG2(:)       = 0.D0
  REFF(:)      = 0.D0
  RCM(:,:)     = 0.D0

  DO i = 1, gNPart
     CLUSTER_INDEX_LIST(i) = i
  END DO

END SUBROUTINE ResetValues


SUBROUTINE WriteOutput
  USE GLOBALS
  IMPLICIT NONE
  INTEGER :: S, I, MAXCLUST
  MAXCLUST    = MAXVAL( NpC )
  WRITE(3,1204) gT, NC , SINGLETS, EPS, THETA, NPC_MEAN, RG2_MEAN, REFF_MEAN, S_MEAN, L_MEAN
  1204 FORMAT( 1X, E15.7E3, 1X, 2(I4), 7(1X, E15.7E3) )

   !WRITE(9,1209) gT, NC , nSINGLETS, nDOUBLETS, nTRIPLETS, nQUADRUPLETS, nPENTUPLETS
  !1209 FORMAT( 1X, E15.7E3, 1X, 1(I4), 5(1X, E15.7E3) )

  
  !DO S = 1, MAXCLUST
   !  WRITE(5,1205)  S, CLUSTDIST( S )
     !WRITE(7,1100)  gT, S, PROB(S)
  !END DO
  
  !WRITE(5,*) '-----'
  1205 FORMAT( 2(1X, I4) )
  
  !WRITE(7,*) '-----'
  1100 FORMAT( 1X, E15.7E3, 1X, I4, 1X, E15.7E3 ) 

END SUBROUTINE WriteOutput

SUBROUTINE CalcEquilValues
  USE GLOBALS
  IMPLICIT NONE
  INTEGER :: I

  !gEQLCOUNT = gEQLCOUNT + 1

  NPC_EQL  = NPC_EQL  + NPC_MEAN
  RG2_EQL  = RG2_EQL  + RG2_MEAN
  REFF_EQL = REFF_EQL + REFF_MEAN

  S_EQL = S_EQL + S_MEAN
  L_EQL = L_EQL + L_MEAN
   
  EPS_EQL   = EPS_EQL   + EPS
  THETA_EQL = THETA_EQL + THETA

  DO I = 1, NC
     WRITE(6,1206)  RG2( I ), NPC( I )
  END DO

  !WRITE(6,*) '-----'
  1206 FORMAT( 1X, E15.7E3, 1X, I4 )

END SUBROUTINE CalcEquilValues

SUBROUTINE RadDist
  USE Globals
  USE VectorFunctions
  IMPLICIT NONE

  INTEGER :: i, j, BIN
  REAL(8) :: RIJM, RIJ2
  REAL(8), DIMENSION(gND) :: RIJ

  DO i = 1, gNPart - 1
     DO j = i + 1, gNPart
        CALL dist(i,j,Rij,RIJ2)
        IF ( RIJ2 .GT. gCutOff2 ) CYCLE
        RIJM = DSQRT(RIJ2)
        BIN  = INT( RIJM / gDR ) + 1
        IF (BIN .LE. gRDBins ) gRD(BIN) = gRD(BIN) + 2.D0
     END DO
  END DO

end subroutine RadDist

SUBROUTINE FinalEquilValues
  USE GLOBALS
  IMPLICIT NONE
  INTEGER :: BIN
  REAL(8) :: rlower, rupper, nideal, const, RBIN
  CONST = gPi43 * gRho

  NPC_EQL  =  NPC_EQL  / DFLOAT(gEQLCOUNT)
  RG2_EQL  =  RG2_EQL  / DFLOAT(gEQLCOUNT)
  REFF_EQL =  REFF_EQL / DFLOAT(gEQLCOUNT)

  S_EQL = S_EQL / DFLOAT(gEQLCOUNT)
  L_EQL = L_EQL / DFLOAT(gEQLCOUNT)

  EPS_EQL   = EPS_EQL   / DFLOAT(gEQLCOUNT)
  THETA_EQL = THETA_EQL / DFLOAT(gEQLCOUNT)
    
  DO BIN = 1, gRDBins
     RLOWER = DFLOAT(BIN-1) * gDR
     RUPPER = RLOWER + gDR
     NIDEAL = CONST * ( RUPPER ** 3 - RLOWER ** 3 )
     RBIN   = (DFLOAT(BIN)) * gDR
     gRD(BIN) = gRD(BIN) / DFLOAT(gEQLCOUNT) / DFLOAT(gNPart) / NIDEAL
     WRITE( 8, * ) RBIN , gRD(BIN)
  END DO
  
  WRITE(4,1203) EPS_EQL, THETA_EQL, NPC_EQL, RG2_EQL, REFF_EQL, S_EQL, L_EQL
  1203 FORMAT(7(1X, E15.7E3) )

END SUBROUTINE FInalEquIlValues


SUBROUTINE EffRad
  USE GLOBALS
  USE VectorFunctions
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: rij2, rij2max
  
  DO K = 1, NC
    rij2max = 0.D0
    DO I = 1, gNPart - 1
       DO J = I + 1, gNPart
          IF ( ( CLUSTER_INDEX_LIST(i) .EQ. K ) .AND. ( CLUSTER_INDEX_LIST(j) .EQ. K ) ) THEN
             rij2 = sqrdist(i,j)
             IF (rij2 .GT. rij2max) rij2max = rij2
          END IF
       END DO
    END DO   
  Reff(K) = ( DSQRT(rij2max) + 2.D0 ) / 2.D0
  END DO

END SUBROUTINE EffRad
