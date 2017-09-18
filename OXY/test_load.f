      PARAMETER(NSP = 13)
      DIMENSION ISPEC(NSP)

      OPEN(UNIT=1,FILE='species_matrix.dat')

      READ(1,201)ISPEC
      WRITE(*, 201)ISPEC

  201 FORMAT(A9)
      STOP
      END
