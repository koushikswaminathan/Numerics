program RungeKutta
  implicit none
  real:: Y1, H, Y2, F1, F2, F1_FN, F2_FN, PI, T
  integer:: I, J
  integer, parameter :: N = 1000, M = 100, L = 1
  real:: K11, K12, K13, K14, K21, K22, K23, K24
  real, dimension (2, N):: Y
  real, parameter :: G = 9.8

  PI = 4.0*ATAN(1.0)

! Time Step

  H = PI/M


! Intial Values

  Y(1,1) = 0.0
  Y(2,1) = 3.0

! Open Output File

  open(1, file = 'pendulum.dat', status = 'old')


! Runge-Kutta Method for Two Differential Equation Functions F1_FN, F2_FN.



  do I = 1, N - 1
    T = H*I
    Y1 = Y(1,I)
    Y2 = Y(2,I)
    K11 = H*F1_FN(Y1,Y2,T)
    K21 = H*F2_FN(Y1,Y2,T)
    K12 = H*F1_FN((Y1+K11/2.0),(Y2+K21/2.0),(T+H/2.0))
    K22 = H*F2_FN((Y1+K11/2.0),(Y2+K21/2.0),(T+H/2.0))
    K13 = H*F1_FN((Y1+K12/2.0),(Y2+K22/2.0),(T+H/2.0))
    K23 = H*F2_FN((Y1+K12/2.0),(Y2+K22/2.0),(T+H/2.0))
    K14 = H*F1_FN((Y1+K13),(Y2+K23),(T+H))
    K24 = H*F2_FN((Y1+K13),(Y2+K23),(T+H))
    Y(1,I+1) = Y(1,I)+(K11+2.0*(K12+K13)+K14)/6.0
    Y(2,I+1) = Y(2,I)+(K21+2.0*(K22+K23)+K24)/6.0

!  Return Y1 back to typical Theta range [-PI, PI]

    Y(1,I+1) = Y(1,I+1)-2.0*PI*NINT(Y(1,I+1)/(2.0*PI))

! Write to File

    write(1,*) I, T, Y(1, I), Y(2, I)

  end do

!  Close File

  close(1)


end program RungeKutta

! Function assignments for Second order Differential equation decomposed in to 2 First order equations.

  function F1_FN (Y1, Y2, T) RESULT (F1)
    implicit none
    real :: Y1, Y2, T, F1
    F1 = Y2
  end function F1_FN

  function F2_FN (Y1, Y2, T) RESULT (F2)
    implicit none
    real, parameter :: G = 9.8, L = 1
    real :: Y1, Y2, T, F2
    F2 = - (G/L)*SIN(Y1)
  end function F2_FN
