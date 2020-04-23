!     Program:         CoxORC.f90
!     Written by:      Xin Zhou
!     Last modified:   Dec 26, 2015
!     Purpose: Cox proportional hazard models with ORC

subroutine der_likelihood_time(mu,beta,gammaobj,tau2, z0, z1, XX, JJ, KK, a, b, &
                        mincomp, maxcomp, GQ, GQX, GQW, derlikelihood, prob)
    implicit none
    ! ---- arg types -----------------------
    integer :: JJ, KK
    double precision :: mu, beta, tau2
    double precision :: gammaobj(JJ)
    ! true values of mu, beta, gammaobj, tau2
    integer :: z0(JJ), z1(JJ)   ! z0: # of subjects with outcome 0; z1: # of subjects with outcome 1.
    integer :: XX(JJ)     ! treatment assignment
    double precision :: a, b     ! integration limits: a - lower, b - upper
    integer :: mincomp(JJ+2), maxcomp(JJ+2)  ! gammaobj(1), ..., gammaobj(JJ), mu, beta
    integer :: GQ   ! number of GQ points
    double precision :: GQX(GQ), GQW(GQ)
    double precision :: derlikelihood(JJ+2)
    double precision :: prob
    ! --------------------------------------
    double precision, parameter :: PI = 3.14159265358979323846   ! pi

    double precision :: likelihoodf_denom, likelihoodf_numer, likelihoodf_denomb2
    double precision :: x
    integer :: i, j, k
    double precision :: ff0, ff1, ff01
    double precision :: ff, ffprob
    double precision :: ff_mu, ff_beta
    double precision :: derlikelihood_mu, derlikelihood_beta, derlikelihood_tau2
    double precision :: ff_gammaobj(JJ-1)
    double precision :: derlikelihood_gammaobj(JJ-1)
    double precision :: temp
    double precision :: eaa, ebb, exx
    double precision :: faeaa, fbebb

  !  open(unit=1, file="debug.txt", status='replace')

    likelihoodf_denom = 0.0d0
    likelihoodf_numer = 0.0d0
    likelihoodf_denomb2 = 0.0d0
    derlikelihood_mu = 0.0d0
    derlikelihood_beta = 0.0d0
    derlikelihood_tau2 = 0.0d0
    derlikelihood_gammaobj = 0.0d0
    
   !     write (unit=1, fmt=*) likelihoodf_denom
        
    prob = 0.0d0
    do i=1,GQ
        x = GQX(i)
        exx = exp(-0.5d0*x*x/tau2)
        
        ff = 1.0d0
        ffprob = 1.0d0
        ff_mu = 0.0d0
        ff_beta = 0.0d0
        do j=1,JJ
            ff1 = mu+beta*XX(j)+gammaobj(j)+x
            ff0 = 1-ff1
            ff = ff*(ff0**z0(j))*(ff1**z1(j))

            ! since GQX are not at limits, we ignore the cases where ff0=0 or ff1=0
            temp = z1(j)/ff1 - z0(j)/ff0
            ff_mu = ff_mu + temp
            ff_beta = ff_beta + temp*XX(j)
            k = j-1
            if (k>0) then
                ff_gammaobj(k) = temp
            end if

            ff01 = ff0*ff1
            ! compute binomial
            ! compute combination number with power of ff0 and ff1
            ! to avoid overflow
            if (z0(j)<z1(j)) then
                ffprob = ffprob*ff1**(z1(j)-z0(j))
                do k=0,(z0(j)-1)
                    ffprob = ffprob*dble(KK-k)/dble(z0(j)-k)*ff01
                end do
            else
                ffprob = ffprob*ff0**(z0(j)-z1(j))
                do k=0,(z1(j)-1)
                    ffprob = ffprob*dble(KK-k)/dble(z1(j)-k)*ff01
                end do
            end if
        end do
        prob = prob + GQW(i)*ffprob*exx
        likelihoodf_denom = likelihoodf_denom + GQW(i)*exx
        likelihoodf_numer = likelihoodf_numer + GQW(i)*ff*exx
        likelihoodf_denomb2 = likelihoodf_denomb2 + GQW(i)*x*x*exx

        derlikelihood_mu = derlikelihood_mu + GQW(i)*ff*ff_mu*exx
        derlikelihood_beta = derlikelihood_beta + GQW(i)*ff*ff_beta*exx
        derlikelihood_gammaobj = derlikelihood_gammaobj + GQW(i)*ff*ff_gammaobj*exx
        derlikelihood_tau2= derlikelihood_tau2 + GQW(i)*ff*x*x*exx
    end do
    


    ! calculate f(a)exp(-0.5*a*a)
    eaa = exp(-0.5d0*a*a/tau2)
    ff = 1.0d0
    do j=1,JJ
        ff1 = mu+beta*XX(j)+gammaobj(j)+a
        ff0 = 1-ff1
        ff = ff*(ff0**z0(j))*(ff1**z1(j))
    end do
    faeaa = ff*eaa
    ! calculate f(b)exp(-0.5*b*b)
    ebb = exp(-0.5d0*b*b/tau2)
    ff = 1.0d0
    do j=1,JJ
        ff1 = mu+beta*XX(j)+gammaobj(j)+b
        ff0 = 1-ff1
        ff = ff*(ff0**z0(j))*(ff1**z1(j))
    end do
    fbebb = ff*ebb

    ! calculate derlikelihood_mu
    derlikelihood_mu = derlikelihood_mu + faeaa*dble(mincomp(JJ+1)) &
                        - fbebb*dble(maxcomp(JJ+1))
    derlikelihood_mu = derlikelihood_mu / likelihoodf_numer - &
            (eaa*dble(mincomp(JJ+1)) - ebb*dble(maxcomp(JJ+1))) / likelihoodf_denom
    ! calculate derlikelihood_beta
    derlikelihood_beta = derlikelihood_beta + faeaa*dble(mincomp(JJ+2)) &
                        - fbebb*dble(maxcomp(JJ+2))
    derlikelihood_beta = derlikelihood_beta / likelihoodf_numer - &
                    (eaa*dble(mincomp(JJ+2)) - ebb*dble(maxcomp(JJ+2))) / likelihoodf_denom
    ! calculate derlikelihood_gammaobj
    do j=2,JJ
        k = j-1
        derlikelihood_gammaobj(k) = derlikelihood_gammaobj(k) + &
                            faeaa*dble(mincomp(j)) - fbebb*dble(maxcomp(j))
        derlikelihood_gammaobj(k) = derlikelihood_gammaobj(k)/likelihoodf_numer &
                        -(eaa*dble(mincomp(j))-ebb*dble(maxcomp(j)))/likelihoodf_denom
    end do

    ! calculate derlikelihood_tau2
    derlikelihood_tau2 = 0.5d0*(derlikelihood_tau2/likelihoodf_numer- &
                    likelihoodf_denomb2/likelihoodf_denom)/tau2/tau2
    prob = prob/likelihoodf_denom

    derlikelihood(1) = derlikelihood_mu
    derlikelihood(2) = derlikelihood_beta
    derlikelihood(3:(JJ+1)) = derlikelihood_gammaobj
    derlikelihood(JJ+2) = derlikelihood_tau2
    
  !  close(unit=1)

end subroutine der_likelihood_time



subroutine der_likelihood_notime(mu, beta, tau2, z00, z01, z10, z11, GQ, GQX, GQW, &
    derlikelihood_mu, derlikelihood_beta, derlikelihood_tau2, prob)
    implicit none
    ! ---- arg types -----------------------
    double precision :: mu, beta, tau2
    ! true values of beta0, beta1, and tau2
    integer :: z00, z01, z10, z11
    integer :: GQ   ! number of GQ points
    double precision :: GQX(GQ), GQW(GQ)
    double precision :: derlikelihood_mu
    double precision :: derlikelihood_beta
    double precision :: derlikelihood_tau2
    double precision :: prob
    ! ---------------------------------------
    double precision :: likelihoodf_denom, likelihoodf_numer, likelihoodf_denomb2

    double precision :: ff00, ff01, ff10, ff11
    double precision :: ff, ff1, ffprob, ff0prob, ff1prob
    double precision :: exx, exx1, exx2, temp
    integer :: z0, z1
    integer :: i, k
    double precision :: x

    double precision, parameter :: PI = 3.14159265358979323846   ! pi

    z0 = z00 + z01
    z1 = z10 + z11
    derlikelihood_mu = 0.0d0
    derlikelihood_beta = 0.0d0
    derlikelihood_tau2 = 0.0d0
    likelihoodf_denom = 0.0d0
    likelihoodf_denomb2 = 0.0d0
    likelihoodf_numer = 0.0d0
    prob = 0.0d0
    do i=1,GQ
        x = GQX(i)
        ff = 1.0d0
        ffprob = 1.0d0
        ff01 = mu+x
        ff00 = 1-ff01
        ff11 = mu+beta+x
        ff10 = 1-ff11
        exx = exp(-0.5d0*x*x/tau2)

        ff = (ff00**z00)*(ff01**z01)*(ff10**z10)*(ff11**z11)
        likelihoodf_numer = likelihoodf_numer + GQW(i) * ff * exx
        likelihoodf_denom = likelihoodf_denom + GQW(i) * exx
        likelihoodf_denomb2 = likelihoodf_denomb2 + GQW(i) * x * x * exx

        ! since GQX are not at limits, we ignore the cases where ff00=0 or ff01=0 or ff10=0 or ff11=0
        ff1 = ff*(z01/ff01-z00/ff00+z11/ff11-z10/ff10)
        derlikelihood_mu = derlikelihood_mu + GQW(i) * ff1 * exx
        ff1 = ff*(z11/ff11-z10/ff10)
        derlikelihood_beta = derlikelihood_beta + GQW(i) * ff1 * exx
        ff1 = ff*x*x
        derlikelihood_tau2 = derlikelihood_tau2 + GQW(i) * ff1 * exx

        ! prob
        ff0prob = ff00*ff01
        ff1prob = ff10*ff11
        ! compute combination number with power of ff00 and ff01
        ! to avoid overflow
        if (z00<z01) then
            ffprob = ffprob*ff01**(z01-z00)
            do k=0,(z00-1)
                ffprob = ffprob*dble(z0-k)/dble(z00-k)*ff0prob
            end do
        else
            ffprob = ffprob*ff00**(z00-z01)
            do k=0,(z01-1)
                ffprob = ffprob*dble(z0-k)/dble(z01-k)*ff0prob
            end do
        end if
        ! compute combination number with power of ff10 and ff11
        ! to avoid overflow
        if (z10<z11) then
            ffprob = ffprob*ff11**(z11-z10)
            do k=0,(z10-1)
                ffprob = ffprob*dble(z1-k)/dble(z10-k)*ff1prob
            end do
        else
            ffprob = ffprob*ff10**(z10-z11)
            do k=0,(z11-1)
                ffprob = ffprob*dble(z1-k)/dble(z11-k)*ff1prob
            end do
        end if
        prob = prob + GQW(i) * ffprob * exx
    enddo

    if (beta>=0) then
        ! we don't consider the case of beta = 0
        exx1 = exp(-0.5d0*mu*mu/tau2)       ! exx on lower limits
        exx2 = exp(-0.5d0*(1-mu-beta)*(1-mu-beta)/tau2)
        if (z01==0) then
            derlikelihood_mu = derlikelihood_mu + ((1-beta)**z10)*(beta**z11)*exx1
        end if
        if (z10==0) then
            temp = ((1-beta)**z01)*(beta**z00)*exx2
            derlikelihood_mu = derlikelihood_mu - temp
            derlikelihood_beta = derlikelihood_beta - temp
        end if
        derlikelihood_mu = derlikelihood_mu / likelihoodf_numer - (exx1-exx2)/likelihoodf_denom
        derlikelihood_beta = derlikelihood_beta / likelihoodf_numer + exx2/likelihoodf_denom
    else
        ! beta < 0
        exx1 = exp(-0.5d0*(mu+beta)*(mu+beta)/tau2)       ! exx on lower limits
        exx2 = exp(-0.5d0*(1-mu)*(1-mu)/tau2)
        if (z00==0) then
            derlikelihood_mu = derlikelihood_mu - ((-beta)**z10)*((1+beta)**z11)*exx2
        end if
        if (z11==0) then
            temp = ((-beta)**z01)*((1+beta)**z00)*exx1
            derlikelihood_mu = derlikelihood_mu + temp
            derlikelihood_beta = derlikelihood_beta + temp
        end if
        derlikelihood_mu = derlikelihood_mu / likelihoodf_numer - (exx1-exx2)/likelihoodf_denom
        derlikelihood_beta = derlikelihood_beta / likelihoodf_numer + exx1/likelihoodf_denom
    end if

    derlikelihood_tau2 = 0.5*(derlikelihood_tau2 / likelihoodf_denom - likelihoodf_denomb2 / likelihoodf_denom)/tau2/tau2
    ! likelihoodf = likelihoodf/sqrt(pi)
    prob = prob / likelihoodf_denom
end subroutine der_likelihood_notime



subroutine vectorsquare(a,n,c)
!============================================================
! Calculate a*a', where a is a vector
! ----------------------------------------------------------
! input ...
! a(n) - vector a
! n    - dimension
! output ...
! c(n,n) - a*a'
!!!!!  next step: modify c to storage compression format
!!!!!  c is a positive definite matrix
!!!!!    stored by rows in lower triangular form 
!!!!!    as a one dimensional array, in the sequence
!===========================================================
implicit none 
integer n
double precision a(n), c(n,n)

integer i, j

do i=1,n-1
    c(i,i) = a(i)*a(i)
    do j=i+1,n
        c(i,j) = a(i)*a(j)
        c(j,i) = c(i,j)
    end do
end do
c(n,n) = a(n)*a(n)
end subroutine vectorsquare



subroutine syminverse(a,c,n)
!============================================================
! Inverse symmetric matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! Both a and c are symmetric
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)

double precision aa(n*(n+1)/2), cc(n*(n+1)/2)
integer i, j, k
integer nullty, ifault

k = 0
do i=1,n
	do j=1,i
		k = k + 1
		aa(k) = a(i,j)
	end do
end do
call syminv(aa, n, cc, nullty, ifault)
k = 0
do i = 1, n
	do j=1,i-1
		k = k + 1
		c(i,j) = cc(k)
		c(j,i) = cc(k)
	end do
	k = k + 1
	c(i,i) = cc(k)
end do
end subroutine syminverse


 subroutine syminv ( a, n, c, nullty, ifault )

!*****************************************************************************80
!
!! SYMINV computes the inverse of a symmetric matrix.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    FORTRAN77 version by Michael Healy
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Michael Healy,
!    Algorithm AS 7:
!    Inversion of a Positive Semi-Definite Symmetric Matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 198-199.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A((N*(N+1))/2), a positive definite matrix stored
!    by rows in lower triangular form as a one dimensional array, in the 
!    sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Output, real ( kind = 8 ) C((N*(N+1))/2), the inverse of A, or generalized
!    inverse if A is singular, stored using the same storage scheme employed
!    for A.  The program is written in such a way that A and U can share 
!    storage.
!
!    Workspace, real ( kind = 8 ) W(N).
!
!    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.  If NULLTY 
!    is zero, the matrix is judged to have full rank.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no error detected.
!    1, N < 1.
!    2, A is not positive semi-definite.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) c((n*(n+1))/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mdiag
  integer ( kind = 4 ) ndiag
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nullty
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x

  ifault = 0

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  nrow = n
!
!  Compute the Cholesky factorization of A.
!  The result is stored in C.
!
  nn = ( n * ( n + 1 ) ) / 2

  call cholesky ( a, n, nn, c, nullty, ifault )

  if ( ifault /= 0 ) then
    return
  end if
!
!  Invert C and form the product (Cinv)' * Cinv, where Cinv is the inverse
!  of C, row by row starting with the last row.
!  IROW = the row number,
!  NDIAG = location of last element in the row.
!
  irow = nrow
  ndiag = nn

  do
!
!  Special case, zero diagonal element.
!
    if ( c(ndiag) == 0.0D+00 ) then

      l = ndiag
      do j = irow, nrow
        c(l) = 0.0D+00
        l = l + j
      end do

    else

      l = ndiag
      do i = irow, nrow
        w(i) = c(l)
        l = l + i
      end do

      icol = nrow
      jcol = nn
      mdiag = nn

      do

        l = jcol

        if ( icol == irow ) then
          x = 1.0D+00 / w(irow)
        else
          x = 0.0D+00
        end if

        k = nrow

        do while ( irow < k )

          x = x - w(k) * c(l)
          k = k - 1
          l = l - 1

          if ( mdiag < l ) then
            l = l - k + 1
          end if

        end do

        c(l) = x / w(irow)

        if ( icol <= irow ) then
          exit
        end if

        mdiag = mdiag - icol
        icol = icol - 1
        jcol = jcol - 1

      end do

    end if

    ndiag = ndiag - irow
    irow = irow - 1

    if ( irow <= 0 ) then
      exit
    end if

  end do

  return
end



subroutine cholesky ( a, n, nn, u, nullty, ifault )

!*****************************************************************************80
!
!! CHOLESKY computes the Cholesky factorization of a PDS matrix.
!
!  Discussion:
!
!    For a positive definite symmetric matrix A, the Cholesky factor U
!    is an upper triangular matrix such that A = U' * U.
!
!    This routine was originally named "CHOL", but that conflicted with
!    a built in MATLAB routine name.
!
!    The missing initialization "II = 0" has been added to the code.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Michael Healy.
!    Modifications by AJ Miller.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Healy,
!    Algorithm AS 6:
!    Triangular decomposition of a symmetric matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 195-197.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A((N*(N+1))/2), a positive definite matrix
!    stored by rows in lower triangular form as a one dimensional array,
!    in the sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, integer ( kind = 4 ) NN, the dimension of A, (N*(N+1))/2.
!
!    Output, real ( kind = 8 ) U((N*(N+1))/2), an upper triangular matrix,
!    stored by columns, which is the Cholesky factor of A.  The program is
!    written in such a way that A and U can share storage.
!
!    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.
!    If NULLTY is zero, the matrix is judged to have full rank.
!
!    Output, integer ( kind = 4 ) IFAULT, an error indicator.
!    0, no error was detected;
!    1, if N < 1;
!    2, if A is not positive semi-definite.
!    3, if NN < (N*(N+1))/2.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) ETA, should be set equal to the smallest positive
!    value such that 1.0 + ETA is calculated as being greater than 1.0 in the
!    accuracy being used.
!
  implicit none

  integer ( kind = 4 ) nn

  real ( kind = 8 ) a(nn)
  real ( kind = 8 ), parameter :: eta = 1.0D-09
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nullty
  real ( kind = 8 ) rsq
  real ( kind = 8 ) u(nn)
  real ( kind = 8 ) w
  real ( kind = 8 ) x

  ifault = 0
  nullty = 0

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  if ( nn < ( n * ( n + 1 ) ) / 2 ) then
    ifault = 3
    return
  end if

  j = 1
  k = 0
  ii = 0
!
!  Factorize column by column, ICOL = column number.
!
  do icol = 1, n

    ii = ii + icol
    x = eta * eta * a(ii)
    l = 0
    kk = 0
!
!  IROW = row number within column ICOL.
!
    do irow = 1, icol

      kk = kk + irow
      k = k + 1
      w = a(k)
      m = j

      do i = 1, irow - 1
        l = l + 1
        w = w - u(l) * u(m)
        m = m + 1
      end do

      l = l + 1

      if ( irow == icol ) then
        exit
      end if

      if ( u(l) /= 0.0D+00 ) then

        u(k) = w / u(l)

      else

        u(k) = 0.0D+00

        if ( abs ( x * a(k) ) < w * w ) then
          ifault = 2
          return
        end if

      end if

    end do
!
!  End of row, estimate relative accuracy of diagonal element.
!
    if ( abs ( w ) <= abs ( eta * a(k) ) ) then

      u(k) = 0.0D+00
      nullty = nullty + 1

    else

      if ( w < 0.0D+00 ) then
        ifault = 2
        return
      end if

      u(k) = sqrt ( w )

    end if

    j = j + icol

  end do

  return
    end


