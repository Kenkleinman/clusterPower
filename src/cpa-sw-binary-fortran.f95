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


      subroutine cdgqf ( nt, kind, alpha, beta, t, wts )

!*********************************************************************72
!
!! CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
!
!  Discussion:
!
!    This routine computes all the knots and weights of a Gauss quadrature
!    formula with a classical weight function with default values for A and B,
!    and only simple knots.
!
!    There are no moments checks and no printing is done.
!
!    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    This FORTRAN77 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer NT, the number of knots.
!
!    Input, integer KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, double precision ALPHA, the value of Alpha, if needed.
!
!    Input, double precision BETA, the value of Beta, if needed.
!
!    Output, double precision T(NT), the knots.
!
!    Output, double precision WTS(NT), the weights.
!
      implicit none

      integer nt

      double precision aj(nt)
      double precision alpha
      double precision beta
      double precision bj(nt)
      integer kind
      double precision t(nt)
      double precision wts(nt)
      double precision zemu

      call parchk ( kind, 2 * nt, alpha, beta )
!
!  Get the Jacobi matrix and zero-th moment.
!
      call class_matrix ( kind, nt, alpha, beta, aj, bj, zemu )
!
!  Compute the knots and weights.
!
      call sgqf ( nt, aj, bj, zemu, t, wts )

      return
      end
      subroutine cgqf ( nt, kind, alpha, beta, a, b, t, wts )

!*********************************************************************72
!
!! CGQF computes knots and weights of a Gauss quadrature formula.
!
!  Discussion:
!
!    The user may specify the interval (A,B).
!
!    Only simple knots are produced.
!
!    Use routine EIQFS to evaluate this quadrature formula.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    This FORTRAN77 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer NT, the number of knots.
!
!    Input, integer KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, double precision ALPHA, the value of Alpha, if needed.
!
!    Input, double precision BETA, the value of Beta, if needed.
!
!    Input, double precision A, B, the interval endpoints, or
!    other parameters.
!
!    Output, double precision T(NT), the knots.
!
!    Output, double precision WTS(NT), the weights.
!
      implicit none

      integer nt

      double precision a
      double precision alpha
      double precision b
      double precision beta
      integer i
      integer kind
      integer mlt(nt)
      integer ndx(nt)
      double precision t(nt)
      double precision wts(nt)
!
!  Compute the Gauss quadrature formula for default values of A and B.
!
      call cdgqf ( nt, kind, alpha, beta, t, wts )
!
!  Prepare to scale the quadrature formula to other weight function with
!  valid A and B.
!
      do i = 1, nt
        mlt(i) = 1
      end do

      do i = 1, nt
        ndx(i) = i
      end do

      call scqf ( nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, &
        a, b )

      return
      end
      subroutine ch_cap ( ch )

!*********************************************************************72
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character CH, the character to capitalize.
!
      implicit none

      character ch
      integer itemp

      itemp = ichar ( ch )

      if ( 97 .le. itemp .and. itemp .le. 122 ) then
        ch = char ( itemp - 32 )
      end if

      return
      end
      function ch_eqi ( c1, c2 )

!*********************************************************************72
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
      implicit none

      character c1
      character c1_cap
      character c2
      character c2_cap
      logical ch_eqi

      c1_cap = c1
      c2_cap = c2

      call ch_cap ( c1_cap )
      call ch_cap ( c2_cap )

      if ( c1_cap .eq. c2_cap ) then
        ch_eqi = .true.
      else
        ch_eqi = .false.
      end if

      return
      end
      subroutine ch_to_digit ( c, digit )

!*********************************************************************72
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     !   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
      implicit none

      character c
      integer digit

      if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

        digit = ichar ( c ) - 48

      else if ( c .eq. ' ' ) then

        digit = 0

      else

        digit = -1

      end if

      return
      end
      subroutine class_matrix ( kind, m, alpha, beta, aj, bj, zemu )

!*********************************************************************72
!
!! CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
!
!  Discussion:
!
!    This routine computes the diagonal AJ and sub-diagonal BJ
!    elements of the order M tridiagonal symmetric Jacobi matrix
!    associated with the polynomials orthogonal with respect to
!    the weight function specified by KIND.
!
!    For weight functions 1-7, M elements are defined in BJ even
!    though only M-1 are needed.  For weight function 8, BJ(M) is
!    set to zero.
!
!    The zero-th moment of the weight function is returned in ZEMU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    This FORTRAN77 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, integer M, the order of the Jacobi matrix.
!
!    Input, double precision ALPHA, the value of Alpha, if needed.
!
!    Input, double precision BETA, the value of Beta, if needed.
!
!    Output, double precision AJ(M), BJ(M), the diagonal and subdiagonal
!    of the Jacobi matrix.
!
!    Output, double precision ZEMU, the zero-th moment.
!
      implicit none

      integer m

      double precision a2b2
      double precision ab
      double precision aba
      double precision abi
      double precision abj
      double precision abti
      double precision aj(m)
      double precision alpha
      double precision apone
      double precision beta
      double precision bj(m)
      integer i
      integer kind
      double precision pi
      parameter ( pi = 3.14159265358979323846264338327950D+00 )
      double precision r8_gamma
      double precision r8_epsilon
      double precision temp
      double precision temp2
      double precision zemu

      temp = r8_epsilon ( )

      call parchk ( kind, 2 * m - 1, alpha, beta )

      temp2 = r8_gamma ( 0.5D+00 )

      if ( 500.0D+00 * temp .lt. abs ( temp2 * temp2 - pi ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CLASS_MATRIX - Fatal error!'
        write ( *, '(a)' ) &
         '  Gamma function does not match machine parameters.'
        stop 1
      end if

      if ( kind .eq. 1 ) then

        ab = 0.0D+00

        zemu = 2.0D+00 / ( ab + 1.0D+00 )

        do i = 1, m
          aj(i) = 0.0D+00
        end do

        do i = 1, m
          abi = i + ab * mod ( i, 2 )
          abj = 2 * i + ab
          bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
        end do

        do i = 1, m
          bj(i) = sqrt ( bj(i) )
        end do

      else if ( kind .eq. 2 ) then

        zemu = pi

        do i = 1, m
          aj(i) = 0.0D+00
        end do

        bj(1) =  sqrt ( 0.5D+00 )
        do i = 2, m
          bj(i) = 0.5D+00
        end do

      else if ( kind .eq. 3 ) then

        ab = alpha * 2.0D+00
        zemu = 2.0D+00**( ab + 1.0D+00 ) &
         * ( r8_gamma ( alpha + 1.0D+00 ) )**2 &
         / r8_gamma ( ab + 2.0D+00 )

        do i = 1, m
          aj(i) = 0.0D+00
        end do

        bj(1) = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )
        do i = 2, m
          bj(i) = i * ( i + ab ) &
           / ( 4.0D+00 * ( i + alpha ) * ( i + alpha ) - 1.0D+00 )
        end do
        do i = 1, m
          bj(i) =  sqrt ( bj(i) )
        end do

      else if ( kind .eq. 4 ) then

        ab = alpha + beta
        abi = 2.0D+00 + ab
        zemu = 2.0D+00**( ab + 1.0D+00 ) * r8_gamma ( alpha + 1.0D+00 ) &
           * r8_gamma ( beta + 1.0D+00 ) / r8_gamma ( abi )
        aj(1) = ( beta - alpha ) / abi
        bj(1) = 4.0D+00 * ( 1.0 + alpha ) * ( 1.0D+00 + beta ) &
         / ( ( abi + 1.0D+00 ) * abi * abi )
        a2b2 = beta * beta - alpha * alpha

        do i = 2, m
          abi = 2.0D+00 * i + ab
          aj(i) = a2b2 / ( ( abi - 2.0D+00 ) * abi )
          abi = abi * abi
          bj(i) = 4.0D+00 * i * ( i + alpha ) * ( i + beta ) &
           * ( i + ab ) / ( ( abi - 1.0D+00 ) * abi )
        end do

        do i = 1, m
          bj(i) =  sqrt ( bj(i) )
        end do

      else if ( kind .eq. 5 ) then

        zemu = r8_gamma ( alpha + 1.0D+00 )

        do i = 1, m
          aj(i) = 2.0D+00 * i - 1.0D+00 + alpha
          bj(i) = i * ( i + alpha )
        end do

        do i = 1, m
          bj(i) = sqrt ( bj(i) )
        end do

      else if ( kind .eq. 6 ) then

        zemu = r8_gamma ( ( alpha + 1.0D+00 ) / 2.0D+00 )

        do i = 1, m
          aj(i) = 0.0D+00
        end do

        do i = 1, m
          bj(i) = ( i + alpha * mod ( i, 2 ) ) / 2.0D+00
        end do

        do i = 1, m
          bj(i) =  sqrt ( bj(i) )
        end do

      else if ( kind .eq. 7 ) then

        ab = alpha
        zemu = 2.0D+00 / ( ab + 1.0D+00 )

        do i = 1, m
          aj(i) = 0.0D+00
        end do

        do i = 1, m
          abi = i + ab * mod ( i, 2 )
          abj = 2 * i + ab
          bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
        end do

        do i = 1, m
          bj(i) =  sqrt ( bj(i) )
        end do

      else if ( kind .eq. 8 ) then

        ab = alpha + beta
        zemu = r8_gamma ( alpha + 1.0D+00 ) &
         * r8_gamma ( - ( ab + 1.0D+00 ) ) &
         / r8_gamma ( - beta )
        apone = alpha + 1.0D+00
        aba = ab * apone
        aj(1) = - apone / ( ab + 2.0D+00 )
        bj(1) = - aj(1) * ( beta + 1.0D+00 ) / ( ab + 2.0D+00 ) &
         / ( ab + 3.0D+00 )
        do i = 2, m
          abti = ab + 2.0D+00 * i
          aj(i) = aba + 2.0D+00 * ( ab + i ) * ( i - 1 )
          aj(i) = - aj(i) / abti / ( abti - 2.0D+00 )
        end do

        do i = 2, m - 1
          abti = ab + 2.0D+00 * i
          bj(i) = i * ( alpha + i ) / ( abti - 1.0D+00 ) * ( beta + i ) &
           / ( abti * abti ) * ( ab + i ) / ( abti + 1.0D+00 )
        end do

        bj(m) = 0.0D+00

        do i = 1, m
          bj(i) =  sqrt ( bj(i) )
        end do

      else if ( kind .eq. 9 ) then

        zemu = pi / 2.0D+00

        do i = 1, m
          aj(i) = 0.0D+00
          bj(i) = 0.5D+00
        end do

      end if

      return
      end
      subroutine get_unit ( iunit )

!*********************************************************************72
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
      implicit none

      integer i
      integer iunit
      logical value

      iunit = 0

      do i = 1, 99

        if ( i .ne. 5 .and. i .ne. 6 .and. i .ne. 9 ) then

          inquire ( unit = i, opened = value, err = 10 )

          if ( .not. value ) then
            iunit = i
            return
          end if

        end if

10      continue

      end do

      return
      end

      subroutine imtqlx ( n, d, e, z )

!*********************************************************************72
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine.
!
!    It has been modified to produce the product Q' * Z, where Z is an input
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.
!    The changes consist (essentially) of applying the orthogonal
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    This FORTRAN77 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, double precision D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, double precision E(N), the subdiagonal entries of the
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, double precision Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
      implicit none

      integer n

      double precision b
      double precision c
      double precision d(n)
      double precision e(n)
      double precision f
      double precision g
      integer i
      integer ii
      integer itn
      parameter ( itn = 30 )
      integer j
      integer k
      integer l
      integer m
      integer mml
      double precision p
      double precision prec
      double precision r
      double precision r8_epsilon
      double precision s
      double precision z(n)

      prec = r8_epsilon ( )

      if ( n .eq. 1 ) then
        return
      end if

      e(n) = 0.0D+00

      do l = 1, n

        j = 0

10      continue

          do m = l, n

            if ( m .eq. n ) then
              go to 20
            end if

            if ( abs ( e(m) ) .le. &
             prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
              go to 20
            end if

          end do

20        continue

          p = d(l)

          if ( m .eq. l ) then
            go to 30
          end if

          if ( itn .le. j ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'IMTQLX - Fatal error!'
            write ( *, '(a)' ) '  Iteration limit exceeded.'
            write ( *, '(a,i8)' ) '  J = ', j
            write ( *, '(a,i8)' ) '  L = ', l
            write ( *, '(a,i8)' ) '  M = ', m
            write ( *, '(a,i8)' ) '  N = ', n
            stop 1
          end if

          j = j + 1
          g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
          r =  sqrt ( g * g + 1.0D+00 )
          g = d(m) - p + e(l) / ( g + sign ( r, g ) )
          s = 1.0D+00
          c = 1.0D+00
          p = 0.0D+00
          mml = m - l

          do ii = 1, mml

            i = m - ii
            f = s * e(i)
            b = c * e(i)

            if ( abs ( g ) .le. abs ( f ) ) then
              c = g / f
              r =  sqrt ( c * c + 1.0D+00 )
              e(i+1) = f * r
              s = 1.0D+00 / r
              c = c * s
            else
              s = f / g
              r =  sqrt ( s * s + 1.0D+00 )
              e(i+1) = g * r
              c = 1.0D+00 / r
              s = s * c
            end if

            g = d(i+1) - p
            r = ( d(i) - g ) * s + 2.0D+00 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
            f = z(i+1)
            z(i+1) = s * z(i) + c * f
            z(i) = c * z(i) - s * f

          end do

          d(l) = d(l) - p
          e(l) = g
          e(m) = 0.0D+00

        go to 10

30      continue

      end do
!
!  Sorting.
!
      do ii = 2, n

        i = ii - 1
        k = i
        p = d(i)

        do j = ii, n
          if ( d(j) .lt. p ) then
            k = j
            p = d(j)
          end if
        end do

        if ( k .ne. i ) then
          d(k) = d(i)
          d(i) = p
          p = z(i)
          z(i) = z(k)
          z(k) = p
        end if

      end do

      return
      end
      subroutine legendre_handle ( order, a, b, x, w )

!*********************************************************************72
!
! LEGENDRE_COMPUTE computes a Legendre quadrature rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [-1,1].
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 September 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ORDER, the order of the rule.
!    1 <= ORDER.
!
!    Output, double precision X(ORDER), the abscissas.
!
!    Output, double precision W(ORDER), the weights.
!
      implicit none

      integer order

      double precision a
      double precision alpha
      double precision b
      double precision beta
      character * ( 255 ) filename
      integer kind
      double precision r(2)
      double precision w(order)
      double precision x(order)

      kind = 1
      alpha = 0.0D+00
      beta = 0.0D+00

      call cgqf ( order, kind, alpha, beta, a, b, x, w )

      return
      end
      subroutine parchk ( kind, m, alpha, beta )

!*********************************************************************72
!
!! PARCHK checks parameters ALPHA and BETA for classical weight functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    This FORTRAN77 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, integer M, the order of the highest moment to
!    be calculated.  This value is only needed when KIND = 8.
!
!    Input, double precision ALPHA, BETA, the parameters, if required
!    by the value of KIND.
!
      implicit none

      double precision alpha
      double precision beta
      integer kind
      integer m
      double precision tmp

      if ( kind .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PARCHK - Fatal error!'
        write ( *, '(a)' ) '  KIND .le. 0.'
        stop 1
      end if
!
!  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
!
      if ( 3 .le. kind .and. kind .le. 8 .and. &
        alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PARCHK - Fatal error!'
        write ( *, '(a)' ) '  3 .le. KIND and ALPHA .le. -1.'
        stop 1
      end if
!
!  Check BETA for Jacobi.
!
      if ( kind .eq. 4 .and. beta .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PARCHK - Fatal error!'
        write ( *, '(a)' ) '  KIND .eq. 4 and BETA .le. -1.0.'
        stop 1
      end if
!
!  Check ALPHA and BETA for rational.
!
      if ( kind .eq. 8 ) then
        tmp = alpha + beta + m + 1.0D+00
        if ( 0.0D+00 .le. tmp .or. tmp .le. beta ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PARCHK - Fatal error!'
          write ( *, '(a)' ) &
           '  KIND .eq. 8 but condition on ALPHA and BETA fails.'
          stop 1
        end if
      end if

      return
      end
      function r8_epsilon ( )

!*********************************************************************72
!
!! R8_EPSILON returns the R8 roundoff unit.
!
!  Discussion:
!
!    The roundoff unit is a number R which is a power of 2 with the
!    property that, to the precision of the computer's arithmetic,
!      1 .lt. 1 + R
!    but
!      1 = ( 1 + R / 2 )
!
!    FORTRAN90 provides the superior library routine
!
!      EPSILON ( X )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision R8_EPSILON, the R8 roundoff unit.
!
      implicit none

      double precision r8_epsilon

      r8_epsilon = 2.220446049250313D-016

      return
      end


      subroutine scqf ( nt, t, mlt, wts, nwts, ndx, swts, st, kind, alpha, beta, a, b )

!*********************************************************************72
!
!! SCQF scales a quadrature formula to a nonstandard interval.
!
!  Discussion:
!
!    The arrays WTS and SWTS may coincide.
!
!    The arrays T and ST may coincide.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    This FORTRAN77 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer NT, the number of knots.
!
!    Input, double precision T(NT), the original knots.
!
!    Input, integer MLT(NT), the multiplicity of the knots.
!
!    Input, double precision WTS(NWTS), the weights.
!
!    Input, integer NWTS, the number of weights.
!
!    Input, integer NDX(NT), used to index the array WTS.
!    For more details see the comments in CAWIQ.
!
!    Output, double precision SWTS(NWTS), the scaled weights.
!
!    Output, double precision ST(NT), the scaled knots.
!
!    Input, integer KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, double precision ALPHA, the value of Alpha, if needed.
!
!    Input, double precision BETA, the value of Beta, if needed.
!
!    Input, double precision A, B, the interval endpoints.
!
      implicit none

      integer nt
      integer nwts

      double precision a
      double precision al
      double precision alpha
      double precision b
      double precision be
      double precision beta
      integer i
      integer k
      integer kind
      integer l
      integer mlt(nt)
      integer ndx(nt)
      double precision p
      double precision r8_epsilon
      double precision shft
      double precision slp
      double precision st(nt)
      double precision swts(nwts)
      double precision t(nt)
      double precision temp
      double precision tmp
      double precision wts(nwts)

      temp = r8_epsilon ( )

      call parchk ( kind, 1, alpha, beta )

      if ( kind .eq. 1 ) then

        al = 0.0D+00
        be = 0.0D+00

        if ( abs ( b - a ) .le. temp ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SCQF - Fatal error!'
          write ( *, '(a)' ) '  |B - A| too small.'
          stop 1
        end if

        shft = ( a + b ) / 2.0D+00
        slp = ( b - a ) / 2.0D+00

      else if ( kind .eq. 2 ) then

        al = -0.5D+00
        be = -0.5D+00

        if ( abs ( b - a ) .le. temp ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SCQF - Fatal error!'
          write ( *, '(a)' ) '  |B - A| too small.'
          stop 1
        end if

        shft = ( a + b ) / 2.0D+00
        slp = ( b - a ) / 2.0D+00

      else if ( kind .eq. 3 ) then

        al = alpha
        be = alpha

        if ( abs ( b - a ) .le. temp ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SCQF - Fatal error!'
          write ( *, '(a)' ) '  |B - A| too small.'
          stop 1
        end if

        shft = ( a + b ) / 2.0D+00
        slp = ( b - a ) / 2.0D+00

      else if ( kind .eq. 4 ) then

        al = alpha
        be = beta

        if ( abs ( b - a ) .le. temp ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SCQF - Fatal error!'
          write ( *, '(a)' ) '  |B - A| too small.'
          stop 1
        end if

        shft = ( a + b ) / 2.0D+00
        slp = ( b - a ) / 2.0D+00

      else if ( kind .eq. 5 ) then

        if ( b .le. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SCQF - Fatal error!'
          write ( *, '(a)' ) '  B .le. 0'
          stop 1
        end if

        shft = a
        slp = 1.0D+00 / b
        al = alpha
        be = 0.0D+00

      else if ( kind .eq. 6 ) then

        if ( b .le. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SCQF - Fatal error!'
          write ( *, '(a)' ) '  B .le. 0.'
          stop 1
        end if

        shft = a
        slp = 1.0D+00 / sqrt ( b )
        al = alpha
        be = 0.0D+00

      else if ( kind .eq. 7 ) then

        al = alpha
        be = 0.0D+00

        if ( abs ( b - a ) .le. temp ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SCQF - Fatal error!'
          write ( *, '(a)' ) '  |B - A| too small.'
          stop 1
        end if

        shft = ( a + b ) / 2.0D+00
        slp = ( b - a ) / 2.0D+00

      else if ( kind .eq. 8 ) then

        if ( a + b .le. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SCQF - Fatal error!'
          write ( *, '(a)' ) '  A + B .le. 0.'
          stop 1
        end if

        shft = a
        slp = a + b
        al = alpha
        be = beta

      else if ( kind .eq. 9 ) then

        al = 0.5D+00
        be = 0.5D+00

        if ( abs ( b - a ) .le. temp ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SCQF - Fatal error!'
          write ( *, '(a)' ) '  |B - A| too small.'
          stop 1
        end if

        shft = ( a + b ) / 2.0D+00
        slp = ( b - a ) / 2.0D+00

      end if

      p = slp ** ( al + be + 1.0D+00 )

      do k = 1, nt

        st(k) = shft + slp * t(k)
        l = abs ( ndx(k) )

        if ( l .ne. 0 ) then
          tmp = p
          do i = l, l + mlt(k) - 1
            swts(i) = wts(i) * tmp
            tmp = tmp * slp
          end do
        end if

      end do

      return
      end
      
      subroutine sgqf ( nt, aj, bj, zemu, t, wts )

!*********************************************************************72
!
!! SGQF computes knots and weights of a Gauss Quadrature formula.
!
!  Discussion:
!
!    This routine computes all the knots and weights of a Gauss quadrature
!    formula with simple knots from the Jacobi matrix and the zero-th
!    moment of the weight function, using the Golub-Welsch technique.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    This FORTRAN77 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer NT, the number of knots.
!
!    Input, double precision AJ(NT), the diagonal of the Jacobi matrix.
!
!    Input/output, double precision BJ(NT), the subdiagonal of the Jacobi
!    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
!
!    Input, double precision ZEMU, the zero-th moment of the weight function.
!
!    Output, double precision T(NT), the knots.
!
!    Output, double precision WTS(NT), the weights.
!
      implicit none

      integer nt

      double precision aj(nt)
      double precision bj(nt)
      integer i
      double precision t(nt)
      double precision wts(nt)
      double precision zemu
!
!  Exit if the zero-th moment is not positive.
!
      if ( zemu .le. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SGQF - Fatal error!'
        write ( *, '(a)' ) '  ZEMU .le. 0.'
        stop 1
      end if
!
!  Set up vectors for IMTQLX.
!
      do i = 1, nt
        t(i) = aj(i)
      end do

      wts(1) = sqrt ( zemu )
      do i = 2, nt
        wts(i) = 0.0D+00
      end do
!
!  Diagonalize the Jacobi matrix.
!
      call imtqlx ( nt, t, bj, wts )

      do i = 1, nt
        wts(i) = wts(i) ** 2
      end do

      return
      end




    function r8_gamma ( x )

!*********************************************************************72
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 .le. X are from reference 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    This FORTRAN77 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
!    Charles Mesztenyi, John Rice, Henry Thatcher, 
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, double precision X, the argument of the function.
!
!    Output, double precision R8_GAMMA, the value of the function.
!
      implicit none

      double precision c(7)
      double precision eps
      double precision fact
      integer i
      integer n
      double precision p(8)
      logical parity
      double precision pi
      double precision q(8)
      double precision r8_gamma
      double precision res
      double precision sqrtpi
      double precision sum
      double precision x
      double precision xbig
      double precision xden
      double precision xinf
      double precision xminin
      double precision xnum
      double precision y
      double precision y1
      double precision ysq
      double precision z
!
!  Mathematical constants
!
      data sqrtpi / 0.9189385332046727417803297D+00 /
      data pi / 3.1415926535897932384626434D+00 /
!
!  Machine dependent parameters
!
      data xbig / 171.624D+00 /
      data xminin / 2.23D-308 /
      data eps / 2.22D-16 /
      data xinf /1.79D+308 /
!
!  Numerator and denominator coefficients for rational minimax
!  approximation over (1,2).
!
      data p / &
      -1.71618513886549492533811d+00, &
       2.47656508055759199108314d+01, &
      -3.79804256470945635097577d+02, &
       6.29331155312818442661052d+02, &
       8.66966202790413211295064d+02, &
      -3.14512729688483675254357d+04, &
      -3.61444134186911729807069d+04, &
       6.64561438202405440627855d+04 /

      data q / &
      -3.08402300119738975254353d+01, &
       3.15350626979604161529144d+02, &
      -1.01515636749021914166146d+03, &
      -3.10777167157231109440444d+03, &
       2.25381184209801510330112d+04, &
       4.75584627752788110767815d+03, &
      -1.34659959864969306392456d+05, &
      -1.15132259675553483497211d+05 /
!
!  Coefficients for minimax approximation over (12, INF).
!
      data c / &
      -1.910444077728D-03, &
       8.4171387781295D-04, &
      -5.952379913043012D-04, &
       7.93650793500350248D-04, &
      -2.777777777777681622553D-03, &
       8.333333333333333331554247D-02, &
       5.7083835261D-03 /

      parity = .false.
      fact = 1.0D+00
      n = 0
      y = x
!
!  Argument is negative.
!
      if ( y .le. 0.0D+00 ) then

        y = - x
        y1 = aint ( y )
        res = y - y1

        if ( res .ne. 0.0D+00 ) then

          if ( y1 .ne. aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
            parity = .true.
          end if

          fact = - pi / sin ( pi * res )
          y = y + 1.0D+00

        else

          res = xinf
          r8_gamma = res
          return

        end if

      end if
!
!  Argument is positive.
!
      if ( y .lt. eps ) then
!
!  Argument < EPS.
!
        if ( xminin .le. y ) then
          res = 1.0D+00 / y
        else
          res = xinf
          r8_gamma = res
          return
        end if

      else if ( y .lt. 12.0D+00 ) then

        y1 = y
!
!  0.0 < argument < 1.0.
!
        if ( y .lt. 1.0D+00 ) then

          z = y
          y = y + 1.0D+00
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
        else

          n = int ( y ) - 1
          y = y - dble ( n )
          z = y - 1.0D+00

        end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
        xnum = 0.0D+00
        xden = 1.0D+00
        do i = 1, 8
          xnum = ( xnum + p(i) ) * z
          xden = xden * z + q(i)
        end do

        res = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0 < argument < 1.0.
!
        if ( y1 .lt. y ) then

          res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
        else if ( y .lt. y1 ) then

          do i = 1, n
            res = res * y
            y = y + 1.0D+00
          end do

        end if

      else
!
!  Evaluate for 12.0 .le. argument.
!
        if ( y .le. xbig ) then

          ysq = y * y
          sum = c(7)
          do i = 1, 6
            sum = sum / ysq + c(i)
          end do
          sum = sum / y - y + sqrtpi
          sum = sum + ( y - 0.5D+00 ) * log ( y )
          res = exp ( sum )

        else

          res = xinf
          r8_gamma = res
          return

        end if

      end if
!
!  Final adjustments and return.
!
      if ( parity ) then
        res = - res
      end if

      if ( fact .ne. 1.0D+00 ) then
        res = fact / res
      end if

      r8_gamma = res

      return
      end
      
      
      
      
      subroutine computeparameter(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
    implicit none
    ! ---- arg types -----------------------
    integer :: JJ
    double precision :: mu, beta, tau2
    double precision :: gamma(JJ)
    double precision :: p11, rho0
    double precision :: p0(JJ)
    ! ---------------------------------------
    ! ---------------------------------------
    integer :: j
    
    mu = p0(1)
    beta = p11 - mu
    tau2 = rho0/(1-rho0)*mu*(1-mu)
    ! gamma
    gamma(1) = 0.0d0
    do j=2,JJ
        gamma(j) = p0(j) - mu
    end do
    end subroutine computeparameter
