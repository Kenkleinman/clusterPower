!     Program:         CoxORC.f90
!     Written by:      Xin Zhou
!     Last modified:   Dec 26, 2015
!     Purpose: Cox proportional hazard models with ORC

subroutine der_likelihood_time(mu,beta,gamma,tau2, z0, z1, XX, JJ, KK, a, b, &
                        mincomp, maxcomp, GQ, GQX, GQW, derlikelihood, prob)
    implicit none
    ! ---- arg types -----------------------
    integer :: JJ, KK
    double precision :: mu, beta, tau2
    double precision :: gamma(JJ)
    ! true values of mu, beta, gamma, tau2
    integer :: z0(JJ), z1(JJ)   ! z0: # of subjects with outcome 0; z1: # of subjects with outcome 1.
    integer :: XX(JJ)     ! treatment assignment
    double precision :: a, b     ! integration limits: a - lower, b - upper
    integer :: mincomp(JJ+2), maxcomp(JJ+2)  ! gamma(1), ..., gamma(JJ), mu, beta
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
    double precision :: ff_gamma(JJ-1)
    double precision :: derlikelihood_gamma(JJ-1)
    double precision :: temp
    double precision :: eaa, ebb, exx
    double precision :: faeaa, fbebb

    likelihoodf_denom = 0.0d0
    likelihoodf_numer = 0.0d0
    likelihoodf_denomb2 = 0.0d0
    derlikelihood_mu = 0.0d0
    derlikelihood_beta = 0.0d0
    derlikelihood_tau2 = 0.0d0
    derlikelihood_gamma = 0.0d0

    prob = 0.0d0
    do i=1,GQ
        x = GQX(i)
        exx = dexp(-0.5d0*x*x/tau2)

        ff = 1.0d0
        ffprob = 1.0d0
        ff_mu = 0.0d0
        ff_beta = 0.0d0
        do j=1,JJ
            ff1 = mu+beta*XX(j)+gamma(j)+x
            ff0 = 1-ff1
            ff = ff*(ff0**z0(j))*(ff1**z1(j))

            ! since GQX are not at limits, we ignore the cases where ff0=0 or ff1=0
            temp = z1(j)/ff1 - z0(j)/ff0
            ff_mu = ff_mu + temp
            ff_beta = ff_beta + temp*XX(j)
            k = j-1
            if (k>0) then
                ff_gamma(k) = temp
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
        derlikelihood_gamma = derlikelihood_gamma + GQW(i)*ff*ff_gamma*exx
        derlikelihood_tau2= derlikelihood_tau2 + GQW(i)*ff*x*x*exx
    enddo

    ! calculate f(a)exp(-0.5*a*a)
    eaa = dexp(-0.5d0*a*a/tau2)
    ff = 1.0d0
    do j=1,JJ
        ff1 = mu+beta*XX(j)+gamma(j)+a
        ff0 = 1-ff1
        ff = ff*(ff0**z0(j))*(ff1**z1(j))
    end do
    faeaa = ff*eaa
    ! calculate f(b)exp(-0.5*b*b)
    ebb = dexp(-0.5d0*b*b/tau2)
    ff = 1.0d0
    do j=1,JJ
        ff1 = mu+beta*XX(j)+gamma(j)+b
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
    ! calculate derlikelihood_gamma
    do j=2,JJ
        k = j-1
        derlikelihood_gamma(k) = derlikelihood_gamma(k) + &
                            faeaa*dble(mincomp(j)) - fbebb*dble(maxcomp(j))
        derlikelihood_gamma(k) = derlikelihood_gamma(k)/likelihoodf_numer &
                        -(eaa*dble(mincomp(j))-ebb*dble(maxcomp(j)))/likelihoodf_denom
    end do

    ! calculate derlikelihood_tau2
    derlikelihood_tau2 = 0.5d0*(derlikelihood_tau2/likelihoodf_numer- &
                    likelihoodf_denomb2/likelihoodf_denom)/tau2/tau2
    prob = prob/likelihoodf_denom

    derlikelihood(1) = derlikelihood_mu
    derlikelihood(2) = derlikelihood_beta
    derlikelihood(3:(JJ+1)) = derlikelihood_gamma
    derlikelihood(JJ+2) = derlikelihood_tau2

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
        exx = dexp(-0.5d0*x*x/tau2)

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
        exx1 = dexp(-0.5d0*mu*mu/tau2)       ! exx on lower limits
        exx2 = dexp(-0.5d0*(1-mu-beta)*(1-mu-beta)/tau2)
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
        exx1 = dexp(-0.5d0*(mu+beta)*(mu+beta)/tau2)       ! exx on lower limits
        exx2 = dexp(-0.5d0*(1-mu)*(1-mu)/tau2)
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

