module rotation

   use star_lib
   use star_def
   use const_def
   use math_lib
   use agn_stars
   use utils_lib

   implicit none

   public
   contains

   subroutine update_J_dist(id, gain, loss, R_bondi)
      real(dp) :: mu_J, var_J, dt, R_bondi, M, R, gain, loss
      real(dp) :: R_acc, h, R_Hill, j_gain_std, j_gain_avg
      real(dp) :: J_crit
      real(dp) :: I_star, lambda, tau, vt
      integer, intent(in) :: id
      integer :: ierr, i
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      M = s%m(1)
      R = s%r(1)
      dt = s%dt

      I_star = 0d0
      do i=1,s%nz
         I_star = I_star + (2d0 / 3d0) * pow2(s%r(i)) * s%dm(i)
      end do
      lambda = M * pow2(R) / I_star

      h = sqrt(2d0) * const_csb / Omega_AGN

      ! Rh = a * (Mstar / 3 MBH)^{1/3}
      ! = a * [(cs^2 R_bondi / 2 G) (1/3) (G / a^3 Omega^2)]^{1/3}
      ! = [(cs^2 R_bondi / Omega^2) (1/6)]^{1/3}
      ! = [(R_bondi h^2) (1/12)]^{1/3}
      R_Hill = pow((1d0 / 12d0) * pow2(h) * R_bondi, 1d0/3d0)
      R_acc = min(R_Hill, R_Bondi)

      vt = const_csb * min(1d0, pow(R_acc / h, 1d0/3d0))
      tau = min(h, R_acc) / vt

      ! Time step quantities
      ! Want to output enough that we can work through all the rotation stuff in post.
      j_gain_avg = min(pow2(R_acc) * Omega_AGN, pow(standard_cgrav * M * R, 0.5d0))
      j_gain_std = vt * vt * tau

      ! Previous step quantities
      mu_J = s%job%extras_rpar(i_mu_J)
      var_J = s%job%extras_rpar(i_var_J)

      ! Compute next step
      mu_J = mu_J + dt * (gain * j_gain_avg - loss * (mu_J / M) * lambda)
      var_J = var_J + dt * (pow2(gain * j_gain_std) / tau - lambda * loss * (var_J / M))

      ! If the mean J exceeds critical, we limit it to critical.

      ! J_crit = I Omega_crit
      ! Omega_crit = Sqrt[G M / R^3]
      ! So
      ! J_crit =  I Sqrt[G M / R^3]
      J_crit = (1d0 / lambda) * pow(standard_cgrav * pow3(M) * R, 0.5d0)

      if (mu_J > J_crit) then
         ! In practice this is what the truncated Gaussian turns into.
         ! We approach critical so quickly that we overshoot dramatically, leading
         ! to numerical problems with a fancier approach.
         mu_J = J_crit
         var_J = 0d0
      end if

      if (is_bad(mu_J)) then
         write(*,*) 'Bad mu_J, halting.'
         write(*,*) 'mu_J',mu_J
         write(*,*) 'var_J', var_J
         write(*,*) 'j_gain_std', j_gain_std
         write(*,*) 'j_gain_avg', j_gain_avg
         write(*,*) 'gain',gain
         write(*,*) 'loss',loss
         write(*,*) 'tau',tau
         write(*,*) 'lambda',lambda
         stop
      end if
                              ! In this case the variance is truly small because almost all of the distribution was truncated.

      s%job%extras_rpar(i_mu_J) = mu_J
      s%job%extras_rpar(i_var_J) = var_J
      s%job%extras_rpar(i_J_crit) = J_crit

   end subroutine update_J_dist

end module rotation