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



   subroutine compute_surface_dj_old(id, gain, loss, R_bondi, dj_dt)
   real(dp), intent(in) :: gain, loss, R_bondi
   real(dp), intent(out) :: dj_dt
   real(dp) :: M, R, h, R_Hill, R_acc, j_acc, j_surf, Omega_surf
   integer, intent(in) :: id
   integer :: ierr
   type(star_info), pointer :: s

   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return

   M = s%m(1)
   R = s%r(1)
   Omega_surf = s%omega(1)

   ! Disk scale height
   h = sqrt(2d0) * const_csb / Omega_AGN

   ! Hill radius: R_H = [R_bondi * h^2 / 12]^{1/3}
   R_Hill = pow((1d0/12d0) * pow2(h) * R_bondi, 1d0/3d0)
   R_acc = min(R_Hill, R_Bondi)

   ! Specific angular momentum of accreted material
   ! Limited by either disk Keplerian value or breakup at stellar surface
   j_acc = min(pow2(R_acc) * Omega_AGN, sqrt(standard_cgrav * M * R))
   ! Specific angular momentum carried away by mass loss (make sure this is not double-counted by MESA's wind routine)
   j_surf = Omega_surf * pow2(R)

   ! Net dJ/dt to apply at surface via other_torque
   dj_dt = abs(gain) * j_acc - abs(loss) * j_surf

   s% job% extras_rpar(i_R_Hill) = R_Hill

end subroutine compute_surface_dj_old



subroutine compute_surface_dj(id, gain, loss, R_bondi, dj_dt)
   real(dp), intent(in) :: gain, loss, R_bondi
   real(dp), intent(out) :: dj_dt
   real(dp) :: M, R, h, R_Hill, R_acc, j_acc, j_surf, Omega_surf
   real(dp) :: v_vcrit_ratio, suppression_factor
   real(dp) :: onset_v_vcrit, suppression_power
   integer, intent(in) :: id
   integer :: ierr
   type(star_info), pointer :: s

   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return

   M = s%m(1)
   R = s%r(1)
   Omega_surf = s%omega(1)

   ! Use MESA's mass-weighted surface v/v_crit average
   v_vcrit_ratio = s% v_div_v_crit_avg_surf

   ! Disk scale height
   h = sqrt(2d0) * const_csb / Omega_AGN

   ! Hill radius and accretion radius
   R_Hill = pow((1d0/12d0) * pow2(h) * R_bondi, 1d0/3d0)
   R_acc = min(R_Hill, R_Bondi)

   ! Specific angular momenta
   j_acc = min(pow2(R_acc) * Omega_AGN, sqrt(standard_cgrav * M * R))
   j_surf = Omega_surf * pow2(R)

   ! Raw torque [g·cm²/s²]
   dj_dt = abs(gain) * j_acc  - abs(loss) * j_surf
   write(*,*) 'Raw dj_dt=', dj_dt, ' gain=', gain, ' loss=', loss, ' j_acc=', j_acc, ' j_surf=', j_surf

   ! Suppress positive torque as we approach critical rotation
   ! Parameters to tune:
   !   onset_v_vcrit: v/v_crit at which suppression begins
   !   suppression_power: higher = sharper cutoff near critical
   onset_v_vcrit = 0.9d0      ! Start suppression at 90% critical
   suppression_power = 4d0    ! Steep falloff

   if (dj_dt > 0d0) then
      if (v_vcrit_ratio >= 1d0) then
         ! At or above critical: no spin-up allowed
         suppression_factor = 0d0
      else if (v_vcrit_ratio > onset_v_vcrit) then
         ! Between onset and critical: smooth but steep suppression
         ! Goes from 1 at onset_v_vcrit to 0 at v_vcrit_ratio = 1
         suppression_factor = pow((1d0 - v_vcrit_ratio) / (1d0 - onset_v_vcrit), suppression_power)
      else
         ! Below onset: no suppression
         suppression_factor = 1d0
      end if
      dj_dt = dj_dt * suppression_factor
   end if

   s% job% extras_rpar(i_R_Hill) = R_Hill
   s% job% extras_rpar(i_v_vcrit_ratio) = v_vcrit_ratio
   s% job% extras_rpar(i_suppression) = suppression_factor

end subroutine compute_surface_dj

subroutine agn_other_torque(id, ierr)
   use star_def
   integer, intent(in) :: id
   integer, intent(out) :: ierr
   type (star_info), pointer :: s
   integer :: k, k_acc
   real(dp) :: gain, loss, R_bondi, dj_dt
   real(dp) :: m_sum, J_surf, scale
   real(dp) :: H_p, r_deposit
   real(dp) :: j_acc, j_surf_specific, dJ_total
   real(dp) :: m_deposit_target
   real(dp) :: torque_check


   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return

   m_deposit_target = 1.0d0 * s% mstar  ! Deposit angular momentum in 100% of stellar mass (tried initially with 1% but led to numerical issues)

   s% extra_jdot(:) = 0d0
   s% extra_omegadot(:) = 0d0

   ! Retrieve gain, loss, and R_bondi from storage
   gain = s% job% extras_rpar(i_gain)
   loss = s% job% extras_rpar(i_loss)
   R_bondi = s% job% extras_rpar(i_R_bondi)

   ! Skip during relaxation
   if (s% doing_relax) return

   call compute_surface_dj(id, gain, loss, R_bondi, dj_dt)

   !write(*,*) 'agn_other_torque: gain=', gain,' loss=',loss,'M_deposit=',m_deposit_target,' dj_dt=',dj_dt  

   ! Store for diagnostics
   s% job% extras_rpar(i_dj_dt) = dj_dt

   ! Nothing to do if no net torque
   if (abs(dj_dt) < 1d-99) then
      s% job% extras_rpar(i_H_p) = 0d0
      s% job% extras_rpar(i_k_acc) = 0d0
      s% job% extras_rpar(i_J_dep) = 0d0
      return
   end if

 ! Determine k_acc and compute total mass in deposition region
   m_sum = 0d0
   k_acc = 1
   do k = 1, s% nz
      m_sum = m_sum + s% dm_bar(k)
      k_acc = k
      if (m_sum >= m_deposit_target) exit
   end do

   ! Recompute m_sum for the actual deposition region
   m_sum = 0d0
   do k = 1, k_acc
      m_sum = m_sum + s% dm_bar(k)
   end do

   ! Deposit torque UNIFORMLY PER UNIT MASS (not per zone!)
   ! This means extra_jdot is the SAME for all zones
   ! extra_jdot [cm²/s²] = dj_dt [g cm²/s²] / m_sum [g]
   do k = 1, k_acc
      s% extra_jdot(k) = dj_dt / m_sum
   end do

   ! Compute diagnostic quantities
   dJ_total = dj_dt * s% dt

   ! Check angular momentum conservation
   torque_check = dot_product(s% extra_jdot(1:k_acc), s% dm_bar(1:k_acc))

   ! Store diagnostics
   s% job% extras_rpar(i_k_acc) = dble(k_acc)
   s% job% extras_rpar(i_J_dep) = dJ_total

   ! Diagnostics output
   j_acc = min(pow2(min(R_bondi, s% job% extras_rpar(i_R_Hill))) * Omega_AGN, &
               sqrt(standard_cgrav * s% m(1) * s% r(1)))
   j_surf_specific = s% omega(1) * pow2(s% r(1))

   write(*,'(A)') '  ===== AGN Torque Diagnostics ====='
   write(*,'(A,ES12.4)') '    gain [g/s]        = ', gain
   write(*,'(A,ES12.4)') '    loss [g/s]        = ', loss
   write(*,'(A,ES12.4)') '    j_acc [cm2/s]     = ', j_acc
   write(*,'(A,ES12.4)') '    j_surf [cm2/s]    = ', j_surf_specific
   write(*,'(A,F8.4)')   '    v/v_crit (surf)   = ', s% job% extras_rpar(i_v_vcrit_ratio)
   write(*,'(A,ES12.4)') '    dJ/dt [g cm2/s2]  = ', dj_dt
   write(*,'(A,ES12.4)') '    dJ this step      = ', dJ_total
   write(*,'(A,I8)')     '    k_acc             = ', k_acc
   write(*,'(A,ES12.4)') '    m_deposit [g]     = ', m_sum
   write(*,'(A,ES12.4)') '    m_deposit/M_star  = ', m_sum / s% mstar
   write(*,'(A,ES12.4)') '    J_star [g cm2/s]  = ', s% total_angular_momentum
   write(*,'(A,ES12.4)') '    dJ / J_star       = ', dJ_total / s% total_angular_momentum
   write(*,'(A,ES12.4)') '    torque_check/dj_dt= ', torque_check / dj_dt
   write(*,'(A)') '  ==================================='


end subroutine agn_other_torque



end module rotation