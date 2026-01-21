module agn_stars

   use star_lib
   use star_def
   use const_def
   use math_lib
   use utils_lib

   implicit none
   public

   ! Radiation-dominated AGN
   real(dp), parameter :: gamma_agn = 4d0 / 3d0
   ! Indices of angular momentum variables in extras_rpar
   integer, parameter :: i_gain = 5
   integer, parameter :: i_loss = 6
   integer, parameter :: i_mu_J = 18
   integer, parameter :: i_var_J = 19
   integer, parameter :: i_J_crit = 20
   integer, parameter :: i_suppression = 8   ! Torque suppression factor for diagnostics
   integer, parameter :: i_v_vcrit_ratio = 9   ! Torque suppression factor for diagnostics 


   ! Index definitions for extras_rpar (add to your existing indices)
   integer, parameter :: i_dj_dt = 1       ! net torque [erg]
   integer, parameter :: i_H_p = 2         ! pressure scale height [cm]
   integer, parameter :: i_k_acc = 3       ! number of zones in deposition region
   integer, parameter :: i_J_dep = 4       ! angular momentum deposited this step [g cm2/s]
   integer, parameter :: i_R_Hill = 15     ! Hill radius [cm]
   integer, parameter :: i_R_Bondi = 16    ! Bondi radius [cm]

   ! Index for extras_lpar
   integer, parameter :: i_verbose_torque = 1  ! toggle verbose output


   ! For time-step control
   real(dp), parameter :: max_delta_ln_m = 3d-3
   real(dp), parameter :: blend_safety = 3d-1
   real(dp), parameter :: burn_safety = 1d2

   real(dp) :: const_rho               ! Ambient density in g/cm^3
   real(dp) :: const_Tb                ! Ambient temperature in K
   real(dp) :: const_csb               ! Ambient sound speed in cm/s, corresponds to const_Tb
   real(dp) :: const_kappa             ! Average opacity of accretion stream (cm^2/s) 
   real(dp) :: atm_blend_time          ! Time (yrs) over which star falls into AGN (used for blending atmospheric boundary condition)
   real(dp) :: mdot_blend_time         ! Time (yrs) over which star falls into AGN (used for blending atmospheric boundary condition)
   real(dp) :: min_Teff                ! Extra term added to Teff^4, phased out during the min_Teff_blend_time.
   real(dp) :: min_Teff_blend_time     ! Time over which to phase out min_Teff
   real(dp) :: Omega_AGN               ! AGN angular frequency
   real(dp) :: cutoff_mass             ! Maximum mass set by gap opening
   real(dp) :: Mass_AGN                ! SMBH mass
   real(dp) :: last_gain               ! mass gain in last time step
   logical :: use_mass_cutoff          ! Use the cutoff mass?
   real(dp) :: atm_blend = 1d-20         ! Current state of the blend
   real(dp) :: mdot_blend = 1d-20       ! Current state of the blend
   real(dp) :: lost_mass(200)          ! Tracks lost mass



   contains


   subroutine my_other_D_mix(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      integer :: j, k, nz
      real(dp) :: Ledd, Leff, R, Dmix, beta, beta_max, beta_trans

      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      nz = s%nz

      if (s%x_logical_ctrl(8)) then
         beta_max = 0d0
         beta_trans = s%x_ctrl(12)
         do k=1,nz
            beta = s%Prad(k) / s%Peos(k)
            Dmix = s%scale_height(k) * pow(max(Lsun,s%L(k)) / (4 * pi * pow2(s%r(1)) * s%rho(k)), 1d0 / 3d0)
            s%D_mix(k) = s%D_mix(k) + Dmix * exp(-min(30d0, beta_trans/beta))
            beta_max = max(beta, beta_max)
         end do

         ! Save variables for history
         s% xtra1_array(3) = beta_max
      else
         Leff = s% L_start(1)
         R = exp(s%lnR_start(1))
         Ledd = eval_Ledd(s, ierr)
         if (s% x_logical_ctrl(7)) then
            Ledd = Ledd * (1d0 - pow2(s%job%extras_rpar(i_mu_J) / s%job%extras_rpar(i_J_crit)))
            Ledd = max(Ledd, Lsun)
         end if
         do k=1,nz
            Dmix = s%scale_height(k) * pow(Leff / (4 * pi * pow2(R) * s%rho(k)), 1d0 / 3d0)
            s%D_mix(k) = s%D_mix(k) + Dmix * tanh(pow(Leff / Ledd, s% x_ctrl(5)))  ! x_ctrl(5) = Power of gamma for increasing Mixing
         end do
         s% xtra1_array(3) = tanh(pow(Leff / Ledd, s% x_ctrl(5))) ! Extra mixing factor (constant through the star) induced by accretion
      end if


   end subroutine my_other_D_mix

   integer function my_other_timestep_limit( &
      id, skip_hard_limit, dt, dt_limit_ratio)
      use const_def, only: dp
      integer, intent(in) :: id
      logical, intent(in) :: skip_hard_limit
      real(dp), intent(in) :: dt
      real(dp), intent(inout) :: dt_limit_ratio
      type (star_info), pointer :: s
      integer :: j
      real(dp) :: eps_nuc, burn_time, mix_time, f_mix, min_burn_time, min_eps_nuc, min_f_mix, diff_integral, mdot
      integer :: ierr

     
      real(dp) :: v_vcrit, rotation_dt_limit
      

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return


      ! Limit timestep as we approach critical rotation
      v_vcrit = s% v_div_v_crit_avg_surf
      
      if (v_vcrit > 0.7d0) then
         ! As we approach critical, force smaller timesteps
         ! This ensures we don't overshoot critical in a single step
         rotation_dt_limit = (1d0 - v_vcrit) * 0.1d0  ! Fraction of remaining "room" to critical
         rotation_dt_limit = max(rotation_dt_limit, 1d-4)  ! Floor to avoid tiny timesteps
         dt_limit_ratio = max(dt_limit_ratio, s%dt / (rotation_dt_limit * s% kh_timescale * s%time_delta_coeff))
         write(*,*) 'Rotation dt limit (v/v_crit):', v_vcrit, ' dt_limit_ratio:', &
                  s%dt / (rotation_dt_limit * s% kh_timescale * s%time_delta_coeff)
      end if


      min_burn_time = 1d100
      diff_integral = 0d0
      do j=2,s%nz
         if (s%m(1) - s%m(j) > last_gain) then
            ! Tracks the diffusion front, assuming a constant source at the surface
            diff_integral = diff_integral + (s%r(1)-s%r(j)) * s%dm(j) / (4d0 * pi * pow2(s%r(j)) * s%rho(j) * s%D_mix(j)) 
            mix_time = diff_integral

            if (mix_time / s%dt < 20d0) then
               f_mix = min(1d0, exp(-mix_time / s%dt))
            else
               f_mix = 1d-50
            end if

            eps_nuc = s%eps_nuc(j)
            burn_time = pow2(s%csound(j)) / (f_mix * eps_nuc) ! (thermal energy / burning energy generation)

            if (burn_time < min_burn_time) then
               min_burn_time = burn_time
               min_eps_nuc = eps_nuc
               min_f_mix = f_mix
            end if
         end if
      end do

      dt_limit_ratio = s%dt / (min_burn_time * burn_safety * s%time_delta_coeff)
      write(*,*) 'Burn dt limit:', s%dt / (min_burn_time * burn_safety * s%time_delta_coeff)

      mdot = max(abs(s%job%extras_rpar(i_gain)), abs(s%job%extras_rpar(i_loss)))
      dt_limit_ratio = max(dt_limit_ratio, mdot * s%dt / (max_delta_ln_m * s%m(1) * s%time_delta_coeff))
      write(*,*) 'Mdot dt limit:', mdot * s%dt / (max_delta_ln_m * s%m(1) * s%time_delta_coeff)

      if (atm_blend < 0.999d0 .and. mdot_blend < 1d-5) then
         dt_limit_ratio = max(dt_limit_ratio, s%dt / (atm_blend_time * secyer * blend_safety * s%time_delta_coeff))
         write(*,*) 'atm blend limit:', s%dt / (atm_blend_time * secyer * blend_safety * s%time_delta_coeff)
      else if (mdot_blend > 1d-5 .and. mdot_blend < 0.99d0) then
         dt_limit_ratio = max(dt_limit_ratio, s%dt / (mdot_blend_time * secyer * blend_safety * s%time_delta_coeff))
         write(*,*) 'mdot blend limit:', s%dt / (mdot_blend_time * secyer * blend_safety * s%time_delta_coeff)
      end if
      write(*,*) 'overall limit', dt_limit_ratio
      my_other_timestep_limit = keep_going
   end function my_other_timestep_limit

   subroutine set_blend(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      real(dp) :: initial_age, mdot_start_time

      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)

      ! 1d5 is because the starting model is 1d5 years old.
      initial_age = 1d5

      atm_blend = 0.5d0 * (1d0 + tanh((s%time / secyer - initial_age - atm_blend_time) / (3d-1 * atm_blend_time))) + 1d-20
      mdot_start_time = 2d0 * atm_blend_time
      mdot_blend = 0.5d0 * (1d0 + tanh((s%time / secyer - initial_age - mdot_start_time) / (3d-1 * mdot_blend_time))) + 1d-20

      write(*,*) 'atm blend, mdot blend:', atm_blend, mdot_blend

   end subroutine set_blend

   real(dp) function get_dtau1(s, ierr)
      type (star_info), pointer :: s
      integer, intent(out) :: ierr
      integer :: k
      real(dp) :: kap, kap_min
      include 'formats'
      ierr = 0
      k = 1
      kap = s% opacity(1)
      get_dtau1 = s% dm(1)*kap/(4*pi*s% rmid(1)*s% rmid(1))
      if (is_bad(get_dtau1)) then
         ierr = -1
         if (.not. s% report_ierr) return
         k = 1
         write(*,2) 'get_dtau1', k, get_dtau1
         write(*,2) 's% dm(1)', k, s% dm(k)
         write(*,2) 's% opacity(1)', k, s% opacity(k)
         write(*,2) 's% rmid(1)', k, s% rmid(k)
         write(*,2) 's% r(1)', k, s% r(k)
         write(*,2) 's% r(2)', 2, s% r(2)
         stop 'get_dtau1'
      end if
   end function get_dtau1

   real(dp) function eval_Ledd(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         !real(dp) :: dqsum, Ledd_sum,  weight
         integer :: k
         !logical, parameter :: dbg = .false.
         !include 'formats'
         ierr = 0
         if (s% cgrav(1) <= 0d0) return
         if (ierr /= 0) return
         eval_Ledd = 3.2d4 * s%m(1) * Lsun / Msun
         if (s% x_logical_ctrl(9)) then
           eval_Ledd = eval_Ledd *0.4d0/(0.2d0*(1.0d0 + s%accretion_species_xa(1)))
         end if

    end function eval_Ledd

end module agn_stars
