! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use agn_stars
      use rotation
      use atmosphere
      use mass_changes

      implicit none



      ! these routines are called by the standard run_star check_model
      contains

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: const_Tb_rad, const_Tb_gas
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         const_rho = s%x_ctrl(1)
         const_csb = s%x_ctrl(2)
         const_kappa = s%x_ctrl(3)
         atm_blend_time = s%x_ctrl(4)
         mdot_blend_time = s%x_ctrl(6)
         Omega_AGN = s%x_ctrl(9)
         cutoff_mass = s%x_ctrl(10)
         Mass_AGN = s%x_ctrl(11)
         use_mass_cutoff = s%x_logical_ctrl(2)
         last_gain = 0d0

         s% job% extras_rpar(i_mu_J) = 0d0
         s% job% extras_rpar(i_var_J) = 0d0
         s% job% extras_rpar(i_J_crit) = 1d0 ! Tiny starting value to avoid NaN.
         s% job% extras_rpar(i_gain) = 0d0
         s% job% extras_rpar(i_loss) = 0d0
         lost_mass = 0d0

         ! Radiation-dominated
         const_Tb_rad = pow((9d0 * pow2(const_csb) * const_rho / (4 * crad)), 0.25d0)

         ! Gas-dominated, assumed non-ionized
         const_Tb_gas = 3*pow2(const_csb)*amu/(5d0 * kerg)

         const_Tb = min(const_Tb_gas, const_Tb_rad)

         write(*,*) 'If radiation dominated, need T=', const_Tb_rad
         write(*,*) 'If gas dominated, need T=', const_Tb_gas
         write(*,*) 'So T=', const_Tb

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

         ! Uncomment these lines if you wish to use the functions in this file,
         ! otherwise we use a null_ version which does nothing.
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

         s% other_surface_PT => other_surface_PT
         s% other_adjust_mdot => other_adjust_mdot
         s% other_D_mix => my_other_D_mix
         s% other_timestep_limit => my_other_timestep_limit
         s% other_torque => agn_other_torque 

         ! Once you have set the function pointers you want,
         ! then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          s% job% warn_run_star_extras = .false.


        ! V_drag for preventing large acceleration of envelope (prevents shocks/pulsations from forming)
        !feel free to turn on in &contols inlist instead of here.
        s% x_ctrl(13) = 0.95 ! min_q_for_drag
        s% x_ctrl(14) = 1d0 ! drag coefficient

        ! re-implemented velocity drag by EbF
        s% other_momentum_implicit => v_drag_other_momentum_implicit
        s% other_energy_implicit => v_drag_other_energy_implicit

      end subroutine extras_controls


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         

         write(*,*) '>>> extras_start_step: rotation_flag =', s% rotation_flag
         write(*,*) '>>> extras_start_step: use_other_torque =', s% use_other_torque
         write(*,*) '>>> extras_start_step: omega(1) =', s% omega(1)

         call set_blend(id, ierr)

      end subroutine extras_startup


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         real(dp) :: blend, const_Tb_rad, const_Tb_gas
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0

         call set_blend(id, ierr)
         const_rho = pow7(atm_blend) * s%x_ctrl(1)
         const_csb = atm_blend * s%x_ctrl(2)

         ! Radiation-dominated
         const_Tb_rad = pow((9d0 * pow2(const_csb) * const_rho / (4 * crad)), 0.25d0)

         ! Gas-dominated, assumed non-ionized
         const_Tb_gas = 3*pow2(const_csb)*amu/(5d0 * kerg)

         const_Tb = min(const_Tb_gas, const_Tb_rad)
         min_Teff = atm_blend * s%x_ctrl(7)

      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

        if (s% log_center_temperature >= 9.50) then
            s% use_other_adjust_mdot = .false.
            s% use_other_momentum_implicit = .true.
            s% use_other_energy_implicit = .true.
        end if
         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 22 + s% species
      end function how_many_extra_history_columns

      subroutine M_outside_LSO(id, ierr, m_outside_a0, m_outside_a1)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp), intent(out) :: m_outside_a0, m_outside_a1
         integer :: j, nz
         real(dp) :: Jtot, Itot, omega, jk, j0, j1

         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         nz = s%nz

         Jtot = s%job%extras_rpar(i_mu_J)
         Itot = 0d0
         do j=1,s%nz
            Itot = Itot + (2d0 / 3d0) * pow2(s%r(j)) * s%dm(j)
         end do
         omega = Jtot / Itot

         m_outside_a0 = 0d0
         m_outside_a1 = 0d0
         do j=nz,1,-1
            j0 = sqrt(6d0)*standard_cgrav*s%m(j)/clight
            j1 = standard_cgrav*s%m(j)/clight
            jk = (2d0/3d0) * pow2(s%r(j)) * omega

            if (jk > j0) then
               m_outside_a0 = m_outside_a0 + s%dm(j)
            end if
            if (jk > j1) then
               m_outside_a1 = m_outside_a1 + s%dm(j)
            end if
         end do

      end subroutine M_outside_LSO


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         use chem_def, only: chem_isos
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         integer :: j, species
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         species = s% species

         names(1) = 'accreted_mass' ! Bondi accretion rate  [ g /s ]
         names(2) = 'lost_mass'     ! Super Eddington mass lost [ g /s ]
         names(3) = 'beta_max'      ! Sets extra mixing
         names(4) = 'L_shock'       ! Accretion Shock Luminosity
         names(5) = 'L_stream'      ! Stream Luminosity
         names(6) = 'R_ph'          ! Photospheric Radius
         names(7) = 'R_bondi'       ! R Bondi
         names(8) = 'mu_J'
         names(9) = 'std_J'
         names(10) = 'mu_J_div_J_crit'
         names(11) = 'J_crit'
         names(12) = 'LdivLedd'
         names(13) = 'R_shock'
         names(14) = 'lnT_shock'
         names(15) = 'Atm_Blend'
         names(16) = 'Mdot_Blend'
         names(17) = 'M_above_a=0_LSO'
         names(18) = 'M_above_a=1_LSO'
         names(19) = 'agn_dj_dt'
         names(20) = 'agn_H_p'
         names(21) = 'agn_k_acc'
         names(22) = 'agn_J_dep'


         do j=1,species
            names(22+j) = chem_isos%name(s%chem_id(j))
         end do

         vals(1) = s% xtra1_array(1)
         vals(2) = s% xtra1_array(2)
         vals(3) = s% xtra1_array(3)
         vals(4) = s% xtra1_array(4)
         vals(5) = s% xtra1_array(5)
         vals(6) = s% xtra1_array(6)
         vals(7) = s% xtra1_array(7)
         vals(8) = s%job%extras_rpar(i_mu_J)
         vals(9) = pow(s%job%extras_rpar(i_var_J), 0.5d0)
         vals(10) = vals(8) / s%job%extras_rpar(i_J_crit)
         vals(11) = s%job%extras_rpar(i_J_crit)
         vals(12) = s% xtra1_array(9)  ! Leff / Ledd (from mass_changes.f90)
         vals(13) = s% xtra1_array(10) ! R Shock (R 'Surface' of MESA model) (from atmosphere.f90)
         vals(14) = s% xtra1_array(11) ! log2 T Shock (T 'Surface' of MESA model) (from atmosphere.f90)
         vals(15) = atm_blend
         vals(16) = mdot_blend
         vals(19) = s% job% extras_rpar(i_dj_dt)
         vals(20) = s% job% extras_rpar(i_H_p)
         vals(21) = s% job% extras_rpar(i_k_acc)
         vals(22) = s% job% extras_rpar(i_J_dep)


         call M_outside_LSO(id, ierr, vals(17), vals(18))

         do j=1,species
            vals(22+j) = lost_mass(j) / Msun
         end do


         !note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.


      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 1
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(1) = 'other_dmix'
         do k=1,nz
            vals(k,1) = s%D_mix(k)
         end do

         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do

      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile,
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve





! extra routines for velocity drag in envelope, reimplemnted by -EbF

subroutine v_drag_other_momentum_implicit(id, ierr)
   use auto_diff
   integer, intent(in) :: id
   integer, intent(out) :: ierr
   type (star_info), pointer :: s
   integer :: k
   real (dp) :: min_q_for_drag, drag_coefficient
   type(auto_diff_real_star_order1) :: drag_force, v_00
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return

    min_q_for_drag = s% x_ctrl(13)
    drag_coefficient = s% x_ctrl(14)
    if (s% v_flag) then
        do k=1,s% nz -1
            if (s% q(k) > min_q_for_drag .and. drag_coefficient > 0) then
               v_00 = wrap_v_00(s,k)
               drag_force = -drag_coefficient*v_00/s% dt
               s% extra_grav(k) = drag_force ! additional drag force
            else
                s% extra_grav(k) = 0 ! no drag
            end if
            !write (*,*) 'k =', k
            !write (*,*) 'drag_force', s% extra_grav(k) %val
        end do
    else
        s% extra_grav(1:s%nz) %val = 0d0
    end if

end subroutine v_drag_other_momentum_implicit

subroutine v_drag_other_energy_implicit(id, ierr)
   use const_def, only: Rsun
   use auto_diff
   integer, intent(in) :: id
   integer, intent(out) :: ierr
   real (dp) :: min_q_for_drag, drag_coefficient
   type(auto_diff_real_star_order1) :: drag_eps, drag_force, v_00, v_p1
   type (star_info), pointer :: s
   integer :: k, nz
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return


    min_q_for_drag = s% x_ctrl(13)
    drag_coefficient = s% x_ctrl(14)
    drag_eps = 0d0

    if (s% v_flag) then
        do k=1, s% nz -1
            drag_eps = 0d0
            !s% extra_heat(k) = 0d0
            if (s% q(k) > min_q_for_drag .and. drag_coefficient > 0) then
               v_00 = wrap_v_00(s,k)
               drag_force = drag_coefficient*v_00/s% dt
               drag_eps = 0.5d0*v_00*drag_force
               ! drag energy for outer half-cell.   the 0.5d0 is for dm/2
            end if
            if (s% q(k+1) > min_q_for_drag .and. drag_coefficient > 0) then
               v_p1 = wrap_v_p1(s,k)
               drag_force = drag_coefficient*v_p1/s% dt
               drag_eps = drag_eps + 0.5d0*v_p1*drag_force
               ! drag energy for inner half-cell.   the 0.5d0 is for dm/2
            end if

            if (s% q(k) > min_q_for_drag .and. drag_coefficient > 0) then
                s% extra_heat(k) = drag_eps
            else
                s% extra_heat(k) = 0d0
            end if
            !write (*,*) 'k =', k
            !write (*,*) 'drag_eps', s% extra_heat(k) %val
        end do
     else
        s% extra_heat(1:s%nz) %val = 0d0
     end if
end subroutine v_drag_other_energy_implicit
      end module run_star_extras
