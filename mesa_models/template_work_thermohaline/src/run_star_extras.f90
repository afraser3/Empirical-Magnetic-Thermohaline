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
      use chem_def !! so it knows what ih2 is !MJoyce 12/1/21
      
      implicit none

      real(dp) :: log_g, log_g_end, log_g_outputs
      real(dp) :: log_g_max = -4

      real(dp) :: start_output_age, end_output_age

      real(dp) :: dense_mesh_coeff, dense_time_coeff, dense_varcontrol_target
      logical :: is_highRes = .false.

      real(dp) :: prev_luminosity, enter_LB, exit_LB 

      
      ! these routines are called by the standard run_star check_model
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         log_g_end = s% x_ctrl(12)
         log_g_outputs = s% x_ctrl(13)


         dense_mesh_coeff = s% x_ctrl(15)
         dense_time_coeff = s% x_ctrl(16)
         dense_varcontrol_target = s% x_ctrl(17)
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

         start_output_age = s% x_ctrl(20)
         end_output_age = s% x_ctrl(21)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

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

      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         enter_LB = 0d0
         exit_LB = 0d0 
      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         real(dp) :: sb_sigma, mu_ideal_gas, R2, R_gas_constant !! pi can be deleted
         type (star_info), pointer :: s
         real(dp) :: power_he_burn, power_c_burn, power_neutrinos, &
         center_h1, center_he4, ocz_top_mass, ocz_bot_mass, &
         ocz_top_radius, ocz_bot_radius!, mass_difference_prev !! no!!
         integer :: nz, j, i, k, k_ocz_top, k_ocz_bot, n_conv_bdy
!         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0

      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         

         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model




      subroutine compute_Ra_Re(id, ierr, v, Ra, Re, Pr, k)
         use const_def

         type (star_info), pointer :: s
         integer, intent(in) :: id, k 
         real(dp), intent(in) :: v
         integer, intent(out) :: ierr
         real(dp), intent(out) :: Ra, Re, Pr


         real(dp) :: alpha, lnLambda, nu
         real(dp) :: dz, dr

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         alpha = 16d0 * boltz_sigma * pow3(s%T(k)) / (3d0 * s%opacity(k) * s%Cp(k) * pow2(s%rho(k)))
         if (s%T(k) > 1.1d5) then
            lnLambda = -17.4d0 + 1.5d0 * log(s%T(k)) - 0.5d0 * log(s%rho(k))
         else
            lnLambda = -12.7d0 + log(s%T(k)) - 0.5d0 * log(s%rho(k))
         end if
         nu = 5.2d-15 * pow(s%T(k),2.5d0) / s%rho(k) / lnLambda


         Pr = nu / alpha
         Re = v * dz / nu
         Ra = Pr * pow2(Re)
      end subroutine compute_Ra_Re


       subroutine compute_Pr(id, ierr, Pr, k)
         use const_def

         type (star_info), pointer :: s
         integer, intent(in) :: id, k 
         integer, intent(out) :: ierr
         real(dp), intent(out) :: Pr


         real(dp) :: alpha, lnLambda, nu
         real(dp) :: dz, dr

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         alpha = 16d0 * boltz_sigma * pow3(s%T(k)) / (3d0 * s%opacity(k) * s%Cp(k) * pow2(s%rho(k)))
         if (s%T(k) > 1.1d5) then
            lnLambda = -17.4d0 + 1.5d0 * log(s%T(k)) - 0.5d0 * log(s%rho(k))
         else
            lnLambda = -12.7d0 + log(s%T(k)) - 0.5d0 * log(s%rho(k))
         end if
         nu = 5.2d-15 * pow(s%T(k),2.5d0) / s%rho(k) / lnLambda


         Pr = nu / alpha
      end subroutine compute_Pr



      subroutine compute_diff_coeffs(id, ierr, kt, kmu, vis, k)
         use chem_def, only: chem_isos, ih1
         use const_def

         type (star_info), pointer :: s
         integer, intent(in) :: id, k
         integer, intent(out) :: ierr
         real(dp), intent(out) :: kt, kmu, vis

         real(dp) :: loglambdah,loglambdacx,loglambdacy,ccx,ccy
         real(dp) :: Bcoeff, Kc, nu_rad, nu_mol
         real(dp) :: chemA, chemZ, acx, acy, qe4
         real(dp), parameter :: sqrt5 = sqrt(5d0)
         integer :: h1, dom_iso


         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         kmu = 1.0
         kt = 1.0
         vis = 1.0 


         h1 = s% net_iso(ih1)

         Kc = 4*crad*clight*s% T(k)*s% T(k)*s% T(k)/(3*s% opacity(k)*s% rho(k))  ! thermal conductivity

         kt = Kc/(s% Cp(k) *s% rho(k))            ! thermal diffusivity (assumes radiatively dominated)

         !write(*,*) "Kc,Kt",kc,kt
         qe4=qe*qe*qe*qe                          ! From const_def 
         
 

         ! Log Lambda for pure H (equation 10 from Proffitt Michaud 93)
         loglambdah = -19.26d0 - 0.5d0*log(s% rho(k)) + 1.5d0*log(s% T(k)) - 0.5d0*log(1d0 + 0.5d0*(1+s% xa(h1,k))) 
         nu_rad = 4d0*crad*s%T(k)*s%T(k)*s%T(k)*s%T(k)/(15d0*clight*s%opacity(k)*s%rho(k)*s%rho(k)) ! radiative viscosity
         nu_mol = 0.406d0*sqrt(amu)*pow(boltzm*s%T(k),2.5d0)/(qe4*loglambdah*s%rho(k)) 
         ! From Spitzer "Physics of Fully Ionized Gases equation 5-54
         ! Assumes pure H. Still trying to work out what it would be for a mixture. 
         vis = nu_mol + nu_rad   ! total viscosity

         ! The following is from Proffitt & Michaud, 1993.
         ! Their constant B (equation 15)
         Bcoeff = (15.d0/16.d0)*sqrt(2.d0*amu/(5*pi))*pow(boltzm,2.5d0)/qe4
         ! Extract what species drives the thermohaline concvection
         dom_iso = s% dominant_iso_for_thermohaline(k)
         chemA = chem_isos%Z_plus_N(dom_iso)
         chemZ = chem_isos%Z(dom_iso)


          if(chemZ.gt.2) then
         ! This is if the driving chemical is NOT He.
            ! Log Lambda for H-dominant chem mixture (equation 10)
            loglambdacx = loglambdah - log(chemz)  
            ! Log Lambda for He-dominant chem mixture (equation 10)
            loglambdacy = loglambdah - log(2.d0*chemz)
            ! Calculation of C_ij coeffs (equation 12)
            ccx = log(exp(1.2d0*loglambdacx)+1.)/1.2d0
            ccy = log(exp(1.2d0*loglambdacy)+1.)/1.2d0
            ! Reduced masses (I had to guess, from Bahcall & Loeb 1990), with H and He
            acx = (1.d0*chemA)/(1.d0+chemA)
            acy = 4*chemA/(4.d0+chemA)
            ! My formula (see notes) based on Proffitt and Michaud 1993
            kmu = 2*Bcoeff*pow(s% T(k),2.5d0)/(sqrt5*s%rho(k)*chemZ*chemZ)/ &
               (s% xa(h1,k)*sqrt(acx)*ccx + (1-s% xa(h1,k))*sqrt(acy)*ccy)

         else
            ! Log Lambda for H-He mixture (equation 10)
            loglambdah = -19.26d0 - log(2d0) - 0.5d0*log(s%rho(k)) + &
               1.5d0*log(s%T(k)) - 0.5d0*log(1d0 + 0.5d0*(1+s% xa(h1,k))) 
            ! Calculation of C_ij coeffs (equation 12)
            ccy = log(exp(1.2d0*loglambdah)+1d0)/1.2d0
            ! My formula (see notes) based on Proffitt and Michaud 1993
            kmu = (Bcoeff*pow(s%T(k),2.5d0)/(s%rho(k)*ccy))*(3+s% xa(h1,k))/((1+s% xa(h1,k))*(3+5*s% xa(h1,k))*(0.7d0+0.3d0*s% xa(h1,k)))
            
         endif

      end subroutine compute_diff_coeffs


    integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         real(dp) :: core_h2, current_h2, surface_d
         integer :: nz  !, phase_of_evolution
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      
      end subroutine data_for_extra_history_columns
      

      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 4
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n), Pr, kt, kmu, vis
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         

         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'

         Pr = 0 

         names(1) = 'Prandtl'
         names(2) = 'kt'
         names(3) = 'kmu'
         names(4) = 'vis'
         do k = 1, nz
            call compute_Pr(id, ierr, Pr, k) 
            vals(k,1) = safe_log10(Pr) 
            call compute_diff_coeffs(id, ierr, kt, kmu, vis, k)
            vals(k,2) = safe_log10(kt)
            vals(k,3) = safe_log10(kmu)
            vals(k,4) = safe_log10(vis)
         end do
         
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
      end subroutine data_for_extra_profile_header_items

      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         real(dp) :: log_surface_gravity
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.
         
         !write(*,*) s% L(1)/lsun, prev_luminosity/lsun, s% L(1)/lsun

         if (( s% L(1)/lsun .gt. 10) .and. &
            ( s% star_age .gt. 5d9) .and. &
             (prev_luminosity .gt. s% L(1)) .and. &
             (enter_LB .eq. 0d0)) then ! Entering luminosity Bump (luminosity decreases)
!         if (( s% L(1)/lsun .gt. 10) .and. (prev_luminosity .gt. s% L(1)) .and. (enter_LB .eq. 0d0)) then ! Entering luminosity Bump (luminosity decreases)

          ! save a profile & update the history
           write(*,*) 'Entering LB'
           s% need_to_update_history_now = .true.
           s% need_to_save_profiles_now = .true.
           enter_LB = s% L(1)


         endif 

         if ((enter_LB .gt. 0d0) .and. (s% L(1) .gt. enter_LB) .and. (exit_LB .eq. 0d0)) then ! Exiting luminosity Bump (same luminosity as entering, but moving up)
            
          ! save a profile & update the history
           write(*,*) 'Exiting LB'
           s% need_to_update_history_now = .true.
           s% need_to_save_profiles_now = .true.
           exit_LB = s% L(1)
           !extras_finish_step = terminate

         endif 

         log_g = safe_log10(s% grav(1))

         if ( (log_g > log_g_max) ) then
             log_g_max = log_g
         endif


         if (.not. is_highRes) then
             write(*,*) 'my value of logg is', log_g, log_g_max, log_g_outputs
             if ( (log_g < log_g_max) .and. (log_g_max > log_g_outputs) \
                  .and. (log_g < log_g_outputs) .and. (s% star_age >= start_output_age) ) then
                 is_highRes = .true.
                 s% write_profiles_flag = .true.
                 s% time_delta_coeff = dense_time_coeff
                 s% mesh_delta_coeff = dense_mesh_coeff
                 s% varcontrol_target = dense_varcontrol_target
             else
                 s% write_profiles_flag = .false.
             endif
         endif

         if ( ( (log_g < log_g_end) .or. (s% star_age >= end_output_age) ) \
                .and. (s% star_age >= 10d6) ) then
            extras_finish_step = terminate
         endif       

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

      end module run_star_extras
      
