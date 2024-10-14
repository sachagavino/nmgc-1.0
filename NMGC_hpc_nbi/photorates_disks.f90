module photorates
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author
!> Sacha Gavino & Valentine Wakelam
!
!> @date 2019
!
! DESCRIPTION: Subroutines to computes the new photorate
!              calculation as in Haeys et al. (2017); Gavino et al. (2021)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    contains

    subroutine compute_molopacity
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! description: computes molecular opacities due to main species
    !              contributing to absorption. Read tables composed
    !              of absorption, dissociation and ionization cross
    !              sections. And uses molecular densities tables. 
    !              --> EQUATION (18) in documentation. VERTICAL
    !    
    ! parameters: -mol_opacity   is \tau_m^V(\nu,r,z)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use utilities
        use global_variables
        implicit none
 
        mol_opacity(:) = 0.D0

        !!!!!! Computation of the opacity due to molecules H2, CO, N2, H2O and a selection of other important molecules
        !  are global variables and the calculation is done at each altitude so that the column density is computed
        ! from the top to the bottom (most dense part).


        mol_opacity(:) = NH2*(cross_sections(:,id_photo_H2,2) + cross_sections(:,id_photo_H2,3) + cross_sections(:,id_photo_H2,4)) ! H2
        mol_opacity(:) = mol_opacity(:) + NCO*(cross_sections(:,id_photo_CO,2) + cross_sections(:,id_photo_CO,3) + &
                    cross_sections(:,id_photo_CO,4)) ! CO
        mol_opacity(:) = mol_opacity(:) + NN2*(cross_sections(:,id_photo_N2,2) + cross_sections(:,id_photo_N2,3) + &
                    cross_sections(:,id_photo_N2,4)) ! N2
        mol_opacity(:) = mol_opacity(:) + NH2O*(cross_sections(:,id_photo_H2O,2) + cross_sections(:,id_photo_H2O,3) + &
                    cross_sections(:,id_photo_H2O,4)) ! H2O
        mol_opacity(:) = mol_opacity(:) + NCH*(cross_sections(:,id_photo_CH,2) + cross_sections(:,id_photo_CH,3) + &
                    cross_sections(:,id_photo_CH,4)) ! CH 
        mol_opacity(:) = mol_opacity(:) + NCH3*(cross_sections(:,id_photo_CH3,2) + cross_sections(:,id_photo_CH3,3) + &
                    cross_sections(:,id_photo_CH3,4)) ! CH3
        mol_opacity(:) = mol_opacity(:) + NCH4*(cross_sections(:,id_photo_CH4,2) + cross_sections(:,id_photo_CH4,3) + &
                    cross_sections(:,id_photo_CH4,4)) ! CH4   
        mol_opacity(:) = mol_opacity(:) + NOH*(cross_sections(:,id_photo_OH,2) + cross_sections(:,id_photo_OH,3) + &
                    cross_sections(:,id_photo_OH,4)) ! OH
        mol_opacity(:) = mol_opacity(:) + NNH*(cross_sections(:,id_photo_NH,2) + cross_sections(:,id_photo_NH,3) + &
                    cross_sections(:,id_photo_NH,4)) ! NH 
        mol_opacity(:) = mol_opacity(:) + NNH2*(cross_sections(:,id_photo_NH2,2) + cross_sections(:,id_photo_NH2,3) + &
                    cross_sections(:,id_photo_NH2,4)) ! NH2  
        mol_opacity(:) = mol_opacity(:) + NNH3*(cross_sections(:,id_photo_NH3,2) + cross_sections(:,id_photo_NH3,3) + &
                    cross_sections(:,id_photo_NH3,4)) ! NH3
        mol_opacity(:) = mol_opacity(:) + NCO2*(cross_sections(:,id_photo_CO2,2) + cross_sections(:,id_photo_CO2,3) + &
                    cross_sections(:,id_photo_CO2,4)) ! CO2
        mol_opacity(:) = mol_opacity(:) + NN2O*(cross_sections(:,id_photo_N2O,2) + cross_sections(:,id_photo_N2O,3) + &
                    cross_sections(:,id_photo_N2O,4)) ! N2O    
        mol_opacity(:) = mol_opacity(:) + NHCN*(cross_sections(:,id_photo_HCN,2) + cross_sections(:,id_photo_HCN,3) + &
                    cross_sections(:,id_photo_HCN,4)) ! HCN
        mol_opacity(:) = mol_opacity(:) + NHNC*(cross_sections(:,id_photo_HNC,2) + cross_sections(:,id_photo_HNC,3) + &
                    cross_sections(:,id_photo_HNC,4)) ! HNC
        mol_opacity(:) = mol_opacity(:) + NCN*(cross_sections(:,id_photo_CN,2) + cross_sections(:,id_photo_CN,3) + &
                    cross_sections(:,id_photo_CN,4)) ! CN
        mol_opacity(:) = mol_opacity(:) + NHCO*(cross_sections(:,id_photo_HCO,2) + cross_sections(:,id_photo_HCO,3) + &
                    cross_sections(:,id_photo_HCO,4)) ! HCO
        mol_opacity(:) = mol_opacity(:) + NH2CO*(cross_sections(:,id_photo_H2CO,2) + cross_sections(:,id_photo_H2CO,3) + &
                    cross_sections(:,id_photo_H2CO,4)) ! H2CO              


  
    end subroutine compute_molopacity




    subroutine compute_local_flux()
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! description: exctract the ISRF and stellar flux spectr from
    !              tables then use them to compute the local flux 
    !              at each coordinates.
    !              --> EQUATION (14) in documentation.
    !
    ! parameters: -Table_local_dust   is I_d(\nu,r,z)
    !             -table_ISRF         is I_{ISRF}(\nu)
    !             -dust_opacity       is \tau_d^V(\nu,r,z)
    !             -Table_Istar        is I_*(\nu,r,z) or I_*^0 (R_0^2/(r^2+z^2))
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        use utilities
        use global_variables
        implicit none

        !local_flux(:) = local_flux_dust(:,x_i)*exp(-(mol_opacity(:)))
        local_flux(:) = local_flux_dust(:,x_i)   !radial compuation

        !UV_FLUX =  sum(local_flux_dust(1:wv_max_factor,1))/sum(flux_isrf_dust(1:wv_max_factor,1)) ! UV factor for the original rates calculation.
        !INT_LOCAL_FLUX =  sum(local_flux(1:wv_max_factor))/sum(flux_isrf_dust(1:wv_max_factor,1)) ! UV factor using molecular opacity.
        INT_LOCAL_FLUX =  sum(local_flux(1:wv_max_factor))/sum(flux_isrf_dust(1:wv_max_factor,x_i)) ! UV factor for radial computation
        !write(*,*) 'INT_LOCAL_FLUX rat6: ', INT_LOCAL_FLUX

    end subroutine compute_local_flux



    subroutine compute_photorates()
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! description: compute the rates (diss and ion) using the 
    !              values of the flux and cross sections. 
    !              --> EQUATION (6) in documentation.
    !
    ! parameters: -cross-section    is sigma_diss(X,\nu)
    !             -local_flux       is I_L(\nu,r,z)
    !             -rate_diss        is K_d(X)   
    !             -rate_ion         is K_i(X)   
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        use utilities
        use global_variables
        implicit none

        !----------------------- VARIABLES -----------------------
        integer :: i, j
        integer :: nb_cross_comments

        real(8) :: d_lambda = 1.D0
        !---------------------------------------------------------


        !------INITIALIZE------    
        nb_cross_comments = 0
        photorates%name = 'empty'
        photorates%k_diss = 0.D0
        photorates%k_ion = 0.D0
        !----------------------

        !------LOOP on name of the molecules in cross-sections folder------
        DO i = 1, nb_line_spec_photorates ! number of molecules with cross-sections
   
            !-----COMPUTE photorates.
            photorates(i)%name = spec_photo(i)
            photorates(i)%k_diss = sum(local_flux(:)*cross_sections(:,i,3))*d_lambda
            photorates(i)%k_ion = sum(local_flux(:)*cross_sections(:,i,4))*d_lambda

        END DO ! end loop on molecule names.
        !-------------------------------------------------
    end subroutine compute_photorates

end module photorates
