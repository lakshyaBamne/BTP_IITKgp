#! for MCW, EXP, IMP, KHI
gfortran Euler_gas_main_2D.f90 initialset_2D.f90 extendcell2D.f90 numerical_flux_Euler2D.f90 SSPRK_Eulergas2D.f90 -o prg.x

#! for 2DR
# gfortran Euler_gas_main_2D.f90 ini_Riemann.f90 extendcell2D.f90 numerical_flux_Euler2D.f90 SSPRK_Eulergas2D.f90 -o prg.x

#! for RT
# gfortran Euler_gas_main_2D.f90 RaleighTaylor/initialset_2DRT.f90 extendcell2D.f90 RaleighTaylor/numerical_flux_RT.f90 RaleighTaylor/SSPRK_EulergasRT.f90 -o prg.x

#! run the program
./prg.x

python plot.py
