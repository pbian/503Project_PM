#!/bin/bash
cd particles_parallel_2D_program
echo "Running parallel particles program"
./particles_parallel_2D
cd ..

cd particles_serial_2D_program
echo "Running serial particles program"
./particles_serial_2D
cd ..
echo "Sys_Cond.txt diff: =====>"
diff particles_serial_2D_program/Sys_Cond.txt particles_parallel_2D_program/Sys_Cond.txt
echo "PS_Phi_Initial.txt diff: =====>"
diff particles_serial_2D_program/PS_Phi_Initial.txt	particles_parallel_2D_program/PS_Phi_Initial.txt
echo "Part_Matrix.txt diff: =====>"
diff particles_serial_2D_program/Part_Matrix.txt particles_parallel_2D_program/Part_Matrix.txt
echo "PS_Phi_Final.txt diff: =====>"
diff particles_serial_2D_program/PS_Phi_Final.txt particles_parallel_2D_program/PS_Phi_Final.txt
echo "PS_ni_Final.txt diff: =====>"
diff particles_serial_2D_program/PS_ni_Final.txt particles_parallel_2D_program/PS_ni_Final.txt
echo "Particle_Count.txt diff: =====>"
diff particles_serial_2D_program/Particle_Count.txt particles_parallel_2D_program/Particle_Count.txt

