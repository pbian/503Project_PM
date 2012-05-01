cd particles_parallel_program
openmpirun --hostfile ../hosts -np 5 particles_parallel
cd ..
cd particles_serial_program
./particles_serial
cd ..
echo "Amatrix diff: =====>"
diff particles_serial_program/Amatrix.txt particles_parallel_program/Amatrix.txt
echo "num_den_profile_final diff: =====>"
diff particles_serial_program/num_den_profile_final.txt particles_parallel_program/num_den_profile_final.txt
echo "part_a_vx diff: =====>"
diff particles_serial_program/part_a_vx.txt particles_parallel_program/part_a_vx.txt
echo "part_b_vx diff: =====>"
diff particles_serial_program/part_b_vx.txt particles_parallel_program/part_b_vx.txt
echo "part_v_hist_final diff: =====>"
diff particles_serial_program/part_v_hist_final.txt particles_parallel_program/part_v_hist_final.txt
echo "velocity_temp_hist diff: =====>"
diff particles_serial_program/velocity_temp_hist.txt particles_parallel_program/velocity_temp_hist.txt
echo "energy diff: =====>"
diff particles_serial_program/energy.txt particles_parallel_program/energy.txt
echo "num_den_profile_initial diff: =====>"
diff particles_serial_program/num_den_profile_initial.txt particles_parallel_program/num_den_profile_initial.txt        
echo "part_a_xpos diff: =====>"
diff particles_serial_program/part_a_xpos.txt particles_parallel_program/part_a_xpos.txt
echo "part_b_xpos diff: =====>"
diff particles_serial_program/part_b_xpos.txt particles_parallel_program/part_b_xpos.txt
echo "part_v_hist_initial diff: =====>"
diff particles_serial_program/part_v_hist_initial.txt particles_parallel_program/part_v_hist_initial.txt
echo "v_th diff: =====>"
diff particles_serial_program/v_th.txt particles_parallel_program/v_th.txt

