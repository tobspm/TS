#[directories]
input_dir  = "/home/birdy/Software/TS_tests/bidule/"
output_dir = "/home/birdy/Software/TS_tests/results/"
log_dir    = "/home/birdy/Software/TS_tests/results/"

#[default]
host_trajectory_file = "host_traj.xyzv"
output_trajectory_file = "result.traj"
output_ephemeris_file = "ephemeris.traj"

#------------- Nima
output_ephemeris_vts_format_file = "ephemeris_vts_format.traj"

#-----------------EOM---

jettison_index = 0
interval = 1./4    # Unit: day (e.g. 1./24 means one step per hour)
steps = 5       # How many steps basing on interval you want? (e.g. 24 steps * 1.24(interval) = 1 day )
show_plot = "False" # True or False
refine_pass = 30    # Maximum number of "fmin_cobyla" function evaluations.

max_deltav = 30
lambert_coarse_delta=1
lambert_coarse_range=50
lambert_fine_delta=1/24
lambert_fine_range=50

