#[directories]
input_dir  = "host-trajectories/"
output_dir = "output-trajectories/"
log_dir    = "/logs/"

#[default]
host_trajectory_file = "GL-01_sun_58122.xyzv"
output_trajectory_file = "result.traj"
output_ephemeris_file = "ephemeris.traj"
jettison_index = 0
interval = 1./4    # Unit: day (e.g. 1./24 means one step per hour)
steps = 990        # How many steps basing on interval you want? (e.g. 24 steps * 1.24(interval) = 1 day )
show_plot = "False" # True or False
refine_pass = 30    # Maximum number of "fmin_cobyla" function evaluations.

max_deltav = 30
lambert_coarse_delta=1
lambert_coarse_range=50
lambert_fine_delta=1/24
lambert_fine_range=50

