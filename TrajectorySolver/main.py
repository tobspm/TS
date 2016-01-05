#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Final modified date: 23/10/2015
# Original author: Lucas Serrano
# Maintenance: Jim, Lin
# Email: f44006076@gmail.com
#
# ==============================================================
# ===========================LICENSE============================
# ==============================================================
# This file is part of DOCKS.
#
# DOCKS is free software: you can redistribute it and/or modify
# it under the terms of the  GNU LESSER GENERAL PUBLIC LICENSE 
# as published by the Free Software Foundation, either version 
# 3 of the License, or any later version.
#
# DOCKS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU LESSER GENERAL PUBLIC LICENSE for more details.
#
# You should have received a copy of the 
# GNU LESSER GENERAL PUBLIC LICENSE along with DOCKS.
# If not, see <http://www.gnu.org/licenses/lgpl.txt>.
#
#
# Update process:
# v1.0 - original version (can not run successfully)
# v2.0 - 1. Modify/update the class name of PyKEP package 
#        2. Import parameters from config file, replace the constant in program with variables
#        3. Add the flag for detemining if the plot will show up or not
#        4. Confirm the data result is similar to the one created from Scilab version
# v2.1 - Change some parameters in order to be implemented by DOCKing system

from __future__ import division
import numpy as np
import PyKEP as pk
import scipy as sp

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot
from scipy.optimize import fmin_cobyla

import extended_classes
import trajectory_io
import io

# Load parameters from config file
import config
# Default parameters of directory setting
input_dir  = "host-trajectories/"
output_dir = "output-trajectories/"
log_dir    = "/logs/"

def get_soi_exit_index(body, dates, positions):
    soi = body.soi
    for index, position in enumerate(positions):
        rel_position = body.get_relative_position(dates[index], position)
        if np.linalg.norm(rel_position) > soi:
            return index
    else:
        raise SOIExitNotFound

def search_min_dv_to_body(body, departure_date, departure_position, host_velocity, min_time_of_flight, time_delta, number):
    time_range = [min_time_of_flight + time_delta*i for i in xrange(number)]
    body_positions = [body.eph(time + departure_date.mjd2000)[0] for time in time_range]
    departure_velocities = [pk.lambert_problem(departure_position, pos, time*pk.DAY2SEC, pk.MU_SUN, False, 0).get_v1()[0] for pos, time in zip(body_positions, time_range)]
    deltaV = [np.linalg.norm(np.array(velocity)-host_velocity) for velocity in departure_velocities]
    index_min = np.array(deltaV).argmin()
    return index_min, departure_velocities[index_min]

def rk4(function, time, values, time_interval, *args):
    semi_interval = time_interval*1./2
    d1 = function(time, values, *args)
    d2 = function(time+semi_interval, values + d1*semi_interval, *args)
    d3 = function(time+semi_interval, values + d2*semi_interval, *args)
    d4 = function(time+time_interval, values + d3*time_interval, *args)

    next_values = values + 1./6*(d1 + 2*d2 + 2*d3 + d4)*time_interval
    return next_values

def function(time, pos_vel, bodies):
    pos = pos_vel[0:3]
    vel = pos_vel[3:6]
    acc = sum(body.get_acceleration(time*pk.SEC2DAY, pos) for body in bodies)
    acc = acc -pk.MU_SUN*pos/(np.linalg.norm(pos)**3)
    return np.append(vel, acc)

def compute_trajectory(start_date, start_position, start_velocity, steps_number, time_delta, bodies):
    pos_vel = np.append(start_position, start_velocity)
    time_delta = time_delta*pk.DAY2SEC
    values = [np.append(start_date.mjd2000, pos_vel)]
    start_date = start_date.mjd2000*pk.DAY2SEC
    time_range = [start_date + time_delta*i for i in xrange(steps_number)]
    for i, time in enumerate(time_range):
        values.append(np.append(time*pk.SEC2DAY, rk4(function, time, values[i][1:7], time_delta, bodies)))
    return values

def minimum_distance_to_body(trajectory, body):
    distances = [np.linalg.norm(body.get_relative_position(datum[0], datum[1:4])) for datum in trajectory]
    min_distance = min(distances)
    return min_distance

def function_to_minimize(velocity, init_date, init_pos, bodies):
    trajectory_values = compute_trajectory(init_date, init_pos, velocity, 6000, config.interval, bodies)
    min_distance = minimum_distance_to_body(trajectory_values, bodies[1])
    print min_distance
    return min_distance

def const(velocity, host_vel):
    dv = np.linalg.norm(velocity - host_vel)
    return 200-dv

if __name__ == '__main__':
	host_dates = []
	host_post = []
	host_vel = []
	with open(input_dir + config.host_trajectory_file) as input_file:
		host_dates, host_pos, host_vel = trajectory_io.parse_trajectory(input_file)
	earth = extended_classes.BasePlanet('earth')
	mars = extended_classes.BasePlanet('mars')
	jupiter = extended_classes.BasePlanet('jupiter')

	soi_exit_index = get_soi_exit_index(earth, host_dates, host_pos)
	min_index, velocity = search_min_dv_to_body(mars, host_dates[soi_exit_index-5], host_pos[soi_exit_index-5], host_vel[soi_exit_index-5], 227-50, 1, 100)
	min_index_2, velocity = search_min_dv_to_body(mars, host_dates[soi_exit_index-5], host_pos[soi_exit_index-5], host_vel[soi_exit_index-5], 227-(min_index-50)-50./24, config.interval, 100)

	bounds = [(axis_vel-20, axis_vel+20) for axis_vel in velocity]
	trajectory_values = compute_trajectory(host_dates[soi_exit_index-5], host_pos[soi_exit_index-5], velocity, config.steps, config.interval, [earth, mars])
	res = fmin_cobyla(function_to_minimize, velocity, const, (host_dates[soi_exit_index-5], host_pos[soi_exit_index-5], [earth, mars]), (host_vel[soi_exit_index-5],), maxfun=config.refine_pass)
	print "initial velocity:"
	print velocity
	print "optimized velocity:"
	print res
	opti_trajectory = compute_trajectory(host_dates[soi_exit_index-5], host_pos[soi_exit_index-5], res, config.steps, config.interval, [earth, mars])

	print "initial DV: %s" % np.linalg.norm(host_vel[soi_exit_index-5]-velocity)
	print "optimized trajectory DV: %s" % np.linalg.norm(host_vel[soi_exit_index-5]-res)

	# Transform the matrix	
	zipped = zip(*trajectory_values)
	zipped2 = zip(*opti_trajectory)

	# Write the trajectory result to specified file
	print '\nwriting trajectory'
	with open(output_dir + config.output_trajectory_file, 'w') as output_file:
		trajectory_io.write_output(output_file, opti_trajectory, earth, mars)
	print "Success for creating the trajectory, the file name is: %s" % config.output_trajectory_file	

	# Write the ephemeris of Jupiter result to specified file ---jim 17092015---
	print '\nwriting ephemeris of Jupiter'
	with open(output_dir + config.output_ephemeris_file, 'w') as output_file:
		trajectory_io.jup_eph_gen(output_file, opti_trajectory, jupiter)
	print "Success for creating the ephemeris of Jupiter, the file name is: %s" % config.output_ephemeris_file
	

	# Plot the result
	if config.show_plot == 'True':	
		ax = trajectory_io.plot_trajectory(host_dates, host_pos, 227-(min_index-50))
		ax.plot(zipped[1], zipped[2], zipped[3], color='blue')
		ax.plot(zipped2[1], zipped2[2], zipped2[3], color='green')
		fig, ax2 = pyplot.subplots()
		ax2.plot([np.linalg.norm(mars.get_relative_position(host_dates[soi_exit_index-5].mjd2000 + i*config.interval, traj[1:4])) for i, traj in enumerate(trajectory_values)], color='blue')
		ax2.plot([np.linalg.norm(mars.get_relative_position(host_dates[soi_exit_index-5].mjd2000 + i*config.interval, traj[1:4])) for i, traj in enumerate(opti_trajectory)], color='green')
		pyplot.show()

	print '\nEnd of program'	
#===End of program===   
