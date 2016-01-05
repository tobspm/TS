#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

import numpy as np
import PyKEP as pk
import scipy as sp
import io
# ---Jim 17092015 ---
import math as m

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot

def parse_trajectory(trajectory_file):
    dates = []
    positions = []
    velocities = []
    for line in trajectory_file.read().splitlines():
        values = [float(element) for element in line.split(' ')]
        dates.append(pk.epoch(values[0], 'jd'))
        positions.append([values[1], values[2], values[3]])
        velocities.append([values[4], values[5], values[6]])
    positions = np.array(positions)*1000
    velocities = np.array(velocities)*1000
    return dates, positions, velocities

def write_output(output_file, trajectory, earth, mars):
    def tab_write(value):
        output_file.write('%s\t' % value)

    for value in trajectory:
        time = value[0]
        pos = value[1:4]
        tab_write(pk.epoch(time).jd)
        tab_write(pos[0]/1000.)  # the position of the CubeSat along the X axis (in km)
        tab_write(pos[1]/1000.)
        tab_write(pos[2]/1000.)
        tab_write(value[4])      # the velocity of the CubeSat along the X axis (in m/s)
        tab_write(value[5])
        tab_write(value[6])
        tab_write(np.linalg.norm(pos)/1000.)  # the radii of the CubeSat from the Sun (in km)
        tab_write(np.linalg.norm(earth.get_relative_position(time, pos))/1000.)  # the radii of the CubeSat from the Earth (in km)
        tab_write(np.linalg.norm(mars.get_relative_position(time, pos))/1000.)   # the radii of the CubeSat from the Mars  (in km)
        tab_write(np.linalg.norm(earth.eph(time))/1000.) 			 # the radii of Earth from the Sun (in km)
        output_file.write('%s' % (np.linalg.norm(mars.eph(time))/1000.) )	 # the radii of Mars from the Sun (in km)
        output_file.write('\n')

def plot_trajectory(dates, positions, flight_duration):
    earth = pk.planet.jpl_lp('earth')
    mars = pk.planet.jpl_lp('mars')

    fig = pyplot.figure()
    ax = fig.gca(projection='3d')
    unzip = zip(*positions)
    pk.orbit_plots.plot_planet(earth, ax=ax, t0=dates[0])
    mars_arrival = pk.epoch(dates[0].mjd2000 + flight_duration)
    pk.orbit_plots.plot_planet(mars, ax=ax, t0=mars_arrival)
    ax.plot(unzip[0], unzip[1], unzip[2], color='red')
    return ax

# --- Jim 17092015 --- function for generating the ephemeris of Jupiter by traj.xyzv 
def jup_eph_gen(output_file, trajectory, jupiter):
	def tab_write(value):
		output_file.write('%s\t' % value)

	def car2sph(car_pos):
		r = np.linalg.norm(car_pos)  # km
		lat = np.arcsin(car_pos[2]/r)*180./np.pi   # degree
		lon = np.arctan2(car_pos[1], car_pos[0])*180./np.pi  # degree
		return np.array([lon,lat,r/1000.])
	
	for value in trajectory:
		time = value[0]
		pos = value[1:4]
		sph_temp = car2sph( jupiter.get_relative_position(time, pos) )
		tab_write(pk.epoch(time).jd)
		tab_write(sph_temp[1]) # latitude 
		tab_write(sph_temp[0]) # longitude  
		tab_write(sph_temp[2]) 
		output_file.write('\n')


