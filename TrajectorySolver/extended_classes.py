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
import PyKEP as pk
import numpy as np

class BasePlanet(pk.planet.jpl_lp):
    def __init__(self, *args, **kwargs):
        super(BasePlanet, self).__init__(*args, **kwargs)
        self.semi_major_axis, self.eccentricity, self.inclination, self.longitude_asc_node, self.arg_periapsis, _ = self.osculating_elements(pk.epoch(0))
        self.soi = self.semi_major_axis*(self.mu_self/self.mu_central_body)**(2./5)

    def get_relative_position(self, date, position):
        pla_pos, _ = np.array(self.eph(date))
        return (position - pla_pos)

    def get_acceleration(self, date, position):
        rel_position = self.get_relative_position(date, position)
        acc = -self.mu_self*rel_position/(np.linalg.norm(rel_position)**3)
        return acc
