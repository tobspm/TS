#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
