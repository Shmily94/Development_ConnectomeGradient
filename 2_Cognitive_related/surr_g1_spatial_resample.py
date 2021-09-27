# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 17:35:08 2021

@author: yeban
"""

from brainsmash.workbench.geo import volume
from brainsmash.mapgen.sampled import Sampled
import scipy.io as sio
coord_file = "voxel_coordinates.txt"
output_dir = ".../Neurosynth_correction"

filenames = volume(coord_file, output_dir)

brain_map = "g1_t.txt"
gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'], resample=True)
surrogate_maps = gen(n=10000)
sio.savemat('surrogate_maps_g1_t_resample.mat',{'surrogate_maps':surrogate_maps})
