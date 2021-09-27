# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 20:26:03 2021

@author: mingr
"""
#import neurosynth as ns
#ns.dataset.download(path='.', unpack=True)

from neurosynth.base.dataset import Dataset
from neurosynth import decode

#dataset = Dataset('database.txt')
#dataset.add_features('features.txt')
#dataset.save('dataset.pkl')
dataset = Dataset.load('dataset.pkl')
decoder = decode.Decoder(dataset)
data = decoder.decode('Reslice_pos_r_correct.nii', save='real_r_grad_pos_all_correct.txt')
data = decoder.decode('Reslice_neg_r_correct.nii', save='real_r_grad_neg_all_correct.txt')

