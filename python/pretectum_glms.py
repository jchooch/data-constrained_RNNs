# pretectum_glms.py

'''
Workflow (model, evaluate, and visualize)
1. Code inputs as in RNN model
2. Try linear-Gaussian
3. Try linear-Poisson
4. Try linear-nonlinear-Poisson
5. Try linear-Poison with spike history
6. Try linear-nonlinear-Poisson with spike history
'''

from datetime import date
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import hankel
from scipy import signal
np.random.seed(2021)

znn_acts = np.genfromtxt('znn_acts.csv', delimiter=',')
stim_course = np.genfromtxt('stimcourse.csv', delimiter=',')
#ds_stim_course =  signal.resample(stim_course, 2000, axis=1)
ds_stim_course = stim_course[:,::2]

print(stim_course.shape)
print(ds_stim_course.shape)
print(np.sum(stim_course))
print(np.sum(ds_stim_course))
print(stim_course)
print(ds_stim_course)

abarr = np.array([[0, 0, 1, 0], [1, 0, 1, 1], [0, 0, 0, 1]])
#ds_abarr = signal.resample(abarr, 2, axis=1)
ds_abarr = abarr[:,::2]
print(abarr)
print(ds_abarr)

integration_window_size = 20     # this becomes the width of the design matrix

print(znn_acts.shape)   # 321 x 2000
print(stim_course.shape)   # 8 x 4001

# X=znn, Y=spikes?, K
# stim_covar = X.T @ X
# spike_trig_ave = X.T @ Y
# k_hat = np.linalg.inv(stim_covar) @ spike_trig_ave

vector = np.linspace(1,100,100)
design_matrix = np.zeros((len(vector) - (integration_window_size - 1), integration_window_size))

for i in range(len(vector) - (integration_window_size - 1)):
    design_matrix[i] = vector[i:i+integration_window_size]

print(design_matrix.shape)
print(design_matrix)

plt.imshow(design_matrix)
plt.show()
# number_of_shifted_rows = length_of_vector - (window_size - 1)
# 81 = 100 - (20 - 1) CHECKS OUT

# FIT N GLMs (ONE FOR EACH NEURON), EACH TAKING IN ALL THE STIMULI

''' ???
# Normalize/clip znn_acts
for i in range(len(znn_acts)):
    znn_acts.iloc[i] = znn_acts.iloc[i] / np.max(np.max(znn_acts))
znn_acts = znn_acts / np.max(np.max(znn_acts)) # normalize activations
print(znn_acts.shape)
for act in znn_acts:
    if act > 0.999: act = 0.999
    if act < -0.999: act = -0.999
    #act = np.min(act, 0.999)
    #act = np.max(act, -0.999)
for i in range(len(znn_acts)):
    znn_acts.iloc[i] = np.amin([znn_acts.iloc[i], 0.999])  # clipping activations above 1
    znn_acts.iloc[i] = np.amax([znn_acts.iloc[i], -0.999])    # clipping activations below -1

# White Noise
WN_amp = 0.0025 # where does this value come from? ("amplitude no more than 1% of baseline according to Andalman")
WN_input = WN_amp * np.random.randn(N, len(model_times)) # used to be inputN
'''