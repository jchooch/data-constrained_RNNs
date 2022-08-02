# testbed.py

import pandas as pd
import datetime
from datetime import date
import time
import numpy as np

### TIME TEST ###

dtData = 0.5
dtModel = 0.25

print(np.min(100, 1000))
znn_acts = pd.read_csv('znn_acts.csv', header=None) # create pandas dataframe of recorded znn calcium activities

#for i in range(len(znn_acts)):
#    znn_acts.iloc[i] = znn_acts.iloc[i] / np.max(np.max(znn_acts))

znn_acts = znn_acts / np.max(np.max(znn_acts))
print(znn_acts)

'''
print(len(znn_acts[0]))
print(len(znn_acts[1]))

R = np.zeros((321, 1000)) # RNN unit activities *** AM I RIGHT IN USING MODEL TIMES HERE? ***
print('len(R[0]) ', len(R[:,0]))

stim_course = pd.read_csv('stimcourse.csv', header=None)
print(stim_course.shape)
print('Number of ZNN neurons: ', znn_acts.shape[0])

data_start_time, real_start_time, model_start_time = np.zeros(3)   
data_end_time = znn_acts.shape[1] # number of timepoints in data (number of activity measurements for each neuron)
real_end_time = data_end_time * dtData  # number of seconds in data (how long were the data collected in real time in seconds)
model_end_time = real_end_time / dtModel # number of model timepoints (model-time is higher resolution! => more integration steps than real or data timepoints)

print('Number of timepoints in data: ', data_end_time)
print('Number of timepoints in model: ', model_end_time)
print('Number of seconds in real world: ', real_end_time)

real_times = np.arange(real_start_time, real_end_time, 1)
data_times = np.arange(real_start_time, real_end_time, dtData)
model_times = np.arange(real_start_time, real_end_time, dtModel)

print('Real timeseries: ', real_times[0:20], '...')
print('Length of real timeseries: ', len(real_times))
print('Data timeseries: ', data_times[0:20], '...')
print('Length of data timeseries: ', len(data_times))
print('Model timeseries: ', model_times[0:20], '...')
print('Length of model timeseries: ', len(model_times))

print(real_times[-1], data_times[-1], model_times[-1])

lastchi2s = np.zeros(len(real_times))
print(lastchi2s.shape)
'''