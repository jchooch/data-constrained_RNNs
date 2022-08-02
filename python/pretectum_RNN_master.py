# Pretectum RNN Master
"""
Terminology: 
    Biological data/vars labelled with "znn". Biological NEURONS, with activities.
    In silico model data/vars labelled with "rnn". Artificial UNITS, with activities.
    Stimuli: we presented moving gratings to the zebrafish. 
        Gratings were presented to each eye independently, in one of four orientations: right, up, left, or down.
        The timecourse of stimuli presented to the fish is one-hot encoded in a matrix loaded from 'stimcourse.csv'.
Hyperparameters(?): nFree, nFreePre, nRunTot, ...
Inputs(?): N, [t], znn_acts, inputs, WN_amp, J0, data_times, ...
Output(?): R, J, W_input, pVars, lastchi2s, chi2, ...
"""

## Setup Libraries/Packages
from datetime import date
import pandas as pd
import numpy as np
np.random.seed(2021)

## Variable Model Settings
test_mode = False
visualize = True
with_stimuli = True
freeze_weights = False # wtf is the freezing about? figure out
nFree = 5 # epochs to freeze learning before final epoch. what does this mean? figure out
nFreePre = 0  # how many epochs from first one that learning is frozen
num_of_inputs = 8

## Load Data
znn_acts = pd.read_csv('znn_acts.csv', header=None) # create pandas dataframe of recorded znn calcium activities
if with_stimuli == True:
    stim_course = pd.read_csv('stimcourse.csv', header=None) # create pandas dataframe of stimuli metadata (which stimuli when). Row order: LR,LU,LL,LD,RR,RU,RL,RD 

## Model Hyperparameters (exogenous to data)
learn = True # what does this do?
nRunTot = 10 # number of steps (what is a step? an epoch?) [maybe this should be called number of runs] #1500
amp_rgc = 0.025 # scale/amplitude of RGC inputs; Andalman et al. recommended 10x larger than noise amplitude
# something about a control parameter here...
g = 1.25   # Peak synaptic conductance constant; Andalman used 1.2-1.3. Elsewhere 1.0-1.5 recommended. (Where?)
dtData = 0.5 # Frame rate for experimental biological data in seconds
dtModel = 0.25 # Integration step rate in seconds. Andalman et al. used 0.25
tau = 2.5 # Neuron membrane time constants; Andalman et al. used 2.5
P0 = 0.1 # What is this? Whit says "free parameter to learning rule for J". Sussillo & Abbott advise 0.01-1

## Model Parameters (endogenous to data)
N = znn_acts.shape[0] # Number of znn neurons (number of rows in znn_acts); also, set N to be a double or floating point or whatever (why? I don't know)
data_start_time, real_start_time, model_start_time = np.zeros(3)   
data_end_time = znn_acts.shape[1] # number of timepoints in data (number of activity measurements for each neuron)
real_end_time = data_end_time * dtData  # number of seconds in data (how long were the data collected in real time in seconds)
model_end_time = real_end_time / dtModel # number of model timepoints (model-time is higher resolution! => more integration steps than real or data timepoints)

# create vectors of timeseries for real world, data, and model
real_times = np.arange(real_start_time, real_end_time, 1)
data_times = np.arange(real_start_time, real_end_time, dtData)
model_times = np.arange(real_start_time, real_end_time, dtModel)

# Normalize znn_acts
for i in range(len(znn_acts)):
    znn_acts.iloc[i] = znn_acts.iloc[i] / np.max(np.max(znn_acts))
    
znn_acts = znn_acts / np.max(np.max(znn_acts)) # normalize activations
print(znn_acts.shape)
for act in znn_acts:
    if act > 0.999: act = 0.999
    if act < -0.999: act = -0.999
    #act = np.min(act, 0.999)
    #act = np.max(act, -0.999)

'''   
for i in range(len(znn_acts)):
    znn_acts.iloc[i] = np.amin([znn_acts.iloc[i], 0.999])  # clipping activations above 1
    znn_acts.iloc[i] = np.amax([znn_acts.iloc[i], -0.999])    # clipping activations below -1
'''

## Train Model (see <visualize> variable in Settings ^^^)
'''
Things to keep track of: loss/error/chi2, updates, pVar, ...?

def fit(nFreePre=0):
    outputname = 'rnn_output_' + str(date.today())
    pass
'''

lastchi2s = np.zeros(len(data_times))
iModelSample = np.zeros(len(data_times)) # this used to use real_times

for i in range(len(data_times)):
    test = np.zeros(len(model_times))
    test[:] = data_times[i]
    iModelSample[i] = np.argmin(np.abs(test - model_times))

#for i in range(len(real_times)):
#    iModelSample[i] = np.argmin(np.abs(model_times - real_times[i]))  # THIS FINDS THE CLOSEST MODEL TIMEPOINT TO EVERY GIVEN DATA TIMEPOINT

stdevData = np.std(znn_acts.to_numpy().flatten()) # why is this computed all the way up here?!

# White Noise
WN_amp = 0.0025 # where does this value come from? ("amplitude no more than 1% of baseline according to Andalman")
WN_input = WN_amp * np.random.randn(N, len(model_times)) # used to be inputN

# Initialize 
J = g * np.random.randn(N, N) / np.sqrt(N) # initialize weight matrix with Gaussian samples scaled by peak conductance constant

#################################### UP TO HERE ####################################

R = np.empty((N, len(model_times))) # RNN unit activities *** AM I RIGHT IN USING MODEL TIMES HERE? ***
AR = np.empty((N + num_of_inputs, len(model_times))) # augmented activity for training input wts *** WHAT IS THIS AND WHY +8? ***
JR = np.empty((N, 1)) # product of weights and unit activity

# Initialize external input weights
W_input = np.random.randn(N, num_of_inputs) / np.sqrt(N)

if learn == True:
    PJ = P0 * np.identity(N + num_of_inputs) # I don't think this should be called "PJ" unless it's related to J, but I'm not sure if it is

print('Data steps: ', len(data_times))
print('Data dt: ', dtData)
print('Model dt: ', dtModel)
print('Number of real seconds: ', len(real_times))
print('Number of data timesteps: ', len(data_times))
print('Number of model timesteps: ', len(model_times))

chi2 = np.zeros(nRunTot) # used to be (1, nRunTot), but I think this is simpler
pVars = np.zeros(nRunTot) # used to be (1, nRunTot), but I think this is simpler

for nRun in range(nRunTot): # Epoch number, out of number of epochs [THIS SHOULD MAYBE BE CALLED A RUN?]
    H = znn_acts[0] # Initialize RNN activities with ZNN activities at first time point
    R[:,0] = np.tanh(H) # Nonlinearly transformed activities
    tLearn = 0 # param for when to update J matrix
    iLearn = 0 # Used to index znn_acts to subtract from model predicted rates in err.
    # [epoch_LR,epoch_LU,epoch_LL,epoch_LD,epoch_RR,epoch_RU,epoch_RL,epoch_RD] = deal(0); # set epochs of external inputs to 0
    input_epochs = np.zeros(num_of_inputs) # this replaced the line above
                                                        
    for tt in range(len(model_times)): # time steps for each epoch; used to be for tt = 2:len(t)-2
        tLearn = tLearn + dtModel # update for each time step THIS IS JUST USED TO COUNT WHETHER THE TIMESTEP IS ONE WHERE AN ERROR CAN BE COMPUTED OR NOT
        R[:, tt] = np.tanh(H) # nonlinear transformation of activities

        for i in range(num_of_inputs):
            if stim_course.iloc[i,tt] == 1: # conditionals for training external input weights when stim is on
                input_epochs[i] == 1
        
        # generate external inputs at each time step - input on if time vectors are 1 at time tt
        inputs = np.empty((num_of_inputs, N, len(model_times)))
        for i in range(inputs.shape[0]):
            for j in range(N):
                inputs[i,j,tt] = amp_rgc * stim_course.iloc[i,tt] * W_input[j,i]

        '''
        I think the stuff with the inputs and input weights here is being done stupidly (above and below).
        Currently, we are creating a 3D tensor of inputs, then converting the same tensor into the weighted inputs with the same dimensions 
        (but the dimensions should be reduced because we should no longer care about the different kinds of input to a neuron, since they should already have been weighted and summed).
        '''

        # Update RNN unit activities: multiply activities by weights, add noise, add inputs *** I DON'T THINK THESE SHOULD BE CALLED "JR"... THEY SHOULD JUST BE CALLED "R" ***
        JR = J @ R[:,tt]
        JR = JR + WN_input[:,tt]
        weighted_inputs = np.sum(inputs, axis=0) # added this line to sum over inputs of different types, to produce one input per neuron per timepoint (I think this makes sense...?)
        JR = JR + weighted_inputs[:,tt] # Is this adding these in the right way? WHAT IS THIS LINE DOING? changed inputs[:,:,tt] to weighted_inputs[:,tt]

        H = H + (dtModel * (-H + JR)) / tau # model prediction of calcium activities at each time step
        
        if tLearn >= dtData: # model updates weights if tLearn exceeds dtData. Since dtData = 2*dt, this happens every other time step.
            tLearn = 0
            print('iLearn: ', iLearn)
            err = JR - znn_acts.iloc[:, iLearn] # As in Andalman. znn_acts has entries every 0.5 s and JR at every 0.25 s. (used to index to iLearn+1)
            meanerr2 = np.mean(np.power(err, 2)) # what is displayed as chi2 when training. Use it to assess model convergence.
            print('meanerr2: ', meanerr2)
            chi2[nRun] = chi2[nRun] + meanerr2 # nRun is epoch number so accumulates meanerr2 every update to weight matrix. Want this to decreaase over training
            print('chi2[nRun]: ', chi2[nRun])
            lastchi2s[iLearn] = meanerr2 # The last meanerr2 for the last weight update during an epoch.
            print('lastchi2s[iLearn]: ', lastchi2s[iLearn])
            if (learn == True) and (nRun <= nRunTot - nFree) and (nRun > nFreePre): # nFree and nFreePre are weight freezing parameters. learn is set to 1 at start of code.
                # augmented Dyn variable. Each trainable external input is added here.
                AR = R[:,tt]            
                for i in range(num_of_inputs):
                    print('AR shape', AR.shape)
                    print('input_epochs shape', input_epochs.shape)
                    print('input_epochs[i] shape', input_epochs[i].shape)
                    print('input_epochs: ', input_epochs)
                    print('input_epochs[i]: ', input_epochs[i])
                    AR = np.stack((AR, input_epochs[i]))

                # compute estimate of inverse cross correlation matrix of network activities, to scale weight update. See Sussillo & Abbott (2009)
                k = PJ @ AR
                rPr = np.transpose(AR) @ k
                c = 1.0 / (1.0 + rPr)
                PJ = PJ - c @ (k @ np.transpose(k))

                # Updating external input weights if they are on
                for i in range(num_of_inputs):
                    if input_epochs[i] == 1:
                        W_input[:,i] = W_input[:,i] - c @ err @ k[i - num_of_inputs] # why is k indexed in this way?
                
                J = J - c @ err @ np.transpose(k[:N]) # update J by err and proportional to inverse cross correlation network rates
                
            iLearn = iLearn + 1 # Set index of znn_acts (formerly Adata) to time of next frame
            input_epochs = np.zeros(len(input_epochs)) # Set epochs of external inputs to 0

    # Summary of model fit - pVar means percentage of variance explained
    rModelSample = np.zeros((N, len(data_times)))
    counter = 0

    for sample_time in iModelSample:
        rModelSample[:, counter] = R[:, int(sample_time)]
        #rModelSample = np.vstack((rModelSample, sample))
        counter += 1

    znn4shape = znn_acts.to_numpy()
    print('znn_acts shape', znn4shape.shape)
    print('rModelSample shape', rModelSample.shape)
    print('R shape', R.shape)
    pVar = 1 - np.power((np.linalg.norm(znn_acts.to_numpy() - rModelSample) / (np.sqrt(N * len(data_times)) * stdevData)), 2)
    pVars[nRun] = pVar
    print('Run: {} \n pVar: {} \n chi2: {}'.format(nRun, pVar, chi2[nRun]))

# WHAT ARE THESE SUMMARY STATS?
varData = np.var(np.reshape(znn_acts, N * len(data_times)), 1)   # transliterated this directly from the old matlab code. not sure if the dims are all right
chi2 = chi2 / (np.sqrt(N * len(data_times)) * varData)  # chi2 = sum_i (observation_i - expectation_i)**2 / expectation_i
lastchi2s = lastchi2s / (np.sqrt(N * len(data_times)) * varData)