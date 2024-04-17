#!/usr/bin/env python

import numpy as np
from obspy.signal.filter import envelope
import matplotlib.pyplot as plt


stalta_tshift_correct = 0.1 #0.1

#------------------------------------------------------------------
TF_run_all = True
# true
input_filename = '/home0/cxd170430/NatureGeo/synthetic/david_lumley/vx_pad.bin'
# false
input_env = '/home0/cxd170430/codes/matlab/vz_test_env.bin'
input_stalta = '/home0/cxd170430/codes/matlab/vz_test_stalta_c.bin'
#------------------------------------------------------------------

#------------------------------------------------------------------
output_bin = True
output_env = '/home0/cxd170430/NatureGeo/synthetic/david_lumley/vx_env.bin'
output_stalta = '/home0/cxd170430/NatureGeo/synthetic/david_lumley/vx_stalta.bin'
#------------------------------------------------------------------

#------------------------------------------------------------------
TF_plot = True
nx = 11218
nt = 3001

# info based on raw input data (unit: sec)
delta = 0.004
stalta_begin = 0.8
stalta_duration = 5.0

# setup STA and LTA time windows
stanpts  = round( 0.1/delta )
ltanpts  = stanpts * 5
#------------------------------------------------------------------


# time info
time_axis = np.linspace(0,(nt-1)*delta,nt)

# params
#                  nt1                    nt2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~raw trace
#    |--------------|                      |--------|
#       LTA length                         STA length
#
# Only [nt1,nt2] part of the raw trace is valid
nt1      = round(stalta_begin/delta)
nt2      = nt1 + round(stalta_duration/delta)
nstalta  = nt2 - nt1 + 1



### Function to calculate STA/LTA ================================
def sta_lta(a, nsta, nlta):
#                   i
# |-------LTA-------|
#      i-ltanpts    |---STA---|
#                    i+stanpts
# length of 'a' includes length of interest and extra LTA plus STA
#
    stalta = np.zeros(len(a)-nlta-nsta)
    for i in range(len(a)-nlta-nsta):
        sta = np.mean(a[i+nlta:i+nlta+nsta])
        lta = np.mean(a[i:i+nlta])
        stalta[i] = sta / lta

    return stalta
### ==============================================================



# read raw data: (nx,nt)
data = np.fromfile(input_filename, dtype=np.float32)
data = data.reshape(nx,nt,order='C')

# main body
if TF_run_all:
    data_env = np.zeros(( nx, nt ))
    data_stalta = np.zeros(( nx, nt ))

    for i in range(nx):
        #print('working on trace #{}/{}'.format(i+1,nx))
        data_env[i,:] = envelope(data[i,:])
        data_stalta[i,nt1:nt2+1] = sta_lta(data_env[i,nt1-ltanpts:nt2+stanpts+1], stanpts, ltanpts)
    
    # correct STA/LTA time shift
    if stalta_tshift_correct != 0:
        itshift = round(abs(stalta_tshift_correct)/delta)
        data_stalta_c = np.zeros((nx, nt))
        if(stalta_tshift_correct > 0):
            data_stalta_c[:,itshift:nt] = data_stalta[:,0:nt-itshift]
        else:
            data_stalta_c[:,0:nt-itshift] = data_stalta[:,itshift:nt]
        
    
    # output binary STA/LTA
    if output_bin:
        data_env.astype(np.float32).tofile(output_env)
        if stalta_tshift_correct == 0:
            data_stalta.astype(np.float32).tofile(output_stalta)
        else:
            data_stalta_c.astype(np.float32).tofile(output_stalta)
        
else:
    data_env = np.fromfile(input_env, dtype=np.float32)
    data_env = data_env.reshape(nx,nt,order='C')
    data_stalta = np.fromfile(input_stalta, dtype=np.float32)
    data_stalta = data_stalta.reshape(nx,nt,order='C')


# plot comparison
if TF_plot:
    itr = 6000
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    ax1.plot(time_axis, data[itr,:], color='b', label='vel wiggle')
    ax1.legend()
    ax1.set_xlim(0, 8)
    ax1.set_ylabel('Amplitude')
    ax1.grid(True)

    ax2.plot(time_axis, data_env[itr,:], color='b', label='envelope')
    ax2.legend()
    ax2.set_ylabel('Amplitude')
    ax2.grid(True)
    
    ax3.plot(time_axis, data_stalta[itr,:], color='b', label='STA/LTA')
    if stalta_tshift_correct != 0 and TF_run_all:
        ax3.plot(time_axis, data_stalta_c[itr,:], color='r', label='corrected STA/LTA')
    ax3.legend()
    ax3.set_ylabel('Amplitude')
    ax3.set_xlabel('Time (s)')
    ax3.grid(True)
    plt.show()

