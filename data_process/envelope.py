import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import hilbert



# input params===============================================================================================
nt = 2001
dt = 0.004
ng = 11218

file_x = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_48.601_67.232_32bit.bin'
file_z = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_48.601_67.232_32bit.bin'
file_x_out = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_48.601_67.232_32bit.bin'
file_z_out = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_48.601_67.232_32bit.bin'

plot_vel_data = False
plot_env_data = False

write_env_data = True
#============================================================================================================



vx = np.fromfile(file_x, dtype=np.float32)
vx = vx.reshape(nt,ng,order='F')

vz = np.fromfile(file_z, dtype=np.float32)
vz = vz.reshape(nt,ng,order='F')

# because written in ‘C’ order
env_x = np.zeros((ng,nt))
env_z = np.zeros((ng,nt))

for i in range(0,ng):
    hilb_vx = hilbert(vx[:,i])
    env_x[i,:] = (vx[:,i]**2 + hilb_vx**2) ** 0.5  # vx[:,i]
    hilb_vz = hilbert(vz[:,i])
    env_z[i,:] = (vz[:,i]**2 + hilb_vz**2) ** 0.5  # vz[:,i]


# always written in ‘C’ order
if write_env_data:
    env_x.astype(np.float32).tofile(file_x_out)
    env_z.astype(np.float32).tofile(file_z_out)


if plot_vel_data:
    t = np.linspace(0,(nt-1)*dt, nt)
    x = np.linspace(1,ng,ng)
    fig = plt.figure(figsize=[12,8])
    ax1 = fig.add_axes([0.06,0.1,0.4,0.75])
    ax2 = fig.add_axes([0.53,0.1,0.4,0.75], sharey=ax1)
    ax1.imshow(vz, cmap="seismic", vmin=-5e-15, vmax=5e-15, extent=(x[0],x[-1],t[-1],t[0]), aspect='auto')
    ax2.imshow(vx, cmap="seismic", vmin=-5e-15, vmax=5e-15, extent=(x[0],x[-1],t[-1],t[0]), aspect='auto')
    ax1.set_title('Z')
    ax2.set_title('R')
    ax1.set_xlabel('Traces in offset order')
    ax1.set_ylabel('Time (s)')
    ax2.set_xlabel('Traces in offset order')
    plt.show()

if plot_env_data:
    t = np.linspace(0,(nt-1)*dt, nt)
    x = np.linspace(1,ng,ng)
    fig = plt.figure(figsize=[12,8])
    #plt.plot(t,env_x[10,:])
    ax1 = fig.add_axes([0.06,0.1,0.4,0.75])
    ax2 = fig.add_axes([0.53,0.1,0.4,0.75], sharey=ax1)
    ax1.imshow(env_z.transpose(), cmap="gray", vmin=0, vmax=3e-14, extent=(x[0],x[-1],t[-1],t[0]), aspect='auto')
    ax2.imshow(env_x.transpose(), cmap="gray", vmin=0, vmax=3e-14, extent=(x[0],x[-1],t[-1],t[0]), aspect='auto')
    ax1.set_title('Z')
    ax2.set_title('R')
    ax1.set_xlabel('Traces in offset order')
    ax1.set_ylabel('Time (s)')
    ax2.set_xlabel('Traces in offset order')
    plt.show()


print('input_x: {}'.format(file_x))
print('input_z: {}'.format(file_z))
print('ouput_x: {}'.format(file_x_out))
print('ouput_z: {}'.format(file_z_out))







