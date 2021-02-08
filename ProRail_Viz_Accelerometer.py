# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 16:40:20 2020

@author: deltares
"""

######################################################################MAIN FUNCTIONS##################################################
  
def ProRail_Viz_load_ACC_Data(ACCFiles,fs,sel_chan):

    from io import  StringIO
    import numpy as np
    
    """
    INPUT:
        
        ACCFiles:*.ASC files to be loaded
        fs:sampling frequency of loaded record
        sel_chan: Channels names to be loaded.
        See description below:
            
            locatie	sensornr	meetrichting	meetkanaal
            1	               V1	v	ch 0
            1	               V2	hl	ch 1
            2	               V3	v	ch 2
            2	               V4	hl	ch 3
            5	               V5	v	ch 4
            5	               V6	hl	ch 5
            6	               V7	v	ch 6
            6	               V9	hl	ch 7
            9	              V10	v	ch 8
            9	              V11	hl	ch 9
            3	              tc1	v	ch 10
            3	              tc1	hl	ch 11
            3	              tc1	he	ch 12
            4	              tc3	v	ch 13
            4	              tc3	hl	ch 14
            4	              tc3	he	ch 15
            7	              tc4	v	ch 16
            7	              tc4	hl	ch 17
            7	              tc4	he	ch 18
            8	              tc5	v	ch 19
            8	              tc5	hl	ch 20
            8	              tc5	he	ch 21
    
    
    
    RETURN:
        
        data_acc: Extracted traces from desired channels
        time_acc: time vector for extractec traces.
    
    """
    # channrl names and indixes used to select calibration factors

    chan_index = np.zeros((len(sel_chan)),dtype=int)
    
    for ix in range(0,len(sel_chan)):
    
        chan_name = sel_chan[ix]
        
        if chan_name=='V1_v':
            index = 0
        elif chan_name=='V2_hl':
            index = 1
        elif chan_name=='V3_v':
            index = 2
        elif chan_name=='V4_hl':
            index = 3
        elif chan_name=='V5_v':
            index = 4
        elif chan_name=='V6_hl':
            index = 5
        elif chan_name=='V7_v':
            index = 6
        elif chan_name=='V9_hl':
            index = 7
        elif chan_name=='V10_v':
            index = 8
        elif chan_name=='V11_hl':
            index = 9
        elif chan_name=='tc1_v':
            index = 10
        elif chan_name=='tc1_hl':
            index = 11
        elif chan_name=='tc1_he':
            index = 12
        elif chan_name=='tc3_v':
            index = 13
        elif chan_name=='tc3_hl':
            index = 14
        elif chan_name=='tc3_he':
            index = 15
        elif chan_name=='tc4_v':
            index = 16
        elif chan_name=='tc4_hl':
            index = 17
        elif chan_name=='tc4_he':
            index = 18
        elif chan_name=='tc5_v':
            index = 19
        elif chan_name=='tc5_hl':
            index = 20
        elif chan_name=='tc5_he':
            index = 21
        
        chan_index[ix]=index
    

    # Calibration factors to convert volts to mg

    cal_fact = np.array([-203.874, -206.249, -207.018,-202.265,-203.190,-202.429,-197.336,-194.420,-198.748,-202.573,-88.917,-93.756,-94.011,-90.921,-88.480,-92.473,-203.046,-94.424,-95.007,-199.900,-203.479,-89.206])


    # Loop to read *.asc files

    f = open(ACCFiles, 'r') # 'r' = read
    lines = f.readlines()
    f.close()
    
    rec_len = 30 # in minutes
    
    time = fs*60*rec_len
    
    data_acc=np.zeros((time,len(chan_index)))
    
    time_acc=np.zeros((time))
    
    for k in range(0,time):
        tmp=np.loadtxt(StringIO(lines[k+1]))
        time_acc[k] = tmp[0]/fs
        data_acc[k,:]=tmp[chan_index+1]


    data_acc = ((data_acc/1000)*cal_fact[chan_index]) # convert volts to mm/s2 (acceleration are divided by 1000 to convert to g's)

    return(data_acc*9.81*1000,time_acc)
    
    


def ProRail_filter(data,fs,fmin,fmax):

    """
    INPUT: 
        
        data: ndim - array
        fs: sanmpling frequency in Hz
        fmin: minimum frequency
        fmax: maximum frequency
        detrend: string set as "True by default"
        filter_signal = string if True filtering is applied, if False signal is not filtered
    
    
    RETURN:    
        tr.data: 1-D array 
    """
    m,n=np.shape(data)
    data_filt = np.zeros((m,n))
    for ik in range(0,n):

        from obspy import Trace
        tr = Trace(header={'station': 'DAS', 'channel': 'Z'})
        tr.data = data[:,ik]
        tr.stats.npts=len(data)
        tr.stats.sampling_rate=fs
        tr.stats.station='TestSite'
        tr.stats.channel='Z'
        tr.stats.network = 'DAS'
        tr.detrend('constant')
        tr.taper(0.05, type='cosine')
        tr.filter("bandpass", freqmin=fmin,freqmax=fmax,corners=4, zerophase=True)
            
        data_filt[:,ik] = tr.data

    return (data_filt)


def acc2vel(acc,fs):     
    """
    INPUT: 
        
        acc : ndim-array e.g. time histories acceleration.
        fs = sampling frequency in Hz
        
    RETURN:
        
        Vel: velocity time history
    
    """
    m,n=np.shape(acc)
    vel = np.zeros((m,n))
    from scipy import integrate
    dt = 1/fs
    for iz in range(0,n):
        Dt = integrate.cumtrapz(acc[:,iz],dx=dt,initial=0) # using Trapezoid method    
        vel[:,iz]=Dt    
    return(vel)



def plot_2C_transducers(data_2c,params):

    
    """
    INPUT:
        
        signal1, signal2 : 1-D arrays
        t_start : initial date and time of selected record 
        t_sel: date and time related to train passing
        t_len : 2*t_len is the time window contnaining the passing train.
        
    RETURN:
            
         Plots the selected signals  
    
    """
    
    t_start = params['t_start']
    t_sel = params['t_sel']
    sel_chan = params['sel_chan']
    plot_all = params['plot_all']
    amplitude_type = params['amplitude_type']
    
    
    import matplotlib 
    import datetime as timedate
    import matplotlib.dates as mdates
    from matplotlib.ticker import AutoMinorLocator
    matplotlib.rc('xtick', labelsize=14) 
    matplotlib.rc('ytick', labelsize=14)
    fig, (ax1,ax2) = plt.subplots(nrows=2,figsize=(15,7))
        
    
    signal1 = data_2c[:,0]
    signal2 = data_2c[:,1]
    
    t_end=t_start + timedate.timedelta(minutes = 30)
    xmin, xmax = mdates.date2num([t_start,t_end])
    
    t_start_0 = t_sel - timedate.timedelta(seconds = t_len) # 30 seconds before train time
    t_end_0 = t_sel + timedate.timedelta(seconds = t_len)   # 30 seconds after train time
    xmin0, xmax0 = mdates.date2num([t_start_0,t_end_0])
    t_range = np.linspace(xmin,xmax,len(signal1))
    ax1.plot(t_range,signal1)
    xlocator = mdates.AutoDateLocator()
    ax1.xaxis.set_major_locator(xlocator)
    ax1.xaxis.set_major_formatter(mdates.ConciseDateFormatter(xlocator,formats=['%Y', '%b', '%d', '%H:%M:%S', '%H:%M:%S', '%S.%f']))
    ax1.set_xlabel('Time [HH:MM:SS]',fontsize=14)
    ax1.set_ylabel(amplitude_type,fontsize=14)
    ax1.set_title(sel_chan[0])
    ax1.set_ylim(-amp,amp)
    
    if plot_all=='true':
        ax1.set_xlim(xmin,xmax)
    elif plot_all=='false':
        ax1.set_xlim(xmin0,xmax0)
    else:
        print("Plot_all can be either true or false")
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.grid(which='minor', color='#CCCCCC', linestyle=':')
    ax1.grid()
    
    ax2.plot(t_range,signal2)
    ax2.xaxis.set_major_locator(xlocator)
    ax2.xaxis.set_major_formatter(mdates.ConciseDateFormatter(xlocator,formats=['%Y', '%b', '%d', '%H:%M:%S', '%H:%M:%S', '%S.%f']))
    ax2.set_xlabel('Time [HH:MM:SS]',fontsize=14)
    ax2.set_ylabel(amplitude_type,fontsize=14)
    ax2.grid()
    ax2.set_title(sel_chan[1])
    ax2.set_ylim(-amp,amp)
    
    
    if plot_all=='true':
        ax2.set_xlim(xmin,xmax)
    elif plot_all=='false':
        ax2.set_xlim(xmin0,xmax0)
    else:
        print("Plot_all can be either true or false")
    ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.grid(which='minor', color='#CCCCCC', linestyle=':')
    fig.tight_layout()
    




def plot_3C_transducers(data_2c,params):

    
    """
    INPUT:
        
        signal1, signal2 : 1-D arrays
        t_start : initial date and time of selected record 
        t_sel: date and time related to train passing
        t_len : 2*t_len is the time window contnaining the passing train.
        
    RETURN:
            
         Plots the selected signals  
    
    """
    
    t_start = params['t_start']
    t_sel = params['t_sel']
    sel_chan = params['sel_chan']
    plot_all = params['plot_all']
    amplitude_type = params['amplitude_type']
    
    
    import matplotlib 
    import datetime as timedate
    import matplotlib.dates as mdates
    from matplotlib.ticker import AutoMinorLocator
    matplotlib.rc('xtick', labelsize=14) 
    matplotlib.rc('ytick', labelsize=14)
    fig, (ax1,ax2,ax3) = plt.subplots(nrows=3,figsize=(15,7))
        
    
    signal1 = data_2c[:,0]
    signal2 = data_2c[:,1]
    signal3 = data_2c[:,2]
    
    
    t_end=t_start + timedate.timedelta(minutes = 30)
    xmin, xmax = mdates.date2num([t_start,t_end])
    
    t_start_0 = t_sel - timedate.timedelta(seconds = t_len) # 30 seconds before train time
    t_end_0 = t_sel + timedate.timedelta(seconds = t_len)   # 30 seconds after train time
    xmin0, xmax0 = mdates.date2num([t_start_0,t_end_0])
    t_range = np.linspace(xmin,xmax,len(signal1))
    ax1.plot(t_range,signal1)
    xlocator = mdates.AutoDateLocator()
    ax1.xaxis.set_major_locator(xlocator)
    ax1.xaxis.set_major_formatter(mdates.ConciseDateFormatter(xlocator,formats=['%Y', '%b', '%d', '%H:%M:%S', '%H:%M:%S', '%S.%f']))
    ax1.set_xlabel('Time [HH:MM:SS]',fontsize=14)
    ax1.set_ylabel(amplitude_type,fontsize=14)
    ax1.set_title(sel_chan[0])
    ax1.set_ylim(-amp,amp)
    
    if plot_all=='true':
        ax1.set_xlim(xmin,xmax)
    elif plot_all=='false':
        ax1.set_xlim(xmin0,xmax0)
    else:
        print("Plot_all can be either true or false")
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.grid(which='minor', color='#CCCCCC', linestyle=':')
    ax1.grid()
    
    ax2.plot(t_range,signal2)
    ax2.xaxis.set_major_locator(xlocator)
    ax2.xaxis.set_major_formatter(mdates.ConciseDateFormatter(xlocator,formats=['%Y', '%b', '%d', '%H:%M:%S', '%H:%M:%S', '%S.%f']))
    ax2.set_xlabel('Time [HH:MM:SS]',fontsize=14)
    ax2.set_ylabel(amplitude_type,fontsize=14)
    ax2.grid()
    ax2.set_title(sel_chan[1])
    ax2.set_ylim(-amp,amp)
    
    
    if plot_all=='true':
        ax2.set_xlim(xmin,xmax)
    elif plot_all=='false':
        ax2.set_xlim(xmin0,xmax0)
    else:
        print("Plot_all can be either true or false")
    ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.grid(which='minor', color='#CCCCCC', linestyle=':')
    fig.tight_layout()
    

    ax3.plot(t_range,signal3)
    ax3.xaxis.set_major_locator(xlocator)
    ax3.xaxis.set_major_formatter(mdates.ConciseDateFormatter(xlocator,formats=['%Y', '%b', '%d', '%H:%M:%S', '%H:%M:%S', '%S.%f']))
    ax3.set_xlabel('Time [HH:MM:SS]',fontsize=14)
    ax3.set_ylabel(amplitude_type,fontsize=14)
    ax3.grid()
    ax3.set_title(sel_chan[1])
    ax3.set_ylim(-amp,amp)
    
    
    if plot_all=='true':
        ax3.set_xlim(xmin,xmax)
    elif plot_all=='false':
        ax3.set_xlim(xmin0,xmax0)
    else:
        print("Plot_all can be either true or false")
    ax3.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax3.grid(which='minor', color='#CCCCCC', linestyle=':')
    fig.tight_layout()

#%% #######################################################END OF MAIN FUNCTIONS#########################################################################################3

if __name__ == '__main__':
    path = r'D:\ProRail_Data_Vizualization' # directory should be replaced by the correct user’s directory
    import glob
    import os
    os.chdir(path)    
    import matplotlib.pylab as plt
    import datetime as timedate
    import numpy as np
    ACCFiles = glob.glob('*.asc')
    
    # Selecting components to be analized
    sel_chan = ['tc1_v','tc1_hl','tc1_he']              # Components names to be processed
    fs      = 1000                                      # Sampling frequency   
    fmin    = 1.0                                       # Minimum frequency limit for filtering
    fmax    = 50.0                                      # Maximum frequency limit for filtering
    
    data_2c_raw,time = ProRail_Viz_load_ACC_Data(ACCFiles[0],fs,sel_chan)  # Loading raw traces
    filt_dat = ProRail_filter(data_2c_raw,fs,fmin,fmax)                    # Filtering extracted traces
    allvel = acc2vel(filt_dat,fs)                                          # Time-integrating filtered data to compute particule velocity
    
    # Dictionary of parameters for processing
    
    amp     = 70                                       # User selected amplitude
    t_len   = 30                                        # time in seconds at the right and left side of train time. 
    t_start = timedate.datetime(2020, 11, 9, 14, 10,4)  # initial record time provided in the culemborg0xxx.seq file
    t_sel   = timedate.datetime(2020, 11, 9, 14, 40,4) # time when tran was passing provided by the user
    plot_all= 'true'                                    # if "true" the whole record length will be ploted. If false the a segment of 2*t_len will be plotted containing the selected trin signal
    amplitude_type = 'mm/s2'                             # units of amplitude e.g. acceleration, velocity or displacement
    params={'fs':fs,'fmax':fmin,'fmin':fmax,'t_len':t_len,'t_start':t_start,'t_sel':t_sel,'plot_all':plot_all,\
            'sel_chan':sel_chan,'amplitude_type':amplitude_type}
    
    # ploting a 3-components transducer
    
    plot_3C_transducers(filt_dat,params)


















