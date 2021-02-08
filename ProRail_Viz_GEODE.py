# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 10:02:29 2021

@author: obandohe
"""

def ProRail_Viz_load_GEODE_Data(GEODEFiles):
    
    
    """
    
    Function for reading seg2 files contaning either combined 4.5 Hz and 1.0 Hz geophones or 4.5 Hz geophones only.
    
    
    INPUT: 
        
       GEODEFiles: Name of seg2 file to be read.
       
    RETURN:
        
        Shotgather_1Hz      : Traces of 12, 1Hz geophones
        offset_1Hz          : Coordinates [m] of 1Hz geophones along the 90 m test site. 
        Shotgather_4_5Hz    : Traces of 48, 4.5Hz geophones
        offset_4_5Hz        : Coordinates [m] of 4.5Hz geophones along the 90 m test site. 
        timestamp           : Time stamp when files was created
        time_range          : Time vector of recorded signals in seconds
    
    
    """
    
    
    from obspy import read  # Obspy modules utilized to read seg2 files
    import numpy as np
    # Search for all tdms files contained in the directory.    
    st=read(GEODEFiles)   # Obspy function to open SEG2 data format
    size=np.shape(st)
    ShotGather_ref=np.zeros((size[1], size[0]))
    for i in range(0, size[0]):
        tr=np.reshape(st[i],(1,size[1]))
        ShotGather_ref[:, i] = tr
    Shotgather_norm=ShotGather_ref[0:15000,:]
    fs = st[0].stats.sampling_rate
    timestamp = st[0].stats.starttime.datetime
    
    # Traces for 1 Hz geophones----------------------------------------------
    index = np.array([0,5,10,15,20,47,45,42,40,35,30,25]) # Indixed according to the field geometry
    index2 = index+48
    All_traces_1Hz = Shotgather_norm
    offset_1Hz = np.array([0,10,20,30,40,45,50,55,60,70,80,90])
    if size[0]>48:
        Shotgather_1Hz = (All_traces_1Hz[:,index2])
    else:    
        Shotgather_1Hz = []
        
    # Traces for 4.5 Hz geophones--------------------------------------------
    Shotgather_4_5Hz=np.zeros((15000,48))
    Shotgather_4_5Hz[:,0:24]=Shotgather_norm[:,24:48]
    Shotgather_4_5Hz[:,24:48]=np.fliplr(Shotgather_norm[:,0:24])
    Shotgather_4_5Hz = np.fliplr(Shotgather_4_5Hz)
    offset_4_5Hz = np.arange(0,48)+23
    time_range = np.arange(0,len(Shotgather_4_5Hz[:,0]))/fs
    
    
    return(Shotgather_1Hz,offset_1Hz,Shotgather_4_5Hz,offset_4_5Hz,timestamp,time_range)
    
  
def ProRail_filter(data,fs,fmin,fmax):

    import numpy as np
    
    """
    INPUT: 
        
        data: ndim - array
        fs: sanmpling frequency in Hz
        fmin: minimum frequency in Hz
        fmax: maximum frequency in Hz
        detrend: string set as "True by default"
        filter_signal = string if True filtering is applied, if False signal is not filtered
    
    
    RETURN:    
        tr.data: 1-D array filtered between [fmin,fmax] range
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
     
    
    
if __name__ == '__main__':
    
    import os
    from segypy import wiggle
    import matplotlib.pylab as plt
    import glob
    import matplotlib 
    matplotlib.rc('xtick', labelsize=14) 
    matplotlib.rc('ytick', labelsize=14)
    path = r'D:\ProRail_Scripts\ProRail_Python_Scripts\sampledata' # directory should be replaced by the correct user’s directory
    os.chdir(path)   
    Allfiles=glob.glob('*.dat')             # Selecting *.dat (seg2) file to be read.
    GEODEFiles = Allfiles[0]
    
    
    #--------------------------------Loading shotgather from 4.5 Hz and 1.0 Hz geophones---------
    
    Shotgather_1Hz,offset_1Hz,Shotgather_4_5Hz,offset_4_5Hz,timestamp,time_range = ProRail_Viz_load_GEODE_Data(GEODEFiles)
    
    fs = 1000           # sanmpling frequency in Hz
    fmin = 1            # minimum frequency in Hz
    fmax = 100          # maximum frequency in Hz
    
    
    # --------------------------------Filtering shotgathers--------------------------------------
    
    Shotgather_1Hz_filtered = ProRail_filter(Shotgather_1Hz,fs,fmin,fmax)
    Shotgather_4_5Hz_filtered = ProRail_filter(Shotgather_4_5Hz,fs,fmin,fmax)
    
    # -------------------------------Displaying shotgathers--------------------------------------
    
    fig, ax = plt.subplots(figsize=(15,8))  # Setting up figure size for ploting 4.5 Hz and 1.0 Hz geophones.
    
    plt.subplot(121)
    wiggle(Shotgather_4_5Hz_filtered,x=offset_4_5Hz,t=time_range)
    plt.title(str(timestamp) + '-- 4.5 Hz geophones',fontsize=14)
    plt.xlabel('Distance [m]',fontsize=14)
    plt.ylabel('Time[s]',fontsize=14)
    
    plt.subplot(122)
    wiggle(Shotgather_1Hz_filtered,x=offset_1Hz,t=time_range)
    plt.title(str(timestamp) + '-- 1.0 Hz geophones',fontsize=14)
    plt.xlabel('Distance [m]',fontsize=14)
    plt.ylabel('Time[s]',fontsize=14)

#%%









    