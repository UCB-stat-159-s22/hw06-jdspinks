import matplotlib as mlb
mlb.use('Agg')
from ligotools import utils
from ligotools import readligo as rl
import matplotlib.mlab as mlab
from scipy.signal import butter, filtfilt
import h5py
import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
import json
import os.path
from os import path

def test_whiten():
    fnjson = "data/BBH_events_v3.json"
    eventname = 'GW150914' 
    events = json.load(open(fnjson,"r"))
    event = events[eventname]
    fs = event['fs']  
    fn_H1 = event['fn_H1'] 
    strain, time, chan_dict = rl.loaddata('data/'+fn_H1, 'H1')
    dt = time[1] - time[0]
    NFFT = 4*fs
    Pxx_H1, freqs = mlab.psd(strain, Fs = fs, NFFT = NFFT)
    psd_H1 = interp1d(freqs, Pxx_H1)
    strain = utils.whiten(strain,psd_H1,dt)
    assert len(strain) == 131072
    
def test_write_wavfile():
    fnjson = "data/BBH_events_v3.json"
    eventname = 'GW150914' 
    events = json.load(open(fnjson,"r"))
    event = events[eventname]
    fs = event['fs'] 
    fband = event['fband'] 
    tevent = event['tevent'] 
    fn_H1 = event['fn_H1'] 
    bb, ab = butter(4, [fband[0]*2./fs, fband[1]*2./fs], btype='band')
    normalization = np.sqrt((fband[1]-fband[0])/(fs/2))
    strain_H1, time_H1, chan_dict = rl.loaddata('data/'+fn_H1, 'H1')
    dt = time_H1[1] - time_H1[0]
    NFFT = 4*fs
    Pxx_H1, freqs = mlab.psd(strain_H1, Fs = fs, NFFT = NFFT)
    psd_H1 = interp1d(freqs, Pxx_H1)
    strain_H1_whiten = utils.whiten(strain_H1,psd_H1,dt)
    strain_H1_whitenbp = filtfilt(bb, ab, strain_H1_whiten) / normalization
    deltat_sound = 2.                     
    indxd = np.where((time_H1 >= tevent-deltat_sound) & (time_H1 < tevent+deltat_sound))
    utils.write_wavfile('audio/'+eventname+"_H1_whitenbp.wav",int(fs), strain_H1_whitenbp[indxd])
    assert path.isfile('audio/GW150914_H1_whitenbp.wav') == True
    
def reqshift():
    fnjson = "data/BBH_events_v3.json"
    eventname = 'GW150914' 
    events = json.load(open(fnjson,"r"))
    event = events[eventname]
    fs = event['fs'] 
    fband = event['fband'] 
    tevent = event['tevent'] 
    fn_H1 = event['fn_H1']
    strain_H1, time_H1, chan_dict = rl.loaddata('data/'+fn_H1, 'H1')
    dt = time_H1[1] - time_H1[0]
    bb, ab = butter(4, [fband[0]*2./fs, fband[1]*2./fs], btype='band')
    normalization = np.sqrt((fband[1]-fband[0])/(fs/2))
    NFFT = 4*fs
    Pxx_H1, freqs = mlab.psd(strain_H1, Fs = fs, NFFT = NFFT)
    psd_H1 = interp1d(freqs, Pxx_H1)
    strain_H1_whiten = whiten(strain_H1,psd_H1,dt)
    strain_H1_whitenbp = filtfilt(bb, ab, strain_H1_whiten) / normalization
    shift = reqshift(strain_H1_whitenbq,fshift=400,sample_rate=4096)
    assert len(strain) == 131072
    assert max(strain) == 1329.5794826006363
    
def test_plot():
    fnjson = "data/BBH_events_v3.json"
    eventname = 'GW150914' 
    events = json.load(open(fnjson,"r"))
    event = events[eventname]  
    tevent = event['tevent'] 
    fn_L1 = event['fn_L1']
    fn_template = event['fn_template'] 
    fs = event['fs'] 
    fband = event['fband']  
    bb, ab = butter(4, [fband[0]*2./fs, fband[1]*2./fs], btype='band')
    normalization = np.sqrt((fband[1]-fband[0])/(fs/2))
    strain_L1, time_L1, chan_dict = rl.loaddata('data/'+fn_L1, 'L1')
    dt = time_L1[1] - time_L1[0]
    data_L1 = strain_L1.copy()
    NFFT = 4*fs
    Pxx_L1, freqs = mlab.psd(strain_L1, Fs = fs, NFFT = NFFT)
    psd_L1 = interp1d(freqs, Pxx_L1)
    strain_L1_whiten = utils.whiten(strain_L1,psd_L1,dt)
    strain_L1_whitenbp = filtfilt(bb, ab, strain_L1_whiten) / normalization

    psd_window = np.blackman(NFFT)
    NOVL = NFFT/2
 
    f_template = h5py.File("data/"+fn_template, "r")
    template_p, template_c = f_template["template"][...]
    template = (template_p + template_c*1.j) 
    dwindow = signal.tukey(template.size, alpha=1./8) 
    template_fft = np.fft.fft(template*dwindow) / fs
    data_psd, freqs_L1 = mlab.psd(data_L1, Fs = fs, NFFT = NFFT, window=psd_window, noverlap=NOVL)
    datafreq = np.fft.fftfreq(template.size)*fs
    df = np.abs(datafreq[1] - datafreq[0])
    power_vec = np.interp(np.abs(datafreq), freqs_L1, data_psd)
    data_fft = np.fft.fft(data_L1*dwindow) / fs
    
    optimal = data_fft * template_fft.conjugate() / power_vec
    optimal_time = 2*np.fft.ifft(optimal)*fs
    sigmasq = 1*(template_fft * template_fft.conjugate() / power_vec).sum() * df
    sigma = np.sqrt(np.abs(sigmasq))
    peaksample = int(data_L1.size / 2)
    SNR_complex = optimal_time/sigma
    SNR = abs(SNR_complex)
    indmax = np.argmax(SNR)
    timemax = time_L1[indmax]
    SNRmax = SNR[indmax]
    
    d_eff = sigma / SNRmax
    horizon = sigma/8
    phase = np.angle(SNR_complex[indmax])
    offset = (indmax-peaksample)
    template_phaseshifted = np.real(template*np.exp(1j*phase))    
    template_rolled = np.roll(template_phaseshifted,offset) / d_eff  
    template_whitened = utils.whiten(template_rolled,interp1d(freqs_L1, data_psd),dt)  
    template_match = filtfilt(bb, ab, template_whitened) / normalization 
    
    pcolor = 'g'
    det = 'L1'
    plottype = 'png'
    eventname = 'GW150914' 
    utils.plot_NSD(SNR, pcolor, det, time_L1, timemax, tevent, plottype, eventname, strain_L1_whitenbp, template_match)
    assert path.isfile('figurs/GW150914_L1_matchtime.png') == True
    assert path.isfile('figurs/GW150914_L1_SNR.png') == True