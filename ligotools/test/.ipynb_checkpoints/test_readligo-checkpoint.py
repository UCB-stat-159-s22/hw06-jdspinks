from ligotools import readligo as rl
import json

def test_ligo1():
    fnjson = "data/BBH_events_v3.json"
    eventname = 'GW150914' 
    events = json.load(open(fnjson,"r"))
    event = events[eventname]
    fn_L1 = event['fn_L1'] 
    strain, time, chan_dict = rl.loaddata('data/'+fn_L1, 'L1')
    assert len(strain) == 131072
    assert len(time) == 131072
    assert len(chan_dict) == 13
    
    
def test_ligo2():
    fnjson = "data/BBH_events_v3.json"
    eventname = 'GW150914' 
    events = json.load(open(fnjson,"r"))
    event = events[eventname]
    fn_H1 = event['fn_H1'] 
    strain, time, chan_dict = rl.loaddata('data/'+fn_H1, 'H1')
    assert len(strain) == 131072
    assert len(time) == 131072
    assert len(chan_dict) == 13
    
def test_ligo3():
    fnjson = "data/BBH_events_v3.json"
    eventname = 'GW150914' 
    events = json.load(open(fnjson,"r"))
    event = events[eventname]
    fn_L1 = event['fn_L1'] 
    strain, time, chan_dict = rl.loaddata('data/'+fn_L1, 'H1')
    DQflag = 'CBC_CAT3'
    segment_list = rl.dq_channel_to_seglist(chan_dict[DQflag])
    assert time[segment_list[0]][0] == 1126259446.0
    
    
def test_ligo4():
    fnjson = "data/BBH_events_v3.json"
    eventname = 'GW150914' 
    events = json.load(open(fnjson,"r"))
    event = events[eventname]
    fn_L1 = event['fn_L1'] 
    strain, time, chan_dict = rl.loaddata('data/'+fn_L1, 'H1')
    segment_list = rl.dq_channel_to_seglist(chan_dict['NO_CBC_HW_INJ'])
    assert time[segment_list[0]][-1] == 1126259477.9997559
    