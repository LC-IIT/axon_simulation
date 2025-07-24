# -*- coding: utf-8 -*-
'for generating Fig 4B data - without slow inactivation'
from neuron import h
from neuron.units import ms, mV, um
h.load_file("stdrun.hoc")
import pandas as pd
import time
import math

h.celsius = 20 #set the global temp to 20deg Celsius (default is squid temp)

AIS = h.Section(name = 'AIS')
axon = h.Section(name = 'axon')
axonJT = h.Section(name = 'axonJT')
daughter1JT = h.Section(name = 'daughter1JT')
daughter2JT = h.Section(name = 'daughter2JT')
daughter1=h.Section(name = 'daughter1')
daughter2=h.Section(name = 'daughter2')
axon.connect(AIS(1),0)
axonJT.connect(axon(1),0)
daughter1JT.connect(axonJT(1),0)
daughter2JT.connect(axonJT(1),0)
daughter1.connect(daughter1JT(1),0)
daughter2.connect(daughter2JT(1),0)

AIS.cm = 1  #specific capacitance
axon.cm = 1 
axonJT.cm = 1
daughter1JT.cm = 1
daughter2JT.cm = 1
daughter1.cm = 1
daughter2.cm = 1
AIS.Ra = 100
axon.Ra = 100
axonJT.Ra = 100
daughter1JT.Ra = 100
daughter2JT.Ra = 100
daughter1.Ra = 100
daughter2.Ra = 100

h.hhy.insert(AIS)
h.hhy.insert(axon)
h.hhy.insert(axonJT)
h.hhy.insert(daughter1JT)
h.hhy.insert(daughter2JT)
h.hhy.insert(daughter1)
h.hhy.insert(daughter2)

AIS.gl_hhy = 0.0003
AIS.el_hhy = -54.3
axon.gl_hhy = 0.0003
axon.el_hhy = -54.3
axonJT.gl_hhy = 0.0003
axonJT.el_hhy = -54.3
daughter1JT.gl_hhy = 0.0003
daughter1JT.el_hhy = -54.3
daughter2JT.gl_hhy = 0.0003
daughter2JT.el_hhy = -54.3
daughter1.gl_hhy = 0.0003
daughter1.el_hhy = -54.3
daughter2.gl_hhy = 0.0003
daughter2.el_hhy = -54.3
AIS.tslow_hhy = 1e6
axon.tslow_hhy = 1e6
axonJT.tslow_hhy = 1e6
daughter1JT.tslow_hhy = 1e6
daughter2JT.tslow_hhy = 1e6
daughter1.tslow_hhy = 1e6
daughter2.tslow_hhy = 1e6


# #add slow inactivation
# axon.smin_hhy = 20 #the minimum time constant of s gate
# axon.a2_hhy = 0.2 
# axonJT.smin_hhy = 20 
# axonJT.a2_hhy = 0.2 
# daughter1JT.smin_hhy = 20 
# daughter1JT.a2_hhy = 0.2 
# daughter2JT.smin_hhy = 20 
# daughter2JT.a2_hhy = 0.2 
# daughter1.smin_hhy = 20 
# daughter1.a2_hhy = 0.2 
# daughter2.smin_hhy = 20 
# daughter2.a2_hhy = 0.2 


# #change the slow inactivation to make it stronger
# axon.a0s_fold_hhy = 0.3 #1 
# axonJT.a0s_fold_hhy = 0.3 #1 
# daughter1JT.a0s_fold_hhy = 0.3 #1 
# daughter2JT.a0s_fold_hhy = 0.3 #1 
# daughter1.a0s_fold_hhy = 0.3 #1 
# daughter2.a0s_fold_hhy = 0.3 #1 


AIS.L = 100 * um
axon.L = 2000 * um
axonJT.L = 200 * um
daughter1JT.L = 200 * um
daughter2JT.L = 200 * um
daughter1.L = 1000 * um
daughter2.L = 1000 * um

#select diameters 
AIS.diam = 1 * um
axon.diam = 1 * um
axonJT.diam = 1 * um
daughter1JT.diam = 0.62996 * um 
daughter1.diam = 0.62996 * um  
daughter2JT.diam = 0.62996 * um 
daughter2.diam = 0.62996 * um 

for sec in [AIS, axon, axonJT, daughter1JT, daughter2JT, daughter1, daughter2]:

 	lambdaf = 1e5*math.sqrt(sec.diam/(4*math.pi*100*sec.Ra*sec.cm))
 	sec.nseg = int((sec.L/(0.1*lambdaf)+0.9)/2)*2 +1
 	print(sec.name(),{lambdaf},{sec.nseg})

gnabar_stand = 0.12
gkbar_stand=0.036
AIS.gkbar_hhy = gkbar_stand
axon.gkbar_hhy = gkbar_stand
axonJT.gkbar_hhy = gkbar_stand
daughter1JT.gkbar_hhy = gkbar_stand
daughter2JT.gkbar_hhy = gkbar_stand
daughter1.gkbar_hhy = gkbar_stand
daughter2.gkbar_hhy = gkbar_stand

h.dt = 0.02 * ms
h.steps_per_ms = 1/h.dt
tstart = 100 
tstop = 120000 
nstim = 1200 
isi_stim =100 

# Loop to apply current clamp at different times
iclamps=[]
for i in range(nstim):
    ic=h.IClamp(AIS(0.1))
    ic.delay = tstart + i * isi_stim  # Apply current clamp at different delay times
    ic.dur = 0.1 # Duration (ms)
    ic.amp = 1.01 # nA  
    iclamps.append(ic)

t = h.Vector().record(h._ref_t)
vAIS = h.Vector().record(AIS(0.5)._ref_v)
vstartax = h.Vector().record(axon(0.007)._ref_v)
vd1 = h.Vector().record(daughter1(0.3)._ref_v)
s_vec = h.Vector().record(daughter1(0.3)._ref_s_hhy)
threshold_crossing_vecs = [h.Vector() for seg in AIS]
threshold_crossing_vecsax = [h.Vector() for seg in axon]
threshold_crossing_vecsaxJT = [h.Vector() for seg in axonJT]
threshold_crossing_vecsd1JT = [h.Vector() for seg in daughter1JT]
threshold_crossing_vecsd1 = [h.Vector() for seg in daughter1]
threshold_crossing_vecsd2JT = [h.Vector() for seg in daughter2JT]
threshold_crossing_vecsd2 = [h.Vector() for seg in daughter2]

for l in [threshold_crossing_vecsax, threshold_crossing_vecsaxJT, threshold_crossing_vecsd1JT, threshold_crossing_vecsd1, threshold_crossing_vecsd2JT, threshold_crossing_vecsd2]:

     threshold_crossing_vecs.extend(l)

monitors = [h.NetCon(seg._ref_v,None, sec = AIS) for seg in AIS]
monitorsax = [h.NetCon(seg._ref_v,None, sec = axon) for seg in axon]
monitorsaxJT = [h.NetCon(seg._ref_v,None, sec = axonJT) for seg in axonJT]
monitorsd1JT = [h.NetCon(seg._ref_v,None, sec = daughter1JT) for seg in daughter1JT]
monitorsd1 = [h.NetCon(seg._ref_v,None, sec = daughter1) for seg in daughter1]
monitorsd2JT = [h.NetCon(seg._ref_v,None, sec = daughter2JT) for seg in daughter2JT]
monitorsd2 = [h.NetCon(seg._ref_v,None, sec = daughter2) for seg in daughter2]

for l in [monitorsax, monitorsaxJT, monitorsd1JT, monitorsd1, monitorsd2JT, monitorsd2]:
        
    monitors.extend(l)

for m in range(len(monitors)):
    monitors[m].threshold = -20 * mV #0 * mV
for monitor, vec in zip(monitors,threshold_crossing_vecs):
    monitor.record(vec)

MonitoredLocations = []
for seg in AIS:
    #print (seg, seg.x * AIS.L - (AIS(0.1).x * AIS.L))
    MonitoredLocations.append({seg, seg.x * AIS.L - (AIS(0.1).x * AIS.L)}) 
for seg in axon:
    #print (seg, seg.x * axon.L + (AIS(0.9).x * AIS.L))
    MonitoredLocations.append({seg, seg.x * axon.L + (AIS(0.9).x * AIS.L)}) 
for seg in axonJT:
    #print(seg, seg.x * axonJT.L + axon.L + (AIS(0.9).x * AIS.L))
    MonitoredLocations.append({seg, seg.x * axonJT.L + axon.L + (AIS(0.9).x * AIS.L)})
for seg in daughter1JT:
    #print(seg, seg.x * daughter1JT.L + axonJT.L + axon.L + (AIS(0.9).x * AIS.L))
    MonitoredLocations.append({seg, seg.x * daughter1JT.L + axonJT.L + axon.L + (AIS(0.9).x * AIS.L)})
for seg in daughter1:
    #print (seg, seg.x * daughter1.L + (AIS(0.9).x * AIS.L + axon.L + axonJT.L + daughter1JT.L))
    MonitoredLocations.append({seg, seg.x * daughter1.L + (AIS(0.9).x * AIS.L + axon.L + axonJT.L + daughter1JT.L)}) 
for seg in daughter2JT:
#     #print(seg, seg.x * daughter2JT.L + axonJT.L + axon.L + (AIS(0.9).x * AIS.L))
    MonitoredLocations.append({seg, seg.x * daughter2JT.L + axonJT.L + axon.L + (AIS(0.9).x * AIS.L)})
for seg in daughter2:
#     #print (seg, seg.x * daughter2.L + (AIS(0.9).x * AIS.L + axon.L + axonJT.L + daughter2JT.L))
    MonitoredLocations.append({seg, seg.x * daughter2.L + (AIS(0.9).x * AIS.L + axon.L + axonJT.L + daughter2JT.L)}) 

   
#if we want to reduce the gnabar
for g in [120]: #range(65,125,5): 

    gnabar_low=g/1000
    
    AIS.gnabar_hhy = gnabar_low
    axon.gnabar_hhy =gnabar_low
    axonJT.gnabar_hhy = gnabar_low
    daughter1JT.gnabar_hhy = gnabar_low
    daughter2JT.gnabar_hhy = gnabar_low
    daughter1.gnabar_hhy = gnabar_low
    daughter2.gnabar_hhy = gnabar_low
    
    
    h.cvode_active(True) # optional. but fixed step will probably do one extra time step
    h.cvode.condition_order(2) # optional. but much more accurate event time evaluation.
    
    start= time.perf_counter()
    h.finitialize(-65 * mV)
    h.continuerun(tstop)
    elapsed_time = time.perf_counter()-start
    
    print(f"This took {elapsed_time/60:.2f} minutes")
     
    SpikeTimesToSave = []
    for vec in range(len(threshold_crossing_vecs)):
        #print("threshold_crossing_vecs:", list(threshold_crossing_vecs[vec]))
        SpikeTimesToSave.append(list(threshold_crossing_vecs[vec]))
        
    pd.DataFrame({"Monitored Location (um)":MonitoredLocations, "Spike Times":SpikeTimesToSave}).to_csv(f"SpikeTimes_gnabar_{gnabar_low}_a{axon.diam}_d1_{daughter1.diam}_d2_{daughter2.diam}_4B_noSI.csv", index = False, sep="\t")
    print(f"gnabar = {gnabar_low} ")
    print("\a")
    
    
    
