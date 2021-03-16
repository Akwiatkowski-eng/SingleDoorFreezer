# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 13:03:11 2021

@author: Adrian
"""
'''
reading files with results
'''
import numpy as np 
from matplotlib import pyplot as plt

presscondoutput = open("HighPressure.txt","r")
#pressure_cond = presscondoutput.readlines()
pressure_cond = []
count = 0
while True:
    count += 1
    pressure_cond_temp = presscondoutput.readline()
    if not pressure_cond_temp:
        break
    if count > 1 and count <962:  #need to figure out how to read number of rows in a text file or save the results in better way
        pressure_cond.append(float(pressure_cond_temp))
    #print (pressure_cond_temp)
    #print (pressure_cond)
    #print (count)
presscondoutput.close()


pressevapoutput = open("LowPressure.txt","r")
pressure_evap = []
count = 0
while True:
    count += 1
    pressure_evap_temp = pressevapoutput.readline()
    if count > 1 and count <962:
        pressure_evap.append(float(pressure_evap_temp))
    if not pressure_evap_temp:
        break
pressevapoutput.close()


powerconsumoutput = open("PowerConsumption.txt","r")
powerconsumption = []
count = 0
while True:
    count += 1
    powerconsumption_temp = powerconsumoutput.readline()
    if count > 1 and count <962:
        powerconsumption.append(float(powerconsumption_temp))
    if not powerconsumption_temp:
        break
powerconsumoutput.close()


coolingcapaoutput = open("CoolingCapacity.txt","r")
coolingcapacity = []
count = 0
while True:
    count += 1
    coolingcapacity_temp = coolingcapaoutput.readline()
    if count > 1 and count <962:
        coolingcapacity.append(float(coolingcapacity_temp))
    if not coolingcapacity_temp:
        break
coolingcapaoutput.close()


heatinleakoutput = open("Heatinleak.txt","r")
heatinleak = []
count = 0
while True:
    count += 1
    heatinleak_temp = heatinleakoutput.readline()
    if count > 1 and count <962:
        heatinleak.append(float(heatinleak_temp))
    if not heatinleak_temp:
        break
heatinleakoutput.close()


tempoutput = open("Temperature.txt","r")
temp = []
count = 0
while True:
    count += 1
    temp_temp = tempoutput.readline()
    if count > 1 and count <962:
        temp.append(float(temp_temp))
    if not temp_temp:
        break
tempoutput.close()


timeoutput = open("Time.txt","r")
time = []
count = 0
while True:
    count += 1
    time_temp = timeoutput.readline()
    if count > 1 and count <962:
        time.append(float(time_temp))
    if not time_temp:
        break
timeoutput.close()

'''Post processing and graphs creation'''

COP = [i / j for i, j in zip(coolingcapacity, powerconsumption)]

Temp = temp
Temp[:] = [i-273 for i in Temp ]

HighPressure = pressure_cond
HighPressure[:] = [i/100000 for i in HighPressure]

LowPressure = pressure_evap
LowPressure[:] = [i/100000 for i in LowPressure]

plt.title("COP vs Time")
plt.xlabel("Time ,s")
plt.ylabel("COP")
plt.plot(time, COP)
plt.show()

plt.title("Temp inside vs Time")
plt.xlabel("Time ,s")
plt.ylabel("Temp, K")
plt.plot(time, Temp)
plt.show()


plt.title("Cooling capacity vs Time")
plt.xlabel("Time ,s")
plt.ylabel("Q, W")
plt.plot(time, coolingcapacity)
plt.show()

plt.title("High pressure vs Time")
plt.xlabel("Time ,s")
plt.ylabel("Pressure, bar")
plt.plot(time, HighPressure)
plt.show()

plt.title("Low pressure vs Time")
plt.xlabel("Time ,s")
plt.ylabel("Pressure, bar")
plt.plot(time, LowPressure)
plt.show()

'''Fancy plot implementation'''
fig, ax = plt.subplots()

axes = [ax, ax.twinx(), ax.twinx()]
fig.subplots_adjust(right = 0.75)
axes[-1].spines['right'].set_position(('axes', 1.2))
axes[-1].set_frame_on(True)
axes[-1].patch.set_visible(False)

axes[0].plot(time,COP, color = 'Green')
axes[0].set_ylabel('COP', color = 'Green')
axes[0].tick_params(axis = 'y', color = 'Green')
axes[1].plot(time,coolingcapacity, color = 'Red')
axes[1].set_ylabel('Q, W', color = 'Red')
axes[1].tick_params(axis = 'y', color = 'Red')
axes[2].plot(time,Temp, color = 'Blue')
axes[2].set_ylabel('Temperature inside, K', color = 'Blue')
axes[2].tick_params(axis = 'y', color = 'Blue')

axes[0].set_xlabel('Time, s')

plt.draw()
plt.show()