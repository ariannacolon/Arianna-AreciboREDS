# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 23:33:27 2022

@author: Arianna Colon Cesani
"""
import matplotlib.pyplot as plt
from scipy.io import readsav
from PyEMD import EMD
import numpy as np
import glob

#Create empty lists to store data
F = [] #List of frequencies
I = [] #List of Stokes I (flux density mean/irradiance/intensity)
trends = [] #List of detrends
DF = [] #List of Dataframes

#Specify the path to access all of the sav files from a specific folder
path = "E:\\Documents\\Arecibo REDS\\New_ReducedData_Gliese436"

#Read all of the .sav files inside the specified folder
Files = glob.glob(path + "\*.sav")

#For loop that iterates through all sav files
for file in Files:
    IDLdata=readsav(file)
    IDLdata_all=IDLdata["data"]
    i=IDLdata_all['I'][0] #Column in the file that corresponds to intensity data
    f=IDLdata_all['f'][0] #Column in the file that corresponds to frequency data
    F.append(f) #appending frequency values to the empty list as an array
    I.append(i) #appending intensity values to the empty list as an array

print('Downloaded the data.')

#CLASSIFYING DATA INTO BOXES FOR L-BAND AND C-BAND
#Creating empty lists that will hold frequency and Irradiance arrays from all of our files, divided into boxes 0-6 for L-band and C-band
freqBox0L=[]
IrrBox0L=[]

freqBox1L=[]
IrrBox1L=[]

freqBox2L=[]
IrrBox2L=[]

freqBox3L=[]
IrrBox3L=[]

freqBox4L=[]
IrrBox4L=[]

freqBox5L=[]
IrrBox5L=[]

freqBox6L=[]
IrrBox6L=[]

freqBox0C=[]
IrrBox0C=[]

freqBox1C=[]
IrrBox1C=[]

freqBox2C=[]
IrrBox2C=[]

freqBox3C=[]
IrrBox3C=[]

freqBox4C=[]
IrrBox4C=[]

freqBox5C=[]
IrrBox5C=[]

freqBox6C=[]
IrrBox6C=[]

#Appending data to frequency and irradiance lists corresponding to each box and band
#The user may change initial values for frequency ranges per box and band depending on their data
for i in range(len(F)):
    for j in range(len(IDLdata_all['f'][0])):
        if F[i][j][0] == 1240.8352:  #Box 0 (L-Band) starts at this frequency
            freqBox0L.append(F[i][j])
            IrrBox0L.append(I[i][j])
        if F[i][j][-1] == 1311.84525:  #Box 1 (L-Band) starts at this frequency
            freqBox1L.append(F[i][j])
            IrrBox1L.append(I[i][j])
        if F[i][j][-1] == 1382.84525: #Box 2 (L-Band) starts at this frequency
            freqBox2L.append(F[i][j])
            IrrBox2L.append(I[i][j])
        if F[i][j][-1] == 1453.84525: #Box 3 (L-Band) starts at this frequency
            freqBox3L.append(F[i][j])
            IrrBox3L.append(I[i][j])
        if F[i][j][-1] == 1524.84525: #Box 4 (L-Band) starts at this frequency
            freqBox4L.append(F[i][j])
            IrrBox4L.append(I[i][j])
        if F[i][j][-1] == 1595.84525: #Box 5 (L-Band) starts at this frequency
            freqBox5L.append(F[i][j])
            IrrBox5L.append(I[i][j])
        if F[i][j][-1] == 1666.84525: #Box 6 (L-Band) starts at this frequency
            freqBox6L.append(F[i][j])
            IrrBox6L.append(I[i][j])
    
        if F[i][j][0] == 4065.6704: #Box 0 (C-Band) starts at this frequency
            freqBox0C.append(F[i][j])
            IrrBox0C.append(I[i][j])
        if F[i][j][0] == 4207.6704: #Box 1 (C-Band) starts at this frequency
            freqBox1C.append(F[i][j])
            IrrBox1C.append(I[i][j])
        if F[i][j][-1] == 4349.6905: #Box 2 (C-Band) starts at this frequency
            freqBox2C.append(F[i][j])
            IrrBox2C.append(I[i][j])
        if F[i][j][-1] == 4491.6905: #Box 3 (C-Band) starts at this frequency
            freqBox3C.append(F[i][j])
            IrrBox3C.append(I[i][j])
        if F[i][j][-1] == 4633.6905: #Box 4 (C-Band) starts at this frequency
            freqBox4C.append(F[i][j])
            IrrBox4C.append(I[i][j])
        if F[i][j][-1] == 4775.6905: #Box 5 (C-Band) starts at this frequency
            freqBox5C.append(F[i][j])
            IrrBox5C.append(I[i][j])
        if F[i][j][-1] == 4917.6905: #Box 6 (C-Band) starts at this frequency
            freqBox6C.append(F[i][j])
            IrrBox6C.append(I[i][j])

print('Divided data into L-Band or C-Band and box number.')

#List of frequencies of spectral lines you wish to observe
#This analysis was performed for 10 possible molecular species
Freq_species=[1612.231, 1665.402, 1667.359, 1720.530, 1728.19000, 4091.78200, 4955.20000, 4829.660, 4176.25700, 1653.10010]
SpeciesNames=['Hydroxyl 1612.231 MHz', 'Hydroxyl 1665.402 MHz', 'Hydroxyl 1667.359 MHz','Hydroxyl 1720.530 MHz', 'Methanol 1728.19000 MHz', 'Methanol 4091.78200 MHz', 'Methanol 4955.20000 MHz', 'Formaldehyde 4829.660 MHz', 'Ammonia 4176.25700 MHz', 'Ammonia 1653.10010 MHz']

#Lists of frequencies and irradiances in a selected range (to later plot), close to each element in Freq_species
fs1=[] #frequencies for species 1
Is1=[] #Irradiances for species 1

fs2=[]
Is2=[]

fs3=[]
Is3=[]

fs4=[]
Is4=[]

fs5=[]
Is5=[]

fs6=[]
Is6=[]

fs7=[]
Is7=[]

fs8=[]
Is8=[]

fs9=[]
Is9=[]

fs10=[]
Is10=[]

print('Normalizing and calculating averages.')
q=2 #select a number to add and subtract to each 'spectral line frequency' to plot the desired range (x-axis)
#Species 1
for i in range(len(freqBox5L)): #iterates through each list of frequencies for box 5, L band.
    list1f=[] 
    list1I=[] 
    for j in range(len(freqBox5L[i])): 
        if freqBox5L[i][j] >= (Freq_species[0] - q) and freqBox5L[i][j]<= (Freq_species[0] + q): #frequencies in desired range
            list1f.append(freqBox5L[i][j]) #add each frequency value in chosen range to a list 
            list1I.append(IrrBox5L[i][j]) #add corresponding intensity value in that range to another list
    fs1.append(list1f) #add the formed list of frequencies to a 2D list for this box and band 
    list1I[:] = [val - np.mean(list1I) for val in list1I] #Normalizing the intensities
    Is1.append(list1I) #add the formed list of intensities to a 2D list for this box and band 
Iprom1 = np.sum(Is1, axis=0)/len(Is1) #Average of intensities
Fprom1 = np.sum(fs1, axis=0)/len(fs1) #Average of frequencies
 
#Species 2
for i in range(len(freqBox5L)):
    list2f = []
    list2I = []
    for j in range(len(freqBox5L[i])):
        if freqBox5L[i][j] >= (Freq_species[1] - q) and freqBox5L[i][j]<= (Freq_species[1] + q): 
            list2f.append(freqBox5L[i][j])
            list2I.append(IrrBox5L[i][j])
    fs2.append(list2f)
    list2I[:] = [val - np.mean(list2I) for val in list2I]
    Is2.append(list2I)
Iprom2 = np.sum(Is2, axis = 0)/len(Is2)
Fprom2 = np.sum(fs2, axis = 0)/len(fs2)

#Species 3
for i in range(len(freqBox5L)): 
    list3f = []
    list3I = []
    for j in range(len(freqBox5L[i])):
        if freqBox5L[i][j] >= (Freq_species[2] - q) and freqBox5L[i][j]<= (Freq_species[2] + q): 
            list3f.append(freqBox5L[i][j])
            list3I.append(IrrBox5L[i][j])
    fs3.append(list3f)
    list3I[:] = [val - np.mean(list3I) for val in list3I]
    Is3.append(list3I)
Iprom3 = np.sum(Is3, axis = 0)/len(Is3)
Fprom3 = np.sum(fs3, axis = 0)/len(fs3)

#Species 4
for i in range(len(freqBox6L)):
    list4f = []
    list4I = []
    for j in range(len(freqBox6L[i])):
        if freqBox6L[i][j] >= (Freq_species[3] - q) and freqBox6L[i][j]<= (Freq_species[3] + q): 
            list4f.append(freqBox6L[i][j])
            list4I.append(IrrBox6L[i][j])
    fs4.append(list4f)
    list4I[:] = [val - np.mean(list4I) for val in list4I]
    Is4.append(list4I)
Iprom4 = np.sum(Is4, axis = 0)/len(Is4)
Fprom4 = np.sum(fs4, axis = 0)/len(fs4)

#Species 5
for i in range(len(freqBox6L)):
    list5f = []
    list5I = []
    for j in range(len(freqBox6L[i])):
        if freqBox6L[i][j] >= (Freq_species[4] - q) and freqBox6L[i][j]<= (Freq_species[4] + q): 
            list5f.append(freqBox6L[i][j])
            list5I.append(IrrBox6L[i][j])
    fs5.append(list5f)
    list5I[:] = [val - np.mean(list5I) for val in list5I]
    Is5.append(list5I)
Iprom5 = np.sum(Is5, axis = 0)/len(Is5)
Fprom5 = np.sum(fs5, axis = 0)/len(fs5)

#Species 6
for i in range(len(freqBox0C)): 
    list6f = []
    list6I = []
    for j in range(len(freqBox0C[i])):
        if freqBox0C[i][j] >= (Freq_species[5] - q) and freqBox0C[i][j]<= (Freq_species[5] + q): 
            list6f.append(freqBox0C[i][j])
            list6I.append(IrrBox0C[i][j])
    fs6.append(list6f)
    list6I[:] = [val - np.mean(list6I) for val in list6I]
    Is6.append(list6I)
Iprom6 = np.sum(Is6, axis = 0)/len(Is6)
Fprom6 = np.sum(fs6, axis = 0)/len(fs6)

#Species 7
for i in range(len(freqBox6C)): 
    list7f = []
    list7I = []
    for j in range(len(freqBox6C[i])):
        if freqBox6C[i][j] >= (Freq_species[6] - q) and freqBox6C[i][j]<= (Freq_species[6] + q): 
            list7f.append(freqBox6C[i][j])
            list7I.append(IrrBox6C[i][j])
    fs7.append(list7f)
    list7I[:] = [val - np.mean(list7I) for val in list7I]
    Is7.append(list7I)
Iprom7 = np.sum(Is7, axis = 0)/len(Is7)
Fprom7 = np.sum(fs7, axis = 0)/len(fs7)

#Species 8
for i in range(len(freqBox5C)): 
    list8f = []
    list8I = []
    for j in range(len(freqBox5C[i])):
        if freqBox5C[i][j] >= (Freq_species[7] - q) and freqBox5C[i][j]<= (Freq_species[7] + q): 
            list8f.append(freqBox5C[i][j])
            list8I.append(IrrBox5C[i][j])
    fs8.append(list8f)
    list8I[:] = [val - np.mean(list8I) for val in list8I]
    Is8.append(list8I)
Iprom8 = np.sum(Is8, axis = 0)/len(Is8)
Fprom8 = np.sum(fs8, axis = 0)/len(fs8)

#Species 9
for i in range(len(freqBox0C)): 
    list9f = []
    list9I = []
    for j in range(len(freqBox0C[i])):
        if freqBox0C[i][j] >= (Freq_species[8] - q) and freqBox0C[i][j]<= (Freq_species[8] + q): 
            list9f.append(freqBox0C[i][j])
            list9I.append(IrrBox0C[i][j])
    fs9.append(list9f)
    list9I[:] = [val - np.mean(list9I) for val in list9I]
    Is9.append(list9I)
Iprom9 = np.sum(Is9, axis = 0)/len(Is9)
Fprom9 = np.sum(fs9, axis = 0)/len(fs9)

#Species 10
for i in range(len(freqBox5L)):
    list10f = []
    list10I = []
    for j in range(len(freqBox5L[i])):
        if freqBox5L[i][j] >= (Freq_species[9] - q) and freqBox5L[i][j]<= (Freq_species[9] + q): 
            list10f.append(freqBox5L[i][j])
            list10I.append(IrrBox5L[i][j])
    fs10.append(list10f)
    list10I[:] = [val - np.mean(list10I) for val in list10I]
    Is10.append(list10I)
Iprom10 = np.sum(Is10, axis = 0)/len(Is10)
Fprom10 = np.sum(fs10, axis = 0)/len(fs10)

'''
#PLOTTING
print('Starting raw-data and de-trended curve plotting process.')

#Defining figure and subplot properties
rows = 5
columns = 2   
c = 1  #plot counter
fig = plt.figure(figsize = (21,21))

plt.subplots_adjust(left = 0.1,
                    bottom = 0.2, 
                    right = 0.558, 
                    top = 0.7, 
                    wspace = 0.2, 
                    hspace = 0.5)

#Plotting subplots (one per molecular species)
for k in range (len(Freq_species)): #k represents Freq_species list index (for SpeciesNames list as well)
    if k == 0 :
        Fprom = Fprom1 #Assigning frequency average list to the first species (Box 5, L Band)
        Iprom = Iprom1 #Assigning intensity average list to the first species (Box 5, L Band)
    if k == 1 :
        Fprom = Fprom2
        Iprom = Iprom2        
    if k == 2 :
        Fprom = Fprom3
        Iprom = Iprom3        
    if k == 3 : 
        Fprom = Fprom4
        Iprom = Iprom4
    if k == 4 : 
        Fprom = Fprom5
        Iprom = Iprom5
    if k == 5 :
        Fprom = Fprom6
        Iprom = Iprom6
    if k == 6 : 
        Fprom = Fprom7
        Iprom = Iprom7
    if k == 7 :    
        Fprom = Fprom8
        Iprom = Iprom8
    if k == 8 :
        Fprom = Fprom9
        Iprom = Iprom9
    if k == 9 :
        Fprom = Fprom10
        Iprom = Iprom10        
        
    plt.subplot(rows, columns, c)
    plt.title('{}'.format(SpeciesNames[k]))
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('I_mean (Jy)')
    plt.plot(Fprom, Iprom)
    
    flatten_lc1, trend_lc1 = flatten(Fprom, Iprom, window_length=int(20), return_trend=True, method='savgol') #detrending
    plt.plot(Fprom, trend_lc1, linewidth=2, color='orange') #plot of detrended curve
    
    
    c = c + 1

#plt.savefig('SpectroscopyAnalysis.png', dpi=170) #Saving image
plt.show()
'''
'''
print('Starting flattened-curve plotting process.')

#Defining figure and subplot properties
rows = 5
columns = 2   
c = 1  #plot counter
fig = plt.figure(figsize = (21,21))

plt.subplots_adjust(left = 0.1,
                    bottom = 0.2, 
                    right = 0.558, 
                    top = 0.7, 
                    wspace = 0.2, 
                    hspace = 0.5)

#Plotting subplots (one per molecular species)
for k in range (len(Freq_species)): #k represents Freq_species list index (for SpeciesNames list as well)
    if k == 0 :
        Fprom = Fprom1 #Assigning frequency average list to the first species (Box 5, L Band)
        Iprom = Iprom1 #Assigning intensity average list to the first species (Box 5, L Band)
    if k == 1 :
        Fprom = Fprom2
        Iprom = Iprom2        
    if k == 2 :
        Fprom = Fprom3
        Iprom = Iprom3        
    if k == 3 : 
        Fprom = Fprom4
        Iprom = Iprom4
    if k == 4 : 
        Fprom = Fprom5
        Iprom = Iprom5
    if k == 5 :
        Fprom = Fprom6
        Iprom = Iprom6
    if k == 6 : 
        Fprom = Fprom7
        Iprom = Iprom7
    if k == 7 :    
        Fprom = Fprom8
        Iprom = Iprom8
    if k == 8 :
        Fprom = Fprom9
        Iprom = Iprom9
    if k == 9 :
        Fprom = Fprom10
        Iprom = Iprom10        
        
    plt.subplot(rows, columns, c)
    plt.title('{}'.format(SpeciesNames[k]))
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('I_mean (Jy)')
    plt.plot(Fprom, Iprom)
    
    flatten_lc1, trend_lc1 = flatten(Fprom, Iprom, window_length=int(30), return_trend=True, method='savgol') #detrending
    plt.plot(Fprom, flatten_lc1, linewidth=2, color='green') #flattened curve    
    
    c = c + 1

#plt.savefig('SpectroscopyAnalysis.png', dpi=170) #Saving image
plt.show()

'''
print('Starting de-trended curve plotting process.')



#Testing for ONE species using EMD 
#Plotting subplots (one per molecular species)
#For species #6 (methanol)
#Empirical Mode Decomposition (EMD) is an iterative procedure which decomposes signal into a set of oscillatory components, called Intrisic Mode Functions (IMFs). 

# Define signal
Fprom = Fprom7 #t
Iprom = Iprom7 #s

# Execute EMD on signal
IMF = EMD().emd(Iprom,Fprom)
N = IMF.shape[0]+1

# Plot results
fig = plt.figure(figsize = (12,20))
plt.subplot(N,1,1)
plt.plot(Fprom, Iprom, 'r')
plt.title("Methanol 4955.2000 MHz")

for n, imf in enumerate(IMF):
    plt.subplot(N,1,n+2)
    plt.plot(Fprom, imf, 'g')
    plt.title("IMF "+str(n+1))
    
newIMF1=[]
newIMF1.append(IMF[-1])
newIMF1.append(IMF[-2])
newIMF1.append(IMF[-3])
det_sum1=np.sum(newIMF1, axis = 0)

plt.figure(figsize = (14,2))
plt.plot(Fprom, Iprom, 'r',label='data')
plt.plot(Fprom,det_sum1,label='detrended data 1',color='blue')
plt.title('Detrended flux for Methanol 4955.2000 MHz')
plt.xlabel('Frequency (MHz)')
plt.ylabel('I_mean (Jy)')
plt.legend()

   
'''
#Detrending manually (subtracting IMFS from the original data; excluded first and last IMF)
print(IMF)
newIMF1=[]
for i in range (1, len(IMF)-1):
    newIMF1.append(IMF[i])
det_sum1=np.sum(newIMF1, axis = 0)# sum of all IMFs except first and last one

#Detrending manually (subtracting IMFS from the original data; excluded first and last IMF)
newIMF2=[]
for i in range (1, len(IMF)-2): #subtracting last two IMFS
    newIMF2.append(IMF[i])
det_sum2=np.sum(newIMF2, axis = 0)# sum of all IMFs except first and last one


plt.figure(figsize = (14,2))
#detrendedI=np.subtract(Iprom,det_sum) #subtracting IMFs from raw data
plt.plot(Fprom,det_sum1,label='detrended data 1',color='blue')
plt.plot(Fprom,det_sum2, label='detrended data 2', color='green')
plt.title('Detrended flux for Methanol 4955.2000 MHz')
plt.xlabel('Frequency (MHz)')
plt.ylabel('I_mean (Jy)')
plt.legend()

# Plot results with detrending
fig = plt.figure(figsize = (12,20))
plt.subplot(N,1,1)
plt.plot(Fprom, Iprom, 'r',label='data')
plt.title("Methanol 4955.2000 MHz")
plt.plot(Fprom,det_sum1, label='detrended data 1', color='blue')
plt.plot(Fprom,det_sum2, label='detrended data 2', color='green')
plt.legend()
'''