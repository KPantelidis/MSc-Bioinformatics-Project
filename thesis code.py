# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 16:34:45 2020

@author: Konstantinos Pantelidis
"""

import numpy
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes,mark_inset
from matplotlib.patches import Rectangle
import seaborn as sns
import os
import shutil
sns.set()

#open the file and create a dataframe
my_list = []
mz_list = []
int_list = []
data = []
min_del = 0
charge_state = 0
with open('peptides_lcms.txt', mode="r") as f:
    my_list = f.read().splitlines()
    d = {}
    for i in my_list:
        (k, v) = i.split()
        d[k] = v   
        data.append((float(k),float(d[k])))
 

columns = ['m/z', 'intensity']
df0 = pd.DataFrame(data, columns = columns)

#sort by intensity 
df = df0.sort_values(['intensity'], ascending = False)


##normalise the intensity
my_list = df['intensity']
from sklearn.preprocessing import minmax_scale
 
def normalise_list_numpy(list_numpy):
    normalised_list = minmax_scale(list_numpy)
    return normalised_list

df['intensity']  = normalise_list_numpy(my_list)

x = df['m/z']
y = df['intensity']

# exploratory scatter plot
plt.scatter(x,y)
plt.xlabel('m/z')
plt.ylabel('intensity')
plt.show()

##### Noise Cut ######
##                  ##
# can be modified by the user
# when this happens, line 243 needs to be modified too

v = df.iloc[:, 1]
noise_cut = 0.015
df1 = df[~(v < noise_cut)]



# exploratory plot with noise cut
fig = plt.figure(figsize = (13,8))
x1 = df1['m/z']
y1 = df1['intensity']
plt.scatter(x1,y1)
plt.xlabel('m/z')
plt.ylabel('intensity')
plt.show()

## create direcotries for output results
os.mkdir('data_results')
os.mkdir('data_results/charge_state_2')
os.mkdir('data_results/charge_state_3')
os.mkdir('data_results/charge_state_4')


####################################################

## Main Algorithm

def find_charge_state(a, b):
    '''
    Parameters
    ----------
    a : Decimal 
        Low limit of the range for a charge state search.
    b : Decimal
        Upper limit of the range for a charge state search.

    Returns
    Labels the charge state for every identified window and creates a plot
    and a .txt file including the m/z and intensity data of each specific window.
    '''
    total = 0
    new = 0
    peaks = [df1['m/z'].iloc[0],]
#finds the first high intensity with distance ~ a,b 
    for i in df1['m/z']:
        if ((df1['m/z'].iloc[0] - i) > a and (df1['m/z'].iloc[0] - i) < b):
            total += 1 
            new = i
            peaks.append(new)
            print(new)
            print(total)
            break
#trying to find if there's any other distance ~ a,b from the new value 
    for i in df1['m/z']:
        if ((new - i) > a and (new - i) < b):
            total += 1 
            new = i
            peaks.append(new)
            print(new)
            print(total)
            break
             
#finds the first high intensity with distance ~ - a, - b 
    for i in df1['m/z']:
        if ((df1['m/z'].iloc[0] - i) < -a and (df1['m/z'].iloc[0] - i) > -b):
            total += 1 
            new = i
            peaks.append(new)
            print(new)
            print(total)
            break
#trying to find if there's any other distance ~ - a, -b from the new value
    for i in df1['m/z']:
        if ((new - i) < -a and (new - i) > -b):
            total += 1 
            new = i
            peaks.append(new)
            print(new)
            print(total)
            break
    peaks.sort()        

#identifying the charge state
    if total >= 2:
        if (a == 0.47 and b == 0.53):
            charge_state = '+2'
            print("the charge state is +2 \n")
        if (a == 0.30 and b == 0.36):
            charge_state = '+3'
            print("the charge state is +3 \n")
        if (a == 0.22 and b == 0.28):
            charge_state = '+4'
            print("the charge state is +4 \n")
    
        print(peaks)
        
        
        #create a list for the relative intensities of the identified peaks
        peaks_intensity = []
        for i in peaks:
            a = df1['intensity'][df1['m/z']==i].to_list()
            peaks_intensity.append(a[0])
        print(peaks_intensity)
        
        
        points = list(zip(peaks,peaks_intensity))
        
        
        # create the values that will help to the window cut
        global min_del
        min_del = peaks[0] - 0.5
        #print(min_del)
        
        global max_del
        max_del = peaks[-1] + 0.5
        #print(max_del)
         
        #######################################
        ### plot
        #####################
        
        ax = plt.subplot()
        ax.set_xlabel('m/z')
        ax.set_ylabel('intensity') 
        ## create legend
        extra = Rectangle((0, 0), 0, 0, fc="w", fill=False, edgecolor='none', linewidth=0)
        ax.legend([extra, ' '], ('Charge State: ' + charge_state,))
        
        ## set the axes
        plt.xlim(min_del - 1 , max_del + 1 )
        plt.ylim(0, (max(peaks_intensity)+ 15*max(peaks_intensity)/100))
        
        ## plot peaks 
        for pt in points:
            plt.plot( [pt[0],pt[0]], [0,pt[1]] )
        
        ## plot the rest/background noise 
        plt.scatter(x,y)
        
        plt.savefig('sample.png')
        plt.show()
        
        
        ## pre work for saving
        df2 = df0[ (df0['m/z'] >= min_del) & (df0['m/z'] <= max_del) ]
        #print(df2['m/z'].size)
        file_start = df2['m/z'].iloc[0]
        file_end = df2['m/z'].iloc[-1]
        
        ######################################
        ### write tabular data in file and ###
        ### move to the right folder       ###        
        f = open('temp_file.txt' , 'w')
        numpy.savetxt(r'data_results/temp_file.txt', df2.values, fmt=['%3.3f', '%.3e'])
        
        if (charge_state == '+2'):
            i=len(os.listdir('data_results/charge_state_2')) + 1
            shutil.move('data_results/temp_file.txt' , 'data_results/charge_state_2/%s_%s.txt' %(file_start, file_end))
            shutil.move('sample.png' , 'data_results/charge_state_2/plot_%s_%s.png'  %(file_start, file_end))
                       
        if (charge_state == '+3'):
            i=len(os.listdir('data_results/charge_state_3')) + 1
            shutil.move('data_results/temp_file.txt' , 'data_results/charge_state_3/%s_%s.txt'  %(file_start, file_end))
            shutil.move('sample.png' , 'data_results/charge_state_3/plot_%s_%s.png'  %(file_start, file_end))

        if (charge_state == '+4'):
            i=len(os.listdir('data_results/charge_state_4')) + 1
            shutil.move('data_results/temp_file.txt' , 'data_results/charge_state_4/%s_%s.txt'  %(file_start, file_end))
            shutil.move('sample.png' , 'data_results/charge_state_4/plot_%s_%s.png'  %(file_start, file_end))

        
        
    else:   #if total<2
        min_del = peaks[0] - 0.5
        max_del = peaks[-1] + 0.5
 
    
 
######## end of function ###################################    
############################################################


# run the main algorithm/function and cut the identified window
#then run the algorithm again with the new highest peak/point of reference        
vector = list(range(0,1000))        
for i in vector:
    find_charge_state(0.22,0.28)
    find_charge_state(0.47,0.53)
    find_charge_state(0.30,0.36)
    
    df1 = df1[~( (df1['m/z'] >= min_del) & (df1['m/z'] <= max_del) )].sort_values(['intensity'], 
                                                           ascending = False)
    
    # need to be changed if the noise cut above is changed by the user
    if df1['intensity'].iloc[0] < 0.016:
        break


############      REPORTS/ANALYSIS   #######################
############################################################

### create a .txt report for the number of different window charges

## count the number of different window charges
c2 = int(len([name for name in os.listdir('data_results/charge_state_2')])/2)
c3 = int(len([name for name in os.listdir('data_results/charge_state_3')])/2)
c4 = int(len([name for name in os.listdir('data_results/charge_state_4')])/2)

p2 = c2*100/(c2+c3+c4)
p3 = c3*100/(c2+c3+c4)
p4 = c4*100/(c2+c3+c4)

### create a report file and add all the information      
f = open('data_results/report_charge_states.txt' , 'w')
s0 = 'REPORT FOR THE SELECTED FILE \nNoise cut at ' + str(noise_cut*100) + '% \n \n'
s1 = 'There were found ' + str(c2+c3+c4) + ' windows in total \nSpecifically : \n \n'
s2 = 'Charge_State +2: ' + str(c2) + '\n'
s3 = 'Charge_State +3: ' + str(c3) + '\n'
s4 = 'Charge_State +4: ' + str(c4) + '\n \n'
s5 = 'with respectful percentages: \n \n'
s6 = 'Charge_State +2: ' + str(p2) + '\n'
s7 = 'Charge_State +3: ' + str(p3) + '\n'
s8 = 'Charge_State +4: ' + str(p4) + '\n'                         

f.write(s0 + s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8)
f.close()

### pie chart showing the percentages of different window charge states 
fig1, ax1 = plt.subplots()
labels = ['Charge State +2', 'Charge State +3', 'Charge State +4']
percentages = [p2, p3, p4]
colours = ['firebrick', 'goldenrod', 'darkolivegreen']
_, _, autopcts = ax1.pie(percentages, labels=labels, colors=colours, autopct='%1.1f%%',shadow = True)
plt.setp(autopcts, **{'color':'white', 'weight':'bold', 'fontsize':12.5})
ax1.set_title('Percentages of Charge States', fontdict={'fontsize': 17})
plt.savefig('data_results/percentages_charge_states.png')
plt.show()
