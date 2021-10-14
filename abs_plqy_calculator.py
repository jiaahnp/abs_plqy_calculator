"""
Absolute PLQY calculator which integrates the excitation and emission portion
of the blank and sample spectra and calculates a PLQY value with error bounds.

@author: Jia-Ahn Pan
"""

import tkinter as tk
import tkinter.filedialog as fd
import matplotlib as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,NavigationToolbar2Tk)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import os
from textwrap import wrap


root = tk.Tk()
root.wm_title("Absolute PLQY Calculator with Error Bounds")

#------Preset Parameters---------
lbound_low = 400
lbound_high = 410
fbound_low = 450  
fbound_high = 570

plt.rcParams['savefig.dpi'] = 500 #resolution of sav

# ------------First frame-Import Files and Display Names-----------------
frame1 = tk.Frame(root, bd="5")

label1 = tk.Label(frame1, text = "Please select blank(s)")
label2 = tk.Label(frame1, text = "Please select sample(s)")

label1.grid(row=0,column=1)
label2.grid(row=1,column=1)

root.blank_files = "Null"
root.sample_files = "Null"

def get_blankfiles():
    root.blank_files = fd.askopenfilenames(title="Select file(s) for blank", filetypes = (("Text files", "*.txt"),("All files","*.*"))) 
    # Display all file names
    blank_string = ""
    for blank_file in root.blank_files: 
        blank_string += str(blank_file) + " \n "
    blank_string = blank_string[:-3] #remove last space
    blank_string = blank_string[0:1000] # max charactesr
    label1.config(text = blank_string)

def get_samplefiles():
    root.sample_files = fd.askopenfilenames(title="Select file(s) for sample", filetypes = (("Text files", "*.txt"),("All files","*.*"))) 
    # Display all file names
    sample_string = ""
    for sample_file in root.sample_files: 
        sample_string += str(sample_file) + " \n "
    sample_string = sample_string[:-3] #remove last space
    sample_string = sample_string[0:1000] # max charactesr
    label2.config(text = sample_string)

button1 = tk.Button(frame1, text="   Blank(s)  ", command=get_blankfiles).grid(row=0, column=0)
button2 = tk.Button(frame1, text="Samples(s)", command=get_samplefiles).grid(row=1, column=0)

frame1.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# ------------Frame 2 - Input calulation bounds-----------------
frame2 = tk.Frame(root, bd="5")

label3 = tk.Label(frame2, text = "Laser low (nm):").grid(row=0,column=0)
entry3 = tk.Entry(frame2, width="5")
entry3.grid (row=0,column=1, padx = "10")
entry3.insert(0, str(lbound_low))

label4 = tk.Label(frame2, text = "Laser high (nm):").grid(row=0,column=2)
entry4 = tk.Entry(frame2,width="5")
entry4.grid (row=0,column=3, padx = "10")
entry4.insert(0, str(lbound_high))

label5 = tk.Label(frame2, text = "PL low (nm):").grid(row=1,column=0)
entry5 = tk.Entry(frame2, width="5")
entry5.grid (row=1,column=1, padx = "10")
entry5.insert(0, str(fbound_low))

label6 = tk.Label(frame2, text = "PL high (nm):").grid(row=1,column=2)
entry6 = tk.Entry(frame2,width="5")
entry6.grid (row=1,column=3, padx = "10")
entry6.insert(0, str(fbound_high))

frame2.pack(side=tk.TOP,fill=tk.BOTH, expand=1)


# ------------Frame 3 - Graph-----------------
fig = Figure(figsize=(6,4), dpi=100, tight_layout=True)

canvas = FigureCanvasTkAgg(fig, master=root)
toolbar = NavigationToolbar2Tk(canvas, root)
toolbar.update()

canvas.get_tk_widget().pack(side=tk.TOP)#, fill=tk.BOTH, expand=0)
    
def add_plot():
    fig.clf()
    canvas.draw()
    canvas.get_default_filename = lambda: (
             os.path.basename(root.sample_files[0])[0:-4] + '.png')
    
    a = fig.add_subplot(111)
    ax = a.axes
    
    lbound_low = int(entry3.get())
    lbound_high = int(entry4.get())
    fbound_low = int(entry5.get())
    fbound_high = int(entry6.get())
    
    # extract only first data
    blank_data_first = np.genfromtxt(root.blank_files[0], skip_header=16, delimiter='\t')
    sample_data_first = np.genfromtxt(root.sample_files[0], skip_header=16, delimiter='\t')
        
    # extract all blank data
    blank_laser_signal = []
    blank_fluor_signal = []    
    for dataname in root.blank_files:
        data = np.genfromtxt(dataname, skip_header=16, delimiter='\t')
        data[:,1] = data[:,1]*np.gradient(data[:,0]) # replace counts with counts*wavelength spacing
        
        laser_signal, fluor_signal = counts_region(data, lbound_low, lbound_high,
                        fbound_low, fbound_high)
        blank_laser_signal.append(laser_signal)
        blank_fluor_signal.append(fluor_signal)
        
    # extract all sample data
    sample_laser_signal = []
    sample_fluor_signal = []    
    for dataname in root.sample_files:
        data = np.genfromtxt(dataname, skip_header=16, delimiter='\t')
        data[:,1] = data[:,1]*np.gradient(data[:,0]) # replace counts with counts*wavelength spacing
        laser_signal, fluor_signal = counts_region(data, lbound_low, lbound_high,
                        fbound_low, fbound_high)
        sample_laser_signal.append(laser_signal)
        sample_fluor_signal.append(fluor_signal)
    
    # Data and stats analysis
    #   Averages
    blank_laser_average = np.average(blank_laser_signal)
    blank_fluor_average = np.average(blank_fluor_signal)
    sample_laser_average = np.average(sample_laser_signal)
    sample_fluor_average = np.average(sample_fluor_signal)    
    # std deviations
    blank_laser_std = np.std(blank_laser_signal)
    blank_fluor_std = np.std(blank_fluor_signal)
    sample_laser_std = np.std(sample_laser_signal)
    sample_fluor_std = np.std(sample_fluor_signal)
    #variances
    blank_laser_var = blank_laser_std**2
    blank_fluor_var = blank_fluor_std**2
    sample_laser_var = sample_laser_std**2
    sample_fluor_var =  sample_fluor_std **2
    
    # Quantum yield average and standard deviation
    qy_average = ((sample_fluor_average-blank_fluor_average)
                 /(blank_laser_average-sample_laser_average))
    qy_std_a = ((sample_fluor_var + blank_fluor_var)
                /((sample_fluor_average-blank_fluor_average))**2)
    qy_std_b = ((blank_laser_var + sample_laser_var)
                /(blank_laser_average-sample_laser_average)**2)
    qy_std = qy_average*((qy_std_a + qy_std_b)**(1/2))
    
    # Display data/stats analysis
    
    a.text(0.3, 0.05, r"Quantum Yield ($\pm \sigma$): "+str(round(qy_average,3))
           + r"$\pm$" + str(round(qy_std,3))
           +"\n\n"
           +r"Blank Laser Signal ($\pm \sigma$): "+str(round(blank_laser_average,3))
           +r"$\pm$" + str(round(blank_laser_std,3))
           +"\n" 
            +r"Sample Laser Signal ($\pm \sigma$): "+str(round(sample_laser_average,3))
           +r"$\pm$" + str(round(sample_laser_std,3))
           +"\n" 
            +r"Blank PL Signal ($\pm \sigma$): "+str(round(blank_fluor_average,3))
           +r"$\pm$" + str(round(blank_fluor_std,3))
           +"\n"
            +r"Sample PL Signal ($\pm \sigma$): "+str(round(sample_fluor_average,3))
           +r"$\pm$" + str(round(sample_fluor_std,3))
           , transform=ax.transAxes
           ,bbox=dict(boxstyle="round",facecolor='white', alpha=0.1))
        
    
    # Graphing
    
    a.plot(blank_data_first[:,0],blank_data_first[:,1], color='black')
    a.plot(sample_data_first[:,0],sample_data_first[:,1], color='red')
    
    ax.set_title("\n".join(wrap(os.path.basename(root.sample_files[0]),60)), fontsize="10")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("PL (a.u.)", color="black")
    ax.set_xlim (375, 750)
    ax.set_ylim(-1)
    
    axins2 = inset_axes(ax, width="100%", height="100%", loc='upper left',
                bbox_to_anchor=(0.25,0.55,0.3,0.4), bbox_transform=ax.transAxes)
    axins2.plot(blank_data_first[:,0],blank_data_first[:,1], color='black')
    axins2.plot(sample_data_first[:,0],sample_data_first[:,1], color='red')
    axins2.set_xlim(lbound_low,lbound_high)
    #a.text(0.25, 0.55,"test",transform=ax.transAxes, fontsize=8)
    axins2.set_title(str(lbound_low)+" < Laser (nm) < "+str(lbound_high), fontsize=8, pad=3)
    axins2.tick_params(axis='both', which='major', labelsize=8)

    axins3 = inset_axes(ax, width="100%", height="100%", loc='upper left',
                bbox_to_anchor=(0.65,0.55,0.3,0.4), bbox_transform=ax.transAxes)
    axins3.plot(blank_data_first[:,0],blank_data_first[:,1], color='black')
    axins3.plot(sample_data_first[:,0],sample_data_first[:,1], color='red')
    axins3.set_xlim(fbound_low,fbound_high)
    axins3.set_title(str(fbound_low)+" < PL (nm) < "+str(fbound_high), fontsize=8, pad=3)
    axins3.tick_params(axis='both', which='major', labelsize=8)
    
    #plot and scale dyamically
    index_low = np.abs(sample_data_first[:,0] - fbound_low).argmin()
    index_high = np.abs(sample_data_first[:,0] - fbound_high).argmin()
    
    ymin_PL = min(blank_data_first[index_low:index_high,1])
    ymax_PL = max(sample_data_first[index_low:index_high,1])
    axins3.set_ylim(ymin_PL, ymax_PL)
    
    canvas.draw()

def counts_region(data, lbound_low, lbound_high, fbound_low, fbound_high):
    # data is array with 2 columns; 1st col: wavelength(nm), 2nd col: Counts
    #lambda_min/max in nm
    laser_signal = 0
    fluor_signal = 0
    
    for i in data:
        if lbound_low < i[0] < lbound_high:
            laser_signal += i[1]*i[0]/1000 #multiply counts by wavelength and wavelength spacing
        elif fbound_low < i[0] < fbound_high:
            fluor_signal += i[1]*i[0]/1000
            
    return laser_signal, fluor_signal
    
def clear_canvas():
    fig.clf()
    canvas.draw()

def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
    


# ------------Frame 4 - Data Buttons-----------------
frame4 = tk.Frame(root, bd="5")

addplotbut = tk.Button(master=root, text="Calculate QY", command=add_plot)
addplotbut.pack()

clearplot = tk.Button(master=root, text="Clear plot", command=clear_canvas)
clearplot.pack()

button = tk.Button(master=root, text="Quit", command=_quit)
button.pack(side=tk.BOTTOM)

frame4.pack(side=tk.BOTTOM,fill=tk.BOTH, expand=1)


# ------------Run-----------------
tk.mainloop()