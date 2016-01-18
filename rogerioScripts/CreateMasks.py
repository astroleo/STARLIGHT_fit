#!/usr/bin/python
import os, glob
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy import interpolate
import re
import time
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import imp
import gc
from PyPopStarsLogo import Logo

from matplotlib import pyplot as plt

class MaskSpec():
    
    def __init__(self, spec,graph,original_mask=np.zeros((1,3))):
        self.graph = graph
        self.spec=spec
        print
        if len(shape(original_mask)) == 2:
                    self.original_mask=original_mask
        else:
                    self.original_mask=reshape(original_mask,(1,3))
                
        self.xso = self.spec[:,0]  
        self.yso = self.spec[:,1]

        if (self.original_mask[0][0]) == 0:
            print "A new mask is beeing created"
            self.xs = []#list(line.get_xdata())
            self.ys = []#list(line.get_ydata())
            self.wei=[]
            self.masked=np.zeros((1,3))            
        else:
            self.xs = []#list(line.get_xdata())
            self.ys = []#list(line.get_ydata())
            self.wei=[]            
            self.masked=self.original_mask
            print "Editting an existing mask"
            for i in arange(0,len(self.original_mask[:,0])):
                self.xs.append(self.original_mask[:,0][i])
                self.xs.append(self.original_mask[:,1][i])
                self.wei.append(self.original_mask[:,2][i])
                if self.original_mask[:,2][i] == 2:
                    cor='green' 
                if self.original_mask[:,2][i] == 0:
                    cor='red' 
                plot(self.xso[(self.xso>=self.xs[-2]) & (self.xso<=self.xs[-1])],self.yso[(self.xso>=self.xs[-2]) & (self.xso<=self.xs[-1])],color=cor,lw='2')
                self.graph.figure.canvas.draw()

        self.cid = graph.figure.canvas.mpl_connect('button_press_event', self)
        self.key= graph.figure.canvas.mpl_connect('key_press_event', self.key_press_callback)        

        self.toremove_x=[]
        self.toremove_y=[]
        self.tmp=[]
        self.temp_remove=[]
    def __call__(self, event):
        if event.button == 1:
            print 'You press button 1 to zomm, to mask use right button'
        if event.button >=2:
            if event.inaxes!=self.graph.axes: return
#            self.bt.append(event.button)
            self.xs.append(event.xdata)
            
            if len(self.xs)%2 == 0:  # only odd numbers 
                if event.button == 3:
                    self.wei.append(0)
                    print 'You press button 3 to mask, maksed regions is', self.xs[-2],self.xs[-1],self.wei[-1]
                    self.tmp=np.column_stack((self.xs[-2],self.xs[-1],self.wei[-1]))
                    if (self.masked[0][0])==0:
                         self.masked=np.row_stack((self.tmp))
                    else:
                        self.masked=np.row_stack((self.masked,self.tmp))
                    plot(self.xso[(self.xso>=self.xs[-2]) & (self.xso<=self.xs[-1])],self.yso[(self.xso>=self.xs[-2]) & (self.xso<=self.xs[-1])],color='red',lw='2')
                    self.graph.figure.canvas.draw()
                if event.button == 2:
                    self.wei.append(2)
                    print 'You press button 2 attributed weight 2 for this region', self.xs[-2],self.xs[-1],self.wei[-1]
                    self.tmp=np.column_stack((self.xs[-2],self.xs[-1],self.wei[-1]))
                    if (self.masked[0][0])==0:
                         self.masked=np.row_stack((self.tmp))
                    else:
                        self.masked=np.row_stack((self.masked,self.tmp))
                    plot(self.xso[(self.xso>=self.xs[-2]) & (self.xso<=self.xs[-1])],self.yso[(self.xso>=self.xs[-2]) & (self.xso<=self.xs[-1])],color='green',lw='2')                   
                    self.graph.figure.canvas.draw()

    
    def key_press_callback(self, event):
        'Used to remove points from the mask'
        if not event.inaxes: return
        if event.key =='d':
            self.toremove_x = event.xdata
            self.toremove_y = event.ydata
            self.temp_remove = self.masked[(self.toremove_x > self.masked[:,0]) & (self.toremove_x<self.masked[:,1])]
            self.masked = self.masked[~((self.toremove_x > self.masked[:,0]) & (self.toremove_x<self.masked[:,1]))]
            plot(self.xso[(self.xso>=self.temp_remove[:,0][0]) & (self.xso<=self.temp_remove[:,1][0])],self.yso[(self.xso>=self.temp_remove[:,0][0]) & (self.xso<=self.temp_remove[:,1][0])],color='blue',lw='2')   
            self.graph.figure.canvas.draw()
            
            print 'You deleted region: ', self.temp_remove[:,0][0],'-', self.temp_remove[:,1][0], ' with weight: ',  int(self.temp_remove[:,2][0])
            
        if event.key == 'q':
            plt.close()
            
            
                #plt.draw()



def starlightMask(filename,maskname='mask.sm',new_mask=True):
    ''' This function uses the class MaskSpec to create mask for the starlight code.
        It alow for edditing existing masks or creating new ones.
        
        **kARGs:
         - filename = ascii file with the spectrum, where the first column is 
           the wavelength and the second is the flux.
         - maskname = the name of the new or to edit mask file.
         - new_mask = True or False. False to create a new.
        ** Interactive Keys
         + Mouse left button (bt1) zomm in button
         + Mouse right button (bt3) to mask region
         + Mouse middle button (bt2) to give weight 2 to the region.
         + Use key 'd' to delete a region from the mask (you have to put the mouse pointer in the region end tipe 'd')
         + Use key 'q' to quit.
         
    example:
    
    starlightMask(filename,maskname='mask.sm',original_mask=True)
    
    '''
    spec=np.loadtxt(filename,usecols=(0,1))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('Zoom= bt 1, Mask = bt3, Weight = bt2, Remove=\"d\" and \"q\" to quit')
    
    graph, = ax.plot(spec[:,0],spec[:,1], color='black')  # empty line

    if not new_mask:
        masked=np.loadtxt(maskname,skiprows=1)
        maskspec = MaskSpec(spec,graph,masked)
    else:
        maskspec = MaskSpec(spec,graph)        
      

    plt.show()
    makemask=maskspec.masked[maskspec.masked[:,0].argsort()]
    s=open(maskname,'w')
    s.write(str(len(makemask[:,0]))+'\n')
    np.savetxt(s,makemask,fmt=('%.2f','%.2f','%.0f'))
    s.close()
    return makemask










