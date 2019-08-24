# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 23:43:01 2018

@author: user
"""

import numpy as np
import math
import pyqtgraph.opengl as gl
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from scipy import linalg
from pyqtgraph.Qt import QtCore, QtGui
from PyQt5 import QtWidgets
from pyqtgraph.opengl.GLGraphicsItem import GLGraphicsItem

EPSILON=0.000001

def f(x,y):
    return ((x-3)*(x-3)/25)+((y-4)*(y-4)/16)-1
def g(x,y):
    return ((x-3)*(x-3)/4)-((y-4)*(y-4)/9)-1
def fx(x,y):
    return (2.0*x-6)/25
def fy(x,y):
    return (1.0*y-4)/8
def gx(x,y):
    return (1.0*x-3)/2
def gy(x,y):
    return -(2.0*y-8)/9
def newton(x,y):
    i=0
    xold=x
    yold=y
    lst = list()
    fold = f(xold, yold)
    gold = g(xold, yold)
    fdx = fx(xold, yold)
    fdy = fy(xold, yold)
    gdx = gx(xold, yold)
    gdy = gy(xold, yold)
    Delta = fdx*gdy - fdy*gdx
    xnew=xold
    ynew=yold
    xold=xnew
    yold=ynew
    lst.append([xnew,ynew,0])

    err = 1.0
    print("     i            xn		              yn               error\n");
    print("------------------------------------------------------------\n")
    while(err>=EPSILON):
        Delta=math.fabs(fx(xold, yold)*gy(xold, yold)-fy(xold, yold)*gx(xold, yold));
        if(Delta==0):
            Delta=EPSILON*10
        if(i==0):
            print("     " ,i,"\t", xnew,"\t",ynew,"\t","---" )
        else:
            print("     " ,i,"\t", xnew,"\t",ynew,"\t", err)
        xold=xnew
        yold=ynew
        h=-((f(xold, yold)*gy(xold, yold)-g(xold, yold)*fy(xold, yold)))/Delta
        k = -((-f(xold, yold)*gx(xold, yold)+g(xold, yold)*fx(xold, yold)))/Delta
        xnew=xold+h
        ynew=yold+k
        err=math.fabs((xnew - xold)*(xnew - xold)+(ynew - yold)*(ynew - yold))**0.5
        xold=xnew
        yold=ynew
        i=i+1
        lst.append([xnew,ynew,0])
    print("------------------------------------------------------------\n")
    print("     " ,i,"\t", xnew,"\t",ynew,"\t", err)
    return i,lst
    
x = float(input('請輸入x： '))
y = float(input('請輸入y： '))
time,lis=newton(x,y)
pos=np.array(lis,'float32')



app = QtGui.QApplication([])
glWidget = gl.GLViewWidget()
glWidget.show()
gridx = gl.GLGridItem()
gridx.rotate(90, 0, 1, 0)        #90度
#gridx.rotate(45, 0, 0, 1)
gridx.translate(0, 0, 0)#-10,0,0
glWidget.addItem(gridx)
gridy = gl.GLGridItem()
gridy.rotate(90, 1, 0, 0)
#gridy.rotate(45, 0, 0, 1)
gridy.translate(0, 0, 0)#0,-10,0
glWidget.addItem(gridy)
gz = gl.GLGridItem()
gz.translate(0, 0, 0)
glWidget.addItem(gz)
#glWidget.opts['distance'] = 24

i=0
connected = np.array([[pos[i],pos[i],pos[i]],[pos[i+1],pos[i+1],pos[i+1]]])
line = gl.GLLinePlotItem(pos=connected, width=3.5, antialias=True)
glWidget.addItem(line)

def update():
    ## update volume colors
    global i
    global time
    i=i+1
    if(i>time-2):
        i=1
    connected = np.array([[pos[i],pos[i],pos[i]],[pos[i+1],pos[i+1],pos[i+1]]])
    line.setData(pos=connected)


t = QtCore.QTimer()
t.timeout.connect(update)
t.start(1000)

if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) :
        QtGui.QApplication.instance().exec_()