'''
Created on 22 mars 2018

@author: Marija
'''

from PyQt5 import QtWidgets # (the example applies equally well to PySide)
import sys
import numpy as np
#import threading
#import time
#from blaze.expr.expressions import symbol
import pyqtgraph as pq

class myWidgetTest(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        #self.settings = QtCore.QSettings('MaxIV', 'GA')
        
        self.gendata()
        self.myLayout()
        
    def resetGraph(self):
        print 'resetting graph'
        self.plotWidget.clear()
    
    def gendata(self):
        self.xdata = np.linspace(0,10,10)
        self.ydata = np.random.rand(10).reshape((10))
        
    def myLayout(self):
        self.layout = QtWidgets.QVBoxLayout(self)
        self.inputLayout = QtWidgets.QHBoxLayout() 
        self.plotLayout = QtWidgets.QHBoxLayout() 
        
        self.layout.addLayout(self.inputLayout)
        self.layout.addLayout(self.plotLayout)
        
        self.resetButton = QtWidgets.QPushButton('reset')
        self.resetButton.clicked.connect(self.resetGraph)
        self.inputLayout.addWidget(self.resetButton)
        
        
        spacerItemV = QtWidgets.QSpacerItem(200, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)        
        self.inputLayout.addSpacerItem(spacerItemV)
        
        self.plotWidget = pq.PlotWidget(useOpenGL=True)
        self.plotWidget.setSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Maximum)
        self.plot1 = self.plotWidget.plot()
        #print len(self.ydata)
        self.plot1.setData(self.ydata)
        #self.plot1.setData(np.random.rand(10,1))
        self.plotLayout.addWidget(self.plotWidget)
                
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    myapp = myWidgetTest()
    myapp.show()
    sys.exit(app.exec_())
