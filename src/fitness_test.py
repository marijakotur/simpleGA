'''
Created on 6 mars 2017

@author: Marija
'''

from PyQt4 import QtGui, QtCore # (the example applies equally well to PySide)
import sys
import numpy as np
#import threading
import time
#from blaze.expr.expressions import symbol
import pyqtgraph as pq
import PyTango as pt

class myWidgetTest(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        #self.settings = QtCore.QSettings('MaxIV', 'GA')
             
        
        self.adName = 'gunlaser/devices/ad5370dac'
        self.adDevice = pt.DeviceProxy(self.adName)
        #self.ad = self.adName.read_attribute('channel0')
        
        #self.cameraName = 'gunlaser/thg/camera'
        self.cameraName = 'gunlaser/cameras/blackfly_test01'
        self.cameraDevice = pt.DeviceProxy(self.cameraName)
        self.cameraDevice.set_timeout_millis(3000)
        #expTime = self.cameraDevice.read_attribute('ExposureTime').value
        #print expTime
        self.cameraDevice.write_attribute('Gain',15)
        self.cameraDevice.write_attribute('ExposureTime',1000)
        self.image = self.cameraDevice.read_attribute('Image').value.astype(np.double)
        
        
        self.ampFactor = 73.0
        self.maxVoltage = 4.2 * self.ampFactor
        #self.voltages = np.linspace(0, self.maxVoltage, 37)
        #self.voltages = np.linspace(0, 36, 37)
        self.noChannels = 40
        self.voltages = [5]*self.noChannels #for 37 electrodes
          
        self.scanTimer = QtCore.QTimer()
        self.scanTimer.timeout.connect(self.updateImage)
          
        self.initGA()       
        self.myLayout()

        self.xData = np.array(range(0, 10))
        self.yData = np.array(range(0, 10)) #np.random.rand(1,10)
        
        #self.measureData()
        
        
    def FWHM(self,):
        pass
    
    def initGA(self):
        a = 1.0
        h = a*np.sqrt(3)/2
        self.electrodeCoordinates = [[-3*a/2,3*h],[-a/2,3*h],[a/2,3*h],[3*a/2,3*h],[-2*a,2*h],[-a,2*h],[0,2*h],[a,2*h],[2*a,2*h],
[-5*a/2,h],[-3*a/2,h],[-a/2,h],[a/2,h],[3*a/2,h],[5*a/2,h],[-3*a,0],[-2*a,0],[-a,0],[0,0],[a,0],[2*a,0],[3*a,0], 
[-5*a/2,-h],[-3*a/2,-h],[-a/2,-h],[a/2,-h],[3*a/2,-h],[5*a/2,-h], 
[-2*a,-2*h],[-a,-2*h],[0,-2*h],[a,-2*h],[2*a,-2*h], 
[-3*a/2,-3*h],[-a/2,-3*h],[a/2,-3*h],[3*a/2,-3*h],[-9*a/2,-3*h],[-7*a/2,-3*h],[-4*a,-2*h]]
        self.electrodeCoordinates = np.array(self.electrodeCoordinates)
        
        self.channelsList = [23,21,0,2,
                             24,22,20,1,3,
                             27,26,25,4,5,6,
                             28,29,30,32,17,18,7,
                             39,38,36,13,15,16,
                             37,34,9,11,14,
                             35,33,10,12,
                             8,19,31]
        self.channelsOuter = [23,21,0,2,24,3,27,6,28,7,39,16,37,14,35,33,10,12]
        self.channelsMiddle = [22,20,1,26,5,29,18,38,15,34,9,11]
        self.channelsInner = [25,4,30,17,36,13]
        self.channelsCenter = [32]
        self.channelsDisconnected = [8,19,31]
              
        self.populationSize = 20
        self.numberOfGenes = 40
        self.population = []   
        
        self.elitismCoef = 0.2
        self.mutationCoef = 0.2
        self.twoPtXoverCoef = 0.6

        
        self.fitness = [None] * self.populationSize
        self.geneLowerBound = 0
        self.geneUpperBound = self.maxVoltage
        self.creepRate = 0.1  
                
        self.fitnessCostWeight = 0.25 / self.populationSize
        self.referencePopulation = [0]*self.numberOfGenes
       


    def sortByFitness(self, population):
        print 'measuring fitness'
        self.fitness = []
        #measure fitness
        image_size_x = self.image.shape[0]
        image_size_y = self.image.shape[1]
        for individual in population:
            s=0
            for chInd, channel in enumerate(self.channelsList):
                chName = 'channel'+ str(self.channelsList[chInd])
                #print chInd, channel
                #if chInd
                self.adDevice.write_attribute(chName,individual[chInd]/self.ampFactor)
                #fitness: central circle white
                #===============================================================
                # for (x,y), pixel in np.ndenumerate(self.image):
                #     if (np.sqrt(x**2+y**2) < l):
                #         s += pixel
                #===============================================================
            #s = np.sum(self.image[image_size_x/2-100:image_size_x/2,image_size_y/2-100:image_size_y/2]) 
            s = np.max(np.sum(self.image,1))
            self.fitness.append(s)
        #sort by fitness          
        j = np.argsort(self.fitness)
        j = np.flipud(j)
        population = population[j]
        population = np.squeeze(population)
        #write the best shape again before updating image
        for chInd, channel in enumerate(self.channelsList):
                chName = 'channel'+ str(self.channelsList[chInd])
                self.adDevice.write_attribute(chName,population[0][chInd]/self.ampFactor)
        self.updateImage()  

        return population                 
   

    def updateImage(self):
        #print 'updating image'
        self.image = self.cameraDevice.read_attribute('Image').value.astype(np.double)
        self.img.setImage(self.image)
        self.view.update()
        #self.cameraWindow.update()
        self.scanTimer.start(100)

#===============================================================================
#     def setVoltages(self):
#                    
#         self.spot_plot_list = []
#         self.plot11 = self.plotWidget1.addPlot()
#         self.scatterPlotItem1 = pq.ScatterPlotItem(size=30, pen=pq.mkPen(None))
#         self.scatterPlotItem1.addPoints(x=self.electrodeCoordinates[:,0], y=self.electrodeCoordinates[:,1])        
#         
#         brushes = []
#         for i in range(0,self.noChannels):
#             chName = 'channel'+ str(self.channelsList[i])
#             self.adDevice.write_attribute(chName,self.voltages[i]/self.ampFactor)
#             if i<=len(self.electrodeCoordinates):
#                 brushes.append(QtGui.QColor(255, self.voltages[i]/self.maxVoltage*255, 0, 220))
# 
#         self.scatterPlotItem1.setBrush(brushes)
#         self.plot11.addItem(self.scatterPlotItem1)
#         self.scatterPlotItem1.sigClicked.connect(self.clicked)
#         
#         self.textItem = pq.TextItem()
#         self.plot11.addItem(self.textItem)
#===============================================================================

    def myLayout(self):
        self.layout = QtGui.QVBoxLayout(self) #the whole window, main layout
        self.inputLayout = QtGui.QHBoxLayout() 
        self.plotLayout = QtGui.QHBoxLayout() 
        self.gridLayout1 = QtGui.QGridLayout()
        self.gridLayout2 = QtGui.QGridLayout()
        self.gridLayout3 = QtGui.QGridLayout()
        self.gridLayout4 = QtGui.QGridLayout()

        #QVBoxLayout1 = QtGui.QHBoxLayout()
        #self.layout1 = self.layout.addLayout(QVBoxLayout1)

        
        self.gridLayout2.addItem(QtGui.QSpacerItem(10, 15,QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding), 3,0)


        self.currentGenLabel = QtGui.QLabel()
        
        self.plotWidget2 = pq.PlotWidget(useOpenGL=True)
        self.plotWidget2.setSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Maximum)
        self.plot21 = self.plotWidget2.plot()
        self.plot22 = self.plotWidget2.plot()

        self.plotWidget1 = pq.PlotWidget()
        #self.plot11 = self.plotWidget1.plot()
        #self.plotWidget1.setSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Maximum)
        self.plot11 = pq.ScatterPlotItem(size=50, pen=pq.mkPen('w'), pxMode=True, symbol='o')
        self.plotWidget1.addItem(self.plot11)
        self.plotWidget1.setAspectLocked(1)

        self.cameraWindow = pq.GraphicsLayoutWidget()
        self.cameraWindow.setMaximumSize(350,350)
        #self.cameraWindow.setSizePolicy(QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding)
        self.view = self.cameraWindow.addViewBox()
        self.view.setAspectLocked(True)
        self.img = pq.ImageItem(border='w')
        self.view.addItem(self.img)
        self.img.setImage(self.image)
        self.updateImage()
                        
        self.plotLayout = QtGui.QHBoxLayout()      
        self.plotLayout.addWidget(self.plotWidget2)
        self.plotWidget1.setVisible(True)
        self.plotWidget1.setMaximumHeight(400)
        self.plotLayout.addWidget(self.plotWidget1)
        self.plotWidget2.setVisible(True)
        self.plotWidget2.setMaximumHeight(400)
        self.plotLayout.addWidget(self.cameraWindow)
        
        self.layout.addLayout(self.inputLayout)
        self.layout.addLayout(self.plotLayout)
        self.inputLayout.addLayout(self.gridLayout1)
        self.inputLayout.addLayout(self.gridLayout2)
        self.inputLayout.addLayout(self.gridLayout3)
        #self.inputLayout.addSpacerItem(QtGui.QSpacerItem(10, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum))
        self.layout.addLayout(self.plotLayout)

        
        
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    myapp = myWidgetTest()
    myapp.show()
    sys.exit(app.exec_())
