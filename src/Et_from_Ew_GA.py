'''
Created on 1 Jul 2018

@author: markot
'''

import warnings
from PyQt5 import QtWidgets, QtCore
#import logging
import sys
import time
#import threading
#import PyTango as pt
import numpy as np
import pyqtgraph as pq
from scipy.signal import tukey
#from scipy.stats import gennorm
from scipy.special import gamma 
from __builtin__ import str
#import scipy as sy
import traceback


if QtCore.QT_VERSION >= 0x50501:
    def excepthook(type_, value, traceback_):
        traceback.print_exception(type_, value, traceback_)
        QtCore.qFatal('')
sys.excepthook = excepthook
    
warnings.resetwarnings()

    
class SpectrometerCamera(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        #super().__init__()
        
        #self.lock = threading.Lock()       
        #self.scanTimer = QtCore.QTimer()
        #self.scanTimer.timeout.connect(self.updateImage)    
        self.w_wl = 1.6
        self.stop_timer = False
        self.shaperPhaseAmpCoef = 0.0
        self.shaperPhaseWidth = 0.5
        self.shaperPhaseCenter = 262.0 
        self.shaperPhaseShapeParam = 2.0
        self.gaussianPhaseCoef = 0.0
        self.ampFilterWidth = 10.0
        self.ampLoss = 1.0
        self.ampMaskOn = False
        self.shaperOn = False          

        #self.absE_wl_mask = np.empty_like(self.wl)
        self.fftsize = 4096*4
        #self.ampMask = np.ones(self.wlsize)
        self.tukey_ind_l = 10
        #self.tukey_ind_r = self.wlsize-10
        self.tukey_param = 0.0
        
        self.maxFitnessTrend = []
        self.avgFitnessTrend = []
        
        #square box
        self.goalTempWidth = 1.0
        self.goalTempShape = np.ones(self.fftsize)  
        self.boxDiff = 1
        #self.defineSpectrum()
        #self.spectralToTemporal()
        
        
        self.w_nu = 1.6
        
        #set up and update widgets
        self.defineNuAxis()        #initial phase parameters
        self.biasPhaseCoef = 0.0 
        self.phi_nu = np.zeros(self.wlsize)
        self.setGAInitParameters()
        self.setup_layout()    
        self.defineSpectralAmp()
        self.setBiasPhase()
        self.calculateBox()
        self.Enu2Et()
        self.updateFFT() 
        
        self.scanTimer = QtCore.QTimer()
        self.scanTimer.timeout.connect(self.scanUpdateAction)

    def setGAInitParameters(self):
        self.populationSize = 20
        self.numberOfGenes = 36
        self.genAmp = 20.0 # allowed gene variation
        self.population = []   
        self.elitismCoef = 0.2   
        self.mutationCoef = 0.3
        self.twoPtXoverCoef = 0.5
        
        self.fitness = [None] * self.populationSize
        self.geneLowerBound = -2.00
        self.geneUpperBound = 2.00
        self.creepRate = 0.1  
                
        self.fitnessCostWeight = 0.25 / self.populationSize
        self.referencePopulation = [0]*self.numberOfGenes
        
    def initGA(self):
        print 'initializina GA'
        #print str(self.popSizeSpinbox.value())
        self.populationSize = self.popSizeSpinbox.value()
        self.numberOfGenes = self.numberOfGenesSpinbox.value()
        self.genAmp = self.geneAmpSpinbox.value()
        
    def generateInitialPopulation(self,popSize,numberOfGenes):
        self.population = []
        noGenes = self.numberOfGenes
        for i in range(0, self.populationSize):
            #print self.genAmp
            genes = ((np.random.rand(noGenes-1)-0.5)*2)*self.genAmp
            print 'genes', genes
            #genes -= genes[0]
            self.population.append(genes)
        self.population = np.array(self.population)
        #print 'population:', self.population

    def newGeneration(self):
        print 'generating the next generation'
        sortedPopulation = np.copy(self.sortByFitness(self.population))
        # print sortedPopulation
        newPopulation = []
        for i in range(0,self.populationSize):
            if i in self.elitismIdx:
                newPopulation.append(sortedPopulation[i])
            elif i in self.mutationIdx:
                #select individual j (from the top half of the population) to mutate 
                j = np.random.randint(0,np.floor(self.populationSize/2))
                #print "individual to mutate is %d" %j
                indToMutate = np.copy(sortedPopulation[j])
                #print indToMutate
                newPopulation.append(self.creepMutation(indToMutate))
            elif i in self.crossOverIdx:
                #loop through all the xOver indices
                #select two individuals to mate
                j1 = np.random.randint(0,np.floor(self.populationSize/2))
                j2 = np.random.randint(0,np.floor(self.populationSize/2))
                # make sure they are different
                while j2 == j1:
                    j2 = np.random.randint(0,np.floor(self.populationSize/2))
                    #print j1, j2
                    # taking only one child from x-over
                #newPopulation.append(self.onePtXover(self.population[j1], self.population[j2]))
                newPopulation.append(np.array(self.onePtXover(self.population[j1], self.population[j2])[1]))
                #self.population = np.array(self.population)
        
        self.population = np.array(newPopulation)
        
        #print type(self.population)
        #print self.population.shape
        #self.population = np.array(self.population)
        
    def scanUpdateAction(self):
        time.sleep(0.0)
                
        if self.generation  == 0:
            #print 'Generating initial population'
            #print 'Current gen: ', self.generation
            self.fitness = 0.0  
            self.calcIdx()
            self.fitnessTrend = []
            self.generateInitialPopulation(self.populationSize, self.numberOfGenes)
        else:     
            #print 'Current gen: ', self.generation
            self.genData = np.append(self.genData, self.generation)
            #print self.generation
            self.newGeneration()
            print 'new gen created'
            self.maxFitnessTrend = np.append(self.maxFitnessTrend,np.max(self.fitness))
            #print self.maxFitnessTrend
            self.avgFitnessTrend = np.append(self.avgFitnessTrend,np.average(self.fitness))
            #print 'avg. fitness', self.avgFitnessTrend
            self.fitnessPlotMax.setData(self.maxFitnessTrend)
            self.fitnessPlotAvg.setData(self.avgFitnessTrend)
            self.fitnessPlotMax.setPen('w')
            self.fitnessPlotAvg.setPen('r')
            self.writeBestShape()
        
        self.generation += 1
        
        #print 'Generations: ', self.genData
        #print 'Fitness trend: ', self.fitnessTrend
        self.scanTimer.start(100)
        #print self.elitismCoef
        
    def startScan(self):
        print 'Scan started'
        self.scanData = None
        self.trendData = None
        #self.running = True
        self.genData = np.array([])
        self.maxFitnessTrend = np.array([])
        self.avgFitnessTrend = np.array([])
        self.generation = 0
        self.initGA()
        self.scanUpdateAction()

    def stopScan(self):
        print 'Scan stopped'
        self.scanData = None
        #self.running = False
        self.scanTimer.stop()

    def calculateBox(self):
        
        #temp_centroid = np.sum(self.Intt*self.t)/np.sum(self.Intt)
        #print 'temp centroid', temp_centroid
        tukey_lim_r = - self.goalTempWidth/2
        tukey_lim_l = self.goalTempWidth/2
        
        #print tukey_lim_l, tukey_lim_r
        self.tukey_ind_r = np.max(np.where(self.t<tukey_lim_l))
        self.tukey_ind_l = np.min(np.where(self.t>tukey_lim_r))
        self.box = np.hstack((np.zeros(self.tukey_ind_l), tukey(self.tukey_ind_r-self.tukey_ind_l,0), np.zeros(self.fftsize-self.tukey_ind_r)))
        print 'sum box', np.sum(self.box)
        #self.box *= np.sum(self.Intt)/np.sum(self.box)        

    def fwhm(self,data):
        data = np.abs(data)
        data /= np.max(data)
        a = [np.diff(np.sign(data-0.5))]
        nonzeros = np.nonzero(a[0])
        fwhm = np.max(nonzeros)-np.min(nonzeros)
        #print np.min(nonzeros), np.max(nonzeros)
        return fwhm   
        
    def measureFitness(self, ind):
        #ind is an array of genes
        #evaluate phase p based on those genes
        #fitness f is the overlap between p and the target function tf
        #print 'measuring fitness'
        #print ind
        
        
        #in the direct basis
        phase_coarse = np.cumsum(ind)-ind   
        self.phi_nu = np.interp(np.linspace(0,self.numberOfGenes-1,len(self.nu)), range(self.numberOfGenes-1), phase_coarse)
        
        #in the polynomial basis
        #self.phi_nu = np.zeros(self.wlsize)
        #n=0
        #for g in ind:
        #    self.phi_nu += g*(self.nu-np.mean(self.nu))**n 
        #    n += 1
            
        
        
        #in the cos basis
        #self.phi_nu = np.zeros(self.wlsize)
        #nu_lowest = 25
        #for ii in range(len(ind)):
        #    self.phi_nu += ind[ii]*np.cos(2*np.pi*ii/nu_lowest * (self.nu-np.mean(self.nu)))
            
        #self.absE_nu = np.exp(-4.0*np.log(2)*(self.nu-np.mean(self.nu))**2/self.w_nu**2)
        self.E_nu = np.sqrt(self.absE_nu) * np.exp(1j*self.phi_nu)
        #calculate Et from the new phase
        self.Enu2Et()
        
        #simple fitness function
        #tf = np.sum(np.abs(self.box*self.Intt))/np.sum(np.abs((1.0-self.box)*self.Intt))
        #inside_box = np.sum(self.box*self.Intt)
        #ind_outside = np.where(self.Intt>self.box)
        #outside_box = np.sum((self.Intt-self.box)[ind_outside])
        #tf = inside_box/outside_box
        #tf = self.fwhm(self.Intt)
        #tf = 
        #print 'ind fitness:', tf
        
        #cross-corelation fitness
        a=np.corrcoef(self.box,self.Intt)
        tf=a[0,1]
        #a=np.correlate(self.box, self.Intt)
        #tf=np.max(a)
        
        #difference
        #tf=-np.sum(np.abs(self.box-self.Intt))
        
        return tf
        #print 'test fitness'       

    def sortByFitness(self, population):
        print 'sorting'
        self.fitness = []
        #measure fitness
        for individual in population:
            #print 'individual', individual
            s = self.measureFitness(individual)
            self.fitness.append(s)
        #sort by fitness
        j = np.argsort(self.fitness)
        j = np.flipud(j)
        population = population[j]
        population = np.squeeze(population)
        self.current_best_individual = population[0]
        #print 'current best', self.current_best_individual
        #write the best shape again before updating image
        print 'sorted by fitness'
        #print 'fitness:', self.fitness
        return population 

    def calcIdx(self):
        # returns the indices of new generation individuals to be created through elitism, mutation and x-over 
        self.elitismIdx = np.arange(0, int(max(1,np.floor(self.populationSize * self.elitismCoef))))
        self.mutationIdx = np.arange(int(max(self.elitismIdx))+1, int(max(self.elitismIdx)+1  + np.floor(self.populationSize * self.mutationCoef)))
        self.crossOverIdx = np.arange(max(self.mutationIdx)+1, self.populationSize)
       
    def onePtXover(self,ind1,ind2):
        cxpoint = np.random.randint(1,self.numberOfGenes) 
        #print cxpoint
        ind1[cxpoint:], ind2[cxpoint:] = ind2[cxpoint:], ind1[cxpoint:]    
        return np.array([ind1, ind2])
        
    def creepMutation(self,individualToMutate):
        #change one gene according to
        #g_new = min( imax, max( imin, g_old + r*s*(imax-imin)
        # r - random number in the [-1,1] range
        # s - creep rate, usually 0.02
        # imin, imax - gene range bounds
        mutationLocation = np.random.randint(0,self.numberOfGenes-1)
        #print "mutation location is %d" %mutationLocation
        #print individualToMutate
        oldGene = individualToMutate[mutationLocation]
        r = np.random.uniform(low=-self.genAmp, high=self.genAmp, size=1)
        newGene = min( self.geneUpperBound, max( self.geneLowerBound, oldGene + r*self.creepRate*(self.geneUpperBound-self.geneLowerBound)))
        individualToMutate[mutationLocation] = newGene
        return individualToMutate    
    
    def setPopSize(self):
        self.populationSize = self.popSizeSpinbox.value()
        self.calcIdx()
    
    def setNumberOfGenes(self):
        self.numberOfGenes = self.numberOfGenesSpinbox.value()
        print 'setting no of genes to ' + str(self.numberOfGenes)
        
    def setGeneAmp(self):
        self.genAmp = self.geneAmpSpinbox.value()
    
    def setElitismCoef(self):
        self.elitismCoef = self.elitismCoefSpinbox.value()
        self.calcIdx()
        print 'Setting elitism to', self.elitismCoef

    def setMutationCoef(self):
        print "setmut"
        self.mutationCoef = self.mutationCoefSpinbox.value()
        self.calcIdx()
                
    def setCrossOverCoef(self):
        self.crossOverCoef = self.crossOverCoefSpinbox.value()
        self.calcIdx()
        
    def gennorm_with_var(self,x,A,beta,alpha,mu):
        try:
            pdf = A*beta/(2*alpha*gamma(1/beta))*np.exp(-(np.abs(x-mu)/alpha)**beta)
            return pdf    
        except ValueError:
            pass    
        
    def defineNuAxis(self):
        # higher to lower wavelengths so the frequencies will be in the right order
        # need for np.interp to work correctly
        print 'defining nu axis'
        self.wl_0 = 262
        self.wlsize = 1024
        self.wl = np.linspace(266, 258, self.wlsize)
        
        #lambda -> nu
        nuI = 3e8/self.wl[0]*1e9/1e12
        nuF = 3e8/self.wl[-1]*1e9/1e12
        self.nu = np.linspace(nuI,nuF,self.wlsize)
        self.d_nu = np.max(np.diff(self.nu))
        self.t=np.linspace(-1/2.0/self.d_nu,1/2.0/self.d_nu,self.fftsize)
        
    def defineSpectralAmp(self):
        #define spectral amp
        self.w_wl = self.spectralFWHMSpinbox.value()
        self.w_nu = 3e5/np.mean(self.wl)**2 * self.w_wl
        self.absE_nu = np.exp(-4.0*np.log(2)*(self.nu-np.mean(self.nu))**2/self.w_nu**2)
        #recalculate the field
        self.E_nu = np.sqrt(self.absE_nu) * np.exp(1j*self.phi_nu)
        
    def setBiasPhase(self):
        #define spectral phase
        print 'changed bias phase'
        self.phi_bias_wl = self.biasPhaseCoef*(self.wl-np.mean(self.wl))**2
        self.phi_bias_nu = self.biasPhaseCoef*(self.nu-np.mean(self.nu))**2
        self.phi_nu = self.phi_bias_nu
        #print self.geomPhaseCoef
        #recalculate the field
        self.E_nu = np.sqrt(self.absE_nu) * np.exp(1j*self.phi_nu)
        
    def genesToPhase(self):
        self.phi_nu = np.zeros(self.wlsize)
        n=0
        for g in ind:
            self.phi_nu += g*(self.nu-np.mean(self.nu))**(n+1) 
            n += 1
        
    def writeBestShape(self):
        #print 'current best', self.current_best_individual
        ind = self.current_best_individual
        
        #phase_coarse = np.cumsum(ind)-ind   
        #self.phi_nu = np.interp(np.linspace(0,self.numberOfGenes-1,len(self.nu)), range(self.numberOfGenes-1), phase_coarse)
        
        self.phi_nu = np.zeros(self.wlsize)
        n=0
        for g in ind:
            self.phi_nu += g*(self.nu-np.mean(self.nu))**(n+1) 
            n += 1
        
        #4self.phi_nu = np.zeros(self.wlsize)
        #nu_lowest = 24
        #for ii in range(len(ind)):
        #    self.phi_nu += ind[ii]*np.cos(2*np.pi*ii/nu_lowest * (self.nu-np.mean(self.nu)))
        
        self.E_nu = np.sqrt(self.absE_nu) * np.exp(1j*self.phi_nu)
        self.Enu2Et()
        self.updateFFT()
        
        self.genePlot.setData(ind)
        
    def Enu2Et(self):
        #pad E_nu with zeros
        pad_length = np.round((self.fftsize - len(self.nu))/2)
        self.Enu_padded = np.pad(self.E_nu,(pad_length,pad_length),'constant')
        #self.Enu_mask_padded = np.pad(self.Enu_mask,(pad_length,pad_length),'constant')

        self.Et = np.fft.fftshift(np.fft.fft(np.fft.fftshift(self.Enu_padded)))
        
        #self.Intt=np.abs(self.Et)**2/np.max(abs(self.Et)**2)
        self.Intt=np.abs(self.Et)**2
        self.Intt *= np.sum(self.box)/np.sum(self.Intt)
        #print 'sum Intt', np.sum(self.Intt)
        #self.Intt_mask = np.abs(self.Et_mask)**2/np.max(abs(self.Et_mask)**2)
        self.phit=np.unwrap(np.angle(self.Et))
        #self.phit_mask = np.unwrap(np.angle(self.Et_mask))       
        
#     def defineSpectrumOld(self):
#         #amp mask
#         if self.ampMaskOn:   
#             tukey_lim_r = self.ampMaskLeftSpinBox.value()
#             self.tukey_ind_r = np.min(np.where(self.wl<tukey_lim_r))
#             tukey_lim_l = self.ampMaskRightSpinBox.value()  
#             self.tukey_ind_l = np.max(np.where(self.wl>tukey_lim_l))
#             #print self.tukey_ind_l, self.tukey_ind_r   
#             self.tukey_param = self.ampMaskTukeyParamSpinBox.value()               
#             self.ampMask = np.hstack((np.zeros(self.tukey_ind_l), tukey(self.tukey_ind_r-self.tukey_ind_l,self.tukey_param), np.zeros(len(self.wl)-self.tukey_ind_r)))
#         else:
#             self.ampMask = np.ones(self.wlsize)
#         
#         self.absE_wl_mask = self.absE_wl*self.ampMask
#         self.ampLoss = np.sum(self.absE_wl_mask)/np.sum(self.absE_wl)
#         geom_phase = self.geomPhaseCoef*(self.wl-self.wl_0)**2
#         if self.shaperOn:
#             self.phi_wl = geom_phase + self.shaper_phase
#         else:
#             self.phi_wl = geom_phase                            
                
    def spectralToTemporalOld(self):
        # lambda -> nu, sample nu uniformly and interpolate
        
        print 'performing FFT'
        #interpolate
        self.absEnu = np.interp(self.nu_uniform, self.nu, self.absE_wl)
        self.absEnu_mask = np.interp(self.nu_uniform, self.nu, self.absE_wl_mask)
        #self.absEnu = np.multiply(self.absEnu,self.ampMask)
        self.phi_nu=np.interp(self.nu_uniform, self.nu, self.phi_wl)
        self.Enu = np.sqrt(self.absEnu) * np.exp(1j*self.phi_nu)
        self.Enu_mask = np.sqrt(self.absEnu_mask) * np.exp(1j*self.phi_nu)

        self.Et_mask = np.fft.fftshift(np.fft.fft(np.fft.fftshift(self.Enu_mask_padded)))
        #Et=np.fft.fftshift(np.fft.fft(abs(self.Enu_padded)))
        #Et=np.fft.fft(abs(self.Enu_padded))
        print 'done with FFT'
    
    def arrayThreshold(self,a):
        #print 'thresholding array'
        indices = np.where(a > np.max(a)/10000.0)
        #indices = np.where(a > 0.0)
        ind_left = np.min(indices)
        ind_right = np.max(indices)
        return ind_left, ind_right
#     
#     def updateElPlotView(self):
#         self.ViewboxYY.setGeometry(self.plotYYItem.vb.sceneBoundingRect())
#         self.ViewboxYY.linkedViewChanged(self.plotYYItem.vb, self.ViewboxYY.XAxis)                
        
    def updateFFT(self):
        print 'updating FFT plots'
        
        self.plotSpectralWidget.plot11.setData(self.wl,self.absE_nu)
        #self.plotSpectralWidget.plot12.setData(self.wl,self.ampMask)
        self.plotSpectralWidget.plot21.setData(self.wl,self.phi_nu)
        #self.plotSpectralWidget.plot22.setData(self.wl,self.shaper_phase)
        
        [ind1, ind2] = self.arrayThreshold(self.Intt)
        self.plotTemporalWidget.plot11.setData(self.t[ind1:ind2],self.Intt[ind1:ind2])
        #self.plotTemporalWidget.plot13.setData(self.t[ind1:ind2],self.Intt_mask[ind1:ind2])
        self.plotTemporalWidget.plot21.setData(self.t[ind1:ind2],self.phit[ind1:ind2]-np.mean(self.phit[ind1:ind2])) 
        self.plotTemporalWidget.plot15.setData(self.t[ind1:ind2],self.box[ind1:ind2])
        
        print 'done updating FFT plots'
        #if self.stop_timer is not True:
        #    self.cameraTimer.start(100)
            
    def closeEvent(self, event):
        self.stop_timer = True
        self.deleteLater()
        print 'stopping'
        
    def updateEverything(self):
        print 'updating everything'
        self.biasPhaseCoef = self.biasPhaseCoefSpinbox.value()
        #print self.biasPhaseCoef
        self.defineSpectralAmp()
        self.setBiasPhase()
        self.Enu2Et()
        self.updateFFT() 
        #print str(self.ampLoss,'%6.2f')
        
    def toggleAmpMask(self):
        if self.ampFilterRadioButton.isChecked():
            self.ampMaskOn = True
        else:
            self.ampMaskOn = False
        self.updateEverything()
        
    def applyShaperPhase(self):
        print 'shaper phase toggled'
        if self.shaperRadioButton.isChecked():
            self.shaperOn = True
        else:
            self.shaperOn = False
        self.defineSpectralAmp()
        self.updateFFT()    
                    
#     def generateRandomShaperPhase(self):
#         noGenes = 20
#         genAmp=3
#         genes = (np.random.rand(noGenes)*2-1)*genAmp
#         genes -= genes[0]
#         phase_coarse = np.cumsum(genes)-genes            
#         phase = np.interp(np.linspace(0,noGenes-1,len(self.nu_uniform)), range(noGenes), phase_coarse)
#         self.updateEverything()
#         return phase
    
    def setup_layout(self):
        
        print 'setting up layout'
        self.layout = QtWidgets.QHBoxLayout(self) #the whole window, main layout
        self.plotLayout = QtWidgets.QVBoxLayout()
        self.controlsLayout = QtWidgets.QVBoxLayout()
        self.GAlayout = QtWidgets.QVBoxLayout()
        self.fitnessLayout = QtWidgets.QVBoxLayout()

        vspaceritem = QtWidgets.QSpacerItem(100,100,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding)


        self.plotSpectralWidget = PlotYYWidget()
        #self.plotSpectralWidget.plot11.setData(self.wl,self.absE_wl)
        #self.plotSpectralWidget.plot12.setData(self.wl,self.ampMask)
        #self.plotSpectralWidget.plot21.setData(self.wl,self.phi_wl)
        #self.plotSpectralWidget.plot22.setData(self.wl,self.shaper_phase)
        self.plotSpectralWidget.show()
        self.plotLayout.addWidget(self.plotSpectralWidget)
        
        self.plotTemporalWidget = PlotYYWidget()
        #[ind1, ind2] = self.arrayThreshold(self.Intt)
        #self.plotTemporalWidget.plot11.setData(self.t[ind1:ind2],self.Intt[ind1:ind2])
        #self.plotTemporalWidget.plot13.setData(self.t[ind1:ind2],self.Intt_mask[ind1:ind2])
        #self.plotTemporalWidget.plot21.setData(self.t[ind1:ind2],self.phit[ind1:ind2]-np.mean(self.phit[ind1:ind2]))
        self.plotLayout.addWidget(self.plotTemporalWidget)       

        self.spectralFWHMSpinbox = QtWidgets.QDoubleSpinBox()
        self.spectralFWHMSpinbox.setDecimals(2)
        self.spectralFWHMSpinbox.setValue(1.6)
        self.spectralFWHMSpinbox.setMinimum(0.1)
        self.spectralFWHMSpinbox.valueChanged.connect(self.updateEverything)
        self.spectralFWHMLabel = QtWidgets.QLabel('spectral width')

           
        self.controlsLayout.addWidget(self.spectralFWHMLabel)
        self.controlsLayout.addWidget(self.spectralFWHMSpinbox)    
        self.controlsLayout.addItem(vspaceritem)
        
        
        #amp mask controls
        self.ampFilterRadioButton = QtWidgets.QRadioButton()
        self.ampFilterRadioButtonLabel = QtWidgets.QLabel()
        self.ampFilterRadioButtonLabel.setText('apply amp mask?')
        self.ampMaskWidthlabel = QtWidgets.QLabel('amp. mask width')
        self.ampMaskWidthSpinBox = QtWidgets.QDoubleSpinBox()
        self.ampMaskWidthSpinBox.valueChanged.connect(self.updateEverything)
        self.ampFilterRadioButton.toggled.connect(self.toggleAmpMask)
        self.ampMaskLeftSpinBox = QtWidgets.QDoubleSpinBox()
        self.ampMaskLeftSpinBox.setMinimum(250)
        self.ampMaskLeftSpinBox.setMaximum(300)
        self.ampMaskLeftSpinBox.setValue(np.min(self.wl)+1)
        self.ampMaskLeftSpinBox.valueChanged.connect(self.updateEverything)
        self.ampMaskLeftLabel = QtWidgets.QLabel('mask left limit')        
        self.ampMaskRightSpinBox = QtWidgets.QDoubleSpinBox()
        self.ampMaskRightSpinBox.setMinimum(250)
        self.ampMaskRightSpinBox.setMaximum(300)
        self.ampMaskRightSpinBox.setValue(np.max(self.wl)-1)
        self.ampMaskRightSpinBox.valueChanged.connect(self.updateEverything)
        self.ampMaskRightLabel = QtWidgets.QLabel('mask right limit')
        self.ampMaskTukeyParamSpinBox = QtWidgets.QDoubleSpinBox()
        self.ampMaskTukeyParamSpinBox.setMinimum(0)
        self.ampMaskTukeyParamSpinBox.setMaximum(1)
        self.ampMaskTukeyParamSpinBox.setValue(0)
        self.ampMaskTukeyParamSpinBox.valueChanged.connect(self.updateEverything)
        self.ampMaskTukeyParamLabel = QtWidgets.QLabel('Tukey param.')
        self.ampMaskLossLabel = QtWidgets.QLabel()
        self.ampMaskLossLabel.setText(str(self.ampLoss))
        self.controlsLayout.addWidget(self.ampFilterRadioButtonLabel)
        self.controlsLayout.addWidget(self.ampFilterRadioButton)
        #self.controlsLayout.addWidget(self.ampMaskWidthlabel)
        #self.controlsLayout.addWidget(self.ampMaskWidthSpinBox)
        self.controlsLayout.addWidget(self.ampMaskLeftLabel)
        self.controlsLayout.addWidget(self.ampMaskLeftSpinBox)
        self.controlsLayout.addWidget(self.ampMaskRightLabel)
        self.controlsLayout.addWidget(self.ampMaskRightSpinBox)
        self.controlsLayout.addWidget(self.ampMaskTukeyParamLabel)
        self.controlsLayout.addWidget(self.ampMaskTukeyParamSpinBox)  
        self.controlsLayout.addWidget(self.ampMaskLossLabel)      
        self.controlsLayout.addItem(vspaceritem)
                
        #GA controls
        self.popSizeSpinboxLabel = QtWidgets.QLabel("Population size")
        self.popSizeSpinbox = QtWidgets.QSpinBox()
        self.popSizeSpinbox.setValue(self.populationSize)
        self.popSizeSpinbox.valueChanged.connect(self.setPopSize)     
        
        self.numberOfGenesSpinboxLabel = QtWidgets.QLabel("Number of genes")
        self.numberOfGenesSpinbox = QtWidgets.QSpinBox()
        self.numberOfGenesSpinbox.setValue(self.numberOfGenes)
        self.numberOfGenesSpinbox.valueChanged.connect(self.setNumberOfGenes)
                
        self.geneAmpSpinboxLabel = QtWidgets.QLabel("Gene amp")
        self.geneAmpSpinbox = QtWidgets.QDoubleSpinBox()
        self.geneAmpSpinbox.setValue(self.genAmp)
        self.geneAmpSpinbox.valueChanged.connect(self.setGeneAmp)       
        
        self.GAlayout.addWidget(self.popSizeSpinboxLabel)
        self.GAlayout.addWidget(self.popSizeSpinbox)
        self.GAlayout.addWidget(self.numberOfGenesSpinboxLabel)
        self.GAlayout.addWidget(self.numberOfGenesSpinbox)
        self.GAlayout.addWidget(self.geneAmpSpinboxLabel)
        self.GAlayout.addWidget(self.geneAmpSpinbox)   
        self.GAlayout.addItem(vspaceritem)     
        #self.GAlayout.addItem(QtWidgets.QSpacerItem(105,QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding))

        self.GAlayout.addWidget(QtWidgets.QLabel("Elitism coef."))
        self.elitismCoefSpinbox = QtWidgets.QDoubleSpinBox()
        self.elitismCoefSpinbox.setValue(self.elitismCoef)
        self.elitismCoefSpinbox.setDecimals(2)
        self.elitismCoefSpinbox.setMaximum(1)
        self.elitismCoefSpinbox.setMinimum(0)
        self.elitismCoefSpinbox.valueChanged.connect(self.setElitismCoef)
        self.GAlayout.addWidget(self.elitismCoefSpinbox)
        
        self.GAlayout.addWidget(QtWidgets.QLabel("Mutation coef."))
        self.mutationCoefSpinbox = QtWidgets.QDoubleSpinBox()
        self.mutationCoefSpinbox.setValue(self.mutationCoef)
        self.mutationCoefSpinbox.setDecimals(2)
        self.mutationCoefSpinbox.setMaximum(1)
        self.mutationCoefSpinbox.setMinimum(0)
        self.mutationCoefSpinbox.valueChanged.connect(self.setMutationCoef)
        self.GAlayout.addWidget(self.mutationCoefSpinbox)
        
        self.GAlayout.addWidget(QtWidgets.QLabel("Crossover coef."))
        self.crossOverCoefSpinbox = QtWidgets.QDoubleSpinBox()
        self.crossOverCoefSpinbox.setValue(self.twoPtXoverCoef)        
        self.crossOverCoefSpinbox.setDecimals(2)
        self.crossOverCoefSpinbox.setMaximum(1)
        self.crossOverCoefSpinbox.setMinimum(0)
        self.crossOverCoefSpinbox.valueChanged.connect(self.setCrossOverCoef)
        self.GAlayout.addWidget(self.crossOverCoefSpinbox)
        self.GAlayout.addItem(vspaceritem)     
        #self.GAlayout.addItem(QtWidgets.QSpacerItem(105,QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding))
 
 
        
        
        
        vspaceritem = QtWidgets.QSpacerItem(100,100,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding)
 
        
        #bias phase
        self.applyBiasPhaseButton = QtWidgets.QCheckBox()
        self.applyBiasPhaseButtonLabel = QtWidgets.QLabel()
        self.applyBiasPhaseButtonLabel.setText('apply bias phase?')
        self.GAlayout.addWidget(self.applyBiasPhaseButtonLabel)
        self.GAlayout.addWidget(self.applyBiasPhaseButton)
        self.biasPhaseCoefSpinbox=QtWidgets.QDoubleSpinBox()
        self.biasPhaseCoefSpinbox.setDecimals(2)
        self.biasPhaseCoefSpinbox.setMinimum(-50)
        self.biasPhaseCoefSpinbox.setMaximum(50)
        self.biasPhaseCoefSpinbox.valueChanged.connect(self.updateEverything)
        self.biasPhaseCoefSpinboxLabel=QtWidgets.QLabel()
        self.GAlayout.addWidget(self.biasPhaseCoefSpinboxLabel)
        self.biasPhaseCoefSpinboxLabel.setText('bias phase quadratic coef')
        self.GAlayout.addWidget(self.biasPhaseCoefSpinbox)
        self.biasPhaseCoefSpinbox.valueChanged.connect(self.updateEverything)
        self.GAlayout.addItem(vspaceritem)     
        

        self.startButton = QtWidgets.QPushButton('Start')
        self.startButton.clicked.connect(self.startScan)
        self.stopButton = QtWidgets.QPushButton('Stop')
        self.stopButton.clicked.connect(self.stopScan)
        self.GAlayout.addWidget(QtWidgets.QLabel("Start scan"))
        self.GAlayout.addWidget(self.startButton)
        self.GAlayout.addWidget(QtWidgets.QLabel("Stop scan"))
        self.GAlayout.addWidget(self.stopButton)
        self.GAlayout.addItem(QtWidgets.QSpacerItem(105,QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding))

        self.currentGenLabel = QtWidgets.QLabel()
        
        #GA plots layout
        self.fitnessPlotWidget = pq.PlotWidget()
        self.fitnessPlotMax = self.fitnessPlotWidget.plot()
        #self.fitnessPlotMax.setPen('w')
        self.fitnessPlotAvg = self.fitnessPlotWidget.plot()
        #self.fitnessPlotAvg.setPen('r')   
        
        self.genePlotWidget = pq.PlotWidget()
        self.genePlot = self.genePlotWidget.plot()
                
        self.fitnessLayout.addWidget(self.fitnessPlotWidget)     
        self.fitnessLayout.addWidget(self.genePlotWidget)
        
        ##
        self.layout.addLayout(self.plotLayout)    
        self.layout.addLayout(self.controlsLayout)  
        self.layout.addLayout(self.GAlayout) 
        self.layout.addLayout(self.fitnessLayout)
                  
#         self.plotWidget3 = pq.PlotWidget(useOpenGL=True)                
#         self.plotWidget3.setSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Maximum)
#         self.plot31 = self.plotWidget3.plot()
#         self.plot31.setData(self.arrayThreshold(self.Intt, self.t)[1],self.arrayThreshold(self.Intt, self.t)[0])
#         self.plot32 = self.plotWidget3.plot()
#         self.plot32.setData(self.t,self.Intt,color='r')
#         self.plot32.setPen('r')
        #self.layout.addWidget(self.plotWidget3) 
        
class PlotYYWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent=parent)
        self.setupLayout()
        
    def setupLayout(self):
        
        self.plotYYWidget = pq.PlotWidget(useOpenGL=True)
        self.plot11 = self.plotYYWidget.plot()
        print dir(self.plot11)
        self.plot11.setPen((10, 200, 70))
        self.plot12 = self.plotYYWidget.plot()
        self.plot12.setPen((255,255,0))
        self.plot13 = self.plotYYWidget.plot()
        plot13pen = pq.mkPen(color=(10, 200, 70), width=1, style=QtCore.Qt.DashLine)
        self.plot13.setPen(plot13pen)
        self.plot14 = self.plotYYWidget.plot()
        plot14pen = pq.mkPen(color=(10, 200, 70), width=1, style=QtCore.Qt.DotLine)
        self.plot14.setPen(plot14pen)
        self.plot15 = self.plotYYWidget.plot()
        plot15pen = pq.mkPen(color=(64, 224, 208), width=1, style=QtCore.Qt.DotLine)
        self.plot15.setPen(plot15pen)

        self.plotYYItem = self.plotYYWidget.plotItem
        self.plotYYItem.setLabels(left='abs(El)')
        self.ViewboxYY = pq.ViewBox()
        self.plotYYItem.showAxis('right')
        self.plotYYItem.scene().addItem(self.ViewboxYY)
        self.plotYYItem.getAxis('right').linkToView(self.ViewboxYY)
        self.ViewboxYY.setXLink(self.plotYYItem)
        self.plotYYItem.getAxis('right').setLabel('Phase / rad')
        self.plot21 = pq.PlotCurveItem()
        self.plot21.setPen((200, 70, 10))
        self.ViewboxYY.addItem(self.plot21)
        self.plot22 = pq.PlotCurveItem()
        self.ViewboxYY.addItem(self.plot22)
        plot22pen = pq.mkPen('r', width=1, style=QtCore.Qt.DashLine)
        self.plot22.setPen(plot22pen)
        self.ViewboxYY.addItem(self.plot22)

        #self.plotYYItem.vb.sigResized.connect(self.updateElPlotView)
        self.ViewboxYY.setGeometry(self.plotYYItem.vb.sceneBoundingRect())
        self.ViewboxYY.linkedViewChanged(self.plotYYItem.vb, self.ViewboxYY.XAxis)
        self.plotYYWidget.setAntialiasing(True)
        self.plotYYWidget.showGrid(True, False)
        self.plotYYWidget.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.update()
        
        self.layout = QtWidgets.QHBoxLayout(self)
        self.layout.addWidget(self.plotYYWidget)
        
class mySlider(QtWidgets.QSlider, QtWidgets.QLabel):
    #doesn't work yet!!!
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent=parent)
        self.setupLayout()
        
    def setupLayout(self):
        self.slider = QtWidgets.QSlider()
        self.slider.setGeometry(10, 10, 10, 10)
        self.slider.setTickInterval(1)
        self.slider.setMaximum(50)
        self.slider.setMinimum(-50)
        #self.slider.setMinimumWidth(40)
        self.slider.setSizePolicy(QtWidgets.QSizePolicy.Minimum,QtWidgets.QSizePolicy.Minimum)
        self.slider.setTickPosition(QtWidgets.QSlider.TicksLeft)
        
        self.label = QtWidgets.QLabel()
        self.label.setText('test label')
        
        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.slider)
        self.layout.addWidget(self.label)
        
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    #app.processEvents()
    myapp = SpectrometerCamera()
    myapp.show()
    #splash.finish(myapp)
    sys.exit(app.exec_())
