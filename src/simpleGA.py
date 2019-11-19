'''
Created on 13 maj 2016

@author: Marija
'''
import numpy as np
import matplotlib.pyplot as plt
#import pyqtgraph as pq

class geneticAlgorithm(object):
    
    def __init__(self):
        
        self.populationSize = 50
        self.numberOfGenes = 6
        self.population = []   
        self.elitismCoef = 0.3   
        self.mutationCoef = 0.1
        self.twoPtXoverCoef = 0.6
        
        self.fitness = [None] * self.populationSize   
        self.geneLowerBound = -2.00
        self.geneUpperBound = 2.00
        self.creepRate = 0.1  
                
        self.fitnessCostWeight = 0.25 / self.populationSize   
        self.referencePopulation = [0]*self.numberOfGenes

    def calcIdx(self):
        # returns the indices of new generation individuals to be created through elitism, mutation and x-over 
        self.elitismIdx = np.arange(0, int(max(1,np.floor(self.populationSize * self.elitismCoef))))
        self.mutationIdx = np.arange(int(max(self.elitismIdx))+1, int(max(self.elitismIdx)+1  + np.floor(self.populationSize * self.mutationCoef)))
        self.crossOverIdx = np.arange(max(self.mutationIdx)+1, self.populationSize)
       
    def newGeneration(self):
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
     
    def onePtXover(self,ind1,ind2):
        cxpoint = np.random.randint(1,self.numberOfGenes) 
        #print cxpoint
        ind1[cxpoint:], ind2[cxpoint:] = ind2[cxpoint:], ind1[cxpoint:]    
        return np.array([ind1, ind2])

    def generateInitialPopulation(self,popSize,numberOfGenes):
        self.population = []
        for i in range(0, self.populationSize):
            self.population.append(np.random.uniform(low=self.geneLowerBound, high=self.geneUpperBound, size=self.numberOfGenes))
        self.population = np.array(self.population)
        print self.population.shape

        
    def sortByFitness(self, pop):        
        
        #=======================================================================
        
        self.fitnesses = np.exp(-(pop[:,0]-1)**2 - (pop[:,1]-1)**2) 
        
        #fitness with cost functional
        fitnessCorrection = []
        
        
        for i in range(0,len(pop)):
            fitnessCorrection.append( np.sum((pop[i,:]-self.referencePopulation)**2) )
        fitnessCorrection = np.array(fitnessCorrection)
            
        self.fitnessesWithCost = self.fitnesses  - self.fitnessCostWeight*fitnessCorrection
        
        #- self.fitnessCostWeight * fitnessCorrection
                
        #print pop
        #self.fitnesses = pop[:,1]
        self.avgFitness = np.average(self.fitnessesWithCost)
        self.maxFitness = np.max(self.fitnessesWithCost) 
        j = np.argsort(self.fitnessesWithCost)
        j = np.flipud(j)
        pop = pop[j]
        return pop  
    
        
    def creepMutation(self,individualToMutate):
        #change one gene according to
        #g_new = min( imax, max( imin, g_old + r*s*(imax-imin)
        # r - random number in the [-1,1] range
        # s - creep rate, usually 0.02
        # imin, imax - gene range bounds
        mutationLocation = np.random.randint(0,self.numberOfGenes)
        #print "mutation location is %d" %mutationLocation
        #print individualToMutate
        oldGene = individualToMutate[mutationLocation]
        r = np.random.uniform(low=-1.0, high=1.0, size=1)
        newGene = min( self.geneUpperBound, max( self.geneLowerBound, oldGene + r*self.creepRate*(self.geneUpperBound-self.geneLowerBound)))
        individualToMutate[mutationLocation] = newGene
        return individualToMutate    
    
    #===========================================================================
    # def plotLayout(self): 
    #     self.fitnessWidget = pq.PlotWidget(useOpenGL=True)
    #     self.fitnessPlot = self.fitnessWidget.plot()           
    #===========================================================================

if __name__ == '__main__':
    ga = geneticAlgorithm()
    ga.calcIdx()
    ga.generateInitialPopulation(ga.populationSize, ga.numberOfGenes)
     
    f = []
    gen = []
    for i in range(1,50):
        #ga.evaluateFitness()
        ga.newGeneration()
        f.append(ga.fitnesses)
        gen.append(i)
            
    #print f
    
    print np.sum(ga.population,axis=0)
    print np.std(ga.population,axis=0)
    print ga.population.sum(axis=0)
    print ga.population.std(axis=0)
    
    plt.subplot(2, 1, 1)
    plt.plot(gen,np.amax(f,1),'b', label='Max fitness')
    plt.plot(gen,np.average(f,1),'g', label='Avg fitness')
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    legend = plt.legend(loc='best', shadow=True, fontsize='x-large')
    
    plt.subplot(2, 1, 2)
    plt.errorbar(range(ga.numberOfGenes), np.mean(ga.population,axis=0), yerr = np.std(ga.population,axis=0))
    
    plt.show()  
    
    
    #===========================================================================
    # a = [1,1,1,2,2,2]
    # b = [3,3,3,4,4,4]
    # print ga.onePtXover(a, b)
    #===========================================================================
    

   
    
    #ga.creepMutation(ga.population[1])
    #print ga.population[1]
   
    #print ga.population
        
    #print ga.onePtXover(ga.population[1], ga.population[2])

#===============================================================================
#      for i in range(1,100):       
#          ga.evaluateFitness()
#          f.append(np.max(ga.fitness))
#          ga.newGeneration()    
#    
#     print np.average(ga.fitness)
# 
#     plt.plot(f)
#     plt.xlabel('Generation')
#     plt.ylabel('Fitness')
#     plt.show()
#     
#===============================================================================

        #=======================================================================
        # cutoffCoef = np.floor(len(self.population)*self.elitismCoef);
        # newPopSize = int(self.populationSize - cutoffCoef)
        # newPop = self.sortByFitness(self.population)
        # #newPop1 = self.generateNewPopulation(self.populationSize - cutoffCoef, self.numberOfGenes)
        # newPop1 = self.generateNewPopulation(newPopSize, self.numberOfGenes)
        # self.population = newPop[0:cutoffCoef,:]
        # self.population = np.vstack([self.population,newPop1])          
        #=======================================================================

        #=======================================================================
        # for i in range(0, self.numberOfGenes):
        #     if i<mutationLocation:
        #         indMutated[i] = indToMutate[i]
        #     else:
        #         indMutated[i] = newGene;
        # return indMutated
        #=======================================================================     
        
    #===========================================================================
    # def generateNewPopulation(self,popSize,numberOfGenes):
    #     newPop = []
    #     for i in range(0, popSize):
    #         newPop.append(np.random.uniform(low=0.0, high=1.0, size=self.numberOfGenes))
    #     newPop = np.array(newPop)  
    #     return newPop      
    #===========================================================================
    
    