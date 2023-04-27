import random
import gene as gene
from operator import attrgetter
import copy

class GeneParameter:

    def __init__(self, poolSize = 128, maxIteration = 20 , crossoverRate = 0.5, ratioOfCrossover = 0.2,
                 mutationRate = 0.7, ratioOfMutation = 0.1, bestGene = None):
        self.poolSize = poolSize
        self.maxIteration = maxIteration
        self.crossoverRate = crossoverRate
        self.ratioOfCrossover = ratioOfCrossover
        self.mutationRate = mutationRate
        self.ratioOfMutation = ratioOfMutation
        self.geneList = []
        self.genePoolList = []
        self.bestGene = gene.Gene()
        self.times = 1
        if bestGene is None:  # all is random
            self.bestGene.f = 1e9
            for i in range(0, self.poolSize):
                g = gene.Gene()
                g.generate()
                self.geneList.append(g)
        else: # 1/4 best's clone, 3/4 random
            self.bestGene = bestGene.clone()
            quarter = self.poolSize / 4
            for i in range(0, int(quarter)):
                g = bestGene.clone()
                self.geneList.append(g)
            for i in range(int(quarter), poolSize):
                g = gene.Gene()
                g.generate()
                self.geneList.append(g)

    def geneIteration(self, input, output):
        for i in range(0, self.maxIteration): # 20 times
            g = self.performGene(input, output)
        return g

    def performGene(self, input, output):
        for i in range(0, self.poolSize):
            self.geneList[i].calculateFitness(input, output)

        self.geneList.sort(key=attrgetter('f'))  # 照適應函數數值進行排序(適應值越小越好)
        if self.geneList[0].f < self.bestGene.f:
            self.bestGene = self.geneList[0].clone()

        if self.times >= self.maxIteration:
            print("bestGene's f = ", self.bestGene.f)
            return self.bestGene

        # reproduction
        self.genePoolList = []
        bestSize = int(self.poolSize / 10)
        for i in range(0, bestSize):
            g = self.geneList[i].clone()
            self.genePoolList.append(g)
        for i in range(bestSize, self.poolSize):
            picked = random.sample(range(self.poolSize), 2)
            g = self.reproduct(self.geneList, picked)
            self.genePoolList.append(g)

        # crossover
        picked = random.sample(range(self.poolSize), self.poolSize) # shuffle
        for i in range(0, self.poolSize, 2):
            if random.uniform(0, 1) < self.crossoverRate:
                self.genePoolList[picked[i]], self.genePoolList[picked[i+1]] = \
                    self.crossover(self.genePoolList[picked[i]], self.genePoolList[picked[i+1]])

        # mutate
        for i in range(0, self.poolSize):
            if random.uniform(0, 1) < self.mutationRate:
                for j in range(0, self.genePoolList[0].DNALength):
                    self.genePoolList[i].DNA[j] = self.mutate(self.genePoolList[i].DNA[j])
            self.genePoolList[i].setrbf()

        self.geneList = copy.deepcopy(self.genePoolList)
        self.times += 1
        return 0

    def reproduct(self, geneList, pickedList):
        g = gene.Gene()
        g.f = 1e9
        for i in pickedList:
            if geneList[i].f < g.f:
                g = geneList[i].clone()
        return g

    def crossover(self, nowx1, nowx2):
        ratio = random.uniform(-1, 1) * self.ratioOfCrossover
        gx1 = gene.Gene()
        gx2 = gene.Gene()
        for i in range(0, len(nowx1.DNA)):
            gx1.DNA[i] = nowx1.DNA[i] + ratio * (nowx1.DNA[i] - nowx2.DNA[i])
            gx2.DNA[i] = nowx2.DNA[i] - ratio * (nowx1.DNA[i] - nowx2.DNA[i])
        nowx1 = gx1
        nowx2 = gx2
        return nowx1, nowx2

    def mutate(self, DNAValue):
        if random.uniform(0, 1) < self.mutationRate:
            ratio = random.uniform(-1, 1) * self.ratioOfMutation
            DNAValue += ratio * DNAValue
        return DNAValue
