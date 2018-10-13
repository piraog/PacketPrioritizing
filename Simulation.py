import numpy as np
import matplotlib.pyplot as plt
from Network import *

class Simulation:
    '''This class runs scheduling algorithms for emergency flows over a given network.
    Args:
            * network (Network): a network we want to run the simulation on
        Kwargs:
            * eFlowDeadLine (int): the deadline of the simulated emergency flow
            * nExp (int): number of simuations we want to run
            * runType (str): specify if we want to run each simulation on every node or just at 1 starting point for the emergency flow
        
        Attributes:
            * schedulability : an array of the schedulability ratio for each algorithm and each experiment
            * nStolenFlows : an array of the average number of flow stolen by each algorithm in each experiment: SFSA, OBSSA
            * nStolenFlowsScheduled : an array of the average number of flow stolen by each algorithm in each experiment when the packet was scheduled: SFSA, OBSSA
            
    '''
    def __init__(self, network, eFlowDeadLine=4, nExp=100, runType='average on all nodes' ):
        self.n = network
        self.nExp = nExp
        self.eFlowDeadLine = eFlowDeadLine
        self.nStolenFlows = np.array([[0 for _ in range(2)] for i in range(nExp)])
        self.nStolenFlowsScheduled = np.array([[0 for _ in range(2)] for i in range(nExp)])
        self.schedulability = np.array([[0 for _ in range(2)] for i in range(nExp)])
        if runType == 'average on all nodes':
            self.runOnAllNodes()
        elif runType == 'simple':
            self.run()
        elif runType == 'average on all starting points':
            self.runOnAllStartpoints()

        self.schedulability = np.around(100*self.schedulability,decimals=0)


    def runOnAllNodes(self):
        '''Runs both algorithms using, in turn, each connected node of the network as starting point. Then average the results.'''                                    
        for i in range(self.nExp):                              
            self.n = Network(self.n.nNodes, self.n.nFlows, self.n.nChannels, self.n.eFork, self.n.schedulingMethod)
            activeNodes = [node for node,w in enumerate(self.n.W) if (w!={} and node!=0)]
            for node in activeNodes:
                SFSA = self.runSFSA(node)
                OBSSA = self.runOBSSA(node)
                self.schedulability[i][0] += SFSA[0]
                self.nStolenFlows[i][0] += SFSA[1]
                self.nStolenFlowsScheduled[i][0] += SFSA[0]*SFSA[1]
                self.schedulability[i][1] += OBSSA[0]
                self.nStolenFlows[i][1] += OBSSA[1]
                self.nStolenFlowsScheduled[i][1] += OBSSA[0]*OBSSA[1]
                # if OBSSA[0]<SFSA[0]:
                #     print(SFSA)
                #     print(OBSSA)
                #     print(self.n.intersections)      
                #     for k,v in enumerate(self.n.W):
                #         print(k, " ", v)

        for i in range(self.nExp):
            if self.schedulability[i][0]!=0:
                self.nStolenFlowsScheduled[i][0] = self.nStolenFlowsScheduled[i][0]/self.schedulability[i][0]
            else:
                self.nStolenFlowsScheduled[i][0] = 0
            if self.schedulability[i][1]!=0:
                self.nStolenFlowsScheduled[i][1] = self.nStolenFlowsScheduled[i][1]/self.schedulability[i][1]
            else:
                self.nStolenFlowsScheduled[i][1] = 0
        self.nStolenFlowsScheduled = np.true_divide(self.nStolenFlowsScheduled.sum(0),(self.nStolenFlowsScheduled!=0).sum(0))
        self.schedulability = np.mean(self.schedulability, 0)
        self.schedulability /= (len(activeNodes))          #*(self.n.globalPeriod-self.eFlowDeadLine)
        self.nStolenFlows = np.mean(self.nStolenFlows, 0)
        self.nStolenFlows /= (len(activeNodes))            #*(self.n.globalPeriod-self.eFlowDeadLine)
    
    def runOnAllStartpoints(self):
        '''Runs both algorithms using, in turn, each node at the begining of a regular flow as starting point. Then average the results.'''                                     
        for i in range(self.nExp):                              
            self.n = Network(self.n.nNodes, self.n.nFlows, self.n.nChannels, self.n.eFork, self.n.schedulingMethod)
            startpoints = [path[0] for path in self.n.H]
            for node in startpoints:
                SFSA = self.runSFSA(node)
                OBSSA = self.runOBSSA(node)
                self.schedulability[i][0] += SFSA[0]
                self.nStolenFlows[i][0] += SFSA[1]
                self.nStolenFlowsScheduled[i][0] += SFSA[0]*SFSA[1]
                self.schedulability[i][1] += OBSSA[0]
                self.nStolenFlows[i][1] += OBSSA[1]
                self.nStolenFlowsScheduled[i][1] += OBSSA[0]*OBSSA[1]
                # if OBSSA[0]<SFSA[0]:
                #     print(SFSA)
                #     print(OBSSA)
                #     print(self.n.intersections)      
                #     for k,v in enumerate(self.n.W):
                #         print(k, " ", v)

        for i in range(self.nExp):
            if self.schedulability[i][0]!=0:
                self.nStolenFlowsScheduled[i][0] = self.nStolenFlowsScheduled[i][0]/self.schedulability[i][0]
            else:
                self.nStolenFlowsScheduled[i][0] = 0
            if self.schedulability[i][1]!=0:
                self.nStolenFlowsScheduled[i][1] = self.nStolenFlowsScheduled[i][1]/self.schedulability[i][1]
            else:
                self.nStolenFlowsScheduled[i][1] = 0
        self.nStolenFlowsScheduled = np.true_divide(self.nStolenFlowsScheduled.sum(0),(self.nStolenFlowsScheduled!=0).sum(0))
        self.schedulability = np.mean(self.schedulability, 0)
        self.schedulability /= (len(startpoints))          #*(self.n.globalPeriod-self.eFlowDeadLine)
        self.nStolenFlows = np.mean(self.nStolenFlows, 0)
        self.nStolenFlows /= (len(startpoints))            #*(self.n.globalPeriod-self.eFlowDeadLine)

    def run(self):
        '''Runs both algorithm on 1 random connected node of the network.'''                        
        for i in range(self.nExp):
            self.n = Network(self.n.nNodes, self.n.nFlows, self.n.nChannels, self.n.eFork, self.n.schedulingMethod)
            origine = random.sample([node for node,w in enumerate(self.n.W) if (w!={} and node!=0)],1)[0]
            # print('-------------------------------------------------')
            SFSA = self.runSFSA(origine)
            OBSSA = self.runOBSSA(origine)
            # print(SFSA)
            # print(OBSSA)
            # for k,v in enumerate(self.n.W):
            #     print(k, " ", v)
            self.schedulability[i][0] += SFSA[0]
            self.nStolenFlows[i][0] += SFSA[1]
            self.schedulability[i][1] += OBSSA[0]
            self.nStolenFlows[i][1] += OBSSA[1]
            

        self.nStolenFlowsScheduled = np.multiply(self.schedulability,self.nStolenFlows)
        self.nStolenFlowsScheduled = np.true_divide(self.nStolenFlowsScheduled.sum(0),(self.nStolenFlowsScheduled!=0).sum(0))
        self.schedulability = np.mean(self.schedulability, 0)
        self.nStolenFlows = np.mean(self.nStolenFlows, 0)

    def runSFSA(self, origine="random", initialTime=0):
        '''Our implementation of SFSA, origine is used to specify the origine node (int) if we want to study a specific transmission. initialTime can be used to offset the departure of the emergency flow (the deadline is relative to this departure time).'''
        if origine == 'random':
            origine = random.sample([node for node,w in enumerate(self.n.W) if (w!={} and node!=0)],1)[0]  #connected origine
        route = [(initialTime, origine)]                     #list of time, node the emergency nodes cross
        D = self.eFlowDeadLine + initialTime
        scheduability = 1
        reachedGateway = 0
        stolenFlows = []
        while scheduability and not reachedGateway:
            #print('route ', route, '  time ', time, '  args ', route[-1][1], ' ', time)
            nextSlot = self.nextSlot(route[-1][1], route[-1][0], D)                 #select the next available slot to live the node
            if nextSlot != False:
                arrivalTime = nextSlot[0]+1
                nextNode = nextSlot[1][0]
                route += [(arrivalTime, nextNode)]
                stolenFlows += [nextSlot[1][1]]
                if nextSlot[1][0]==0:
                    reachedGateway = 1
            else:
                #print('last node ', self.n.W[route[-1][1]])
                scheduability = 0
        return (scheduability, len(set(stolenFlows)), route)            #return the schedulability, 
                                                                        #nb of stolen flows and route followed


    def nextSlot(self, node, time, deadLine):                 #select the next flow exiting a node, if it is within deadlines
        slots = self.n.W[node]
        availableSlots = [slot for slot in slots.items() if (deadLine-slot[0]-1)>=0 and slot[0]>=time]
        if len(availableSlots)>0:
            nextSlot = min(availableSlots[:])
        else:
            return False
        return nextSlot


    def runOBSSA(self, origine = 'random', initialTime=0):
        '''Our implementation of OBSSA, origine is used to specify the origine node (int) if we want to study a specific transmission. initialTime can be used to offset the departure of the emergency flow (the deadline is relative to this departure time).'''
        if origine == 'random':
            origine = random.sample([node for node,w in enumerate(self.n.W) if (w!={} and node!=0)],1)[0]  #connected origine
        route = [(initialTime, origine)]                     #list of time, node the emergency nodes cross
        D = self.eFlowDeadLine + initialTime
        schedulability = 1
        reachedGateway = 0
        stolenFlows = []
        while schedulability and not reachedGateway:
            if stolenFlows!=[] and self.n.A[route[-1][1]][stolenFlows[-1]][2]+1<self.eFlowDeadLine: #schedulable on the current flow
                reachedGateway = 1
            elif self.goodAvailableStolenFlow(route[-1][0], route[-1][1], stolenFlows)!=[]: #schedulable on a previously stolen flow
                reachedGateway = 1
            elif  self.optimalSlot(route[-1][0], route[-1][1])!='false':    #finding the optimal path, if it exists
                flow = self.optimalSlot(route[-1][0], route[-1][1])         #adding the optimal slot to our route
                arrivalTime = self.n.A[route[-1][1]][flow][0]+1
                destination = self.n.A[route[-1][1]][flow][1]
                route += [(arrivalTime,destination)]
                stolenFlows.append(flow)
                if destination == 0:
                    reachedGateway=1
            else:
                schedulability = 0
        return (schedulability,len(set(stolenFlows)),route)

    def goodAvailableStolenFlow(self, time, node, stolenFlows):
        flows = list(set(stolenFlows)&set([flow for flow,info in self.n.A[node].items() if (info[0]>=time and self.eFlowDeadLine>=info[2])]))
        if flows!=[]:
            flow = flows[0]
            return [flow]
        return []

    def optimalSlot(self,time,node):
        criterion = []
        for flow,info in self.n.A[node].items():
            N = info[0] - time
            if N>=0 and info[0]<self.eFlowDeadLine:
                a = len([1 for i in self.n.intersections[flow][1] if i[0]<node])
                criterion.append((flow,N/(a+1),N))
        criterion.sort(key=lambda t:t[2])
        if criterion!=[]:
            optimalFlow = min(criterion,key=lambda t:t[1])[0]
            return optimalFlow
        return 'false'
                
            


if __name__ == '__main__':



    '''comparison of OBSSA and SFSA, 1000 exp, averaged for each exp on every edge node'''
    # network = Network(300, 15, nChannels=2, eFork=5, schedulingMethod='FCFS')
    # simulation = Simulation(network, nExp = 100, eFlowDeadLine=8,runType='average on all starting points')
    # print(simulation.schedulability, simulation.nStolenFlows, simulation.nStolenFlowsScheduled)
    '''We have the same results in average for randomly generated networks (the cases where OBSSA presents a real advantage are 
    extremely rare, even with a high forking'''


    nNodes = [10*i for i in range(2,11)]
    schedulability1 = []
    stolenFlows1 = []
    stolenFlowsScheduled1 = []
    for n in nNodes:
        network = Network(n, 15, nChannels=2, eFork=5, schedulingMethod='RoundRobin')
        simulation = Simulation(network, nExp = 500, eFlowDeadLine=6)
        schedulability1 += [[simulation.schedulability[0],simulation.schedulability[1]]]
        stolenFlows1 += [[simulation.nStolenFlows[0],simulation.nStolenFlows[1]]]
        stolenFlowsScheduled1 += [[simulation.nStolenFlowsScheduled[0],simulation.nStolenFlowsScheduled[1]]]
    print(schedulability1)
    print(stolenFlowsScheduled1)
    schedulability2 = []
    stolenFlows2 = []
    stolenFlowsScheduled2 = []
    # for n in nNodes:
    #     network = Network(n, 15, nChannels=2, eFork=5, schedulingMethod='RoundRobin')
    #     simulation = Simulation(network, nExp = 1000, eFlowDeadLine=i)
    #     schedulability2 += [[simulation.schedulability[0],simulation.schedulability[1]]]
    #     stolenFlows2 += [[simulation.nStolenFlows[0],simulation.nStolenFlows[1]]]
    #     stolenFlowsScheduled2 += [[simulation.nStolenFlowsScheduled[0],simulation.nStolenFlowsScheduled[1]]]
    # print(schedulability2)

    plt.plot(nNodes, [i[0] for i in schedulability1])
    plt.plot(nNodes, [i[1] for i in schedulability1])
    plt.legend(['SFSA RR','OBSSA RR'])
    plt.title('schedulability versus number of nodes for 15 flows, 2 channels, D=6')
    plt.xlabel('Number of nodes')
    plt.ylabel('Schedulability Ratio')
    plt.show()
    plt.plot(nNodes, [i[0] for i in stolenFlows1])
    plt.plot(nNodes, [i[1] for i in stolenFlows1])
    plt.legend(['SFSA RR','OBSSA RR'])
    plt.title('stolen flows versus number of nodes for 15 flows, 2 channels, D=6')
    plt.xlabel('Number of nodes')
    plt.ylabel('Stolen Flows')
    plt.show()
    plt.plot(nNodes, [i[0]*100.0 for i in stolenFlowsScheduled1])
    plt.plot(nNodes, [i[1]*100.0 for i in stolenFlowsScheduled1])
    plt.legend(['SFSA RR','OBSSA RR'])
    plt.title('stolen flows versus number of nodes for 15 flows, 2 channels, D=6')
    plt.xlabel('Number of nodes')
    plt.ylabel('Stolen Flows for scheduled packets')
    plt.show()
