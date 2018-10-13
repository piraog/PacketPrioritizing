import numpy as np
import random

"""
.. module:: Network
.. moduleauthor:: Marc Belicard


"""

class Network:
    '''This class defines a network, and creates one according to given parameters. 
    It will later be used as an input for our simulation.
        Args:
            * nNodes (int): The number of nodes in the network
            * nFlows (int): The number of flows we generate (going from one node to the gateway, they define the static routing/scheduling for regular packets)

        Kwargs:
            * nChannels (int): the number of channels the nodes can use (basicaly the number of transmissions that can occur at the same time between 2 pairs of nodes)
            * eFork (int): the expectancy of the Poisson law used to determine the number of connections exiting a node
            * schedulingMethod (str): the method used to schedule the regular flows (FCFS or Round Robin)
        Attributes:
            * connectionTable : nNode*nNode map of all existing connections between nodes
            * H : list of each flow path  
            * W : list of dictionaries, for each node:  {time of leaving : destination, flowID, arrival}
            * A : list of dictionaries, for each node:  {flowID : time of leaving, destination, arrival}
            * flows : list of the flows parameters, for each flow: [number of hops, arrival time, deadline]
            * transmissionView : list of every transmission occuring at each time step
            * globalPeriod : period of the entire schedule
            * utilisationRate : utilisation rate of the network (ration of occupied transmission slots)
            * connectedNodes : number of nodes actually part of a flow
    '''


    def __init__(self, nNodes, nFlows, nChannels = 2, eFork = 1, schedulingMethod='FCFS'):
        self.nNodes = nNodes
        self.nFlows = nFlows
        self.eFork = eFork
        self.nChannels = nChannels
        self.schedulingMethod = schedulingMethod
        self.H = [[] for _ in range(nFlows)]
        self.connectionTable = [[0 for i in range(nNodes)] for j in range(nNodes)]
        self.W = [{} for _ in range(nNodes)]
        self.flows = [[] for flow in self.H]
        self.addConnections()
        self.addFlows()
        if schedulingMethod == 'FCFS':
            self.scheduleTransmissions()
        else:
            self.scheduleTransmissions2()
        self.listFlows()
        # print(self.H)
        # print(self.W)
        # print(self.globalPeriod)
        self.visualizeTransmissions()
        self.updateW()
        self.createA()
        self.intersectionIdentification()
        self.connectedNodes = len([i for i in self.W if i!={}])

    def addConnections(self):
        '''For each node, randomly add output connections to other nodes (a node can only have an output connection to a lesser node).
        The number of connections is defined by a Poisson law of parameter *eFork*.'''
        for node in range(1, self.nNodes):
            self.addChildren(node)

    def addChildren(self, i):
        nChildren = int(np.random.poisson(self.eFork))
        for child in range(0,max(1, nChildren)):
            self.connectionTable[i][np.random.randint(0,i)] = 1     #change 0 by max(0,i-4) to stay in the same area

    def addFlows(self):
        '''Randomly generate *nFlow* Flows from a random starting node to the gateway, every next node of the flow is randomly selected amongst the possible outputs of the current node.'''
        sourceNodes = random.sample(range(1, self.nNodes), self.nFlows)
        for (flow, sourceNode) in enumerate(sourceNodes):
            self.createFlow(flow, sourceNode)

    def createFlow(self, flow, sourceNode):
        self.H[flow].append(sourceNode)
        node = sourceNode
        while node != 0:
            connexions = [i for (i,j) in enumerate(self.connectionTable[node]) if j==1]
            selectedConnexion = random.sample(connexions, 1)[0]
            self.H[flow].append(selectedConnexion)
            node = selectedConnexion

    def scheduleTransmissions(self):
        '''Schedule the transmissions in the network, knowing the connections of each flow. We do so in a FCFS manner.
        There can be *nChannel* transmissions at the same time (between 2 distinct pair of emitter/receiver.
        We achieve a high utilization rate with this technique, and the average transmission time for a flow is short.'''        
        occupiedNodes = [[0,[]] for _ in range(self.nFlows*self.nNodes)]
        globalPeriod = 0
        for (nFlow, flow) in enumerate(self.H):       #FCFS, farest node first
            i = 0
            for (origine, destination) in zip(flow[:(len(flow)-1)], flow[1:]):
                not_scheduled = 1
                while not_scheduled:
                    if (origine not in occupiedNodes[i][1]) and (destination not in occupiedNodes[i][1])\
                    and occupiedNodes[i][0]<self.nChannels:
                        occupiedNodes[i][1] += [origine, destination]
                        occupiedNodes[i][0] += 1
                        i += 1
                        self.W[origine][i] = [destination, nFlow]
                        if destination==0:                              #store the arrival time of each flow in W[0]
                            self.W[0][nFlow] = i+1
                        globalPeriod = max(globalPeriod, i)
                        not_scheduled = 0
                    else:
                        i += 1
            self.flows[nFlow] += [i]
        self.globalPeriod = globalPeriod           #repete the operation after all flows are scheduled
                                                    #to try to send more packets before globalPeriod
                                                    #decribed by utilization rate (other method called last)

    def scheduleTransmissions2(self):
        '''Schedule the transmissions in the network, knowing the connections of each flow. We do so in a Round Robin manner.
        There can be *nChannel* transmissions at the same time (between 2 distinct pair of emitter/receiver.
        We achieve a high utilization rate with this technique, and the average transmission time for a flow is short.''' 
        occupiedNodes = [[0,[]] for _ in range(self.nFlows*self.nNodes)]
        globalPeriod = 0
        flag = True
        step = [0 for _ in self.H]
        times = [0 for _ in self.H]
        flow = 0
        while flag:
            if len(self.H[flow])>step[flow]+1:
                origine = self.H[flow][step[flow]]
                destination = self.H[flow][step[flow]+1]
                flag2=True
                i=times[flow]
                while flag2:
                    if (origine not in occupiedNodes[i][1]) and (destination not in occupiedNodes[i][1])\
                    and occupiedNodes[i][0]<self.nChannels:
                        occupiedNodes[i][1] += [origine, destination]
                        occupiedNodes[i][0] += 1
                        self.W[origine][i] = [destination, flow]
                        if destination==0:                              #store the arrival time of each flow in W[0]
                            self.W[0][flow] = i+1
                            globalPeriod = max(globalPeriod, i+1)
                            self.flows[flow] += [i+1]
                        step[flow] += 1
                        times[flow] = i+1
                        flag2=False
                    else:
                        i+=1
            flow+=1
            flow = flow%(len(self.H))
            if step == [len(i)-1 for i in self.H]:
                flag = False
        self.globalPeriod = globalPeriod
            #ordre niqué




    def visualizeTransmissions(self):   
        '''Create a list of every transmissions at each time step (maximum *nChannel* transmissions) for visualization purposes.'''
        tab = [[] for _ in range(self.globalPeriod)]
        for (origine, transmissions) in enumerate(self.W[1:]):
            for time, destination in transmissions.items():
                tab[time-1] += [[origine+1, destination[0]]]
        self.transmissionView = tab
        self.utilizationRate = sum([len(slot) for slot in self.transmissionView])/(self.globalPeriod*self.nChannels)


    def listFlows(self):
        '''Creates a list containing each flows parameters (as describes in the article).
        We create a random deadline to each flow in [current arrival time, global period of the schedule].'''
        self.flows = [[len(flow)-1] + self.flows[nFlow] for nFlow, flow in enumerate(self.H)]

        for flow in self.flows:
            flow += [np.random.randint(flow[1], self.globalPeriod+1)]

    def updateW(self):              #add the arrival time of each transmission
        for node in self.W[1:]:
            for time, info in node.items():
                info.append(self.W[0][info[1]])

    def createA(self):
        '''Creates a list of dictionnaries, for each node we have a shuffled view pf W (similar to the A described in the article, but with a different structure)'''

        self.A = [{info[1] :[time, info[0], info[2]] for time, info in node.items()} for n,node in enumerate(self.W) if n!=0]
        self.A = [0] + self.A

    def intersectionIdentification(self):
        '''List of the intersection parameters for each flows, not used in the simulation, but described ni the article.'''
        self.intersections = [[0,[]] for _ in range(self.nFlows)]
        for i,flowi in enumerate(self.H):
            for j,flowj in enumerate(self.H):
                if i!=j:
                    for k,step in enumerate(flowi[:-2]):
                        if (step in flowj) and (flowi[k+1] == flowj[flowj.index(step)+1]):
                            self.intersections[i][0] += 1
                            self.intersections[i][1] += [(step,j)]
        for flow in self.intersections:
            flow[1].sort(key=lambda tup: tup[0],reverse=True)                  #sorting according to the position on the flow


if __name__ == '__main__':

    network = Network(nNodes=50,nFlows=15,nChannels=2,eFork=5,schedulingMethod='FCFS')

    print('Characteristics of a network with the parameters: 50 nodes, 15 flows, 2 channels, in average 5 outputs per node, and using round robin for the regular schedule: \n')
    print('Connections: ',np.array(network.connectionTable))      #list of available connections between the nodes
    print('\nH ',network.H)        #list of the flows (each flow is a list of nodes)
    print('\nW ',network.W)    #list of dictionaries  {time of leaving : destination, flow, arrival} for each node
    print('\nA ', network.A)
    print('\nintersections: ',network.intersections)
    print('\nGlobal period: ',network.globalPeriod)     #period of the longest flow
    print('\nTransmissions: ',network.transmissionView)      #array of he transmissions for every timeslot
    print('\nFlow parameters: ',network.flows)       #list of the flows parameters
    print('Utilization rate of the network: ',network.utilizationRate)
    print('Average number of hops of a flow: ', np.average(np.array([len(flow) for flow in network.H]))) 
    print('Average period of a flow: ',np.average(np.array([time for flow, time in network.W[0].items()])))
    print('Number of nodes used: ',len([i for i in network.W if i!={}]))
    print('\n ^^^ Characteristics of a network with the parameters: 50 nodes, 15 flows, 2 channels, in average 5 outputs per node, and using round robin for the regular schedule\n')
    
