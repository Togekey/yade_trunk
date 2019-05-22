# domain decomposition : Orthogonal Recursive Bisection Algorithm : both serial and parallel available 
# deepak kunhappan : deepak.kunhappan@3sr-grenoble.fr, deepak.kn1990@gmail.com
# #decomp fibers?!
from __future__ import division 
from math import log10, log
from mpi4py import MPI 
from yade.wrapper import * 
from math import log, ceil 

# init mpis (if not initialized) 
#TODO: 
#add driver functions for parallel and serial. 
#some more functions to be added in serail decomposition -> setting of subdomains directly.  
# splitScene function friendly.  
# Testing. 


comm = MPI.COMM_WORLD 
rank = comm.Get_rank() 
threads = comm.Get_size() 
numLevels = int (ceil (log (numthreads))) + 1 
TAG_DATA=999   # tag for sending/recvng 



class decompBodiesParallel: 
	
	def __init__(self, ndim): 
		self.numLevels = numLevels 
		self.numThreads = numThreads 
		self.ndim = ndim 
		self.splitAxisInit = 0 
		self.levelCount = None 
	  
	  
	def scatterData(self, in_data): 
		sendData = [];
		prcData = in_data 
		for nlevel in range(self.numLevels): 
		  self.levelCount = nlevel 
		  toRecv, sourceRank = self.ifScatterRecv(nlevel, rank)
		  if toRecv: 
			prcData = comm.recv(source=sourceRank,tag=TAG_DATA)
		  toSend,destRank = self.ifScatterSend(nlevel,rank,self.numThreads)
		  if toSend: 
			prcData, sendData = self.processData(prcData,nlevel)
			comm.send(sendData, dest=destRank, tag=TAG_DATA)
			
		return prcData 
	      
	
	
	def ifScatterRecv(self, level, inRank): 
		myRank = inRank+1 
		levelPw2 = 1 << level; levelPw2_1 = 1 << (level+1) 
		if (myRank >= levelPw2) and (myRank < levelPw2_1): 
			sourceRank = myRank - levelPw2
			return True, sourceRank
		else:
			return False, -1 
		      
	
	def ifScatterSend(self, inRank, numThreads): 
		myRank = inRank-1
		levelPw2 =  1 << level 
		if myRank < levelPw2: 
			destRank = myRank + levelPw2 
			if destRank > numThreads: 
				return False, -1
			return True, destRank+1 
		      
		else: 
			return False, -1 
		      
	
	def gatherData(self, inData): 
		if rank == 0: return None
		recvData = [] 
		nlevel = self.numLevels-1
		while nlevel >= 0: 
			toSend, destRank = self.ifGatherSend(nlevel, rank)
			if toSend: 
				comm.send(inData,dest=destRank,tag=TAG_DATA)
			toRecv, sourceRank = self.ifGatherRecv(nlevel, rank, self.numThreads)
			if toRecv:
				recvData = comm.recv(source=SourceRank,tag=TAG_DATA)
				inData += recvData
			nlevel = nlevel-1
		return inData


	def ifGatherSend(self,level,inRank): 
		myRank = inRank -1 
		levelPw2 = 1 << level; levelPw2_1 = 1<<(level+1)
		
		if level==0 and myRank > 1: 
			return False, -1 
		
		elif (myRank >= levelPw2) and (myRank < levelPW2_1): 
			destRank = myRank-levelPw2
			return True, destRank+1 
			
		else:
			return False, -1
	
	def ifGatherRecv(self, level, inRank, numThreads):
		myRank = inRank-1 
		levelPw2 = 1 << level
		if level == 0 and myRank > 0: 
			return False, -1 
		elif (myRank < levelPw2): 
			sourceRank = myRank + levelPw2 
			if sourceRank >= numThreads:
				return False, -1 
		else:
			return False, -1
		      
		      
	def processData(self, inData, level): 
		if self.levelCount == 0: splitAxis = splitAxisInit 
		splitAxis = level%self.ndim 
		sortedData  = self.sortData(inData, splitAxis)
		splitPt = len(sortedData)//2 
		inData = sortedData[:splitPt]; sendData = sortedData[splitPt+1:]
		sendData.append(sortedData[splitPt])
		return inData, sendData
	
	
	def processBodies(self,bodyList, level):
		posList = []; 
		for b in bodyList:
			posList.append((b.state.pos, b.id))
		in1, in2 = self.processData(posList,level)
		ids1 = [x[1] for x in in1]; ids2 = [x[1] for x in in2]
		#bodyList = [O.bodies[idx] for idx=x[1] for x in in1]
		bodyList = [O.bodies[idx] for idx in ids1]
		partList = [O.bodies[idx] for idx in ids2]
		return bodyList, partList 
	
	def sortData(self, inData, splitAxis): 
		return sorted(inData, key=lambda x: x[0][splitAxis])
	      
	      

class decompBodiesSerial: 
  

	class treeNode:
		def __init__(self, data): 
			self.data = data 
			self.left = None 
			self.right = None 
			self.level = None 
			self.splitPoint = None 

	class partTree:
		
		def __init__(self, npart, ndim, data):
			self.root = treenode(data)
			self.npart = npart 
			self.ndim = ndim 
			self.numLevel = int (log10(self.npart)/log10(2))
			
			
		def buildTree(self): 
			self._buildTree(self.root, depth=0)
			
			
		def _buildTree(self,curNode,depth): 
			if curNode.level == self.numLevel
				for x in curNode.data: 
					O.bodies[x[1].subdomain = colScale[self.ncount-1]]
				self.count += 1 
				return None 
			splitAxis = None 
			if depth==0: splitAxis = splitAxisInit 
			splitAxis = depth%self.ndim 
			sortedData = sorted(curNode.data, key = lambda x:x[0][splitAxis])
			n = len(sortedData)//2 
			
			
			if curNode.left == None: 
				curNode.left = treeNode[sortedData[:n]]
				curNode.left.splitPoint = sortedData[n]
				curNode.left.level = depth+1
				curNode.left.data.append(curNode.left.splitPoint)
				self._buildTree(curNode.left, depth+1)
			
			if curNode.right == None: 
				curNode.right = treeNode(sortData[n+1:])
				curNode.right.splitPoint = sortedData[n] 
				curNode.right.level = depth+1
				curNode.right.data.append(curNode.right.splitPoint)
				self._buildTree(curNode.right, depth+1)
				
				
	def __init__(self, numSubdomains): 
			pass 