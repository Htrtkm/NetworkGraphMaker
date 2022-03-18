#!/usr/bin/env python

import sys
import re
import math

#information about plasmid scaffold to analyze with Hi-C ~~~~~~~~~~~~~~~~~~~~~
class hgtTargets:
	def __init__(self,hgtScaffoldList):
		self.targetList = set()
		with open(hgtScaffoldList,"r")as data:
			for line in data:
				line_trim = line.strip("\n")
				self.targetList.add(line_trim)

#informaiton about bins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class bins:
	def __init__(self,GTDB_result):
		self.bin_allInfo = {}
		with open(GTDB_result,"r")as data:
			for line in data:
				line_trim = line.strip("\n")
				table = line_trim.split("\t")
				Lineage = lineage(table)
				self.bin_allInfo.setdefault(Lineage.binName,Lineage)

class lineage:
	def __init__(self,lineage_table):
		self.binName = lineage_table[0]
		fullLineage = lineage_table[1].split(";")
		if len(fullLineage) == 1:
			return
		self.Domain = fullLineage[0][3:]
		self.Phylum = fullLineage[1][3:]
		self.Class = fullLineage[2][3:]
		self.Order = fullLineage[3][3:]
		self.Family = fullLineage[4][3:]
		self.Genus = fullLineage[5][3:]
		self.Species = fullLineage[6][3:]

class checkM:
	def __init__(self,checkM_result):
		self.binName_stats = {}
		with open(checkM_result,"r")as data:
			for line in data:
				line_trim = line.strip("\n")
				table = line_trim.split("\t")
				if float(table[8]) < 10 and float(table[7]) > 50:
					self.binName_stats.setdefault("bin" + table[0],[float(table[7]),float(table[8])])

#information about mapping result ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class mappingResult:
	def __init__(self,samFile):
		self.HiCLink_weight = {}
		self.readNumber = 0
		with open(samFile,"r")as data:
			preLine =  []
			for line in data:
				self.readNumber += 1
				line_trim = line.strip("\n")
				table = line_trim.split("\t")
				if preLine == []:
					preLine = table
					continue
				if preLine[0] == table[0]:
					if int(preLine[4]) == 60 and int(table[4]) == 60 and preLine[2] != table[2] and (preLine[2] in HGTtargets.targetList or table[2] in HGTtargets.targetList):
						HiCLink = hicLink(preLine,table)
						if HiCLink.edgeName in self.HiCLink_weight:
							self.HiCLink_weight[HiCLink.edgeName].LinkCount += 1
						self.HiCLink_weight.setdefault(HiCLink.edgeName,HiCLink)
					continue
				if preLine[0] != table[0]:
					preLine = table
					continue
	def normalize_readcount(self,readcountFile):
		self.scaffold_readcount = {}
		with open(readcountFile,"r")as data:
			for line in data:
				line_trim = line.strip("\n")
				table = line_trim.split("\t")
				readcount = int(table[1])
				if readcount == 0:
					readcount = 1
				self.scaffold_readcount.setdefault(table[0],readcount)
		for k,v in self.HiCLink_weight.items():
			normalizeNumber =  self.scaffold_readcount[v.preScaffoldName] * self.scaffold_readcount[v.postScaffoldName]
			loggedWeight = math.log10(v.LinkCount/(normalizeNumber*self.readNumber))
			if loggedWeight < -14.84 or v.LinkCount <= 5:#default->10
				continue
			v.normalizeWeight(self.scaffold_readcount,self.readNumber)
	def extractStrongLink(self,gmlFile,HGTtargets,CheckM):
		self.node_id = {}
		nodeID = 0
		alreadyMentionedNode = set()
		for k,v in self.HiCLink_weight.items():
			if v.preScaffoldName in HGTtargets.targetList and v.postScaffoldName not in HGTtargets.targetList:
				if v.preScaffoldName not in alreadyMentionedNode:
					self.node_id.setdefault(v.preScaffoldName,nodeID)
					alreadyMentionedNode.add(v.preScaffoldName)
					nodeID += 1
				binName = v.postScaffoldName.split("_")[3]
				if binName not in alreadyMentionedNode and binName in CheckM.binName_stats:
					self.node_id.setdefault(binName,nodeID)
					alreadyMentionedNode.add(binName)
					nodeID += 1
			if v.preScaffoldName not in HGTtargets.targetList and v.postScaffoldName in HGTtargets.targetList:
				binName = v.preScaffoldName.split("_")[3]
				if binName not in alreadyMentionedNode and binName in CheckM.binName_stats:
					self.node_id.setdefault(binName,nodeID)
					alreadyMentionedNode.add(binName)
					nodeID += 1
				if v.postScaffoldName not in alreadyMentionedNode:
					self.node_id.setdefault(v.postScaffoldName,nodeID)
					alreadyMentionedNode.add(v.postScaffoldName)
					nodeID += 1
		sourceTarget_normedWeight = {}
		for k,v in self.HiCLink_weight.items():
			if v.preScaffoldName in HGTtargets.targetList and v.postScaffoldName not in HGTtargets.targetList and v.isValid == 1:
				binName = v.postScaffoldName.split("_")[3]
				sourceTarget = v.preScaffoldName + "+" + binName
				if binName in CheckM.binName_stats:
					if sourceTarget in sourceTarget_normedWeight:
						sourceTarget_normedWeight[sourceTarget] += v.normalizedWeight
					sourceTarget_normedWeight.setdefault(sourceTarget,v.normalizedWeight)
			if v.preScaffoldName not in HGTtargets.targetList and v.postScaffoldName in HGTtargets.targetList and v.isValid == 1:
				binName = v.preScaffoldName.split("_")[3]
				sourceTarget = v.postScaffoldName + "+" + binName
				if binName in CheckM.binName_stats:
					if sourceTarget in sourceTarget_normedWeight:
						sourceTarget_normedWeight[sourceTarget] += v.normalizedWeight
					sourceTarget_normedWeight.setdefault(sourceTarget,v.normalizedWeight)
		with open(gmlFile,"w")as f:
			f.write("Creator \"Takumi Hattori\"\ngraph\n[\n")
			nodes = ""
			alreadyMentionedScaffold = set()
			binId = 0
			for k,v in self.node_id.items():
				node = "\tnode\n\t[\n\t\tid " + str(v) + "\n\t\tbin \"" + k + "\"\n"
				if k in Bins.bin_allInfo:
					node += "\t\tLabel" + " \"" + k + "\"\n"
					node += "\t\tplasmid \"No\"\n"
					node += "\t\tDomain" + " \"" + Bins.bin_allInfo[k].Domain + "\"\n"
					node += "\t\tPhylum" + " \"" + Bins.bin_allInfo[k].Phylum + "\"\n"
					node += "\t\tClass" + " \"" + Bins.bin_allInfo[k].Class + "\"\n"
					node += "\t\tOrder" + " \"" + Bins.bin_allInfo[k].Order + "\"\n"
					node += "\t\tFamily" + " \"" + Bins.bin_allInfo[k].Family + "\"\n"
					node += "\t\tGenus" + " \"" + Bins.bin_allInfo[k].Genus + "\"\n"
					node += "\t\tSpecies" + " \"" + Bins.bin_allInfo[k].Species + "\"\n"
					node += "\t]\n"
					nodes += node
				if k not in Bins.bin_allInfo and k in HGTtargets.targetList:
					node += "\t\tLabel" + " \"" + k + "\"\n"
					node += "\t\tplasmid \"Yes\"\n"
					node += "\t\tDomain" + " \"Plasmid\"\n"
					node += "\t\tPhylum" + " \"Plasmid\"\n"
					node += "\t\tClass" + " \"Plasmid\"\n"
					node += "\t\tOrder" + " \"Plasmid\"\n"
					node += "\t\tFamily" + " \"Plasmid\"\n"
					node += "\t\tGenus" + " \"Plasmid\"\n"
					node += "\t\tSpecies" + " \"Plasmid\"\n"
					node += "\t]\n"
					nodes += node
				if k not in Bins.bin_allInfo and k not in HGTtargets.targetList:
					node += "\t\tLabel" + " \"" + k + "\"\n"
					node += "\t\tplasmid \"No\"\n"
					node += "\t\tDomain" + " \"unclassified\"\n"
					node += "\t\tPhylum" + " \"unclassified\"\n"
					node += "\t\tClass" + " \"unclassified\"\n"
					node += "\t\tOrder" + " \"unclassified\"\n"
					node += "\t\tFamily" + " \"unclassified\"\n"
					node += "\t\tGenus" + " \"unclassified\"\n"
					node += "\t\tSpecies" + " \"unclassified\"\n"
					node += "\t]\n"
					nodes += node
			f.write(nodes)
			nodes = ""
			edges = ""
			alreadyMentionedEdge = set()
			for k,v in sourceTarget_normedWeight.items():
				source = k.split("+")[0]
				target = k.split("+")[1]
				if source in HGTtargets.targetList and target in HGTtargets.targetList:
					continue
				edges += "\tedge\n\t[\n\t\tsource " + str(self.node_id[source]) + "\n\t\ttarget " + str(self.node_id[target]) + "\n\t\tweight " + str(v) + "\n\t]\n"
			f.write(edges)
			f.write("]")
			edges = ""
			
			nodeWithEdge = {}
			for k,v in sourceTarget_normedWeight.items():
				source = k.split("+")[0]
				target = k.split("+")[1]
				if source in nodeWithEdge:
					nodeWithEdge[source] += 1
				nodeWithEdge.setdefault(source,1)
				if target in nodeWithEdge:
					nodeWithEdge[target] += 1
				nodeWithEdge.setdefault(target,1)
			nodeWithEdge_sorted = sorted(nodeWithEdge.items(),key = lambda x:-x[1])
			for i in range(len(nodeWithEdge_sorted)):
				print(nodeWithEdge_sorted[i][0] + "\t" + str(nodeWithEdge_sorted[i][1]))

class hicLink:
	def __init__(self,preLine,currentLine):
		self.isValid = 0
		preSCFnumber = int(preLine[2].split("_")[0][8:])
		postSCFnumber = int(currentLine[2].split("_")[0][8:])
		if preSCFnumber < postSCFnumber:
			self.setInformation(preLine,currentLine)
		else:
			self.setInformation(currentLine,preLine)
	def setInformation(self,former,latter):
		self.edgeName = former[2] + "+" + latter[2]
		self.preScaffoldName = former[2]
		self.postScaffoldName = latter[2]
		self.LinkCount = 1
	def normalizeWeight(self,scaffold_readcount,allReadNumber):
		self.isValid = 1
		self.normalizedWeight = int((self.LinkCount/(scaffold_readcount[self.preScaffoldName]*scaffold_readcount[self.postScaffoldName]*allReadNumber))*(10**15))

if __name__ == "__main__":
	args = sys.argv
	if len(args) != 7:
		print("usage:", args[0], "scaffold_list.txt GTDB-Tk.tsv CheckM.tsv read_map.sam read_count.tsv output.gml", file=sys.stderr)
		sys.exit(1)

	HGTtargets = hgtTargets(args[1])
	Bins = bins(args[2])
	CheckM = checkM(args[3])
	MappingResult = mappingResult(args[4])
	MappingResult.normalize_readcount(args[5])
	MappingResult.extractStrongLink(args[6], HGTtargets, CheckM)
"""
args[1] -> scaffold list
args[2] -> GTDB-Tk result
args[3] -> CheckM result
args[4] -> sam file
args[5] -> readcount 
args[6] -> output
"""
