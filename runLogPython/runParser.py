import run14List
import runData

for runI in run14List.runNumbers :
	#runI = run14List.runNumbers[ 0 ]
	file = open( "./run14/" + str( runI  ) + ".xml", "r")	
	page = file.read()
	runData.printData( page )
	file.close()
