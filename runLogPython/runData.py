import xml.etree.cElementTree as ET

def printData( page ) :
	xml = ET.fromstring( page )
	
	runNumber = int( xml.attrib[ 'number'] )
	#print "Run Number: " + str( runNumber )
	
	for events in xml.iter( "events" ) :
		totalEvents = int( events.text )
	#print "Total Events: " + str( totalEvents )

	# rhic summary
	for fill_id in xml.iter( "fill_id" ) :
		fillID =  fill_id.text 
	
	#print "Fill ID: " + fillID 

	#for detectors in xml.iter( "detectors" ) :
	detectors = xml[ 0 ].find( "detectors").text 
	
	#print "Detectors : " + detectors
	
	#if "tof" in detectors :
 	#	print "TOF ACTIVE"

 	#if "tpx" in detectors :
 	#	print "TPX ACTIVE"
 	vpdEvents = 0
 	vpdEff = 0
	triggers = xml[ 2 ]
	for trigger in triggers.iter( 'trigger' ) : 
		if "VPD_mb" == trigger.get( 'name' ) :
			vpdEvents = int( trigger.get( 'nevents' ) )
			#print "VPD Events: " + str( vpdEvents )
			vpdEff = 0
			if totalEvents != 0 :
				vpdEff = (float(vpdEvents) / float(totalEvents)) * 100
			
				
			#print "Eff: " + str( vpdEff )

	print "%8i\t%8i\t%8s\t%8i\t%8f" % (runNumber, totalEvents, fillID, vpdEvents, vpdEff ) 