from xml.dom import minidom
import runPageLister
import run14List


for runI in run14List.runNumbers :

	file = open( "./run14/" + str( runI  ) + ".xml", "w")
	file.write( runPageLister.getRunXML( runI ) )
	file.close()


#xml = minidom.parseString( page )
#xml = minidom.parse( urlopen( "http://online.star.bnl.gov/RunLog/index.php?r=15047086&method=xml") )
#triggers = xml.getElementsByTagName("trigger")
#print len( triggers )
