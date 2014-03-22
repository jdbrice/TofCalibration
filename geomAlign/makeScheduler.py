import sys, os

print str( sys.argv )

if len( sys.argv ) < 4:
	print "makeScheduler.py filename outDir inputDir (filesPerProcess filesPerHour)"
	sys.exit(0)


filename = sys.argv[ 1 ]
outDir = sys.argv[ 2 ]
inputDir = sys.argv[ 3 ]

if len( sys.argv ) >= 5 :
	fpp = sys.argv[ 4 ]
else :
	fpp = 8
if len( sys.argv ) >= 6 :
	fph = sys.argv[ 5 ]
else :
	fph = 2


if not os.path.exists( outDir ) :
	os.mkdir( outDir )
else :
	print "Data directory already exists\n"
	sys.exit(1)
if not os.path.exists( outDir + "/log" ) :
	os.mkdir( outDir + "/log" )
else :
	print "Data directory already exists\n"
	sys.exit(1)
if not os.path.exists( outDir + "/output" ) :
	os.mkdir( outDir + "/output" )
else :
	print "Data directory already exists\n"
	sys.exit(1)
if not os.path.exists( outDir + "/generator" ) :
	os.mkdir( outDir + "/generator" )
else :
	print "Data directory already exists\n"
	sys.exit(1)

f = open( "schedule_" + filename, 'w')

f.write( "<?xml version=\"1.0\" encoding=\"utf-8\" ?>\n" )
f.write( "<job maxFilesPerProcess=\"" + str(fpp) + "\" filesPerHour=\"" + str(fph) + "\">\n" )
f.write( "\t<command>\n")
f.write( "\t\tstarver SL14a\n \t\tpwd\n" )
f.write( "\t\troot4star -b -q " + filename +  ".C\(\"$FILELIST\",\"$JOBID.ntuple.root\"\) &gt;&amp; ${JOBID}.log\n" )
f.write( "\t</command>\n\n")
f.write( "\t<input URL=\"file:" + inputDir + "\"/>\n\n")
f.write( "\t<stdout URL=\"file:" + outDir + "/log/$JOBID.out\" />\n\n" )
f.write( "\t<stderr URL=\"file:" + outDir + "/log/$JOBID.err\" />\n\n" )
f.write( "\t<output fromScratch=\"*\" toURL=\"file:" + outDir + "/output/\" />\n\n" );

f.write( "\t<SandBox>\n \
\t\t<Package>\n \
\t\t\t<File>file:./doMuDstEvents.C</File> \n \
\t\t\t<File>file:./.sl64_gcc447/</File> \n \
\t\t</Package> \n \
\t</SandBox>\n")

f.write( "\t<Generator> \n \
\t\t<Location>" + outDir + "/generator</Location> \n \
\t</Generator> \n")

f.write( "</job>" )

f.close()





