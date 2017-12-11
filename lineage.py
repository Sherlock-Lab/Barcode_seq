import docopt
import LineageTracking
import pooling_fileDA

if __name__ == "__main__":
	arguments = docopt.docopt("\n".join(open("lineage.txt").readlines()), version="0.1")
	if arguments["align"]:
		LineageTracking.align(arguments)
	if arguments["barcodes"]:
		LineageTracking.barcodes()
	if arguments["demultiplexL"]:
		LineageTracking.demultiplex(arguments)
	if arguments["BC1cluster"]:
		LineageTracking.BC1cluster()
	if arguments["BCtables"]:
		LineageTracking.BCtables()
	if arguments["bartender"]:
		LineageTracking.bartender()
	if arguments["trajectories"]:
		pooling_fileDA.trajectories()