// Version: $Id: lamodel.cpp 172 2014-02-12 10:06:07Z gk $
/* 
 
 
lamodel is a network simulator that simulates memory engram formation 
in a population consisting of excitatory and inhibitory neurons with independent
dendritic subunits. 

This is the entry point for the network simulator. 
This file parses the command line parameters 

Basic command line options for the simulator:

Example: ./lamodel -P 2 -T 1400 [-G] [-L] [-s datadir-name] [-S random-seed]

-P <patterns>: number of memories  to enoode
-T <minutes>: Interval between memories (minutes)
-G: Global-only protein synthesis 
-L: Dendritic-only protein synthesis 
-s <dirname>: Name of the directory to store the data (inside the data/ folder)
-S <random-seed>: Intitialize the random generator

-o <param-name>=<param-value>: Set various simulation parameters:

-o nlTypePV=[0,1,2,3]   : Set basket cell dendrite nonlinearity type. 0=supra, 1=sub, 2=linear, 3=mixed supra/sub
-o nlTypeSOM=[0,1,2,3]  : Set SOM+ cell dendrite nonlinearity type. 0=supra, 1=sub, 2=linear, 3=mixed supra/sub
-o INClustered=[0,1]    : Set whether IN synapses  should target all dendrites randomly (dispersed) or only 33% of them (clustered)


*/


#include "constructs.h"
#include <iostream>
#include <cstring>
#include <string>
#include <unistd.h>
#include <getopt.h>


// Parse command line  arguments
int main( int argc, char* argv[])
{
	
	int c;

	// Set defaults
	int nneurons = 500;
	int nbranches = 20;
	int ninputs = 400;
	int nperinput = 1;
	int npatterns = ninputs;
	int nonesperpattern = 1;
	int interstim = 60;
	int rseed = 1980;
	char* suffix = NULL;
	bool storeData = false;
	bool disableCreb = false;
	int patternsOverlapping = -1;

	LANetwork net; 
	net.enablePruning = true; // default

	while ((c = getopt(argc, argv, "M:N:H:B:I:i:P:p:T:S:s:d:w:O:g:l:b:c:o:t:xnLDRJCGhU"))!= -1)
	{
		switch (c)
		{
			case '?':
			case 'h':
			cout <<
			"Basic command line options for the simulator:" << endl << 
			"Example: ./lamodel -P 2 -T 1400 [-G] [-L] [-s datadir-name] [-S random-seed]" << endl << 
			"-P <patterns>: number of memories  to enoode"<< endl << 
			"-T <minutes>: Interval between memories (minutes)"<< endl << 
			"-G: Global-only protein synthesis "<< endl << 
			"-L: Dendritic-only protein synthesis "<< endl << 
			"-s <dirname>: Name of the directory to store the data (inside the data/ folder)"<< endl << 
			"-S <random-seed>: Intitialize the random generator"<< endl << 
			""<< endl << 
			"-o <param-name>=<param-value>: Set various simulation parameters:"<< endl << 
			""<< endl << 
			"-o nlTypePV=[0,1,2,3]   : Set basket cell dendrite nonlinearity type. 0=supra, 1=sub, 2=linear, 3=mixed supra/sub"<< endl << 
			"-o nlTypeSOM=[0,1,2,3]  : Set SOM+ cell dendrite nonlinearity type. 0=supra, 1=sub, 2=linear, 3=mixed supra/sub" << endl << 
			"-o INClustered=[0,1]    : Set whether IN synapses  should target all dendrites randomly (dispersed) or only 33% of them (clustered)" << endl;
			return 1;
			break;

			case 'B': nbranches = atoi(optarg); break;
			//case 'I': ninputs = atoi(optarg); break;
			//case 'i': nperinput = atoi(optarg); break;
			case 'N': nneurons = atoi(optarg); break;
			case 'P': npatterns = atoi(optarg); break;
			case 'p': nonesperpattern = atoi(optarg); break;
			case 'T': interstim = atoi(optarg); break;
			case 'S': rseed = ( atoi(optarg)); break;
			case 's': suffix = strdup(optarg); break;
			case 'd': patternsOverlapping = atoi(optarg); break;

			case 'x': storeData = true; break;
			case 'n': disableCreb = true; break;
			case 'w': net.isWeakMem.push_back(atoi(optarg)-1); break;

			case 'L': net.localProteins = true; break;
			case 'G': net.globalProteins = true; break;
			case 'D': net.debugMode = true; break;
			case 'R': net.repeatedLearning = true; break;
			case 'J': net.pretraining = true; break;
			case 'C': net.altConnectivity = true; break;
			case 'O': net.branchOverlap = atof(optarg); break;
			case 'H': net.homeostasisTime = atof(optarg); break;


			case 'o': 
				char* o = strstr(optarg, "=");
				if (o)
				{
					*o = '\0';
					char* val = o+1;

					if (!strcmp(optarg, "connectivityParam")) net.connectivityParam = atof(val); 
					else if (!strcmp(optarg,  "BSPTimeParam")) net.BSPTimeParam = atof(val); 
					else if (!strcmp(optarg,  "homeostasisTimeParam")) net.homeostasisTimeParam = atof(val); 
					else if (!strcmp(optarg,  "CREBTimeParam")) net.CREBTimeParam = atof(val); 
					else if (!strcmp(optarg,  "inhibitionParam")) net.inhibitionParam = atof(val); 
					else if (!strcmp(optarg,  "globalPRPThresh")) net.globalPRPThresh = atof(val); 
					else if (!strcmp(optarg,  "localPRPThresh")) net.localPRPThresh = atof(val); 
					else if (!strcmp(optarg,  "dendSpikeThresh")) net.dendSpikeThresh = atof(val); 
					else if (!strcmp(optarg,  "initWeight")) net.initWeight*= atof(val); 
					else if (!strcmp(optarg,  "maxWeight")) net.maxWeight*= atof(val); 
					else if (!strcmp(optarg,  "stimDurationParam")) net.stimDurationParam = atof(val); 
					else if (!strcmp(optarg,  "nNeuronsParam")) nneurons *= atof(val); 
					else if (!strcmp(optarg,  "nBranchesParam")) nbranches *= atof(val); 
					else if (!strcmp(optarg,  "nBranchesTurnover")) net.nBranchesTurnover = atoi(val); 
					else if (!strcmp(optarg,  "INClustered")) net.INClustered = atoi(val); 
					else if (!strcmp(optarg,  "nlTypePV")) net.nlTypePV = atoi(val); 
					else if (!strcmp(optarg,  "nlTypeSOM")) net.nlTypeSOM = atoi(val); 

					printf("Param name='%s' value='%f'\n", optarg, atof(val));
				}
			break;
		}
	}

	printf("Params:\n net.nBranchesTurnover=%d\n net.connectivityParam=%f\n net.BSPTimeParam=%f\n net.homeostasisTimeParam=%f\n net.CREBTimeParam=%f\n net.inhibitionParam=%f\n net.globalPRPThresh=%f\n net.localPRPThresh=%f\n net.dendSpikeThresh=%f\n net.initWeight=%f\n net.maxWeight=%f\n net.nlTypePv = %d \nnet.nlTypeSOM=%d INClustered=%d\n-------\n" , 

	net.nBranchesTurnover, net.connectivityParam , net.BSPTimeParam , net.homeostasisTimeParam , net.CREBTimeParam , net.inhibitionParam , net.globalPRPThresh , net.localPRPThresh , net.dendSpikeThresh , net.initWeight, net.maxWeight, net.nlTypePV, net.nlTypeSOM, net.INClustered );

	LANetwork::SetRandomSeed(rseed);
	net.disableCreb = disableCreb;

	ninputs = nonesperpattern * npatterns;
	printf("\nNinputs=%d, nPerInput=%d, patterns=%d\n", ninputs, nperinput, npatterns);

	// Create the network 
	net.CreateFearNet(nneurons, nbranches, ninputs, nperinput);

	char buf[512];
	if (suffix)
		sprintf(buf, "./data/%s", suffix );
	else
		sprintf(buf, "./data/N%d.B%d.I%d.i%d.P%d.p%d.T%d.S%d.w%d_%s", nneurons, nbranches, ninputs, nperinput, npatterns, nonesperpattern, interstim, rseed, (int)net.isWeakMem.size(),  suffix ? suffix : "");

	cout << "Data directory is "<< buf <<  endl;
	net.SetDataDir( buf);

	if (net.pretraining)
	{
		char buf2[1024];
		sprintf(buf2, "%s/%s", buf, "pre-syn.dat");
		net.SaveSynapseState(buf2);
	}

	cout << "Running main simulation..."<< endl;
	net.RunStoreTest(npatterns, nonesperpattern, interstim, 0, patternsOverlapping);

	cout << "Storing data files ..."<< endl;
	cout<<buf << endl;
	net.StoreDataFiles( storeData);

	printf("Done!\n");


	return 0;
}



