//
//    Version: $Id: constructs.h 170 2014-01-29 14:00:04Z gk $
//
/* 
 
lamodel is a network simulator that simulates memory engram formation 
in a population consisting of excitatory and inhibitory neurons with independent
dendritic subunits. 


This file contains definitions of data structures.

*/
#ifndef  __CONSTRUCTS_H__
#define  __CONSTRUCTS_H__



#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <iostream>

#include <algorithm>
#include <vector>
#include <fstream>
#include <map>


using namespace std;


// leakage current
#define param_E_L  0.0 
// resting potential
#define param_V_rest   param_E_L 



// Running mode
enum  {
	RUN_TRAINING = 1,
	RUN_PRE = 2,
	RUN_TEST = 3,
};



//forward  declarations
struct LANeuron;
struct LABranch;



// Random number in [0,1]
inline double rgen()
{
	return double(std::rand())/double(RAND_MAX);
}



// Synapse  data
struct LASynapse
{
	int sid;    	// synapse id
	bool isPlastic; // Whether it undergoes plasticity or not
	float weight, iltp, eltp, pos, calcium, stag, trigger, ltag, stdpTag ;

	float location; // Position in branch , [0,1]

	vector<float> weightHistory;
	vector<float> tagHistory;

	LANeuron* source_nrn;
	LANeuron* target_nrn;
	LABranch* target_branch;

	LASynapse()
	{
		Reset();
	}

	~LASynapse() 
	{ 
		Reset();
	}

	void Reset()
	{
		source_nrn = target_nrn  =0;
		target_branch = 0;
		weight = 0.9999;
		isPlastic = false;
		iltp = ltag = stag = eltp = calcium =0.0;
		pos = rgen();
		sid = -1;
		stdpTag =0;
	}
};





// Dendritic nonlinearity types for interneurons
enum dend_conds {
	DEND_SUPRA = 0,
	DEND_SUB,
	DEND_LINEAR,
	DEND_MIXED,
};



// A neuronal branch , tracking voltage and plasticity related vars
struct LABranch
{
	int bid, branch_spikes ;
	LANeuron* neuron;


	vector<LASynapse*> synapses;
	float depol;
	float protein, proteinRate, strength, dreset, strengthTag, dspike, totcalc;
	int dspikeT;
	vector<pair<float,float> > prpTransients; // Timestamp of latest PRP transients on this branch

	float turnoverRate; // Synapse turnover rate for this branch
	float turnoverPoint;

	int nlType; // Nonlinearity type of this branch

	vector<float> branchStrengthHistory;	// For graphs 
	vector<float> branchProteinHistory;	// For graphs 
	vector<float> branchVoltageHistory; 		// For graphs
	vector<float> branchSpikesHistory; 		// For graphs
	vector<float> branchCalciumHistory;


	LABranch() 
	{
		Reset();
	}

	~LABranch() 
	{ 
		Reset();
	}


	void Reset()
	{
		turnoverRate= 0.;
		turnoverPoint = rgen();
		depol =0.0;
		bid = -1;
		bid = -1;
		proteinRate = protein = 0.;
		strength = 1.0;
		strengthTag =0.0;
		branch_spikes = 0;
		dspike = dreset =0.0;
		totcalc =0.0;
		dspikeT =0;
		nlType = DEND_SUPRA;
	}
};



struct LANetwork;

/* Tracks somatic potential, somatic protein availability, spiking and backpropagating */
struct LANeuron
{
	int nid, input_id, input_type;
	char type;

	float pos_x,pos_y,pos_z; // For 3d graphs
	float glx, gly;

	float V, w, crebLevel, protein, proteinRate, synScaling, totcalc;

	bool isSpiking;
	float wadapt, actvar, vreset, vspike;
	float stdp_x, branch_scaling;
	float synapticWeightsInitialSum;
	int lastSpikeT;
	float exc_cur;
	float inh_cur;
	LANetwork* network;

	vector<pair<float,float> > prpTransients;

	vector<float> voltageHistory; 
	vector<int>   spikeTimings; 

	vector<float>   proteinHistory; 
	vector<float>   crebHistory; 

	int total_spikes, dend_spikes;

	vector<LABranch*> branches;
	vector<LASynapse*> outgoing;
	vector<LASynapse*> incoming;

	LANeuron()
	{
		Reset();
	}

	~LANeuron() 
	{
		Reset();
	}

	void Reset()
	{
		type = ' ';
		w = crebLevel = protein =  wadapt = actvar = vreset= 0.0;
		nid = -1;
		V = 0.;

		pos_x = rgen();
		pos_y = rgen();
		pos_z = rgen();
		stdp_x=0.0;
		lastSpikeT =0;

		protein = proteinRate = 0.0;
		totcalc =0.0;

		synScaling = total_spikes =0;
		branch_scaling =0;
		input_id = -1;
		vspike =0.;
		isSpiking = false;
		dend_spikes =0;
		exc_cur =inh_cur= 0;
	}

};



/* Artificial spike generator neurons, used to provide input stimulation to the network */
struct LAInput: public LANeuron
{
	int* spikeTimes;
	int curSpike;
	int nextSpikeT;
	int totalSpikes;
	int groupIdx; // id of this neuron in the group representing a memory


	LAInput()
	{
		LANeuron();
		spikeTimes = 0;
		groupIdx =-1;
	}

	~LAInput()
	{
		Reset();
	}

	void Reset() 
	{ 
		curSpike = -1;
		delete[] spikeTimes;
		spikeTimes = 0;
		totalSpikes =0;
		nextSpikeT = -1;
	}

	

	int Program(int tstart, int duration, float freq, float randomness)
	{
		Reset();
		int total = round((float(duration) * freq)/1000.0);

		if (!total) return 0;

		float period = 1000.0/freq;
		spikeTimes = new int[total];
		for (int i =0; i < total; i++)
		{
			spikeTimes[i]  =  tstart + period*i + rgen()*randomness*period;
		}

		this->curSpike =0;
		this->nextSpikeT = this->spikeTimes[this->curSpike];
		this->totalSpikes = total;
		return total;
	}



	int CopyShuffled(LAInput &other, float randomness)
	{
		Reset();
		this->spikeTimes = new int[other.totalSpikes];
		for (int i=0; i < other.totalSpikes; i++)
		{
			this->spikeTimes[i] = other.spikeTimes[i] + (rgen() - 0.5)*randomness;
		}

		this->totalSpikes = other.totalSpikes;
		return this->totalSpikes;
	}



};



/* used by wxWidgets */
struct Arr2D {
	float* data;
	int nx, ny;

	Arr2D(int nx, int ny) 
	{ 
		this->data = new float[nx*ny];
		this->nx = nx;
		this->ny = ny;
	}

	float& at(int x, int y)
	{
		return data[x*nx+y];
	}

}; 


/* used by wxWidgets */ 
class LAWindow;




/* iterator shortcuts */
typedef vector<LANeuron*> 			nrn_list;
typedef vector<vector<LANeuron*> >::iterator 	input_iter;
typedef vector<LANeuron*>::iterator 		nrn_iter;
typedef vector<LABranch*>::iterator 		branch_iter;
typedef vector<LASynapse*>::iterator 		syn_iter;
typedef vector< pair<float, float> >::iterator 	pair_iter;



/* Global structure to hold network configuration */ 
class LANetwork
{

	public:

	vector<LASynapse*> synapses; /* List of all synapses */
	vector<LANeuron*> neurons;
	vector<LABranch*> branches;



	// neuron lists
	vector<LANeuron*> pyr_list;
	vector<LANeuron*> in_pv;
	vector<LANeuron*> in_som;
	vector<LANeuron*> bg_list;
	vector<LANeuron*> noise_inputs;

	vector< vector<LANeuron*> > inputs_cs;  /*  conditioned inputs */
	vector< vector<LANeuron*> > inputs_us; /* unconditioned inputs */

	vector<LANeuron* >    inputs_binary;

	vector< vector<int> > spikesPerPattern;

	vector<float> dbgNeuron;

	vector< vector<int> > spikeTimings;  // stores time of spikes during stimulation  only!
	vector< vector<int> > spikesPerStim;  // stores time of spikes during stimulation  only!
 
	vector< vector<float> > nrnVoltages; // ditto for voltages
	map< pair<int, int>, double> distances; // holds euclidean between neurons if needed

	static int RSEED;

	int runningMode, runningPatternNo;
	float localPRPThresh, globalPRPThresh;
	float homeostasisTimeParam; // Time that it takes for synaptic scaling to be applied
	float BSPTimeParam; // Time that it takes for BSP to  be applied
	float CREBTimeParam; // Time that it takes for CREB to fall
	float connectivityParam, inhibitionParam, stimDurationParam; // Multipliers for doing sensitivity analysis 

	int nlTypePV, nlTypeSOM;


	int weakMemId;
	int nBranchesTurnover;
	vector<int> isWeakMem;

	pthread_mutex_t synapses_mutex;


	ofstream mfile, vfile, sumweightsFile;

	vector< vector<int> > patterns;
	float homeostasisTime;

	int synapsesCounter;


	LAWindow* wx_window;

	FILE* spikesFile; // Save spiking info 'ere 

	Arr2D* voltageData;

	int     n_neurons, 
		n_branches_per_neuron,
		n_inputs,
		n_neurons_per_input,
		Tstimulation, // Total stimulated time
		T; // Simulation clock time 
	bool enablePlasticity; /* Is plasticity enabled? */
	bool isInterstim;
	bool disableCreb;
	bool debugMode;

	bool isRecall;

	bool localProteins, repeatedLearning, pretraining, altConnectivity, globalProteins;

	char* conditionsString;

	string datadir ;

	float branchOverlap;

	float initWeight, maxWeight, dendSpikeThresh;
	bool enablePruning, isPruningSynapses, INClustered;

	LANetwork()
	{
		enablePruning = isPruningSynapses = false;
		synapsesCounter =0;
		T  = Tstimulation =0;
		isRecall = false;
		wx_window = NULL;
		n_neurons = n_branches_per_neuron = n_inputs = n_neurons_per_input =0;
		enablePlasticity = true;
		isInterstim = false;
		nBranchesTurnover=0;
		spikesFile = NULL;
		runningMode = RUN_TRAINING;
		runningPatternNo = 0;
		weakMemId = -1;
		repeatedLearning = globalProteins =  localProteins = false;
		pthread_mutex_init(&this->synapses_mutex, NULL);
		pretraining = false;
		debugMode = false;
		altConnectivity=false;
		conditionsString = NULL;
		datadir = "./";
		branchOverlap = -1.0; 
		disableCreb = false;
		homeostasisTime = 24.0;

		localPRPThresh=1.0;
		globalPRPThresh=1.0; 
		homeostasisTimeParam = 1.0;
		BSPTimeParam = 1.0;
		CREBTimeParam = 1.0;
		inhibitionParam = 1.0;
		connectivityParam = 1.0;
		stimDurationParam = 1.0;
		dendSpikeThresh = 1.0;
		initWeight = 0.3;
		maxWeight = 1.0;


		// Default is PV dendrites = supra, SOM dendrites = sublinear
		nlTypePV = DEND_SUPRA;
		nlTypeSOM = DEND_SUB;

		INClustered=0;
	}




	static void SetRandomSeed( int seed)
	{
		LANetwork::RSEED = seed+80; // XXX
		std::srand(seed+80);
	}



	~LANetwork() 
	{
		Cleanup();
	}


	/* Perform synapse turnover (delete and add new synapses ) */
	void DoTurnover(float durationSecs );

	/* Set up the network, neurons and connections */
	void CreateFearNet(int, int , int, int);


	void Cleanup(void)
	{
		ResetSpikeCounters();

		for (nrn_iter ni = this->neurons.begin(); ni != this->neurons.end(); ni++)
			delete (*ni);

		for (branch_iter ni = this->branches.begin(); ni != this->branches.end(); ni++)
			delete (*ni);

		for (syn_iter ni = this->synapses.begin(); ni != this->synapses.end(); ni++)
			delete (*ni);

		this->neurons.clear();
		this->branches.clear();
		this->synapses.clear();

		if (this->spikesFile)
			fclose(this->spikesFile);
	}


;
	/* Connect two sets of neurons , can specify minimum  / maximum allowed distances between pairs of neurons */
	int ConnectNeurons(vector<LANeuron*> fromList, vector<LANeuron*> toList, bool IsClustered, float toDistance, int nNeuronPairs, int nSynapsesPerNeuron,float, bool, bool, float );


	/* Add a synapse to the network */
	void AddSynapse(LANeuron* a, LABranch* br, float weight, bool isPlastic);

	/* Connect two lists of neurons, randomly */
	int ConnectInputs(vector<LANeuron*> fromList, vector<LANeuron*> toList, int nSynapses);


	/* Delete a number of input synapses */
	int PurgeInputSynapses( int totalToRemove, float);

	/* Add a number of input synapses */
	int CreateInputSynapses( int totalToAdd);

	/* Create a set of neurons and append them to specified list */
	void CreateNeurons(int number, int branches_per_neuron, char type, vector<LANeuron*>* appendTo, int, int );

	/* Calculate distances between  neurons */
	void CalculateDistances();



	/* Simulate detailed voltage dynamics (1msec per step) */
	void StimDynamics(int duration);

	/* Simulate inter-stimulus dynamics (protein / creb level / synapse weights changes only */
	void Interstim(int duration);


	void ResetSpikeCounters(void)
	{
		for (nrn_iter i = neurons.begin(); i != neurons.end(); i++)
		{
			(*i)->total_spikes =0;
			(*i)->dend_spikes =0;
			for (branch_iter b = (*i)->branches.begin(); b != (*i)->branches.end(); ++b)
			{
				(*b)->branch_spikes = 0;
			}
		}
	}

	void ResetCrebLevels() 
	{
		for (nrn_iter i = neurons.begin(); i != neurons.end(); i++)
			(*i)->crebLevel =0.0;
	}

	/* Perform main simulation */
	void RunStoreTest(int  np, int app, int dm, int test=0, int a=-1);

	/* Stimulate using a single pattern */
	void RunPattern(vector<int>&, float, float, int , int );

	void StoreDataFiles( bool);

	void RecordInitialWeightSums()
	{
		for (nrn_iter i = neurons.begin(); i != neurons.end(); i++)
		{
			LANeuron* nrn = *i;
			nrn->synapticWeightsInitialSum =0.;
			
			for (branch_iter b = (*i)->branches.begin(); b != (*i)->branches.end(); ++b)
				for (syn_iter si = (*b)->synapses.begin(); si != (*b)->synapses.end(); ++si)
				{
					nrn->synapticWeightsInitialSum += (*si)->weight;
				}
		}
	}


	void SaveSpikeCounters(ofstream& ratesdat)
	{
		for (nrn_iter na = neurons.begin(); na != neurons.end(); ++na)
		{
			LANeuron* nrn = *na;
			ratesdat << nrn->total_spikes << " ";
		}
		ratesdat << endl;
	}

	
	void SaveSnapshot(char*);

	bool HasCondition(char* cond)
	{
		if (this->conditionsString && strstr(this->conditionsString, cond))
			return true;
		return false;
	}

	bool SaveSynapseState(char* filename)
	{
		ofstream synstatedat(filename);

		for (syn_iter si =this->synapses.begin(); si != this->synapses.end(); ++si)
		{
			LASynapse* s = *si;
			if (s->isPlastic)
			{
			synstatedat << s->sid<<" " 
					<< s->target_branch->bid<<" "
					<< s->target_nrn->nid  << " "
					<< s->source_nrn->nid <<" " 
					<< s->source_nrn->input_id<< " "
					<< s->target_branch->strength  << " "
					<< s->weight << " " 
					<<endl;


			}
		}
		return true;
	}

	void ReportSumWeights()
	{

		float consolidated[this->n_inputs] ;

		for (int i=0; i < this->n_inputs; i++) consolidated[i] =0;

		if (!this->sumweightsFile.is_open())
		{
			this->sumweightsFile.open((this->datadir + "/sum-weights.txt").c_str(), std::ofstream::out );
		}

		int totalPot =0;
		for (syn_iter si = this->synapses.begin(); si != this->synapses.end(); ++si)
		{
			LASynapse* s =*si;
			if ( s->source_nrn->input_id >=0 && s->target_nrn->type == 'P')
			{
				if (s->weight > 0.7)
					totalPot ++;
				consolidated[s->source_nrn->input_id] += s->weight;
			}
		}
		
		for (int i=0; i < this->n_inputs; i++)
		{
			//cout << " ["<< i <<"]-> "<< consolidated[i]<< endl;

			this->sumweightsFile << consolidated[i] <<  " " ;
			cout << consolidated[i] <<  " " ;
		}
		this->sumweightsFile << endl;
		cout << endl;

		cout <<"Total psyn: "<<totalPot<<endl;
	}

	void SetDataDir(string dir)
	{
		this->datadir = dir;
		mkdir(this->datadir.c_str(), 0755);

	}

	void PrintSynapsesSnapshot(string outfile)
	{
		ofstream  fout(outfile.c_str());

		for (nrn_iter ni = this->pyr_list.begin(); ni != this->pyr_list.end(); ni++)
		{
			LANeuron*  nrn = *ni;
			for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
			{
				LABranch* b = *bi;
				for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse* s =*si;
					LAInput* src = (LAInput*)s->source_nrn;
					if (src->input_id >=0) 
					{
						fout  << src->input_id << " " << src->groupIdx << " " << s->target_branch->bid << " " << s->target_nrn->nid << " " << s->weight << " "  << endl;
					}
					
				}
			}
		}
	}


	void SaveCalcs()
	{
		ofstream ooo("./data/calc.dat");
		for (syn_iter si = this->synapses.begin(); si!= this->synapses.end(); si++)
		{
			LASynapse* s = *si;
			if (s->source_nrn->input_id ==0)
			{
			//	ooo << s->calcium << endl;
			}

		}
		ooo.close();
	}
	
	void RunTests();
};





#endif
