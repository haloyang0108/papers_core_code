#include <utility>
#include <random>
#include <algorithm>
#include <iomanip>
#include <fstream> 
#include <string>
#define SIZE 500
#define c 1
#define lower_bound 0
using namespace std;

typedef struct Particle
{
	int strategy;
	int paststrategy;
	double payoff;
	int state[4]; 
	int paststate[4];
	double reputation;
	double pastreputation;
	int tNeighbor[2];	
	int bNeighbor[2];		
	int lNeighbor[2];		
	int rNeighbor[2];		
}Node;



int transmissionfunc(Node arrNetworkSL[SIZE][SIZE], Node* swarmNode, double omega, double omega_2, double upper_bound)
{
	int st;
	double q;
	Node* surrounding[4] = {
		&arrNetworkSL[swarmNode->lNeighbor[0]][swarmNode->lNeighbor[1]],
		&arrNetworkSL[swarmNode->rNeighbor[0]][swarmNode->rNeighbor[1]],
		&arrNetworkSL[swarmNode->tNeighbor[0]][swarmNode->tNeighbor[1]],
		&arrNetworkSL[swarmNode->bNeighbor[0]][swarmNode->bNeighbor[1]]
	};

	double swarmReputation = swarmNode->pastreputation;

	for (int i = 0; i < 4; i++)
	{
		if (i == 1 || i == 3)
		{
			continue;
		}
		double neighborReputation = surrounding[i]->pastreputation;
		double R_min = min(swarmReputation, neighborReputation);
		double R_max = max(swarmReputation, neighborReputation);

		double R_minRatio = (double)R_min / upper_bound;
		double R_maxRatio = (double)R_max / upper_bound;
		if (swarmNode->paststate[i] == 1)
		{
			double R_weighted = omega * R_maxRatio + (1 - omega) * R_minRatio;
			q = R_weighted;

			if (q == 0.0)
			{
				swarmNode->state[i] = 2;
			}
			else
			{
				swarmNode->state[i] = ((rand() / (double)RAND_MAX) <= q) ? 1 : 2;
			}
		}
		else if (swarmNode->paststate[i] == 2)
		{
			double R_weighted = omega_2 * R_maxRatio + (1 - omega_2) * R_minRatio;
			q = R_weighted;

			if (q == 0.0)
			{
				swarmNode->state[i] = 2;
			}
			else
			{
				swarmNode->state[i] = ((rand() / (double)RAND_MAX) <= q) ? 1 : 2;
			}
		}
		int oppositeStateIndex = (i < 2) ? 1 - i : 5 - i;
		surrounding[i]->state[oppositeStateIndex] = swarmNode->state[i];
	}
	return 0;
}

double arrPayoffHGPD(Node arrNetworkSL[SIZE][SIZE], Node* swarmNode, double b_1, double b_2, double upper_bound)
{
	int i, j;
	int strategy, g[4] = { 0 };
	double pi[4] = { 0.0 }, payoffCum;
	Node* centerNode;
	Node* groupNode[5], * colonyNode[5];

	strategy = swarmNode->strategy;
	groupNode[0] = swarmNode;
	groupNode[1] = &arrNetworkSL[swarmNode->lNeighbor[0]][swarmNode->lNeighbor[1]];
	groupNode[2] = &arrNetworkSL[swarmNode->rNeighbor[0]][swarmNode->rNeighbor[1]];
	groupNode[3] = &arrNetworkSL[swarmNode->tNeighbor[0]][swarmNode->tNeighbor[1]];
	groupNode[4] = &arrNetworkSL[swarmNode->bNeighbor[0]][swarmNode->bNeighbor[1]];

	centerNode = groupNode[0];
	colonyNode[0] = centerNode;
	colonyNode[1] = &arrNetworkSL[centerNode->lNeighbor[0]][centerNode->lNeighbor[1]];
	colonyNode[2] = &arrNetworkSL[centerNode->rNeighbor[0]][centerNode->rNeighbor[1]];
	colonyNode[3] = &arrNetworkSL[centerNode->tNeighbor[0]][centerNode->tNeighbor[1]];
	colonyNode[4] = &arrNetworkSL[centerNode->bNeighbor[0]][centerNode->bNeighbor[1]];

	g[0] = colonyNode[1]->strategy;
	g[1] = colonyNode[2]->strategy;
	g[2] = colonyNode[3]->strategy;
	g[3] = colonyNode[4]->strategy; 

	if (strategy == 1 && g[0] == 1)
	{
		pi[0] = (colonyNode[0]->state[0] == 1) ? (b_1 - c) : (b_2 - c);

	}
	else if (strategy == 1 && g[0] == 0)
	{
		pi[0] = -c;
	}
	else if (strategy == 0 && g[0] == 1)
	{
		pi[0] = (colonyNode[0]->state[0] == 1) ? b_1 : b_2;
	}
	else if (strategy == 0 && g[0] == 0)
	{
		pi[0] = 0;
	}
	
	if (strategy == 1 && g[1] == 1)
	{
		pi[1] = (colonyNode[0]->state[1] == 1) ? (b_1 - c) : (b_2 - c);

	}
	else if (strategy == 1 && g[1] == 0)
	{
		pi[1] = -c;
	}
	else if (strategy == 0 && g[1] == 1)
	{
		pi[1] = (colonyNode[0]->state[1] == 1) ? b_1 : b_2;
	}
	else if (strategy == 0 && g[1] == 0)
	{
		pi[1] = 0;
	}
	
	if (strategy == 1 && g[2] == 1)
	{
		pi[2] = (colonyNode[0]->state[2] == 1) ? (b_1 - c) : (b_2 - c);

	}
	else if (strategy == 1 && g[2] == 0)
	{
		pi[2] = -c;
	}
	else if (strategy == 0 && g[2] == 1)
	{
		pi[2] = (colonyNode[0]->state[2] == 1) ? b_1 : b_2;
	}
	else if (strategy == 0 && g[2] == 0)
	{
		pi[2] = 0;
	}

	if (strategy == 1 && g[3] == 1)
	{
		pi[3] = (colonyNode[0]->state[3] == 1) ? (b_1 - c) : (b_2 - c);

	}
	else if (strategy == 1 && g[3] == 0)
	{
		pi[3] = -c;
	}
	else if (strategy == 0 && g[3] == 1)
	{
		pi[3] = (colonyNode[0]->state[3] == 1) ? b_1 : b_2;
	}
	else if (strategy == 0 && g[3] == 0)
	{
		pi[3] = 0;
	}
	return payoffCum = pi[0] + pi[1] + pi[2] + pi[3];
}
