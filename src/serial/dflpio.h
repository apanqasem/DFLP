#include<string>

using namespace std;

void ReadConf(string confFile, int &size, int &years, int **permute, int** uRC, 
	      string &fileprefix);


void PrintConfigInfo(int size, int years, int *permute, int *relocCost, string fileprefix);
void PrintPermutations(int size, int n, int **permutes);

void printPath(int node, int previous[], int **permutations, int size, int fact, int years);
void printOptimalPath(int distance[], int previous[], int **permutations, int fact,
		      int size, int years);

