#include<iostream>
#include<fstream>
#include<cstdlib>

#include<dflpio.h>

void ReadConf(string confFile, int &size, int &years, int **permute, int** uRC, 
	      string &fileprefix) {

  ifstream input;
  input.open(confFile.c_str(), ios::in);
  if (!input) {
    cerr << "Error opening file, exiting ...";
    exit(0);
  }

  if (!(input >> size))
    cerr << "Error reading file contents";
  if (!(input >> years))
    cerr << "Error reading file contents";

   (*permute) = new int[size];
   (*uRC) = new int[size];
  
  for (int i = 0; i < size; i++)
    if (!(input >> (*permute)[i])) {
      cerr << "Error reading file contents";
      exit(1);
    }

  for (int i = 0; i < size; i++)
    if (!(input >> (*uRC)[i])) {
      cerr << "Error reading file contents";
      exit(1);
    }

  if (!(input >> fileprefix)) {
    cerr << "Error reading file contents";
    exit(1);
  }

  input.close();
  return;
}

void PrintConfigInfo(int size, int years, int *permute, int *uRC, string fileprefix) {
  cout << "Problem size: " << size << endl;
  cout << "Number of years: " << years << endl;
  cout << "Initial permuation: ";
  for (int i = 0; i < size; i++)
    cout <<  permute[i] << " "; 
  cout << endl;
  cout << "Relocation cost: ";
  for (int i = 0; i < size; i++)
    cout <<  uRC[i] << " "; 
  cout << endl;

  cout << "Flow data prefix: " << fileprefix << endl;
  return;
}

void PrintPermutations(int size, int n, int **permutes) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < size; j++)
      cout << permutes[i][j] << " ";
    cout << endl;
  }
  return;
}


/* 
 * Recursively print a path in the layout
 */
void printPath(int node, int previous[], int **permutations, int size, int fact, int years) {
  if (node == 0)
    return;

  if (previous[node] == -1) {
    cerr << "No path from source to " << (node + 1) << endl;
    return; 
  }
  
  printPath(previous[node], previous, permutations, size, fact, years);
  if (node != fact * years + 1) {
    int j = (node % fact) - 1;
    for (int i = 0; i < size; i++)
      cout << permutations[j][i] << " ";
    cout << node - 1 << endl;
  }

  return;
}

/*
 * Prints optimal path; make use of printPath() 
 */
void printOptimalPath(int distance[], int previous[], int **permutations, int fact,
		      int size, int years) {
  for (int i = fact * years + 1; i < fact * years + 2; i++) {
      cout << "Optimal path: " << endl;
      printPath(i, previous, permutations, size, fact, years);
      cout << endl << "Cost: " << distance[fact * years + 1] << endl;
  }
}
