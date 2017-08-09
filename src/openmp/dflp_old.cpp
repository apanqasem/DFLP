#include<iostream>
#include<fstream>
#include<vector>
#include<cstdio>
#include<algorithm>
#include<stdio.h>
#include<cstdlib>
#include <sstream>
#include<string>
#include<time.h>
#include<functional>
#include<iomanip>
#include<omp.h>
#define INFINITY 999999999

using namespace std;
//Finding the total number of permutations possible
int findFactorial(int n) {
	int fact = 1;
	if (n == 0 || n == 1)
		return fact;
	else
		return n * findFactorial(n - 1);
}

//Build the exact permutations in a particular order if the factorial value is less than or equal to the limit  
void buildPermutations(int size, int *a, int fact, int **permutations) {
	std::sort(a, a + size);
	int i = -1;
	do {
		i += 1;
		if (i < fact)
			for (int j = 0; j < size; j++) {
				permutations[i][j] = a[j];
			}
	} while (std::next_permutation(a, a + size));
}

//Build the random permutations equal to the limit if the factorial value is greater than the limit
void buildingPermutations(int size, int *a, int fact, int **permutations) {
	int j;
	for (int k = 0; k < fact; k++) {
		for (int i = size - 1; i > 0; i--) {
			j = rand() % (i + 1);
			int temp = a[i];
			a[i] = a[j];
			a[j] = temp;
			for (int l = 0; l < size; l++) {
				permutations[k][l] = a[l];
			}
		}
	}
}

//Sub function of finding out the material handling cost
float calculating(float *h_flows, float *h_dist, int *init_sol, int nsize,
		int row) {
	float calcost;
	for (int i = 0; i < row; i++) {
		calcost = 0;
		#pragma omp parallel for reduction(+:calcost)
		for (int j = 0; j < nsize - 1; j++) {
			for (int k = j + 1; k <= nsize - 1; k++) {
				calcost = calcost
						+ (h_flows[(init_sol[(i * nsize) + j] - 1) * nsize
						           + (init_sol[(i * nsize) + k] - 1)])
						           * h_dist[j * nsize + k];
			}
		}
		#pragma omp parallel for reduction(+:calcost)
		for (int k = 1; k < nsize; k++) {
			for (int l = 0; l < k; l++) {
				calcost = calcost
						+ h_flows[(init_sol[(i * nsize) + k] - 1) * nsize
						          + (init_sol[(i * nsize) + l] - 1)]
						          * h_dist[k * nsize + l];
			}
		}
	}
	return calcost;
}

//finding out the material handling cost
void computeMHC(int **permutations, int fact, int *costarray, int years,
		int size, string filename) {
	int count = 0;
	int arraySizeX, arraySizeY, a = 0, b = 0;
	float num;
	ifstream input;
	for (int i = 1; i <= years; i++) {
		a = 0;
		b = 0;
		std::ostringstream ss;
		ss << filename << i << ".txt";
		input.open(ss.str().c_str());
		if (!input)
			cout << "error opening file";
		if (input.eof())
			cerr << "Error in reading file contents";
		int nsize = size;
		arraySizeX = 2 * nsize;
		arraySizeY = nsize;
		float ** array;
		array = (float**) malloc(arraySizeX * sizeof(float*));
		for (int i1 = 0; i1 < arraySizeX; i1++)
			array[i1] = (float*) malloc(arraySizeY * sizeof(float));
		for (int row1 = 0; row1 < (arraySizeX); row1++) {
			for (int col1 = 0; col1 < arraySizeY; col1++) {
				array[row1][col1] = 0;
			}
		}
		int size_A = nsize * nsize;
		int mem_size_A = sizeof(int) * size_A;
		float *h_dist = (float *) malloc(mem_size_A);
		float *h_flows = (float *) malloc(mem_size_A);
		for (int i2 = 0; i2 < nsize; i2++) {
			for (int j = 0; j < nsize; j++) {
				h_dist[i2 * nsize + j] = 0;
				h_flows[i2 * nsize + j] = 0;
			}
		}
		while (!input.eof()) {
			input >> num;
			if (b == nsize) {
				a++;
				b = 0;
			}
			if (a != (nsize * 2) && b != nsize) {
				array[a][b] = num;
				b++;
			}
		}
		input.close();
		for (int row = 0; row < nsize; row++) {
			for (int col = 0; col < nsize; col++) {
				h_flows[row * nsize + col] = array[row][col];
			}
		}
		int irow = 0;
		for (int row = nsize; row < nsize * 2; row++) {
			int icol = 0;
			for (int col = 0; col < nsize; col++) {
				h_dist[irow * nsize + icol] = array[row][col];
				icol++;
			}
			irow++;
		}
		int calcost;
		int row = 1;
		int *supply = new int[size];
		for (int in = 0; in < fact; in++) {
			for (int j = 0; j < size; j++) {
				supply[j] = permutations[in][j];
			}
			count += 1;
			calcost = calculating(h_flows, h_dist, supply, nsize, row);
			costarray[count] = calcost;
		}
	}
}

//Computing the relocation cost
void computeRC(int **permutations, int fact, int size, int years,
		int *relocationcostarray, int **costing) {
	int cost = 0,i,j,k,l;
	for(int i=0;i<fact+1;i++){
		for(int j=0;j<fact+1;j++){
			if(i==j)
			costing[i][j]=0;
		}
	}

	#pragma omp parallel for
	for (i = 1; i < fact + 1; i++) {
	  for (j = 1; j < fact + 1; j++) {
	    if(i<j){
	      for (k = 0; k < size; k++) {
		if (permutations[i - 1][k] != permutations[j - 1][k]) {
		  l = permutations[i - 1][k];
		  cost += relocationcostarray[l - 1];
		}
	      }
	      costing[i][j] = costing[j][i] = cost;
	      cost = 0;
	    }
	  }
	}

/*
	#pragma omp parallel for
	for (i = 1; i < fact + 1; i++) {
		for (j = 1; j < fact + 1; j++) {
			if (costing[i][j] == 1) {
				for (k = 0; k < size; k++) {
					if (permutations[i - 1][k] != permutations[j - 1][k]) {
						l = permutations[i - 1][k];
						cost = cost + relocationcostarray[l - 1];
					}
				}
			}
			costing[i][j] = cost;
			cost = 0;
		}
	}*/
}

//Find out the total cost to get to the last node
void calculateDistance(bool mark[], int distance[], int fact, int previous[],
		int **permutations, int years, int size, int *costarr, int **costing) {
	int *sample = new int[fact];
	for (int i = 0; i < fact; i++)
		sample[i] = 0;
	int i;
	int closestUnmarkedNode, closestUnmarkedNode1;
	int count = 0;
	for (i = 0; i < fact * years + 1; i++) {
		int minDistance = INFINITY;
		//for (i = 0; i < fact * years + 1; i++) {
			if ((minDistance >= distance[i])) {
				minDistance = distance[i];
				if (i <= fact)
					closestUnmarkedNode = i;
				else if (i > fact) {
					if (i % fact == 0)
						closestUnmarkedNode = fact;
					else
						closestUnmarkedNode = i % fact;
				}
				closestUnmarkedNode1 = i;
			}
		//}
		//cout << closestUnmarkedNode1 << endl;
		mark[closestUnmarkedNode1] = true;
		if (closestUnmarkedNode1 == 0) {
			for (int i = 1; i < fact + 1; i++) {
				//if ((mat[0].costing[closestUnmarkedNode][i] >= 0)) {
				if (distance[i] > distance[closestUnmarkedNode1]
				                           //+ mat[0].costing[closestUnmarkedNode][i]
				                           + costarr[i]) {
					distance[i] = distance[closestUnmarkedNode1]
					                       //+ mat[0].costing[closestUnmarkedNode][i]
					                       + costarr[i];
					previous[i] = closestUnmarkedNode1;
					//}
				}
			}
		}

		if (closestUnmarkedNode1 >= 1 && closestUnmarkedNode1 <= (years-1) * fact) {
			for (int j = 0; j < (years - 1); j++) {
				if (closestUnmarkedNode1 >= j * fact + 1
						&& closestUnmarkedNode1 <= (j + 1) * fact) {
					for (int i = (j + 1) * fact + 1; i < (j + 2) * fact + 1;
							i++) {
						if (i % fact != 0) {
							if (costing[closestUnmarkedNode][i % fact] >= 0)
								if (distance[i]
								             > distance[closestUnmarkedNode1]
								                        + costing[closestUnmarkedNode][i
								                                                       % fact] + costarr[i]) {
									distance[i] = distance[closestUnmarkedNode1]
									                       + costing[closestUnmarkedNode][i
									                                                      % fact] + costarr[i];
									previous[i] = closestUnmarkedNode1;
								}
						} else if (i % fact == 0) {
							if (costing[closestUnmarkedNode][i % fact] >= 0)
								if (distance[i]
								             > distance[closestUnmarkedNode1]
								                        + costing[closestUnmarkedNode][fact]
								                                                       + costarr[i]) {
									distance[i] = distance[closestUnmarkedNode1]
									                       + costing[closestUnmarkedNode][fact]
									                                                      + costarr[i];
									previous[i] = closestUnmarkedNode1;
								}
						}
					}
				}
			}
		}
		
		if (closestUnmarkedNode1 > (years - 1) * fact) {
			//if ((mat[5].costing[closestUnmarkedNode][fact + 1] >= 0)) {
			if (distance[years * fact + 1] > distance[closestUnmarkedNode1]
			                                      //+ mat[5].costing[closestUnmarkedNode][fact + 1]
			                                      + costarr[years * fact + 1]) {
				distance[years * fact + 1] = distance[closestUnmarkedNode1]
				                                  //+ mat[5].costing[closestUnmarkedNode][fact + 1]
				                                  + costarr[years * fact + 1];
				previous[years * fact + 1] = closestUnmarkedNode1;
				//}
			}
		}
	}
}
//Printing the optimum path
void printPath(int node, int previous[], int **permutations, int size, int fact,
		int source,int years) {
	if (node == source) {
		return;
	} else if (previous[node] == -1)
		cout << "No path from source to " << (node + 1) << endl;
	else if (node == fact * years + 1)
		printPath(previous[node], previous, permutations, size, fact, source,years);
	else {
		printPath(previous[node], previous, permutations, size, fact, source,years);
		for (int i = 0; i < size; i++)
			cout << permutations[(node % fact) - 1][i] << " ";
		cout << node - 1 << endl;
	}
}

//Printing the optimum path
void output(int distance[], int previous[], int **permutations, int fact,
		int size, int years, int source) {
	for (int i = fact * years + 1; i < fact * years + 2; i++) {
		if (i == source)
			cout << (source + 1) << " " << source + 1;
		else {
			cout << "optimum path is " << endl;
			printPath(i, previous, permutations, size, fact, source,years);
			cout << " \noptimum cost is " << distance[fact * years + 1] << endl;
		}
	}
}

//Printing the optimum path
void optimumPath(int **permutations, int fact, int size, int years,
		int *costarray, int **costing, int source) {
	int *previous = new int[fact * years + 2];
	int *distance = new int[fact * years + 2];
	bool *mark = new bool[fact * years + 2];
	for (int i = 0; i < fact * years + 2; i++) {
		mark[i] = false;
		previous[i] = -1;
		distance[i] = INFINITY;
	}
	distance[0] = 0;
	calculateDistance(mark, distance, fact, previous, permutations, years, size,
			costarray, costing);
	output(distance, previous, permutations, fact, size, years, source);
}

int main(int argc, char *argv[]) {
	cout << setprecision(4) << fixed;
	clock_t start, end;
	start = clock();
	ifstream input;
	int source = 0;
	int number = *(argv[2]);
	srand(number);
	cout<<number<<endl;
	int size, fact, years, limit = 85000, var;
	string filename;
	//std::ostringstream ss;
	//ss << argv[1] << i << ".txt";
	cout<<argv[1]<<" ";
	input.open(argv[1],ios::in);
	if (!input)
		cout << "error opening file";
	if (input.eof())
		cerr << "Error in reading file contents";
	//input.open("inputs.txt");
	if (input.fail())
		cout << "Error opening the file";
	else {
		input >> size;
		input >> years;
	}
	int a[size];
	int *relocationcostarray = new int[size];
	if (input.eof())
		cout << "Error in reading file contents";
	else
		for (int i = 0; i < size; i++)
			input >> a[i];
	if (input.eof())
		cout << "Error in reading file contents";
	else
		for (int i = 0; i < size; i++)
			input >> relocationcostarray[i];
	if (input.eof())
		cout << "Error in reading file contents";
	else
		input >> filename;
	input.close();
	fact = findFactorial(size);
	var = fact;
	if (fact > limit)
		fact = limit;
	int **costing = new int*[fact + 2];
	for (int j = 0; j < fact + 2; j++) {
		costing[j] = new int[fact + 2];
	}
	for (int i = 0; i < fact + 2; i++) {
		for (int j = 0; j < fact + 2; j++)
			costing[i][j] = 1;
	}
	for (int i = 0; i < fact + 2; i++) {
		for (int j = 0; j < fact + 2; j++) {
			if (i == 0 || i == fact + 1)
				costing[i][j] = -1;
			if (j == 0 || j == fact + 1)
				costing[i][j] = -1;
		}
	}
	int **permutations = new int*[fact];
	for (int i = 0; i < fact; i++)
		permutations[i] = new int[size];
	if (fact == var)
		buildPermutations(size, a, fact, permutations);
	else
		buildingPermutations(size, a, fact, permutations);
	cout << fact << endl;
	int *costarray = new int[fact * years + 2];
	costarray[0] = costarray[fact * years + 1] = 0;
#pragma omp sections 
{
#pragma omp section 
{
  computeMHC(permutations, fact, costarray, years, size, filename);
 }
#pragma omp section 
{
  computeRC(permutations, fact, size, years, relocationcostarray, costing);
 }	
}
	optimumPath(permutations, fact, size, years, costarray, costing, source);
	end = clock();
	cout << "Time taken ";
	float total = (end - start);
	cout << total / CLOCKS_PER_SEC << endl;
	return 0;
}
