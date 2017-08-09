/* 
 * Sequential version of DFLP solution 
 * 
 * @author: Chandra Kolla
 * @date:
 *
 * Revisions 
 * 
 * 1. Major revision 
 *       @author: Apan Qasem
 *       @date: 06/24/15
 *       major re-write of code  
 *       Significant increase in performance for limit <= 40K
 */

#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<string>
#include<cstdlib>

#include<algorithm> // for fill()
#include<dflpio.h>  // dflp I/O routines 

#ifdef PROFILE 
#include<time.h>
#endif

#define INFINITY 999999999

using namespace std;

struct arrays{
  int *finalmhc;
  int **fpermutations;
};

struct relocation{
  int **rc;
};

/* 
 * compute factorial recursively 
 * needed for determining total number of permutations 
 * if problem size is n, then permuations = n!
 *
 */ 
unsigned long factorial(int n) {
  if (n == 0 || n == 1)
    return 1;
  else
    return n * factorial(n - 1);
}

void InitPermuteRCost(int ***prcost, int n) {

  int **lprcost = new int*[n + 2];
  for (int j = 0; j < n + 2; j++)
    lprcost[j] = new int[n + 2];

  fill(&lprcost[0][0], &lprcost[0][0] + sizeof(lprcost), 1);
  for(int i = 0; i < n + 2; i++) {
    lprcost[i][0] = -1;
    lprcost[i][n + 1] = -1;
    lprcost[0][i] = -1;
    lprcost[n + 1][i] = -1;
    lprcost[i][i] = 0;
  }
  (*prcost) = lprcost;
  return;
}
/* 
 * Generate all permutations for a given problem size 
 * Only invoked when number of permutations is < LIMIT 
 */ 
void BuildAllPermutations(int size, int *initPermute, int n, int **permutations) {
  std::sort(initPermute, initPermute + size);
  int i = 0;
  do {
    if (i < n) {
      for (int j = 0; j < size; j++)
	permutations[i][j] = initPermute[j];
      i++;
    }
    else 
      break;

  } while (std::next_permutation(initPermute, initPermute + size));
}

/*
 * Generate a partial set of random permutations of all possible permutations for the problem 
 * size. Number of permutations generated is based on the global constant LIMIT 
 */
void BuildPartialPermutations(int size, int *initPermute, int n, int **permutations) {
  int j;
  for (int k = 0; k < n; k++) {
    for (int i = size - 1; i > 0; i--) {
      j = rand() % (i + 1);
      int temp = initPermute[i];
      initPermute[i] = initPermute[j];
      initPermute[j] = temp;
      for (int l = 0; l < size; l++) {
  	permutations[k][l] = initPermute[l];
      }
    }
  }
}

/* 
 * Compute material handling cost between two permutations 
 */  
float calcMHC(float *h_flows, float *h_dist, int *init_sol, int nsize,
		  int row) {
  float calcost;
  for (int i = 0; i < row; i++) {
    calcost = 0;
    //#pragma omp parallel for reduction(+:calcost)
    for (int j = 0; j < nsize - 1; j++) {
      for (int k = j + 1; k <= nsize - 1; k++) {
	calcost = calcost
	  + (h_flows[(init_sol[(i * nsize) + j] - 1) * nsize
		     + (init_sol[(i * nsize) + k] - 1)])
	  * h_dist[j * nsize + k];
      }
    }
    //#pragma omp parallel for reduction(+:calcost)
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
    if (!input) {
      cerr << "Error opening file " << ss.str().c_str() << ". Exiting ..." << endl;
      exit(0);
    } 
    if (input.eof()) {
      cerr << "Error in reading contents from " << ss.str().c_str() << ". Exiting ..." << endl;
      exit(0);
    }      

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
      calcost = calcMHC(h_flows, h_dist, supply, nsize, row);
      costarray[count] = calcost;
    }
  }
}

void swap(int& a, int& b, int **permutations, int x, int y, int size, int n) {
  /* swap cost */
  int temp = a;
  a = b;
  b = temp;

  /* swap permutations */
  int *temp1 = new int[size];
  for (int i = 0; i < size; i++)
    temp1[i] = permutations[x][i];
  for (int i = 0; i < size; i++)
    permutations[x][i] = permutations[y][i];
  for (int i = 0; i < size; i++)
    permutations[y][i] = temp1[i];
}

int pivot(int a[], int first, int last, int** permutations, int size, int fact) {
  int p = first;
  int rare = first;
  int pivotElement = a[first];

  for (int i = first + 1; i <= last; i++) {
    // If you want to sort the list in the other order, change "<=" to ">"                                                
    if (a[i] <= pivotElement) {
      rare++;
      p++;
      swap(a[i], a[p], permutations, i, rare, size, fact);
    }
  }

  swap(a[p], a[first], permutations, rare, first, size, fact);

  return p;

}

void sortCost(int a[], int first, int last,int **permutations, int size,
              int fact, int top[],int temp) {

  int pivotElement;
  if (first < last) {
    pivotElement = pivot(a, first, last, permutations, size, fact);
    sortCost(a, first, pivotElement - 1, permutations, size, fact, top, temp);
    sortCost(a, pivotElement + 1, last, permutations, size, fact, top, temp);
  }
  for (int i = 0; i < temp; i++) {
    top[i] = a[i];
#ifdef DEBUG
    cout<<a[i]<<endl;
#endif
  }
}



/*                                                                                                                        
 * computing relocation costs for all pairs of permutations                                                               
 */

void computeRC_chandra(int **permutations1, int **permutations2, int n, int size, int *uRC, int **pRC) {
  int i, j, k;
#pragma omp parallel for private(i,j,k)
  for (i = 0; i < n; i++) {
    int cost = 0;
    for (j = 0; j < n; j++) {
      for (k = 0; k < size; k++) {
	if (permutations1[i][k] != permutations2[j][k]) {
	  cost += uRC[permutations1[i][k] - 1];
	}
      }
      pRC[i][j]=cost;
      cost = 0;
    }
  }
}

/* 
 * computing relocation costs for all pairs of permutations 
 */
void computeRC(int **permutations, int n, int size, int *uRC, int **pRC) {
  
  int i, j, k;
#pragma omp parallel for num_threads(16) private(i,j,k)
  for (i = 1; i < n + 1; i++) {
    int cost = 0;
    for (j = i + 1; j < n + 1; j++) {
      for (k = 0; k < size; k++) {
	if (permutations[i - 1][k] != permutations[j - 1][k]) {
	  cost += uRC[permutations[i - 1][k] - 1];
	}
      }
      pRC[i][j] = pRC[j][i] = cost;
      cost = 0;
    }
  }
}



void calculateDistance_chandra(bool mark[], int distance[], int fact, int previous[],
		       arrays *a, int years, relocation *r) {

int x;
int closestUnmarkedNode, closestUnmarkedNode1;
for (int i = 0; i < fact * years + 1; i++) {
  int minDistance = INFINITY;
  //for (i = 0; i < fact * years + 1; i++) {                                                                      
  if ((minDistance >= distance[i])) {
    minDistance = distance[i];
    if (i <= fact){
      closestUnmarkedNode = i;
      x=1;
    }
    else if (i > fact) {
      if (i % fact == 0){
	x=i/fact;
	closestUnmarkedNode = fact;
      }
      else{
	closestUnmarkedNode = i % fact;
	x=i/fact+1;
      }
    }
    closestUnmarkedNode1 = i;
  }
  //cout<<endl<<x<<endl;                                                                                          
  //}                                                                                                             
  //cout<<endl<<closestUnmarkedNode<<endl;                                                                        
  //cout << closestUnmarkedNode1 << endl;                                                                         
  mark[closestUnmarkedNode1] = true;
  if (closestUnmarkedNode1 == 0) {
    for (int i = 1; i < fact + 1; i++) {
      //if ((mat[0].costing[closestUnmarkedNode][i] >= 0)) {                                                      
      if (distance[i] > distance[closestUnmarkedNode]
	  //+ mat[0].costing[closestUnmarkedNode][i]                                                              
	  + a[x].finalmhc[i-1]) {
                distance[i] = distance[closestUnmarkedNode]
                  //+ mat[0].costing[closestUnmarkedNode][i]                                                              
                  + a[x].finalmhc[i-1];
                previous[i] = closestUnmarkedNode;
                //}                                                                                                       
      }
    }
  }
  //cout<<x<<endl;                                                                                                
  if (closestUnmarkedNode1 >= 1 && closestUnmarkedNode1 <= (years-1) * fact) {
    for (int j = 0; j < (years - 1); j++) {
      if (closestUnmarkedNode1 >= j * fact + 1
	  && closestUnmarkedNode1 <= (j + 1) * fact) {
	for (int i = (j + 1) * fact + 1; i < (j + 2) * fact + 1;
	     i++) {
	  if (i % fact != 0) {
	    //if (r[x-1].rc[closestUnmarkedNode][i % fact-1] >= 0)                                                
	    if (distance[i]
                        > distance[closestUnmarkedNode1]
                        + r[x].rc[closestUnmarkedNode-1][i
                                                         % fact-1] + a[x+1].finalmhc[(i-1)%fact]) {

                      distance[i] = distance[closestUnmarkedNode1]
                        + r[x].rc[closestUnmarkedNode-1][i
                                                         % fact-1] + a[x+1].finalmhc[(i-1)%fact];
                      previous[i] = closestUnmarkedNode1;
		      /*cout<<closestUnmarkedNode<<endl;                                                            
			cout<<endl<<i<<" "<<distance[i]<<endl;                                                      
			cout<<endl<<"distance[node1] "<<distance[closestUnmarkedNode1]<<" ";                        
			cout<<endl<<"rc "<<r[x].rc[closestUnmarkedNode-1][i% fact-1]<<endl;                         
			cout<<endl<<"finalmhc "<<a[x+1].finalmhc[i-1]<<endl;*/
	    }
	  }
	  else if (i % fact == 0) {
	    //if (r[x-1].rc[closestUnmarkedNode-1][i % fact-1] >= 0)                                        
	    if (distance[i]
                              > distance[closestUnmarkedNode1]
                              + r[x].rc[closestUnmarkedNode-1][fact-1]
		+ a[x+1].finalmhc[fact-1]) {
                            distance[i] = distance[closestUnmarkedNode1]
                              + r[x].rc[closestUnmarkedNode-1][fact-1]
                              + a[x+1].finalmhc[fact-1];
                            previous[i] = closestUnmarkedNode1;
                            /*cout<<"closest "<<closestUnmarkedNode<<endl;                                                
                              cout<<endl<<i<<" "<<distance[i]<<endl;                                                      
                              cout<<endl<<"distance[node1] "<<distance[closestUnmarkedNode1]<<" ";                        
                              cout<<endl<<"rc "<<r[x].rc[closestUnmarkedNode-1][fact-1]<<endl;                            
                              cout<<endl<<"finalmhc "<<a[x+1].finalmhc[fact-1]<<endl;*/
	    }
	  }
	}
      }
    }
    //cout<<x<<endl;                                                                                        
  }

  if (closestUnmarkedNode1 > (years - 1) * fact) {
    //if ((mat[5].costing[closestUnmarkedNode][fact + 1] >= 0)) {                                           
    if (distance[years * fact + 1] > distance[closestUnmarkedNode1]){
      //+ mat[5].costing[closestUnmarkedNode][fact + 1]                                                     
      // + a[x][years * fact + 1]) {                                                                        
      distance[years * fact + 1] = distance[closestUnmarkedNode1];
      //+ mat[5].costing[closestUnmarkedNode][fact + 1]                                                     
      // + costarr[years * fact + 1];                                                                       
      previous[years * fact + 1] = closestUnmarkedNode1;
      //}                                                                                                   
    }
  }
 }
}


//Find out the total cost to get to the last node
void calculateDistance(bool mark[], int distance[], int fact, int previous[],
		       int **permutations, int years, int *mhc, int **rc) {

  int i = 0, j = 0;

  int cunRC, cunDist;
  int iter;
  int allPerms = fact * years + 1;

  for (i = 0; i < allPerms; i++) {

    if (INFINITY >= distance[i]) {
      iter = i % fact;
      if (i <= fact)
	cunRC = i;
      else {
	if (iter == 0)
	  cunRC = fact;
	else
	  cunRC = iter;
      }
      cunDist = i;
    }

    mark[cunDist] = true;

    if (cunDist == 0) {
      for (j = 1; j < fact + 1; j++) {
	int dist = distance[cunDist] + mhc[j];
	if (distance[j] > dist) {
	  distance[j] = dist;
	  previous[j] = cunDist;
	}
      }
    }

    if (cunDist >= 1 && cunDist <= (years - 1) * fact) {
      for (int j = 0; j < (years - 1); j++) {
	if (cunDist >= j * fact + 1 && cunDist <= (j + 1) * fact) {
	  int lb = (j + 1) * fact + 1;
#pragma omp parallel for num_threads(16) private(iter)
	  for (int i = lb; i < lb + fact; i++) {
	    iter = i % fact;
	    if (rc[cunRC][iter] >= 0) {
	      if (iter == 0) 
		iter = fact;
	      int dist = distance[cunDist] + rc[cunRC][iter] + mhc[i];
	      if (distance[i] > dist) {
		distance[i] = dist;
		previous[i] = cunDist;
	      }
	    }
	  }
	}
      }
    }
		
    if (cunDist > (years - 1) * fact) {
      if (distance[allPerms] > distance[cunDist] + mhc[allPerms]) {
	distance[allPerms] = distance[cunDist] + mhc[allPerms];
	previous[allPerms] = cunDist;
      }
    }
  }
}

/*
 * Determine optimal path in multi-year layout 
 */
void optimumPath(int **permutations, int fact, int size, int years,
		 int *mhc, int **rc) {

  int *previous = new int[fact * years + 2];
  int *distance = new int[fact * years + 2];
  bool *mark = new bool[fact * years + 2];
  for (int i = 0; i < fact * years + 2; i++) {
    mark[i] = false;
    previous[i] = -1;
    distance[i] = INFINITY;
  }
  distance[0] = 0;
  calculateDistance(mark, distance, fact, previous, permutations, years, mhc, rc);
  printOptimalPath(distance, previous, permutations, fact, size, years);
}

int main(int argc, char *argv[]) {

  if (argc != 4) {
    cout << "usage: " << endl;
    cout << "\t./dflp conf_file rand_seed limit" << endl;
    exit(1);
  }

#ifdef DEBUG
  cout <<"config file name: "<< argv[1] << endl; 
  cout <<"random seed: "<< argv[2] << endl; 
  cout <<"max permutations: "<< argv[3] << endl; 
#endif

#ifdef PROFILE
  clock_t start, end;
  start = clock();
#endif

  string filename = argv[1];
  unsigned int seed = atoi(argv[2]);
  unsigned long limit = atoi(argv[3]);

  if (seed == 0)
    srand(time(0));
  else
    srand(seed);

  int size, years;
  int *rcost;
  int *init;
  string fileprefix;
  ReadConf(filename, size, years, &init, &rcost, fileprefix);

#ifdef DEBUG
  PrintConfigInfo(size, years, init, rcost, fileprefix);
#endif 

  bool partial = false;
  unsigned long n = factorial(size);
  if (n > limit) {
    n = limit;
    partial = true;
  }

  // initialize relocation costs  
  int **prcost;
  InitPermuteRCost(&prcost, n);

  // buld permutations 
  int **permutations = new int*[n];
  for (int i = 0; i < n; i++)
    permutations[i] = new int[size];
  if (partial)
    BuildPartialPermutations(size, init, n, permutations);
  else
    BuildAllPermutations(size, init, n, permutations);

  // compute material handling cost
  int *mhc = new int[n * years + 2];
  mhc[0] = mhc[n * years + 1] = 0;
  computeMHC(permutations, n, mhc, years, size, fileprefix);

  // compute relocation cost
  computeRC(permutations, n, size, rcost, prcost);

  // Chandra's code multi-year relocation cost 

  // allocation                                                                                                           

  relocation *r=new relocation[years-1];
  for(int i=1;i<years;i++){
    r[i].rc=new int*[n];
    for(int j=0;j<n;j++)
      r[i].rc[j]=new int[n];
  }

  for(int i = 1; i < years; i++){
    computeRC_chandra(permutations, permutations, n, size, rcost, prcost);
    for(int j=0;j<n;j++){
      for(int k=0;k<n;k++){
        r[i].rc[j][k]=prcost[j][k];
      }
    }
  }


  // derive optimal path 
  optimumPath(permutations, n, size, years, mhc, prcost);

#ifdef PROFILE
  end = clock();
  cout << "Time taken ";
  float total = (end - start);
  cout << setprecision(2) << fixed;
  cout << total / CLOCKS_PER_SEC << endl;
#endif
  return 0;
}
