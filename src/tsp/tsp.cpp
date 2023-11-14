#include "tsp.h"

//---------------------------------------------------------------------------
//Constructor
TSP::TSP(std::vector<TSP_Point> aPointList) : points(aPointList) {
    n = points.size();

    graph = new double *[n];
    for (int i = 0; i < n; i++) {
        graph[i] = new double[n];
        for (int j = 0; j < n; j++) {
            graph[i][j] = 0.;
        }
    }

    adjlist = new std::vector<int>[n];
}

//---------------------------------------------------------------------------
//Destructor
TSP::~TSP() {
    for (int i = 0; i < n; i++) {
        delete[] graph[i];
    }

    delete[] graph;
    delete[] adjlist;
}

//---------------------------------------------------------------------------
void TSP::Solve() {
    fillMatrix();
    findMST();

    perfectMatching();

    // Loop through each index and find shortest path
    double best = std::numeric_limits<double>::max();
    int bestIndex;
    for (long t = 0; t < n; t++) {
        double result = findBestPath(t);
        if (result < best) {
            bestIndex = t;
            best = result;
        }
    }

    //Create path for best tour
    euler_tour(bestIndex, circuit);
    make_hamiltonian(circuit, pathLength);
}

//---------------------------------------------------------------------------

std::vector<int> TSP::getResult() {
    return circuit;
}

//---------------------------------------------------------------------------
//Special care for integers

//double TSP::get_distance(struct TSP_Point c1, struct TSP_Point c2) {
//    double dx = std::pow((double) (c1.x - c2.x), 2);
//    double dy = std::pow((double) (c1.y - c2.y), 2);
//    double dz = std::pow((double) (c1.z - c2.z), 2);
//
//    return std::floor(std::sqrt(dx + dy + dz) + 0.5);
//}

//---------------------------------------------------------------------------

double TSP::get_distance(struct TSP_Point c1, struct TSP_Point c2) {
    double dx = std::pow(c1.x - c2.x, 2);
    double dy = std::pow(c1.y - c2.y, 2);
    double dz = std::pow(c1.z - c2.z, 2);

    return std::sqrt(dx + dy + dz);
}

//---------------------------------------------------------------------------

void TSP::fillMatrix() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            graph[i][j] = graph[j][i] = get_distance(points[i], points[j]);
        }
    }
}

/******************************************************************************
  This function uses Prim's algorithm to determine a minimum spanning tree on
    the graph
******************************************************************************/

void TSP::findMST() {
    double *key = new double[n];
    bool *included = new bool[n];
    int *parent = new int[n];

    for (int i = 0; i < n; i++) {

        // set each key to infinity
        key[i] = std::numeric_limits<double>::max();

        // node node yet included in MST
        included[i] = false;

    }

    // root of MST has distance of 0 and no parent
    key[0] = 0.;
    parent[0] = -1;

    for (int i = 0; i < n - 1; i++) {

        // find closes vertex not already in tree
        int k = getMinIndex(key, included);

        // set included to true for this vertex
        included[k] = true;

        // examine each unexamined vertex adjacent to most recently added
        for (int j = 0; j < n; j++) {

            // node exists, is unexamined, and graph[k][j] less than previous
            // key for u
            if (graph[k][j] && included[j] == false && graph[k][j] < key[j]) {

                // update parent
                parent[j] = k;

                // update key
                key[j] = graph[k][j];
            }
        }
    }

    delete[] key;
    delete[] included;

    // construct a tree by forming adjacency matrices
    for (int i = 0; i < n; i++) {
        int j = parent[i];
        if (j != -1) {
            adjlist[i].push_back(j);
            adjlist[j].push_back(i);
        }
    }

    delete[] parent;
}

/******************************************************************************
  find the index of the closest unexamined node
******************************************************************************/

int TSP::getMinIndex(double key[], bool mst[]) {
    // initialize min and min_index
    double min = std::numeric_limits<double>::max();
    int min_index;

    // iterate through each vertex
    for (int i = 0; i < n; i++) {

        // if vertex hasn't been visited and has a smaller key than min
        if (mst[i] == false && key[i] < min) {

            // reassign min and min_index to the values from this node
            min = key[i];
            min_index = i;
        }
    }
    return min_index;
}

/******************************************************************************
  find all vertices of odd degree in the MST. Store them in an subgraph O
******************************************************************************/

void TSP::findOdds() {
    for (int i = 0; i < n; i++) {

        // if degree of vertex i is odd
        if ((adjlist[i].size() % 2) != 0) {

            // push vertex to odds, which is a representation of subgraph O
            odds.push_back(i);
        }
    }
}

//---------------------------------------------------------------------------

void TSP::perfectMatching() {
    /************************************************************************************
     find a perfect matching M in the subgraph O using greedy algorithm but not minimum
    *************************************************************************************/
    int closest;
    double length; //int d;
    std::vector<int>::iterator tmp, first;

    // Find nodes with odd degrees in double to get subgraph O
    findOdds();

    // for each odd node
    while (!odds.empty()) {
        first = odds.begin();
        std::vector<int>::iterator it = odds.begin() + 1;
        std::vector<int>::iterator end = odds.end();
        length = std::numeric_limits<int>::max();
        for (; it != end; ++it) {
            // if this node is closer than the current closest, update closest and length
            if (graph[*first][*it] < length) {
                length = graph[*first][*it];
                closest = *it;
                tmp = it;
            }
        } // two nodes are matched, end of list reached
        adjlist[*first].push_back(closest);
        adjlist[closest].push_back(*first);
        odds.erase(tmp);
        odds.erase(first);
    }
}

//---------------------------------------------------------------------------
//find an euler circuit

void TSP::euler_tour(int start, std::vector<int> &path) {
    //Create copy of adj. list
    std::vector<std::vector<int> > tempList;

    for (int i = 0; i < n; i++) {
        tempList.push_back(adjlist[i]);
    }

    std::stack<int> stack;
    int pos = start;
    path.push_back(start);
    while (!stack.empty() || tempList[pos].size() > 0) {
        //Current node has no neighbors
        if (tempList[pos].empty()) {
            //add to circuit
            path.push_back(pos);
            //remove last vertex from stack and set it to current
            pos = stack.top();
            stack.pop();
        }
            //If current node has neighbors
        else {
            //Add vertex to stack
            stack.push(pos);
            //Take a neighbor
            int neighbor = tempList[pos].back();
            //Remove edge between neighbor and current vertex
            tempList[pos].pop_back();
            for (unsigned int i = 0; i < tempList[neighbor].size(); i++) {
                if (tempList[neighbor][i] == pos) {
                    tempList[neighbor].erase(tempList[neighbor].begin() + i);
                }
            }
            //Set neighbor as current vertex
            pos = neighbor;
        }
    }
    path.push_back(pos);
}

//---------------------------------------------------------------------------
//Make euler tour Hamiltonian

void TSP::make_hamiltonian(std::vector<int> &path, double &pathCost) {
    //remove visited nodes from Euler tour
    bool *visited = new bool[n];

    for (int i = 0; i < n; i++) {
        visited[i] = 0;
    }

    pathCost = 0.;

    int root = path.front();
    std::vector<int>::iterator cur = path.begin();
    std::vector<int>::iterator iter = path.begin() + 1;
    visited[root] = 1;

    //iterate through circuit
    while (iter != path.end()) {
        if (!visited[*iter]) {
            pathCost += graph[*cur][*iter];
            cur = iter;
            visited[*cur] = 1;
            iter = cur + 1;
        } else {
            iter = path.erase(iter);
        }
    }

    delete[] visited;

    //Add distance to root
    pathCost += graph[*cur][*iter];
}

//---------------------------------------------------------------------------

double TSP::findBestPath(int start) {
    std::vector<int> path;
    euler_tour(start, path);

    double length;
    make_hamiltonian(path, length);

    return length;
}

//---------------------------------------------------------------------------

void TSP::printResult() {
    for (std::vector<int>::iterator it = circuit.begin(); it != circuit.end(); ++it) {
        std::cout << *it << std::endl;
    }
}

//---------------------------------------------------------------------------

void TSP::printPath() {
    std::cout << std::endl;
    for (std::vector<int>::iterator it = circuit.begin(); it != circuit.end() - 1; ++it) {
        std::cout << *it << " to " << *(it + 1) << " ";
        std::cout << graph[*it][*(it + 1)] << std::endl;
    }
    std::cout << *(circuit.end() - 1) << " to " << circuit.front();
    std::cout << "\nLength: " << pathLength << std::endl << std::endl;
}

//---------------------------------------------------------------------------

void TSP::printAdjList() {
    for (int i = 0; i < n; i++) {
        std::cout << i << ": "; //print which vertex's edge list follows
        for (unsigned int j = 0; j < adjlist[i].size(); j++) {
            std::cout << adjlist[i][j] << " "; //print each item in edge list
        }
        std::cout << std::endl;
    }
}
//---------------------------------------------------------------------------