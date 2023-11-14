#ifndef TSP_H
#define TSP_H

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <stack>
#include <string>
#include <cstdio>
#include <vector>
#include <limits>

using namespace std;

struct point_3d {
    double x;
    double y;
    double z;
};

//---------------------------------------------------------------------------
class TSP {
public:


protected:
    // List of odd nodes
    std::vector<int> odds;

    //Adjacency list
    std::vector<int> *adjlist;

    // Number of points
    int n;

    //Shortest path length
    double pathLength;

    //euler circuit
    std::vector<int> circuit;

    // n x n, pairwise distances between points
    double **graph;

    // Point list
    std::vector<point_3d> points;

    // -
    void findOdds();

    //Find perfect matching
    void perfectMatching();

    //Find Euler tour
    void euler_tour(int start, std::vector<int> &path);

    //Find Hamiltonian path
    void make_hamiltonian(std::vector<int> &path, double &pathCost);

    double get_distance(struct point_3d c1, struct point_3d c2);

    //T get_distance(struct point_3d c1, struct point_3d c2);

    void findMST();

    int getMinIndex(const double key[], const bool mst[]);

    void fillMatrix();

    double findBestPath(int start);

public:
    TSP(std::vector<point_3d> aPointList);

    ~TSP();

    void solve();

    void printResult();

    void printPath();

    void printAdjList();

    std::vector<int> getResult();
};

#endif
