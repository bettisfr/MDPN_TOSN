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
#include <stdio.h>
#include <vector>
#include <limits>

struct TSP_Point {
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
    std::vector<TSP_Point> points;

    // -
    void findOdds();

    //Find perfect matching
    void perfectMatching();

    //Find Euler tour
    void euler_tour(int start, std::vector<int> &path);

    //Find Hamiltonian path
    void make_hamiltonian(std::vector<int> &path, double &pathCost);

    double get_distance(struct TSP_Point c1, struct TSP_Point c2);

    //T get_distance(struct TSP_Point c1, struct TSP_Point c2);

    void findMST();

    int getMinIndex(double key[], bool mst[]);

    void fillMatrix();

    double findBestPath(int start);

public:
    TSP(std::vector<TSP_Point> aPointList);

    ~TSP();

    void Solve();

    void printResult();

    void printPath();

    void printAdjList();

    std::vector<int> getResult();
};

#endif
