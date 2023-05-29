#ifndef __lez9__
#define __lez9__

#include <armadillo>
#include <iostream>
#include <cmath>
#include "random.h"

using namespace std;
using namespace arma;

// Check function with efficiency O(n)
void Check(vec v)
{
    vec appo(v.n_elem); // elements are initialized to zero automatically
    for (int i = 0; i < (int)v.n_elem; i++)
    {
        appo[v[i]]++;
    }

    // If v follows bounds every element of appo should have been incremented once and only once
    for (int i = 0; i < (int)v.n_elem; i++)
    {
        if (appo(i) == 0)
        {
            cout << "City number " << i + 2 << " does not appear!" << endl;
        }
        if (appo(i) > 1)
        {
            cout << "City number " << i + 2 << " appears more than once!" << endl;
        }
    }
}

// Swaps the a and b elements
void Swap(vec &v, int a, int b)
{
    double mute = v[a];
    v[a] = v[b];
    v[b] = mute;
}

void Mutate1(vec &v, Random &rand)
{
    if (v.n_elem < 2)
    {
        cout << "Error: mutated vector is too small!" << endl;
        return;
    }

    int p1 = rand.Discrete(0, v.n_elem - 1);
    int p2 = rand.Discrete(0, v.n_elem - 1);

    Swap(v, p1, p2);
}

void Mutate2(vec &v, Random &rand)
{
    int vsize = v.n_elem;
    if (vsize < 2)
    {
        cout << "Error: mutated vector is too small!" << endl;
        return;
    }

    int N = rand.Discrete(1, vsize - 1); // Number of elements to shift forward
    // cout << "N = " << N << endl;
    int m = rand.Discrete(1, vsize - 1); // Number of places of which they are shifted
    m = 3;
    // cout << "m = " << m << endl;
    int p1 = rand.Discrete(0, vsize - N); // Position of first numbr of antiguous group
    p1 = 3;
    // cout << "p1 = " << p1 << endl;

    vec vnew = v;
    // Shifts the elements one step at a time m times
    for (int i = 0; i < m; i++)
    {

        // Shifts forward the elements one by one
        for (int j = N - 1; j >= 0; j--)
        {

            // Determines the lement to shift forward
            int moved = p1 + j + i;
            if (moved > vsize - 1)
                moved -= vsize;

            // Checks if it is shifting forward the last element and shifts
            if (moved == vsize - 1)
            {
                Swap(v, moved, 0);
                // cout << "Swapping " << moved << " and " << 0 << endl;
            }
            else
            {
                Swap(v, moved, moved + 1);
                // cout << "Swapping " << moved << " and " << moved+1 << endl;
            }
        }
    }
}

void Mutate3(vec &v, Random &rand)
{
    if (v.n_elem < 2)
    {
        cout << "Error: mutated vector is too small!" << endl;
        return;
    }

    int N = rand.Discrete(1, v.n_elem / 2); // Size of two groups
    // cout << "N = " << N << endl;
    int p1 = rand.Discrete(0, v.n_elem - 2 * N); // Position of first value of first group
    // cout << "p1 = " << p1 << endl;
    int p2 = rand.Discrete(p1 + N, v.n_elem - N); // Position of first value of second group
    // cout << "p2 = " << p2 << endl;
    vec v1 = v.subvec(p1, p1 + N - 1);
    vec v2 = v.subvec(p2, p2 + N - 1);
    for (int i = 0; i < N; i++)
    {
        v[p1 + i] = v2[i];
        v[p2 + i] = v1[i];
    }
}

void Mutate4(vec &v, Random &rand)
{
    if (v.n_elem < 2)
    {
        cout << "Error: mutated vector is too small!" << endl;
        return;
    }

    int N = rand.Discrete(1, v.n_elem);
    // cout << "N = " << N << endl;
    int p1 = rand.Discrete(0, v.n_elem - 1 - N);
    // cout << "p1 = " << p1 << endl;

    vec v1 = arma::reverse(v.subvec(p1, p1 + N - 1));
    for (int i = 0; i < N; i++)
    {
        v[p1 + i] = v1[i];
    }
}

void Mutate( vec &v, double p1, double p2, double p3, double p4, Random &rand ) {
    double ran = rand.Rannyu();
    if( ran < p1 ) {
        Mutate1(v, rand);
    } else if( ran < p1 + p2 ) {
        Mutate2(v, rand);
    } else if( ran < p1 + p2 + p3 ) {
        Mutate3( v, rand);
    }else if( ran < p1 + p2 + p3 + p4 ) {
        Mutate4(v, rand);
    }

}

void Crossover(vec &v1, vec &v2, Random &rand)
{
    int vsize = v1.n_elem;
    if ((int)v2.n_elem != vsize)
    {
        cout << "Error: crossed vectors are different sizes!" << endl;
        return;
    }

    int p1 = rand.Discrete(1, vsize - 1); // Position of first number to be swapped
    // cout << "p1 = " << p1 << endl;
    int ncut = vsize - p1; // Number of elements that are cut

    vec v1old = v1;
    vec v2old = v2;

    int count1 = 0;
    int count2 = 0;
    // Cycles over all the elements
    for (int i = 0; i < vsize; i++)
    {
        // Checks if they correspond to any of the cut away ones
        for (int j = 0; j < ncut; j++)
        {

            if (v1old[i] == v2old[p1 + j])
            {
                // Puts the element in the first free position
                v2[p1 + count2] = v1old[i];
                count2++;
            }

            if (v2old[i] == v1old[p1 + j])
            {
                // Puts the element in the first free position
                v1[p1 + count1] = v2old[i];
                count1++;
            }
        }
    }
}

int Selection(vec lenghts, Random &rand)
{
    // Fitness is inversely proportional to lenght
    vec fitness = pow( lenghts, -1);
    double ftot = arma::sum(fitness);
    double r = rand.Rannyu(0, ftot);
    double fdym = 0;
    for (int i = 0; i < (int)fitness.n_elem; i++)
    {
        fdym += fitness[i];
        if (r < fdym)
            return i;
    }
    cout << "Error with selection!" << endl;
    return datum::nan;
}

// Generates a point located on the unitary circumeference with center (1, 1)
rowvec CityCircle(Random &rand)
{
    rowvec city(2);
    double theta = rand.Rannyu(0, 2 * M_PI);
    city(0) = sin(theta);
    city(1) = cos(theta);
    return city;
}

// Genarates a point located on the square with unitary side and center (1/2, 1/2)
rowvec CitySquare(Random &rand)
{
    rowvec city(2);
    city[0] = rand.Rannyu();
    city[1] = rand.Rannyu();
    return city;
}

// Generates the distance between cities matrix
mat Distances(const mat &Cities)
{
    int n = Cities.n_rows;

    mat result(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            result(i, j) = arma::norm(Cities.row(i) - Cities.row(j), 2);
        }
    }

    return result;
}

// Generates a path
vec Path(int n_cities, Random &rand)
{
    vec result = arma::regspace(0, n_cities - 2);
    for (int i = 0; i < 10 * n_cities; i++)
    {
        Mutate1(result, rand);
        Mutate2(result, rand);
        Mutate3(result, rand);
        Mutate4(result, rand);
    }
    Check(result);
    return result;
}

// Generates the population
mat Population(int n_cities, int n_paths, Random &rand)
{
    mat result = Path(n_cities, rand);
    for (int i = 1; i < n_paths; i++)
    {
        result = arma::join_horiz(result, Path(n_cities, rand));
    }
    return result;
}

// Calculates the lenght of a path
double lenght(vec v, const mat &distance)
{
    int vlen = v.n_elem;
    if (vlen + 1 != (int)distance.n_rows)
    {
        cout << "Error: path and distance have different number of cities!" << endl;
    }
    double result = 0;
    result += distance(v[0], vlen) + distance(v[vlen - 1], vlen);
    for (int i = 0; i < vlen - 1; i++)
    {
        result += distance(v[i], v[i + 1]);
    }
    return result;
}

// Calculates the lenght of path vector
vec lenghts(mat population, const mat &distance)
{
    int n_paths = population.n_cols;
    vec result(n_paths);
    for (int i = 0; i < n_paths; i++)
    {
        result(i) = lenght(population.col(i), distance);
    }
    return result;
}

// Sorts the population from shortest to longest
void LSort(mat &population, const mat &distance)
{

    // Create an index vector to store the sorting order
    arma::uvec sortedIndices = arma::sort_index(lenghts(population, distance));

    // Rearrange the columns of the matrix based on the sorting order
    population = population.cols(sortedIndices);
}

// Returns a slected path
vec Selection(mat population, const mat &distance, Random &rand)
{
    int index = Selection(lenghts(population, distance), rand);
    return population.col(index);
}

// Elitary GA, New Generation Evolution
void NGE(mat &population, mat &distance, Random &rand)
{
    // Sorts best to worst
    LSort(population, distance);
    int n_paths = population.n_cols;
    mat oldgen = population;

    // Keeps the best 20%, replaces the rest with offsprings
    for ( int i = n_paths / 5; i < n_paths; i+=2 ) {

        // Selects, crosses and possibly mutates two paths
        vec v1 = Selection( oldgen, distance, rand );
        vec v2 = Selection( oldgen, distance, rand );
        Crossover(v1, v2, rand);
        Mutate( v1, 0.02, 0.02, 0.02, 0.02, rand );
        Mutate( v2, 0.02, 0.02, 0.02, 0.02, rand );
        Check(v1);
        Check(v2);

        // Substitutes the new paths in the new population
        population.col(i) = v1;
        if( i+1 < n_paths ) population.col(i+1) = v2;
    }

}

double Avlen( mat pop, const mat &dist) {
    vec half = lenghts(pop, dist);
    half.set_size(half.n_elem/2);
    return arma::mean( half );
}

// Positions of cities ordered in the best path
mat BestPath( mat pop, const mat &dist, mat cities ){
    LSort(pop, dist);

    mat result = cities.row(cities.n_rows-1);
    for( int i = 0; i < (int) pop.n_rows; i++ ){
        int index = pop(i, 0);
        result = arma::join_vert( result, cities.row(index) );
    }
    return result;
}

#endif