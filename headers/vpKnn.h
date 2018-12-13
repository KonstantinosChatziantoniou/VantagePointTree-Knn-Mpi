#ifndef vpKNN_H
#define vpKNN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SMALLFLOATTHRESSHOLD 0.0000001


typedef struct{
    int k;
    int maxLevel;
    int dimensions;
    float *points;
    float* medians;
    float* nbrs;
    float* vp;
    float tau;
    int nbrsLength;

}knnVars;



knnVars* _initForKnnSearch(int k, int maxLevel , int dimensions , float* points , float* medians , int vpIndex);

void _finishKnnSearch(knnVars* self);
void startSearchKnn(knnVars* self, int level , int index);
float findDistance(knnVars* self,float* point , float* vp);
void addNeighbour(knnVars* self , float* point ,float dist);

float* getNeigbours(knnVars* self);

#endif