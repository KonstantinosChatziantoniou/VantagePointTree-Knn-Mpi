#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../headers/vpKnn.h"
#include <math.h>
#include <limits.h>




knnVars* _initForKnnSearch(int k, int maxLevel , int dimensions , float* points , float* medians , int vpIndex){
    knnVars* v =(knnVars*)malloc(sizeof(knnVars));
    v->k = k; 
    v->maxLevel = maxLevel;
    v->points = points;
    v->medians = medians;
    v->dimensions = dimensions;
    v->nbrs = (float*)malloc((dimensions+1)*sizeof(int)*k);
    v->vp = (float*)malloc(sizeof(float) * dimensions);
    memcpy(v->vp , &points[vpIndex*dimensions] , sizeof(float)*dimensions);
    v->nbrsLength = 0;
    v->tau = 100000000;
    return v;
}

void _finishKnnSearch(knnVars* self){
    free(self->vp);
    free(self->nbrs);
    free(self);
}

void startSearchKnn(knnVars* self, int level , int index){
    if(level > self->maxLevel){
        return;
    }

    float dist = findDistance(self ,&self->points[index*self->dimensions] , self->vp);
    //printf("dist tau %f   %f\n",dist , self->tau);
    if(dist > SMALLFLOATTHRESSHOLD){
        if(dist < self->tau){
            addNeighbour(self, &self->points[index*self->dimensions] , dist);
        }        
    }
    //printf("dist median %f %f\n",dist , self->medians[index]);
    if(dist <= self->medians[index]){
        startSearchKnn(self , level+1 , index + 1);
        if(dist + self->tau >= self->medians[index]){
            startSearchKnn(self , level+1 , index + pow(2,self->maxLevel - level));
        }
    }else{
        startSearchKnn(self , level+1 , index + pow(2,self->maxLevel - level));
        if(dist - self->tau <= self->medians[index]){
            startSearchKnn(self , level+1 , index + 1);
        }
    }
    

}

float findDistance(knnVars* self,float* point , float* vp){
    float dist = 0;
    for(int j = 0; j < self->dimensions; j++){
        dist += pow(point[j] - vp[j] , 2);
    }
    dist = sqrt(dist);
    return dist;
}

void addNeighbour(knnVars* self , float* point ,float dist){
    float max = self->nbrs[self->dimensions];
    int indexMax = 0;
    //printf("--------adding ---------\n");
    //----- if the neighbours list is not full add anyway -----//
    if(self->nbrsLength < self->k){

            //printf("--------adding --------- %f\n",dist);
        memcpy(&self->nbrs[self->nbrsLength*(self->dimensions+1)] , point , sizeof(float)*self->dimensions);
        self->nbrs[self->nbrsLength*(self->dimensions+1) + self->dimensions] = dist;
        self->nbrsLength++;
        //return;
    }else{
        //----- if it is full, swap the given point with the furthest neighbour -----//
        max = self->nbrs[self->dimensions];
        indexMax = 0;
        for(int i = 1; i < self->nbrsLength; i++){
            if(self->nbrs[i*(self->dimensions+1) + self->dimensions] > max){
                indexMax = i;
                max = self->nbrs[i*(self->dimensions+1) + self->dimensions];
            }
        }
        //---- set Tau equal to max ----//
        if(dist < max){
            memcpy(&self->nbrs[indexMax*(self->dimensions+1)] , point , sizeof(float)*self->dimensions);
            self->nbrs[indexMax*(self->dimensions+1) + self->dimensions] = dist;
        }
    }
    if(self->k == self->nbrsLength){
        max = self->nbrs[self->dimensions];
        indexMax = 0;
        for(int i = 0; i < self->nbrsLength; i++){
           // printf("finding max %f\n",self->nbrs[i*(self->dimensions) + self->dimensions]);
            if(self->nbrs[i*(self->dimensions+1) + self->dimensions] > max){
                indexMax = i;
                max = self->nbrs[i*(self->dimensions+1) + self->dimensions];
            }
        }
    self->tau = max;
    }
}

float* getNeigbours(knnVars* self){

        //printf("tau = %f \n",self->tau);
        printf("vp");
    for(int j = 0; j < self->dimensions; j++){
            printf("x%d %f\t",j,self->vp[j]);
        }
        printf("\n");
    for(int i = 0; i < 3; i++){
        printf("neighbour %d: \t",i);
        for(int j = 0; j < self->dimensions; j++){
            printf("x%d %f\t",j,self->nbrs[i*(self->dimensions+1)+j]);
        }
        printf("distance = %f\n",self->nbrs[i*(self->dimensions+1)+self->dimensions]);
    }
    return self->nbrs;
}



