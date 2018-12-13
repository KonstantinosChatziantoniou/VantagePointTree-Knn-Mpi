#ifndef MEDIAN_FUNCS_H
#define MEDIAN_FUNCS_H

void _set_dimensions_for_indexing(int dimensions);

int _index(int x ,int y);

void generatePoints(float* arr , int size , int seed);

void printPoints(float* arr , int size , int rank);

void calculateDistances(float* points , int size ,float* vp, float* numberPart);

void swap_values(float *array,int x,int y,float* points);

void swap_points(float** array, int x , int y, int dimensions);

void removeElement(int *array, int *size, int element);

void validationST(float median,int size,float *numberPart, int* out_countMin, int* out_countMax , int* out_countEq);

void partition (float *array, int elements, float pivot, float **arraysmall, float **arraybig, int *endsmall, int *endbig,float* points);

float selection(float *array,int number,float* points);


#endif