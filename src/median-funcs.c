#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>


int _dimensions;


void generatePoints(float* arr , int size , int seed){
    srand(seed);
    for(int i = 0; i < size; i++){
       for(int j = 0;  j < _dimensions; j++){
            arr[i*(_dimensions) + j] = ((float)rand()/(float)RAND_MAX)*10; 
        }
    }
}

void printPoints(float* arr , int size , int rank){
    for(int i = 0 ; i < size ; i++){
        printf("Point of rank: %d -->",rank);
        for(int j = 0; j < _dimensions; j++ ){
            printf("x%d = %f \t",j,arr[i*_dimensions + j]);

        }
        printf("\n");
    }
}

void calculateDistances(float* points , int size ,float* vp, float* numberPart){
    for(int i = 0; i < size; i++){
        numberPart[i] = 0;
        for(int j = 0; j < _dimensions; j++){
            numberPart[i] += pow(vp[j] - points[i*_dimensions + j], 2);
        }

    numberPart[i] = sqrt(numberPart[i]);;
    }
}
void _set_dimensions_for_indexing(int dimensions){

    _dimensions = dimensions; 
}

int _index(int x ,int y){
    return x*(_dimensions) + y;
}

void swap_values(float *array,int x,int y,float* points)
{
    float temp;
    temp=array[x];
    array[x]=array[y];
    array[y]=temp;
    float* temp2 = (float*)malloc(sizeof(float)*_dimensions);
    memcpy(temp2 , &points[x*_dimensions] ,sizeof(float)*_dimensions);
    memcpy(&points[x*_dimensions] , &points[y*_dimensions] , sizeof(float)*_dimensions); 
    memcpy(&points[y*_dimensions] ,temp2 , sizeof(float)*_dimensions);
    free(temp2);
}


void swap_points(float** array, int x , int y, int dimensions){
    float* temp =(float*)malloc(dimensions*sizeof(float));
    memcpy(temp , array[x],dimensions*sizeof(float));
    memcpy(array[x],array[y],dimensions*sizeof(float));
    memcpy(array[y],temp,dimensions*sizeof(float));
    free(temp);
}


/***Kills processes that have no values left in their arrays****/
void removeElement(int *array, int *size, int element)
{
    int i;
    int flag=0;
    for(i=0;i<*size;i++)
    {
        if(flag==1)
            array[i]=array[i+1];
        if(array[i]==element&& flag==0)
        {
            array[i]=array[i+1];
            flag=1;
        }
    }
    *size=*size-1;
}

//void removePoint(float* array , int *size , int index)
/***Validates the stability of the operation (Single Threaded)****/
void validationST(float median,int size,float *numberPart, int* out_countMin, int* out_countMax , int* out_countEq)
{
	int countMin=0;
    int countMax=0;
    int countEq=0;
    int i;
    for(i=0;i<size;i++)
    {
        if(numberPart[i]>median)
            countMax++;
        else if(numberPart[i]<median)
            countMin++;
        else
            countEq++;
    }
  /*if((countMax<=size/2)&&(countMin<=size/2))  //Checks if both the lower and higher values occupy less than 50% of the total array.
        printf("VALIDATION PASSED!\n");
    else
        printf("VALIDATION FAILED!\n");

	printf("Values greater than median: %d\n",countMax);
        printf("Values equal to median: %d\n",countEq);
        printf("Values lower than median: %d\n",countMin);
*/
    *out_countEq = countEq;
    *out_countMax = countMax;
    *out_countMin = countMin;
}
/*========================FIND MEDIAN FUNCTIONS====================================
 * ================================================================================
 * ================================================================================
*/


/****Partitions the Array into larger and smaller than the pivot values****/
void partition (float *array, int elements, float pivot, float **arraysmall, float **arraybig, int *endsmall, int *endbig, float* points)
{
    int right=elements-1;
    int left=0;
    int pos;
    if(elements==1)
    {
        if(pivot>array[0])
        {
            *endsmall=1;  //One value in the small part
            *endbig=0;   //Zero on the big one
            *arraysmall=array;   //There is no big array therefore NULL value
            *arraybig=NULL;
        }
        else if(pivot<=array[0])
        {
            *endsmall=0;    //The exact opposite of the above actions.
            *endbig=1;
            *arraysmall=NULL;
            *arraybig=array;
        }
    }
    else if(elements>1)
    {
        while(left<right)
        {
            while(array[left]<pivot)
            {
                left++;
                if(left>=elements)
                {
                    break;
                }
            }
            while(array[right]>=pivot)
            {
                right--;
                if(right<0)
                {
                    break;
                }
            }
            if(left<right)
            {
                swap_values(array,left,right,points);
            }
        }
        pos=right;
        if(pos<0)                   //Arrange the arrays so that they are split into two smaller ones.
        {                               //One containing the small ones. And one the big ones.
            *arraysmall=NULL;   
            //printf("NUUUUUUUUUUUUUULLLLLLLLLLLLLLL = %d\n",pos);        //However these arrays are virtual meaning that we only save the pointer values of the beging and end
        }                               //of the "real" one.
        else
        {
            *arraysmall=array;
        }
        *endsmall=pos+1;
        *arraybig=&array[pos+1];
        *endbig=elements-pos-1;
    }
}


/***==============================================***/
/***==============================================***/
/***=============SERIAL SELECTION==============***/
/***==============================================***/
/***==============================================***/

float selection(float *array,int number,float* points)
{
    //printf("SELECTION: length %d\n",number);
    //for(int i = 0; i < number; i++){
    //    printf("SELECTION: DISTSANCES \t %f\n",array[i]);
    //}
    //for(int i = 0; i < number; i++){
    //    printf("SELECTION: POINTS \t %d",i);
    //    for(int j = 0; j <_dimensions; j++){
    //        printf("%f\t",points[i*_dimensions+j]);
    //    }
    //    printf("\n");
    //}
    float *arraybig;
    float *arraysmall;
    int endsmall=0;
    int endbig=0;
    float *arraytobeused;
    int i;
    int counter=0;
    int k;
    float pivot;
    float median;
    k=(int)number/2+1;
    arraytobeused=array;
    for(;;)
    {
        pivot=arraytobeused[rand() % number];

        //printf("SELECTION pivot = %f\n",pivot);
        partition(arraytobeused,number,pivot,&arraysmall,&arraybig,&endsmall,&endbig,points);
        if(endbig>k)
        {
            number=endbig;
            arraytobeused=arraybig;
            for(i=0;i<endbig;i++)
            {
                if(pivot==arraybig[i])
                    counter++;
                else
                    break;
            }
            if(counter==endbig)
            {
                median=arraybig[0];
                break;
            }
            else
                counter=0;
            //end of count equals
        }
        else if(endbig<k)
        {
            number=endsmall;
            arraytobeused=arraysmall;
            k=k-endbig;
        }
        else
        {
            median=pivot;
            break;
        }
    }
    //printf("SELECTION MEDIA = %f\n",median);
    return median;
}

