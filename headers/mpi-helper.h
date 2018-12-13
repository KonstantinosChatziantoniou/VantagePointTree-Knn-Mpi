#ifndef MPI_HELPER_H
#define MPI_HELPER_H

#include <mpi.h>


void slavePart(int processId,int partLength,float *numberPart,int size,MPI_Comm group, float* points);

float masterPart(int noProcesses,int processId,int size,int partLength,float *numberPart,MPI_Comm group,float* points);

void validation(float median,int partLength,int size,float *numberPart,int processId,MPI_Comm group);

void sendLengths(int size,int noProcesses , MPI_Comm group);


#endif