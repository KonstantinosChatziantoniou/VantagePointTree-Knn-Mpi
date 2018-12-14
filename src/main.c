#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include  <time.h>
#include <sys/time.h>

#include "../headers/median-funcs.h"
#include "../headers/mpi-helper.h"
#include "../headers/vpKnn.h"

void changeMpiTreeStructure(float* oldTree , float* newTree , int levels , int c_lvl , int parentOld , int parentNew ,int demensions);
void STVantagePointTree(float* pointArr ,float* numberPart , int length , int dimensions);

void PrintToCSV(float* points , float* medians , int length , int dimensions , int rank);

void PrintMpiTreeToCsv(float* tree , int length , int dimensions);

void assignMpiLeavesToProcs(int level , int maxLevelplus ,  int index , int* counter , int* procs);


void printPointsToCsv(float* points ,int length , int dimensions , float *medians , int rank , char* args);


void checkKnnGlobally(int currentPos, int trlen , int* Status, int* nodeAdded, int* nodeChecked ,float* vp , float* mpiTree , int dimensions, knnVars* obj);
  
int isLeaf(int pos, int trlen);
void checkIfChildrenChecked(int pos, int* nodeChecked);
int main(int argc , char **argv){

    int processId,noProcesses;
    int size, partLength;
    //-------------- For Time Measuring -------------------------//
    FILE* fileTime = fopen("time.csv","w");
    struct timeval startVal, endVal;
    struct timezone tz;
    //-------------- Every process uses these vars --------------//
    float median;
    float *pointArr;
    float *vp;
    float *mpiTreeSaver;
    int mpiTreeCounter;
    float *numberPart;
    float* points_to_receive;
    float *mediansTree;

    size=atoi(argv[1]);
    size = pow(2 , size) - 1;
    int dimensions = atoi(argv[2]);
    //int randomGenOrReadFile = atoi(argv[3]);
    _set_dimensions_for_indexing(dimensions);

    //printPointsToCsv(points_to_receive , -1 , -1 ,points_to_receive , processId , "w");

    //------------------------ Start MPI -------------------------//
    MPI_Init (&argc, &argv);	/* starts MPI */
    MPI_Comm_rank (MPI_COMM_WORLD, &processId);	/* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &noProcesses);	/* get number of processes */
   // printf("Problem size = %d, Number of processes = %d\n",size , noProcesses);



    //------- Single Threaded if only 1 process is available ------//
    if(noProcesses == 1){
        partLength = size;
        pointArr = (float*)malloc(partLength*dimensions*sizeof(float));
        /*if(randomGenOrReadFile == 0){
           
        }else{
            File data = fopen("../data/data.csv" , "r");
            for(int i = 0; i < size; i++){
                for(int j = 0; j < dimensions; j++ ){
                    fscanf(data , "%f",&pointArr[i*dimensions + j]);
                }
            }
        }*/
        generatePoints(pointArr , partLength , processId + 1);
        mediansTree = (float*)malloc(sizeof(float)*partLength);
        STVantagePointTree(pointArr , mediansTree , partLength , dimensions);
        PrintToCSV(pointArr , mediansTree , partLength , dimensions , processId);
        //-------------------------- find knn ---------------------------//
        knnVars*  obj;
        obj = _initForKnnSearch(3 , log2(partLength + 1) - 1 , dimensions , pointArr , mediansTree , 0);
        startSearchKnn(obj , 0 , 0);

        float* nbrs = getNeigbours(obj);
        
        for(int i = 0; i < 3; i++){
           // printf("neighbour %d: \t",i);
            for(int j = 0; j < dimensions; j++){
               // printf("x%d %f\t",j,nbrs[i*(dimensions+1)+j]);
            }
          //  printf("distance = %f\n",nbrs[i*(dimensions+1)+dimensions]);
        }
            _finishKnnSearch(obj);
            free(pointArr);
            free(mediansTree);
            MPI_Finalize();
            return 0;
    }


     //----- Generate Random Points -----//
    if(processId == 0){
        //----- Calculate the partLength og each process and send it -----//
        if(size%noProcesses==0)
            partLength=(size/noProcesses);
        else
            partLength=(size/noProcesses)+1;
        sendLengths(size,noProcesses,MPI_COMM_WORLD);   
    }else
    {
        MPI_Recv(&partLength,1,MPI_INT,0,1,MPI_COMM_WORLD,NULL);        
    }

    mpiTreeSaver = (float*)malloc((log2(noProcesses))*(1+dimensions)*sizeof(float)); //TREE_LEN_TAG
    pointArr = (float*)malloc(partLength*dimensions*sizeof(float));
    generatePoints(pointArr , partLength , processId + 1);
    //printPoints(pointArr , partLength , processId);

     //-----to split processes to work on different sets-----//
    MPI_Comm group_comm;
    int group_rank = processId;
    int group_size = noProcesses;
    int group_point_size;
    //int level = noProcesses;
    //ITER_TAG
    gettimeofday(&startVal, &tz);
    for(int level = noProcesses; level > 1; level = level/2){

         //----- Split the communicators -----//
        //-- each subcommunicator will act on a different set --//
        int color = processId / level;
        MPI_Comm_split(MPI_COMM_WORLD, color, processId, &group_comm);
        MPI_Comm_rank(group_comm, &group_rank);
        MPI_Comm_size(group_comm, &group_size);
        int length_received;


        //---- Calculate the size ----//
        MPI_Allreduce(&partLength , &group_point_size , 1 , MPI_INT , MPI_SUM ,  group_comm);


        vp = (float*)malloc(dimensions*sizeof(float));
        //----- Choose Vantage Point and send it to slaves -----//
        if(group_rank == 0){
            memcpy(vp , &pointArr[(partLength-1)*dimensions] , dimensions*sizeof(float));
           // printf("Vantage Point%d ",processId);
            for(int j = 0; j < dimensions; j++){
               // printf("x%d %f\t",j,vp[j]);
            }
         //   printf("\n");
            partLength--;
            group_point_size--;
            MPI_Bcast(vp , dimensions , MPI_FLOAT , 0 , group_comm);            
        }else{
            group_point_size--;
            MPI_Bcast(vp , dimensions , MPI_FLOAT , 0 , group_comm);
        }
        //----- add vp to tree -----//
        memcpy(&mpiTreeSaver[mpiTreeCounter*(dimensions + 1)] , vp ,dimensions*sizeof(float));
        mpiTreeCounter++;
        numberPart = (float*)malloc(sizeof(float)*partLength);
        calculateDistances(pointArr , partLength , vp , numberPart);
        printPointsToCsv(pointArr , partLength , dimensions ,numberPart , processId , "a");

         //-------------------- Find the median ---------------------//
        if(group_rank == 0)
        {
            median = masterPart(group_size , group_rank , group_point_size , partLength , numberPart , group_comm , pointArr );
            MPI_Bcast(&median , 1 , MPI_FLOAT , 0 , group_comm);
           // printf("-------------> MEDIAN = %f\n",median);
        }else{
            slavePart(group_rank , partLength , numberPart , group_point_size , group_comm , pointArr);
            MPI_Bcast(&median , 1 , MPI_FLOAT , 0  , group_comm);
        }
            calculateDistances(pointArr , partLength , vp , numberPart);

        //----- add median to tree -----//
        mpiTreeSaver[mpiTreeCounter*(dimensions+1) - 1] = median;

        //---- Separate in 2 arrays, smaller and larger than median -----//
        float* arr_big;
        float* arr_small;
        int len_arr_big;
        int len_arr_small;
        partition(numberPart , partLength , median , &arr_small , &arr_big , &len_arr_small , &len_arr_big , pointArr);

        //printPoints(pointArr , partLength , group_rank);


        if(group_rank == 0){
             //----- Create an array to store hi and lo sizes from slaves -----//
            int** hiloarr = (int**)malloc(group_size*sizeof(int*));
            hiloarr[0] = (int*)malloc(5*sizeof(int));
            hiloarr[0][0] = len_arr_small;
            hiloarr[0][1] = len_arr_big;
            int countMax , countMin , countEq;
            validationST(median , partLength , numberPart , &countMin , &countMax , &countEq);
            hiloarr[0][2] = countEq;
            hiloarr[0][3] = countMin;
            hiloarr[0][4] = countMax;
            //------- Balance high an low arrs so they have the same length ---------//
            if(countEq == 1){
                for(int i = len_arr_small; i < partLength;i++){
                    if(numberPart[i] == median){
                        if(i == len_arr_small){
                            break;
                        }else{
                            swap_values(numberPart , len_arr_small , i , pointArr);
                        }
                    }
                    
                }
                len_arr_small++;
                len_arr_big--;
                hiloarr[0][0] = len_arr_small;
                hiloarr[0][1] = len_arr_big;
                arr_big = &numberPart[len_arr_small];
            }

            //-------- Receive high an low lengths from slaves --------//
            for(int j = 1; j < group_size; j++)
            {
                hiloarr[j] = (int*)malloc(5*sizeof(int));
                MPI_Recv(hiloarr[j],5,MPI_INT,j,0,group_comm,NULL);
            }/*
            printf("---------->rank: %d low: %d high: %d  eq: %d  cmin: %d   cmax: %d\n",processId,hiloarr[0][0] , hiloarr[0][1] , hiloarr[0][2], hiloarr[0][3], hiloarr[0][4]);
            for(int j = 1; j < group_size; j++)
            {
                printf("---------->rank: %d low: %d high: %d  eq: %d cmin: %d   cmax: %d\n",j,hiloarr[j][0] , hiloarr[j][1] , hiloarr[j][2], hiloarr[j][3], hiloarr[j][4]);
            }*/

            //----- Calculate how many points each slave should receive from others -----//
            int tmpSum = 0;
            int tmpMod = 0;
            for(int i = 0; i < group_size; i++){
                tmpSum += hiloarr[i][0];
            }
            
            tmpMod = tmpSum%(group_size/2);
            tmpSum = tmpSum/(group_size/2);
            int* lengths_to_receieve = (int*)malloc(sizeof(int)*group_size);
            for(int i = 0; i < group_size/2; i++){
                lengths_to_receieve[i] = tmpSum - hiloarr[i][0];
                if(tmpMod > 0){
                    lengths_to_receieve[i]++;
                    tmpMod--;
                }
            }
            length_received = lengths_to_receieve[0];
            tmpSum = 0;
            tmpMod = 0;
            for(int i = 0; i < group_size; i++){
                tmpSum += hiloarr[i][1];
            }

            tmpMod = tmpSum%(group_size/2);
            tmpSum = tmpSum/(group_size/2);
            for(int i = group_size/2; i < group_size; i++){
                lengths_to_receieve[i] = tmpSum - hiloarr[i][1];
                if(tmpMod > 0){
                    lengths_to_receieve[i]++;
                    tmpMod--;
                }
            }
            for(int i = 0; i < group_size; i++){
              //  printf("length to receive:%d   %d\n",i , lengths_to_receieve[i]);
            }

            //------ Send the lengths of the arr to expect for malloc -----//
            for(int i = 1; i < group_size; i++){
                MPI_Send(&lengths_to_receieve[i] , 1, MPI_INT , i , 0 , group_comm);
            }
            points_to_receive = (float*)malloc(sizeof(float)*dimensions*lengths_to_receieve[0]);
            for(int i = 0; i < group_size; i++){
                //printf("lengths to receieve %d %d\n",i,lengths_to_receieve[i]);
                
            }

            //----- Send instructions for sending and receiving points -----//
            int message[3];
            int jc = group_size/2;
            int tempHaveReceived = 0;
            for(int i = 0; i < group_size/2; i++){
                while(lengths_to_receieve[i] > 0){
                    if(lengths_to_receieve[i] >= hiloarr[jc][0]){
                        lengths_to_receieve[i] -= hiloarr[jc][0];
                        message[0] = 0; //0 -> send
                        message[1] = i; //to
                        message[2] = hiloarr[jc][0];
                        hiloarr[jc][0] = 0;
                        MPI_Send(message , 3 , MPI_INT, jc , 1 , group_comm);
                        if(i == 0){
                           // printf("started1 receiving rank %d points %d from %d\n",group_rank , message[2] , jc);
                            MPI_Recv(&points_to_receive[tempHaveReceived] , message[2]*dimensions , MPI_FLOAT , jc , 2, group_comm , NULL);
                           // printf("done1 receiving rank %d points %d from %d\n",group_rank , message[2] , jc);
                            tempHaveReceived += message[2]*dimensions;
                        }else{
                            message[0] = 1; //1 -> receive
                            message[1] = jc; //from
                            MPI_Send(message , 3 , MPI_INT , i , 1 , group_comm);
                          //  printf("asdasd22222222222222asdasdasdasdas\n");
                        }
                        jc++;
                    }else{
                        message[2] = lengths_to_receieve[i];
                        message[0] = 0;
                        message[1] = i;
                        hiloarr[jc][0] -= lengths_to_receieve[i];
                        lengths_to_receieve[i] = 0;
                        MPI_Send(message , 3 , MPI_INT , jc , 1 , group_comm);
                        if(i == 0){
                          //  printf("started2 receiving rank %d points %d from %d\n",group_rank , message[2] , jc);
                            MPI_Recv(&points_to_receive[tempHaveReceived] , message[2]*dimensions , MPI_FLOAT , jc , 2 , group_comm , NULL);
                            tempHaveReceived += message[2]*dimensions;
                          //  printf("done2 receiving rank %d points %d from %d\n",group_rank , message[2] , jc);
                        }else{
                            message[0] = 1;
                            message[1] = jc;
                            MPI_Send(message , 3 , MPI_INT , i , 1 , group_comm);
                           // printf("asdasdasdasdasdasdas\n");
                        }

                    }
                }
            }
            //-- signal that first transfer is over --//
            message[0] = -1;
            for(int i = 1; i < group_size; i++){
                MPI_Send(message , 3,  MPI_INT , i , 1 , group_comm);
            }
            //-- second part of transfer --//
            int tempHaveSent = 0;
            jc = 0;
            for(int i = group_size/2; i < group_size; i++){
                while(lengths_to_receieve[i] > 0){
                 //   printf("DEBUGGGGGGGGGGGGGGG i %d , length , %d , jc %d , lng %d\n",i,lengths_to_receieve[i],jc,hiloarr[jc][1]);
                    if(lengths_to_receieve[i] >= hiloarr[jc][1]){
                        message[0] = 1;
                        message[1] = jc;
                        message[2] = hiloarr[jc][1];
                        lengths_to_receieve[i] -= hiloarr[jc][1];
                        hiloarr[jc][1] = 0;
                        MPI_Send(message , 3 , MPI_INT , i , 1 , group_comm);
                        if(jc == 0){
                         //   printf("started22 sending rank %d points %d to %d\n",group_rank , message[2] , i);
                            MPI_Send(&pointArr[len_arr_small*dimensions +  tempHaveSent] , message[2]*dimensions, MPI_FLOAT , i , 2 , group_comm);
                         //   printf("done22 sending rank %d points %d to %d\n",group_rank , message[2] , i);
                            tempHaveSent += message[2]*dimensions;
                        }else{
                            message[0] = 0;
                            message[1] = i;
                            MPI_Send(message , 3 , MPI_INT , jc , 1 , group_comm);
                        }
                        jc++;
                    }else{
                        message[0] = 1;
                        message[1] = jc;
                        message[2] = lengths_to_receieve[i];
                        hiloarr[jc][1] -= lengths_to_receieve[i];
                        MPI_Send(message , 3 , MPI_INT , i , 1 , group_comm);
                        if(jc == 0){
                         //   printf("started22b sending rank %d points %d to %d\n",group_rank , message[2] , i);
                            MPI_Send(&pointArr[len_arr_small*dimensions +  tempHaveSent] , message[2]*dimensions , MPI_FLOAT , i , 2 , group_comm);
                          //  printf("done22b sending rank %d points %d to %d\n",group_rank , message[2] , i);
                            tempHaveSent += message[2]*dimensions;
                        }else{
                            message[0] = 0;
                            message[1] = i;
                            MPI_Send(message , 3 , MPI_INT , jc , 1 , group_comm);
                        }

                        lengths_to_receieve[i] = 0;
                    }
                }
            }
            //-- signal that second transfer is over --//
            message[0] = -1;
            for(int i = 1; i < group_size; i++){
                MPI_Send(message , 3,  MPI_INT , i , 1 , group_comm);
            }
            for(int i = 0; i < group_size; i++){
                free(hiloarr[i]);
            }
            free(hiloarr);
            free(lengths_to_receieve);

    
        }else{
            int hl[5];
            hl[0] = len_arr_small;
            hl[1] = len_arr_big;
            int countMax , countMin , countEq;
            validationST(median , partLength , numberPart , &countMin , &countMax , &countEq);
            hl[2] = countEq;
            hl[3] = countMin;
            hl[4] = countMax;
            //----- Balance high and low arrays to have the same length -----//
            if(countEq == 1){
                for(int i = len_arr_small; i < partLength;i++){
                    if(numberPart[i] == median){
                        if(i == len_arr_small){
                            break;
                        }else{
                            swap_values(numberPart , len_arr_small , i , pointArr);
                        }
                    }
                    
                }
                len_arr_small++;
                len_arr_big--;
                hl[0] = len_arr_small;
                hl[1] = len_arr_big;
                arr_big = &numberPart[len_arr_small];
            }
            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
            MPI_Send(hl , 5 , MPI_INT , 0 , 0 , group_comm);

            //----- Get the length to malloc for swaps -----//
            int length_to_malloc;
            MPI_Recv(&length_to_malloc , 1 , MPI_INT , 0 , 0 , group_comm , NULL);
            length_received = length_to_malloc;
            points_to_receive = (float*)malloc(sizeof(float)*dimensions*length_to_malloc);
            
            //----- Transfering points -----//
            int message[3];
            int tempHaveReceived = 0;
            int tempHaveSent = 0;
            while(1){
                //printf("dddddddddddddddccccccchhhhhhhhiiiiiiiiiiillllllllllldddddddddd \n");
                MPI_Recv(message , 3, MPI_INT , 0 , 1 , group_comm , NULL);
                if(message[0] == -1){
                    break;
                }else if(message[0] == 0){
                  //  printf("started1 sending rank %d points %d to %d\n",group_rank , message[2] , message[1]);
                    MPI_Send(&pointArr[tempHaveSent] , message[2]*dimensions , MPI_FLOAT , message[1] , 2 , group_comm);
                  //  printf("done1 sending rank %d points %d to %d\n",group_rank , message[2] , message[1]);
                    tempHaveSent += message[2]*dimensions;
                }else if(message[0] == 1){
                 //   printf("started1 receiving rank %d points %d from %d\n",group_rank , message[2] , message[1]);
                    MPI_Recv(&points_to_receive[tempHaveReceived] , message[2]*dimensions , MPI_FLOAT , message[1] , 2 , group_comm , NULL );
                   // printf("started1 receiving rank %d points %d from %d\n",group_rank , message[2] , message[1]);
                    tempHaveReceived += message[2]*dimensions;
                }
            }
            //-- second part of tranfer --//
            tempHaveReceived = 0;
            tempHaveSent = 0;
            while(1){
                MPI_Recv(message , 3, MPI_INT , 0 , 1 , group_comm , NULL);
                if(message[0] == -1){
                    break;
                }else if(message[0] == 0){
                   // printf("started22 sending rank %d points %d to %d\n",group_rank , message[2]*dimensions , message[1]);
                    MPI_Send(&pointArr[len_arr_small*dimensions +  tempHaveSent] , message[2]*dimensions , MPI_FLOAT , message[1] , 2 , group_comm);
                   // printf("done22 sending rank %d points %d to %d\n",group_rank , message[2]*dimensions , message[1]);
                    tempHaveSent += message[2]*dimensions;
                }else{
                 //   printf("started22 getting rank %d points %d from %d\n",group_rank , message[2] , message[1]);
                    MPI_Recv(&points_to_receive[tempHaveReceived] , message[2]*dimensions , MPI_FLOAT , message[1] , 2 , group_comm , NULL );
                   // printf("done22 getting rank %d points %d from %d\n",group_rank , message[2] , message[1]);
                    tempHaveReceived += message[2]*dimensions;
                }
            }

        }


       // printPoints(pointArr , partLength , 0);
       /* for(int i = 0; i<length_received; i++){
            printf("grouprank:%d \t",group_rank);
            for(int j = 0; j < dimensions; j++){
                printf("x%d %f\t",j,points_to_receive[i*dimensions+j]);
            }
            printf("\n");
        }*/
        //,medians[i]
        group_point_size = group_point_size/2;
        //------------- Merge the recieved array with the points ---------------------//
        //----- The first half of processes keeps the points nearer than median ------//
        //----- The second keeps the larger ------------------------------------------//
        if(group_rank < group_size/2){
            //printf("============================== length received %d %d rank %d\n",length_received,len_arr_big,group_rank);
            if(len_arr_big == length_received){
                memcpy(&pointArr[len_arr_small*dimensions],points_to_receive , sizeof(float)*length_received*dimensions);
            }else{
                float* tempsave = (float*)malloc(sizeof(float)* len_arr_small*dimensions);
                memcpy(tempsave , &pointArr[0], sizeof(float) * len_arr_small * dimensions);
                free(pointArr);
                pointArr = (float*)malloc(sizeof(float)*dimensions*(len_arr_small+length_received));
                memcpy(pointArr , tempsave , sizeof(float)*len_arr_small*dimensions);
               //pointArr = (float*)realloc(pointArr , sizeof(float)*dimensions*(len_arr_small+length_received));
               memcpy(&pointArr[len_arr_small*dimensions],points_to_receive , sizeof(float)*length_received*dimensions);
               partLength = len_arr_small + length_received;
               free(tempsave);
            }
        }else{
            //printf("================================================== %d %d  rank %d\n",len_arr_small , length_received,group_rank);
            if(len_arr_small == length_received){
                memcpy(pointArr,points_to_receive , sizeof(float)*length_received*dimensions);
            }else{
                float* tempsave = (float*)malloc(sizeof(float)* len_arr_big*dimensions);
                memcpy(tempsave , &pointArr[len_arr_small*dimensions], sizeof(float) * len_arr_big * dimensions);
                //pointArr = (float*)realloc(pointArr , sizeof(float)*dimensions*(len_arr_big+length_received));
                free(pointArr);
                pointArr = (float*)malloc(sizeof(float)*dimensions*(len_arr_big+length_received));
                memcpy(pointArr , tempsave , sizeof(float)*len_arr_big*dimensions);
                memcpy(&pointArr[len_arr_big*dimensions],points_to_receive , sizeof(float)*length_received*dimensions);
                partLength = len_arr_big + length_received;
                free(tempsave);
            }
        }
       // printPoints(pointArr,partLength,group_rank);

        //------- CHECK PLEASE if points are separated correctly ----------//
        free(numberPart);
        numberPart = (float*)malloc(partLength*sizeof(float));
        calculateDistances(pointArr , partLength , vp ,numberPart);
        printPointsToCsv(pointArr , partLength , dimensions , numberPart , processId , "a");
        if(group_rank < group_size/2){
            for(int i = 0; i < partLength; i++){
                float dist = 0;
                for(int k = 0; k < dimensions; k++){
                    dist += pow(vp[k] - pointArr[i*dimensions + k] , 2);
                }
                dist = sqrt(dist); //SQUARED_TAG
                //printf("dist%d = %f\n",group_rank,dist);
                if(dist > median){
                    //printf("YOU SUCK! rank %d\n",group_rank);
                    //for(int k = 0; k < dimensions; k++){
                        //printf("problemx%d-%d %f \t",i,group_rank,pointArr[i*dimensions + k]);
                   // }
                   // printf("\n");
                }
            }
        }else{
            for(int i = 0; i < partLength; i++){
                float dist = 0;
                for(int k = 0; k < dimensions; k++){
                    dist += pow(vp[k] - pointArr[i*dimensions + k] , 2);
                }
                dist = sqrt(dist); //SQUARED_TAG

                //printf("dist%d = %f\n",group_rank,dist);
                if(dist < median){
                   // printf("YOU SUCK! rank %d\n",group_rank);
                    //for(int k = 0; k < dimensions; k++){
                       //printf("problemx%d-%d %f \t",i,group_rank,pointArr[i*dimensions + k]);
                    //}
                   // printf("\n");
                }
            }
        }
        
           
        
        free(points_to_receive);
        free(numberPart);
        free(vp);
    

    



    }//END_ITER_TAG

    MPI_Barrier(MPI_COMM_WORLD);
    if(processId == 0){
    gettimeofday(&endVal , &tz);
    fprintf(fileTime,"vp time %ld s , %ld us\n", endVal.tv_sec -startVal.tv_sec  ,endVal.tv_usec - startVal.tv_usec );
    gettimeofday(&startVal,&tz);
    }
 
    //-------------- ROOT recieves mpitrees from even processes to merge it --------------//
    int treeLength = (noProcesses-1);
    float* finalMpiTree = (float*)malloc(sizeof(float)*(dimensions+1)*(noProcesses-1));
    int trlen = log2(noProcesses);

    if(processId == 0){

        float **allMPITrees = (float**)malloc(sizeof(float*)*noProcesses/2);
        for(int i = 0; i < noProcesses; i=i+2){
            if(i == 0){
                allMPITrees[0] = mpiTreeSaver;
            }else{

        //printf("restructuring tree -------------------- %d\n",i);
                allMPITrees[i/2] = (float*)malloc(log2(noProcesses)*(1+dimensions)*sizeof(float));
                MPI_Recv(allMPITrees[i/2] , log2(noProcesses)*(1+dimensions) , MPI_FLOAT , i , 4 , MPI_COMM_WORLD , NULL); //POSSIBLEERRORTAG
            }
        }

        //printf("restructuring tree --------------------\n");

        int tmpcnt = 0;
        float* mpiTree2 = (float*)malloc(sizeof(float)*(dimensions+1)*(noProcesses-1));
        for(int i = 0; i < trlen; i++){
            for(int j = 0; j < noProcesses/2; j = j + ceil((double)(noProcesses/2)/pow(2,i))){
                //printf("%d %d limit -> %d %d , tmp %d\n",i , j,trlen , noProcesses/2 , tmpcnt);
                memcpy(&mpiTree2[tmpcnt*(dimensions+1)] , &allMPITrees[j][i*(dimensions+1)] , (dimensions+1)*sizeof(float));
                tmpcnt++;
                //printf("end\n");
            }
        }

       /* printf("restructuring tree --------------------\n");
        for(int i = 0; i < treeLength; i++){
            printf("index %d\t",i);
            for(int j = 0; j < dimensions; j++){
                printf("%f\t",mpiTree2[i*(dimensions+1) + j]);
            }
            printf("median %f\n",mpiTree2[i*(dimensions+1) + dimensions]);
        }*/
            //----- change tree structure -----//

        //-------- CHANGE TREE STRUCTURE TO low_pos = root_pos + 1 ----------------------------//
        //--------------------------------- high_pos = root_pos + 2^(maxLevel - currentLevel) -//
        
        //printf("restructuring tree --------------------\n");
        //changeMpiTreeStructure(mpiTree2 , finalMpiTree , trlen , 0 , 0 , 0,dimensions);
        //printf("restructuring tree --------------------\n");
       /* for(int i = 0; i < treeLength; i++){
            printf("index %d\t",i);
            for(int j = 0; j < dimensions; j++){

            printf("aasdasdasdasdasdasdasdasdasdasd = %d\n",i*(dimensions+1)+j);
                printf("%f\t",finalMpiTree[i*(dimensions+1) + j]);
            }
            printf("median %f\n",finalMpiTree[i*(dimensions+1) + dimensions]);
        }
            */
        //INSTEAD
        memcpy(finalMpiTree , mpiTree2 ,sizeof(float)*(dimensions+1)*(noProcesses-1));
        PrintMpiTreeToCsv(finalMpiTree , noProcesses-1 , dimensions);


        for(int i = 0; i < noProcesses; i++){
            MPI_Send(finalMpiTree ,(dimensions+1)*(noProcesses-1) , MPI_FLOAT , i , 5 , MPI_COMM_WORLD );
        }
        free(mpiTree2);
        free(allMPITrees);





    }else{
        //------ Processes with even ID send their mpitree to the root ---------------//
        //--------- even IDs have the same tree as their previous even ---------------//
        if(processId%2 == 0)
            MPI_Send(mpiTreeSaver ,log2(noProcesses)*(1+dimensions) , MPI_FLOAT , 0 , 4 , MPI_COMM_WORLD);

        MPI_Recv(finalMpiTree , (dimensions+1)*(noProcesses-1) , MPI_FLOAT , 0 , 5 , MPI_COMM_WORLD , NULL);
    }

    //---------------------- SINGLE THREAD VPTree CREATION --------------------------------//
    mediansTree = (float*)malloc(sizeof(float)*partLength);
    STVantagePointTree(pointArr , mediansTree , partLength , dimensions);
    PrintToCSV(pointArr , mediansTree , partLength , dimensions , processId);
    /*for(int i = 0; i < partLength; i++){
        printf("ndx%d\t",i);
        for(int j = 0; j < dimensions; j++){
            printf("x%d=%f\t",j,pointArr[i*dimensions+j]);
        }
        printf("median = %f\n",mediansTree[i]);
    }*/
    //-------------------------- SINGLE THREAD KNN SEARCH ---------------------------------//
    char rankchar[5] ;
    sprintf(rankchar , "%d",processId);
    char name[20]  = "nbrs";
    strcat(name , rankchar);
    strcat(name , ".csv");
    FILE* file  = fopen(name , "w");
    int k = 3;
    knnVars*  obj;
    float* nbrsSave = (float*)malloc(sizeof(float)*(dimensions+1)*k*partLength);
    for(int i = 0; i < partLength; i++){
       obj = _initForKnnSearch(k , log2(partLength + 1) - 1 , dimensions , pointArr , mediansTree , i);       
        startSearchKnn(obj , 0 , 0);
        float* nbrs = getNeigbours(obj);
        memcpy(&nbrsSave[i*(dimensions+1)*k],nbrs,sizeof(float)*k*(dimensions+1)); 
        int* nodeChecked = (int*)malloc(sizeof(int)*(treeLength+noProcesses)); 
        int* nodeAdded = (int*)malloc(sizeof(int)*(treeLength+noProcesses));
       // printf("DEBUG4-%d----%d\n",treeLength+noProcesses,treeLength+processId);
        for(int j = 0; j < treeLength + noProcesses; j++){
            nodeChecked[j] = 0;
            if(j == processId+treeLength){
                //printf("did it\n");
                nodeChecked[j] = 2;
            }
            if(j == treeLength + 2){
               //printf("tried it\n");
                //nodeChecked[j] = 2;
                //printf("did it\n");
            }
        }

        //----- Print negihbours to file ----///
        //nodeChecked[treeLength+processId] = 1;
        int Status66 = -1;
        checkKnnGlobally(treeLength + processId , treeLength , &Status66  ,nodeAdded , nodeChecked , obj->vp , finalMpiTree , dimensions , obj);
      //  printf("STAAAAAAAAATTTTTTTTUUUUUUUUUSSSSSSSSSSSS =-%d- %d\n",Status66,processId);
        for(int j = 0; j < dimensions; j++){
            fprintf(file , "%f,",obj->vp[j]);
        }fprintf(file,"%f\n",0.0000);
        for(int l = 0; l < k; l++){
            for(int j = 0; j < dimensions; j++){
                fprintf(file,"%f,",obj->nbrs[l*(dimensions+1)+j]);
            }
            fprintf(file,"%f\n",obj->nbrs[l*(dimensions + 1)+ dimensions]);
        }
        _finishKnnSearch(obj);
        free(nodeChecked); 
        free(nodeAdded);
    }
    fclose(file);

    MPI_Barrier(MPI_COMM_WORLD);
    if(processId == 0){
    gettimeofday(&startVal , &tz);
    fprintf(fileTime ,"vp time %ld s , %ld us\n", startVal.tv_sec-endVal.tv_sec   , startVal.tv_usec -endVal.tv_usec );
    }
    fclose(fileTime);
    
    free(nbrsSave);
    free(finalMpiTree); 
    free(mediansTree);
    free(mpiTreeSaver);
    free(pointArr);

    MPI_Finalize();
    return 0;
}


void printPointsToCsv(float* points ,int length , int dimensions , float* medians , int rank , char* args){
    char rankchar[5] ;
    sprintf(rankchar , "%d",rank);
    char name[20]  = "POINTS";
    strcat(name , rankchar);
    strcat(name , ".csv");
    FILE* file  = fopen(name , args);
    for(int i = 0; i < length; i++){
        for(int j = 0; j < dimensions; j++){
            fprintf(file , "%f,",points[i*(dimensions)+j]);
        }
        fprintf(file , "%f\n",medians[i]);
    }
    fprintf(file ,"\n");
    fclose(file);

}
void changeMpiTreeStructure(float* oldTree , float* newTree , int levels , int c_lvl , int parentOld , int parentNew ,int dimensions){
    
    if(c_lvl > levels-1 ){
        return;
    }
    //printf("restructuring  %d %d clvl %d levels %d\n", parentOld,parentNew,c_lvl,levels );
    memcpy(&newTree[parentNew*(dimensions+1)] , &oldTree[parentOld*(dimensions+1)] , sizeof(float)*(dimensions+1));
    changeMpiTreeStructure(oldTree , newTree , levels , c_lvl+1 , parentOld*2 + 1 , parentNew + 1, dimensions);
    changeMpiTreeStructure(oldTree , newTree , levels , c_lvl+1 , parentOld*2 + 2, parentNew + pow(2 , levels - c_lvl - 1) , dimensions);


}

void STVantagePointTree(float* pointArr ,float* numberPart , int length , int dimensions){
    if(length == 1){
        return;
    }

    
    //----- Get Vantage Point -----//
    float* vp = (float*)malloc(sizeof(float)*(dimensions));
    memcpy(vp , &pointArr[0] ,sizeof(float)*(dimensions));
    //printf("VP = \t");
    //for(int j = 0; j < dimensions; j++){
    //    printf("%f\t",vp[j]);
   // }
    /*printf("\n");
    for(int i = 0; i < length; i++){
        for(int j = 0; j < dimensions; j++){
            printf("%f\t",pointArr[i*dimensions+j]);
        }
        printf("distance %f\n",numberPart[i]);
    }*/
    
    length = length - 1;
    float* arrToUse = &pointArr[dimensions];
    //float* numberPart = (float*)malloc(sizeof(float)*length*(dimensions+1));
    calculateDistances(arrToUse , length , vp , &numberPart[1]);
    //printf("\n");
    for(int i = 0; i < length; i++){
        for(int j = 0; j < dimensions; j++){
           // printf("%f\t",arrToUse[i*dimensions+j]);
        }
      //  printf("distance %f\n",numberPart[i]);
    }
    float median = selection(&numberPart[1] , length , arrToUse);

    //printf("MEDIAN %f\n",median);
    memcpy(&numberPart[0] , &median, sizeof(float));
    float* numPartToUse = &numberPart[1];
    calculateDistances(arrToUse, length , vp , numPartToUse);
    
    float* arr_big;
    float* arr_small;
    int len_arr_big;
    int len_arr_small;
    partition(numPartToUse , length , median , &arr_small , &arr_big , &len_arr_small , &len_arr_big , arrToUse);
    //---- swap the median to the low part ----//
    //int countMax , countMin , countEq;
    //validationST(median , length , numberPart , &countMin , &countMax , &countEq);
    //if()
    //printf("LLLLLLLLEEEEEEEENGGGGGGGGTHHHHSSS %d big %d median %f\n",len_arr_small , len_arr_big,median);
    for(int i = len_arr_small; i < length;i++){
        if(numPartToUse[i] == median){
            if(i == len_arr_small){
                break;
            }else{
                swap_values(numPartToUse , len_arr_small , i , arrToUse);
            }
        }
        
    }
    if(len_arr_small == 0){
        arr_small = numPartToUse;
    }
    //printf("len arr small %d        len arr big %d\n",len_arr_small , len_arr_big);
    len_arr_small++;
    len_arr_big--;
    for(int i = 0; i < len_arr_small; i++){
       // printf("arrsmall: %f\n",arr_small[i]);

    }

    for(int i = 0; i < len_arr_big; i++){
       // printf("arrbig: %f\n",arr_big[i]);
    }
    //printf("LLLLLLLLEEFIIXXXXXEEEDDDTHHHHSSS %d big %d\n",len_arr_small , len_arr_big);
    arr_big = &numPartToUse[len_arr_small];
    STVantagePointTree(arrToUse , arr_small , len_arr_small , dimensions);
    STVantagePointTree(&arrToUse[len_arr_small*dimensions] , arr_big , len_arr_big , dimensions);

    free(vp);

}

void PrintToCSV(float* points , float* medians , int length , int dimensions , int rank)
{
    char rankchar[5] ;
    sprintf(rankchar , "%d",rank);
    char name[20]  = "tree";
    strcat(name , rankchar);
    strcat(name , ".csv");
    FILE* file  = fopen(name , "w");
    for(int i = 0; i < length; i++){
        for(int j = 0; j < dimensions; j++){
            fprintf(file , "%f,",points[i*(dimensions)+j]);
        }
        fprintf(file , "%f\n",medians[i]);
    }
    fclose(file);
}

void PrintMpiTreeToCsv(float* tree , int length , int dimensions){
    FILE* file  = fopen("mpiTree.csv" , "w");
    for(int i = 0; i < length; i++){
        for(int j = 0; j < dimensions+1; j++){
            fprintf(file , "%f",tree[i*(dimensions+1)+j]);
            if(j == dimensions){
                fprintf(file , "\n");
            }else{
                fprintf(file , ",");
            }
        }
        
    }
    fclose(file);
}

void PrintNeigboursListToCsv(){

}

void assignMpiLeavesToProcs(int level , int maxLevelplus ,  int index , int* counter , int* procs){
    if(level == maxLevelplus){
        procs[*counter] = index;
        //printf("COUNTERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR %d\n",*counter);
        *counter = *counter + 1;
        return;
    }
    assignMpiLeavesToProcs(level + 1 , maxLevelplus , index+1 , counter , procs);
    assignMpiLeavesToProcs(level + 1 , maxLevelplus , index + pow(2 , maxLevelplus - level),counter,procs);

}

void checkKnnGlobally(int currentPos, int trlen , int* Status66, int* nodeAdded, int* nodeChecked ,float* vp , float* mpiTree , int dimensions, knnVars* obj){
    if(isLeaf(currentPos,trlen) == 0){
        checkIfChildrenChecked(currentPos, nodeChecked);
        if(nodeAdded[currentPos] == 1){
            nodeAdded[currentPos] = 1;
            float dist = 0;
            for(int i = 0; i < dimensions; i++){
                dist += pow(vp[i] - mpiTree[currentPos*(dimensions+1)+i],2);
            }
            dist = sqrt(dist);
            addNeighbour(obj,&mpiTree[currentPos*(dimensions+1)],dist);
        }
    }
    if(currentPos == 0 && nodeChecked[0] > 0){
        *Status66 = 666;
        return;
    }
    if(nodeChecked[currentPos] > 0){

        if(currentPos%2 == 0)
            checkKnnGlobally((currentPos-2)/2 , trlen , Status66 , nodeAdded , nodeChecked , vp  ,  mpiTree ,  dimensions , obj);
        else    
            checkKnnGlobally((currentPos-1)/2, trlen , Status66 , nodeAdded , nodeChecked, vp  ,  mpiTree ,  dimensions , obj);
    }else{
        if(isLeaf(currentPos , trlen)){

            // printf("STUFF %d----------%d\n",(currentPos-1)/2,currentPos);
            *Status66 = currentPos - trlen;
            nodeChecked[currentPos] = 3;
        }else{
            float dist = 0;
            float median = mpiTree[currentPos*(dimensions+1)+dimensions];
            for(int i = 0; i < dimensions; i++){
                dist += pow(vp[i] - mpiTree[currentPos*(dimensions+1)+i],2);
            }
            dist = sqrt(dist);
            if(dist <= median){
                if(nodeChecked[2*currentPos+1] == 0){
                    checkKnnGlobally(2*currentPos+1, trlen , Status66 , nodeAdded , nodeChecked, vp  ,  mpiTree ,  dimensions,obj);
                }else{
                    if(dist + obj->tau >= median){
                        checkKnnGlobally(2*currentPos+2, trlen , Status66 , nodeAdded , nodeChecked, vp  ,  mpiTree ,  dimensions,obj);     
                    }else{
                        *Status66 = 666;
                    }               
                }
            }else{
                 if(nodeChecked[2*currentPos+2] == 0){
                    checkKnnGlobally(2*currentPos+2, trlen , Status66 , nodeAdded , nodeChecked, vp  ,  mpiTree ,  dimensions,obj);
                }else{
                    if(dist - obj->tau <= median){
                        checkKnnGlobally(2*currentPos+1, trlen , Status66 , nodeAdded , nodeChecked, vp  ,  mpiTree ,  dimensions,obj);
                    }else{
                        *Status66 = 666;
                    }
                }
            }
        }
    }
}

int isLeaf(int pos, int trlen){
    if(pos < trlen){
        return 0;
    }
    return 1;
}

void checkIfChildrenChecked(int pos, int* nodeChecked){
    if(nodeChecked[2*pos+1] >0 && nodeChecked[2*pos + 2] >0){
        nodeChecked[pos] = 4;
        return;
    }
    nodeChecked[pos] = 0;
}