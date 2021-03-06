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
    //struct timezone tz;
    //-------------- Every process uses these vars --------------//
    float median;
    float *pointArr;
    float *vp;
    float *mpiTreeSaver;
    int mpiTreeCounter = 0;
    float *numberPart;
    float* points_to_receive;
    float *mediansTree;

    size=atoi(argv[1]);
    size = pow(2 , size) - 1;
    int dimensions = atoi(argv[2]);
    int input_k_for_knn = atoi(argv[3]);
    //int randomGenOrReadFile = atoi(argv[3]);
    _set_dimensions_for_indexing(dimensions);

    //printPointsToCsv(points_to_receive , -1 , -1 ,points_to_receive , processId , "w");

    //------------------------ Start MPI -------------------------//
    MPI_Status mpistatus;
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
        MPI_Recv(&partLength,1,MPI_INT,0,1,MPI_COMM_WORLD,&mpistatus);        
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
    gettimeofday(&startVal, NULL);
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

            for(int j = 1; j < group_size; j++)
            {
                hiloarr[j] = (int*)malloc(5*sizeof(int));
                MPI_Recv(hiloarr[j],5,MPI_INT,j,0,group_comm,&mpistatus);
            }  
            for(int i = 0; i < group_size; i++){
                for(int j = 0; j < 5; j++){
                    printf("%d\t",hiloarr[i][j]);
                }
                printf("\n");
            }          
            int sumSmall = 0;
            int sumBig = 0;
            for(int i = 0; i < group_size; i++){
                sumSmall += hiloarr[i][0];
                sumBig += hiloarr[i][1];

            }

            int* howManySwaps = (int*)malloc(sizeof(int)*group_size);
            for(int i = 0; i < group_size; i++){
                howManySwaps[i] = 0;
            }
            int procCounter = 0;
            for(int i = 0; i < group_size; i++){
                
                while(sumSmall != sumBig && hiloarr[i][2] > 0){
                    hiloarr[i][2]--;
                    sumSmall++;
                    sumBig--;
                    howManySwaps[i]++;

                }

                if(sumSmall == sumBig) break;
            }
            printf("SMAAAAAALLLL %d     BIIIIIIIIGGGGGG %d\n",sumSmall , sumBig);
            for(int i = 1; i < group_size; i++){
                MPI_Send(&howManySwaps[i],1,MPI_INT , i , 9 ,  group_comm);

            }

            if(howManySwaps[0] > 0){
                for(int i = len_arr_small; i < partLength;i++){
                    if(numberPart[i] == median){
                        
                        swap_values(numberPart , len_arr_small , i , pointArr);
                        howManySwaps[0]--;
                        len_arr_small++;
                        len_arr_big--;
                        
                    }
                    if(howManySwaps[0] == 0){
                        break;
                    }
                    
                }
                hiloarr[0][0] = len_arr_small;
                hiloarr[0][1] = len_arr_big;
                arr_big = &numberPart[len_arr_small];
            }

            for(int i = 1; i < group_size; i++){
                hiloarr[i][0] += howManySwaps[i];
                hiloarr[i][1]-= howManySwaps[i];
            }
            //------- Balance high an low arrs so they have the same length ---------//
            free(howManySwaps);
             for(int i = 0; i < group_size; i++){
                for(int j = 0; j < 5; j++){
                    printf("|%d|\t",hiloarr[i][j]);
                }
                printf("\n");
            }          

            //-------- Receive high an low lengths from slaves --------//
            /*
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
               printf("length to receive:%d   %d\n",i , lengths_to_receieve[i]);
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
                            MPI_Recv(&points_to_receive[tempHaveReceived] , message[2]*dimensions , MPI_FLOAT , jc , 2, group_comm , &mpistatus);
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
                            MPI_Recv(&points_to_receive[tempHaveReceived] , message[2]*dimensions , MPI_FLOAT , jc , 2 , group_comm , &mpistatus);
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
            MPI_Send(hl , 5 , MPI_INT , 0 , 0 , group_comm);
            int howManySwaps;
            MPI_Recv(&howManySwaps , 1 , MPI_INT , 0 , 9 , group_comm , &mpistatus);
            //----- Balance high and low arrays to have the same length -----//

            if(howManySwaps > 0){
                for(int i = len_arr_small; i < partLength;i++){
                    if(numberPart[i] == median){
                        
                        swap_values(numberPart , len_arr_small , i , pointArr);
                        howManySwaps--;
                        len_arr_small++;
                        len_arr_big--;
                        
                    }
                    if(howManySwaps == 0){
                        break;
                    }
                    
                }
                
            }
            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
            printf("rank %d small %d big  %d\n",processId , len_arr_small , len_arr_big);
            hl[0] = len_arr_small;
                hl[1] = len_arr_big;
                arr_big = &numberPart[len_arr_small];

            //----- Get the length to malloc for swaps -----//
            int length_to_malloc;
            MPI_Recv(&length_to_malloc , 1 , MPI_INT , 0 , 0 , group_comm , &mpistatus);
            length_received = length_to_malloc;
            points_to_receive = (float*)malloc(sizeof(float)*dimensions*length_to_malloc);
            
            //----- Transfering points -----//
            int message[3];
            int tempHaveReceived = 0;
            int tempHaveSent = 0;
            while(1){
                //printf("dddddddddddddddccccccchhhhhhhhiiiiiiiiiiillllllllllldddddddddd \n");
                MPI_Recv(message , 3, MPI_INT , 0 , 1 , group_comm , &mpistatus);
                if(message[0] == -1){
                    break;
                }else if(message[0] == 0){
                  //  printf("started1 sending rank %d points %d to %d\n",group_rank , message[2] , message[1]);
                    MPI_Send(&pointArr[tempHaveSent] , message[2]*dimensions , MPI_FLOAT , message[1] , 2 , group_comm);
                  //  printf("done1 sending rank %d points %d to %d\n",group_rank , message[2] , message[1]);
                    tempHaveSent += message[2]*dimensions;
                }else if(message[0] == 1){
                 //   printf("started1 receiving rank %d points %d from %d\n",group_rank , message[2] , message[1]);
                    MPI_Recv(&points_to_receive[tempHaveReceived] , message[2]*dimensions , MPI_FLOAT , message[1] , 2 , group_comm , &mpistatus );
                   // printf("started1 receiving rank %d points %d from %d\n",group_rank , message[2] , message[1]);
                    tempHaveReceived += message[2]*dimensions;
                }
            }
            //-- second part of tranfer --//
            tempHaveReceived = 0;
            tempHaveSent = 0;
            while(1){
                MPI_Recv(message , 3, MPI_INT , 0 , 1 , group_comm , &mpistatus);
                if(message[0] == -1){
                    break;
                }else if(message[0] == 0){
                   // printf("started22 sending rank %d points %d to %d\n",group_rank , message[2]*dimensions , message[1]);
                    MPI_Send(&pointArr[len_arr_small*dimensions +  tempHaveSent] , message[2]*dimensions , MPI_FLOAT , message[1] , 2 , group_comm);
                   // printf("done22 sending rank %d points %d to %d\n",group_rank , message[2]*dimensions , message[1]);
                    tempHaveSent += message[2]*dimensions;
                }else{
                 //   printf("started22 getting rank %d points %d from %d\n",group_rank , message[2] , message[1]);
                    MPI_Recv(&points_to_receive[tempHaveReceived] , message[2]*dimensions , MPI_FLOAT , message[1] , 2 , group_comm , &mpistatus );
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
                   printf("YOU SUCK! rank %d\n",group_rank);
                    //for(int k = 0; k < dimensions; k++){
                       //printf("problemx%d-%d %f \t",i,group_rank,pointArr[i*dimensions + k]);
                    //}
                   // printf("\n");
                }
            }

       
        }
         printf("PARTLENGTH %d %d\n",processId , partLength);
           
        
        free(points_to_receive);
        free(numberPart);
        free(vp);
    

    



    }//END_ITER_TAG

    /*MPI_Barrier(MPI_COMM_WORLD);
    if(processId == 0){
    gettimeofday(&endVal , NULL);
    fprintf(fileTime,"vp time %ld s , %ld us\n", endVal.tv_sec -startVal.tv_sec  ,endVal.tv_usec - startVal.tv_usec );
    gettimeofday(&startVal,NULL);
    }*/
    printf("ENDDDDDDDD PARTLENGTH %d %d\n",processId , partLength);
    
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
                MPI_Recv(allMPITrees[i/2] , log2(noProcesses)*(1+dimensions) , MPI_FLOAT , i , 4 , MPI_COMM_WORLD , &mpistatus); //POSSIBLEERRORTAG
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

        memcpy(finalMpiTree , mpiTree2 ,sizeof(float)*(dimensions+1)*(noProcesses-1));
        PrintMpiTreeToCsv(finalMpiTree , noProcesses-1 , dimensions);


        for(int i = 1; i < noProcesses; i++){
            MPI_Send(finalMpiTree ,(dimensions+1)*(noProcesses-1) , MPI_FLOAT , i , 5 , MPI_COMM_WORLD );
        }
        free(mpiTree2);
        free(allMPITrees);





    }else{
        //------ Processes with even ID send their mpitree to the root ---------------//
        //--------- even IDs have the same tree as their previous even ---------------//
        if(processId%2 == 0)
            MPI_Send(mpiTreeSaver ,log2(noProcesses)*(1+dimensions) , MPI_FLOAT , 0 , 4 , MPI_COMM_WORLD);

        MPI_Recv(finalMpiTree , (dimensions+1)*(noProcesses-1) , MPI_FLOAT , 0 , 5 , MPI_COMM_WORLD , &mpistatus);
    }

    //---------------------- SINGLE THREAD VPTree CREATION --------------------------------//
    mediansTree = (float*)malloc(sizeof(float)*partLength);
    printf("started making tree rank %d\n",processId);
    STVantagePointTree(pointArr , mediansTree , partLength , dimensions);
    printf("ended making tree rank %d\n",processId);

    MPI_Barrier(MPI_COMM_WORLD);
    if(processId == 0){
    gettimeofday(&endVal , NULL);
    fprintf(fileTime,"vp time %ld s , %ld us\n", endVal.tv_sec -startVal.tv_sec  ,endVal.tv_usec - startVal.tv_usec );

   
    gettimeofday(&startVal,NULL);
    }

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







    int k_for_Knn = input_k_for_knn;
    knnVars*  obj;
    float* nbrsSave = (float*)malloc(sizeof(float)*(dimensions+1)*k_for_Knn*partLength);
    for(int i = 0; i < partLength; i++){
       obj = _initForKnnSearch(k_for_Knn , log2(partLength + 1) - 1 , dimensions , pointArr , mediansTree , i);       
        startSearchKnn(obj , 0 , 0);
        float* nbrs = getNeigbours(obj);
        memcpy(&nbrsSave[i*(dimensions+1)*k_for_Knn],nbrs,sizeof(float)*k_for_Knn*(dimensions+1)); 
        _finishKnnSearch(obj);
        
    }


    int* statusArray = (int*) malloc(sizeof(int)*partLength);
    int* nodeChecked = (int*)malloc(sizeof(int)*(treeLength+noProcesses)*partLength); 
    int* nodeAdded = (int*)malloc(sizeof(int)*(treeLength+noProcesses)*partLength);
    int flagForSharedNbrs = 1;

    //--- initialize node info ---//
    for(int i = 0; i < partLength; i++){
        for(int j = 0; j < treeLength + noProcesses; j++){
            nodeChecked[i*(treeLength+noProcesses) + j] = 0;
            nodeAdded[i*(treeLength+noProcesses) + j] = 0;
            if(j == processId+treeLength){
                nodeChecked[i*(treeLength+noProcesses) + j] = 2;
            }
        }
        statusArray[i] = processId;
    }
    int debugLoopCounter = 0;
    while(flagForSharedNbrs == 1){
        //----------- CHECK the whole tree --------------//
        for(int i = 0; i < partLength; i++){
            if(statusArray[i] != 666){
                obj = _initForKnnSearch(k_for_Knn , log2(partLength + 1) - 1 , dimensions , pointArr , mediansTree , i);    
                //--- load previous neighbours
                memcpy(obj->nbrs, &nbrsSave[i*(dimensions+1)],sizeof(float)*k_for_Knn*(dimensions+1)); 
                obj->nbrsLength = k_for_Knn;
                
                //--- search the tree
                int status66 = -1;
                checkKnnGlobally(treeLength+statusArray[i], treeLength, &status66, &nodeAdded[i*(treeLength+noProcesses)], &nodeChecked[i*(treeLength+noProcesses)], obj->vp, finalMpiTree , dimensions, obj);
                statusArray[i] = status66;
                //printf("STATUS %d of point %d of proc %d\n",status66, i, processId);
                 memcpy( &nbrsSave[i*(dimensions+1)],obj->nbrs,sizeof(float)*k_for_Knn*(dimensions+1));
                _finishKnnSearch(obj);
            }
        }

        int* counterForPointsToSend = (int*) malloc(sizeof(int)*noProcesses);
        int* integralCounterForPointsToSent = (int*) malloc(sizeof(int)*noProcesses);
        int* counterToStackPointsTosend = (int*)malloc(sizeof(int)*noProcesses);

        //--- init counters
        for(int i = 0; i < noProcesses; i++){
            counterForPointsToSend[i] = 0;
            integralCounterForPointsToSent[i] = 0;
            counterToStackPointsTosend[i] = 0;
        }
        int rand666flag = 0;
        for(int i = 0; i < partLength; i++){
            if(statusArray[i] != 666){
                rand666flag = 1;
                counterForPointsToSend[statusArray[i]]++;
            }
        }

        int rand666sum = 0;
        MPI_Allreduce(&rand666flag ,  &rand666sum , 1 , MPI_INT , MPI_SUM  , MPI_COMM_WORLD);
        if(rand666sum == 0) break;



        int totalPointsToSent = 0;
        for(int i = 0; i < noProcesses; i++){
            totalPointsToSent += counterForPointsToSend[i];
        
            if(i > 0){
                for(int j = 0; j < i; j++){
                    integralCounterForPointsToSent[i] += counterForPointsToSend[j];
                }
            }
        }
        //-- malloc and point based on inits
        float* pointsToSentForKnn = (float*)malloc(sizeof(float)*totalPointsToSent*(dimensions));
        float** pointerForPointsToSend = (float**)malloc(sizeof(float*)*noProcesses);
        for(int i = 0; i < noProcesses; i++){
            pointerForPointsToSend[i] = &pointsToSentForKnn[integralCounterForPointsToSent[i]*dimensions];
        }
        for(int i = 0; i < partLength; i++){
            if(statusArray[i] != 666){
            // printf("status arr %d , conuter of stat %d\n",statusArray[i],counterToStackPointsTosend[statusArray[i]]*dimensions);
                memcpy( &pointerForPointsToSend[statusArray[i]][counterToStackPointsTosend[statusArray[i]]*dimensions],&pointArr[i*dimensions],sizeof(float)*dimensions);
                counterToStackPointsTosend[statusArray[i]]++;
        }
    }

        //-- receive and send points
        int* numberOfPointsToReceive = (int*)malloc(noProcesses*sizeof(int));
        float** holderForArrPointsToRecv = (float**)malloc(sizeof(float*)*noProcesses);
        for(int i = 0; i < noProcesses; i++){
            if(processId == i){
                for(int j = 0; j < noProcesses; j++){
                    if(j != processId){
                        MPI_Send(&counterForPointsToSend[j] , 1 , MPI_INT , j , 6 , MPI_COMM_WORLD);
                        printf("rank %d sending length %d\n",processId , counterForPointsToSend[j]);
                        MPI_Send(pointerForPointsToSend[j] , counterForPointsToSend[j]*dimensions , MPI_FLOAT , j , 7  , MPI_COMM_WORLD);
                    }
                }
            }else{
                MPI_Recv(&numberOfPointsToReceive[i],1,MPI_INT , i , 6 , MPI_COMM_WORLD , &mpistatus);
                printf("rank %d to recv %d\n",processId , numberOfPointsToReceive[i]);
                holderForArrPointsToRecv[i] = (float*)malloc(numberOfPointsToReceive[i]*sizeof(float)*dimensions);
                MPI_Recv(holderForArrPointsToRecv[i],numberOfPointsToReceive[i]*dimensions , MPI_FLOAT , i , 7 , MPI_COMM_WORLD , &mpistatus);

            }
        }
        numberOfPointsToReceive[processId] = 0;
       /* MPI_Barrier(MPI_COMM_WORLD);
        printf("MEEEEEEEEEEEEEINNNNNNNNNNNNN %d %d\n",processId,numberOfPointsToReceive[processId]);
        for(int i = 0; i < noProcesses; i++){
            for(int j = 0; j < numberOfPointsToReceive[i]; j++){
                printf(">>>>>>>>>>>>>>>>>rank %d from %d [%d]\t",processId,i,numberOfPointsToReceive[i]);
                for(int d = 0; d < dimensions; d++){
                    printf("xy%d %f\t",d,holderForArrPointsToRecv[i][j*dimensions + d]);
                }
            }printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);*/
        //------------------------ Find Their Knn To send to origin -------------------------//
        
        float** holderForNbrsToSendBack = (float**)malloc(sizeof(float*)*noProcesses);
        
        for(int rnk = 0; rnk < noProcesses; rnk++){
            if(processId != rnk){
                holderForNbrsToSendBack[rnk] = (float*)malloc(sizeof(float)*(dimensions+1)*k_for_Knn*numberOfPointsToReceive[rnk]);
                //printf("MAALLLLOOOOCCCCCCCCCCC %d %d    of %d\n",rnk,(dimensions+1)*k_for_Knn*numberOfPointsToReceive[rnk],processId);
                for(int p = 0; p < numberOfPointsToReceive[rnk]; p++){
                    if(processId == rnk){
                      //  printf("ALLLLLLLLLLLLL ABOOOOOOOOOOOOOOOOORT\n");
                    }
                    knnVars*  obj;
                    obj = _initForKnnSearch(k_for_Knn , log2(partLength + 1) - 1 , dimensions , pointArr , mediansTree , 0);
                    memcpy(obj->vp , &holderForArrPointsToRecv[rnk][p*dimensions] , dimensions*sizeof(float));
                    startSearchKnn(obj , 0 , 0);
                    float* nbrs_temp = getNeigbours(obj);
                    float* testDrive = (float*)malloc(sizeof(float)*k_for_Knn*(dimensions+1));
                    memcpy(testDrive,nbrs_temp,sizeof(float)*(dimensions+1)*k_for_Knn);
                    memcpy(&holderForNbrsToSendBack[ rnk ][ p * k_for_Knn * (dimensions+1) ] , nbrs_temp , sizeof(float) * (dimensions+1) * k_for_Knn);

                /*  for(int i = 0; i < k_for_Knn; i++){
                        printf("Points%d from %d to %d\t",i,processId,rnk);
                        for(int d = 0; d < dimensions; d++){
                            printf("x%d %f\t",d,nbrs_temp[i*(dimensions+1) + d]);
                            //printf("x%d %f\t",d,testDrive[i*(dimensions+1) + d]);
                            printf("x%d %f\t",p*(dimensions+1)*k_for_Knn + i*(dimensions+1) + d,holderForNbrsToSendBack[rnk][p*(dimensions+1)*k_for_Knn + i*(dimensions+1) + d]);
                        }printf("\n");
                    }printf("------------------------------\n"); 
                    */

                    _finishKnnSearch(obj);
                }
            }
        }

       // for(int i = 0; i < noProcesses; i++){
//if(processId != i){
      //          for(int j = 0; j < numberOfPointsToReceive[i]*(dimensions+1)*k_for_Knn; j++){
       //             printf("rnk %d ---%d/%f\n",processId,j,holderForNbrsToSendBack[i][j]);
         //       }
         //   }
      //  }

        float* nbrsToRecv = (float*)malloc(sizeof(float)*totalPointsToSent*(dimensions+1)*k_for_Knn);
        float** pointerForNbrsToRecv = (float**)malloc(sizeof(float*)*noProcesses);
        for(int i = 0; i < noProcesses; i++){
            pointerForNbrsToRecv[i] = &nbrsToRecv[integralCounterForPointsToSent[i]*(dimensions+1)*k_for_Knn];
            //printf("---------------------------------------------rank %d from %d integralCounter %d \n",processId , i,integralCounterForPointsToSent[i]);
        }

        for(int i = 0; i < noProcesses; i++){
            if(processId == i){
                for(int j = 0; j < noProcesses; j++){
                    if(processId != j){
                        printf("rank %d sending %d nbrs sets to %d\n",processId,numberOfPointsToReceive[j],j);
                        MPI_Send(holderForNbrsToSendBack[j], k_for_Knn*(dimensions + 1)*numberOfPointsToReceive[j], MPI_FLOAT,j,4,MPI_COMM_WORLD );
                    }
                }
            }else{

                    printf("rank %d receiving %d nbrs sets from %d\n",processId,counterForPointsToSend[i],i);
                    MPI_Recv(pointerForNbrsToRecv[i], k_for_Knn*(dimensions + 1)*counterForPointsToSend[i],MPI_FLOAT,i,4,MPI_COMM_WORLD,&mpistatus);
            }
        }


        /*
        for(int i = 0; i < totalPointsToSent; i++){
            
            for(int k = 0; k < k_for_Knn; k++){
                printf("rank %d -> nbr%d -->\t",processId,k);
                for(int j = 0; j < dimensions; j++){
                    printf("x%d %f\t",j,nbrsToRecv[i*(dimensions+1)*k_for_Knn + k*(dimensions+1) + j]);
                }
            }   printf("\n");
        }*/




        int* counterNbrsChecked = (int*)malloc(sizeof(int) * noProcesses);
        for(int i = 0; i < noProcesses; i++){
            counterNbrsChecked[i] = 0;
        }
        for(int p = 0; p < partLength; p++){
            if(statusArray[p] != 666){
                int fromRank = statusArray[p];
                float* nbrsToUse = &pointerForNbrsToRecv[fromRank][counterNbrsChecked[fromRank]*(dimensions+1)*k_for_Knn];
                knnVars*  obj;
                obj = _initForKnnSearch(k_for_Knn , log2(partLength + 1) - 1 , dimensions , pointArr , mediansTree , p);
                memcpy(obj->nbrs, &nbrsSave[p*(dimensions+1)*k_for_Knn],sizeof(float)*k_for_Knn*(dimensions+1)); 
                obj->nbrsLength = k_for_Knn;
                for(int kn = 0; kn < k_for_Knn; kn++){
                    addNeighbour(obj, &nbrsToUse[(dimensions+1)*kn] , nbrsToUse[(dimensions+1)*kn+dimensions]);
                }   
                memcpy( &nbrsSave[p*(dimensions+1)*k_for_Knn],obj->nbrs,sizeof(float)*k_for_Knn*(dimensions+1)); 
                counterNbrsChecked[fromRank] ++ ;//= k_for_Knn*(dimensions+1);
            }
        }
        free(counterForPointsToSend);
        free(integralCounterForPointsToSent);
        free(counterToStackPointsTosend);

        free(pointsToSentForKnn);
        free(pointerForPointsToSend);

        free(numberOfPointsToReceive);
        free(holderForArrPointsToRecv);

        for(int i = 0; i < noProcesses; i++){
            if(processId != i)
                free(holderForNbrsToSendBack[i]);
            }
        free(holderForNbrsToSendBack);

        free(nbrsToRecv);
        free(pointerForNbrsToRecv);

        //flagForSharedNbrs = 0;
        //printf("reached end %d\n",processId);
        printf("DEBUG ITER NO %d , of rank %d\n",++debugLoopCounter , processId);
    }





    free(nodeChecked); 
    free(nodeAdded);

    free(statusArray);


    //----------- KNN FINISH ------------//
    MPI_Barrier(MPI_COMM_WORLD);
    if(processId == 0){
    gettimeofday(&startVal , NULL);
    fprintf(fileTime ,"vp time %ld s , %ld us\n", startVal.tv_sec-endVal.tv_sec   , startVal.tv_usec -endVal.tv_usec );
    }
    fclose(fileTime);
    
    for(int i = 0; i < partLength; i++){
        for(int j = 0; j < dimensions; j++){
            fprintf(file , "%f,",pointArr[i*dimensions+j]);
        }fprintf(file,"%f\n",0.0000);
        for(int l = 0; l < k_for_Knn; l++){
            for(int j = 0; j < dimensions; j++){
                fprintf(file,"%f,",nbrsSave[i*(dimensions+1)*k_for_Knn+l*(dimensions+1)+j]);
            }
            fprintf(file,"%f\n",nbrsSave[i*(dimensions+1)*k_for_Knn+l*(dimensions + 1)+ dimensions]);
        }
    }

    

    fclose(file);
    free(nbrsSave);
    free(finalMpiTree); 
    free(mediansTree);
    free(mpiTreeSaver);
    free(pointArr);
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name(hostname,&len);
    printf("rank %d/%d running on %s\n",processId ,noProcesses , hostname);
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
    if(length == 2){
        printf("GOTCHA2\n");
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
    int v1,v2,v3;
    validationST(median, length , numberPart , &v1,&v2,&v3);

    //printf("MEDIAN %f\n",median);
    memcpy(&numberPart[0] , &median, sizeof(float));
    float* numPartToUse = &numberPart[1];
    calculateDistances(arrToUse, length , vp , numPartToUse);
    
    float* arr_big;
    float* arr_small;
    int len_arr_big;
    int len_arr_small;
    partition(numPartToUse , length , median , &arr_small , &arr_big , &len_arr_small , &len_arr_big , arrToUse);
     if(len_arr_small == 0){
        arr_small = numPartToUse;
    }
    //---- swap the median to the low part ----//
    //int countMax , countMin , countEq;
    //validationST(median , length , numberPart , &countMin , &countMax , &countEq);
    //if()
    //printf("LLLLLLLLEEEEEEEENGGGGGGGGTHHHHSSS %d big %d median %f\n",len_arr_small , len_arr_big,median);
    for(int i = len_arr_small; i < length;i++){
        if(numPartToUse[i] == median){
            swap_values(numPartToUse , len_arr_small , i , arrToUse);
            len_arr_small++;
            len_arr_big--;
            
        }
        if(len_arr_small == len_arr_big){
            break;
        }
        
    }
   
    //printf("len arr small %d        len arr big %d\n",len_arr_small , len_arr_big);
    //len_arr_small++;
    //len_arr_big--;
    if(len_arr_small != len_arr_big){
        printf("gotchaaaaaaaaaaa3 small %d   , big %d\n",len_arr_small , len_arr_big);
    }
    for(int i = 0; i < len_arr_small; i++){
       // printf("arrsmall: %f\n",arr_small[i]);

    }

    for(int i = 0; i < len_arr_big; i++){
       // printf("arrbig: %f\n",arr_big[i]);
    }
   // printf("LLLLLLLLEEFIIXXXXXEEEDDDTHHHHSSS %d big %d\n",len_arr_small , len_arr_big);
    if(len_arr_big != len_arr_small){
        printf("ERRRRRRRRRROOOOOOOOOORRRRRRRRRR   small %d big %d\n",len_arr_small , len_arr_big);
    }
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
   // printf("DEBUG current poss %d\n",currentPos);
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