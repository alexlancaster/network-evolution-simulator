#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>
#include <errno.h>

#include "random.h"
#include "lib.h"
#include "netsim.h"

typedef int (*compfn)(const void*, const void*);

int intcmp(const void *a, const void *b)
{
    return *(int *)a - *(int *)b;
}

int compare(struct AllTFBindingSites  *elem1, struct AllTFBindingSites *elem2){
    if(elem1->leftEdgePos < elem2->leftEdgePos)
        return -1;
    else if(elem1->leftEdgePos > elem2->leftEdgePos)
         return 1;
    else
         return 0;
}

 struct Coltype {
    int colnum;
    float *kval;
};
    
struct Ttype {
  int row;
  struct Coltype *col;
  int colCount;
};  

int convertToDecimal(int *bits, int TFBS){
     int n = TFBS-1;
     int count =0;
     int record[TFBS];
     int rec=0;
     while( n>=0 ){
         if(bits[n]==0){
            count++;
         }else{
             record[rec]=(TFBS-1)-n;
             rec++;
         }
         n--;
     }
     int i;
     int decimal=0;
     for(i=0; i<rec;i++){
        decimal+= pow(2,record[i]);
     }
     return decimal;
}              

int isHindered(int bindSite, int *bits, int *startPos){
    int s = bindSite - 1;
    int count=0;
    //printf("%d\n", s);
    //printf("%d\n", startPos[bindSite]);
    int check = startPos[bindSite]-HIND_LENGTH;
    if(check<0){
        check=0;}
    //printf("%d\n", check);
    while(s>=0 && startPos[s] <= startPos[bindSite] && startPos[s] >= check){
       if(bits[s] == 1){
          return 1;
       } else {
          s--;
          count++;
       }
    }
    //printf("%d\n",count);
    return 0;
  }
  
void configure(int bindSite, int *bits, int *numStates, int *statesArray, int TFBS, int *startPos){

     if(bindSite<TFBS-1){
        bits[bindSite] = 0;
        configure(bindSite+1,bits,numStates,statesArray,TFBS,startPos);
        if(!isHindered(bindSite, bits,startPos)){
           bits[bindSite] = 1;
           configure(bindSite+1,bits,numStates,statesArray,TFBS,startPos);
        }
     } else {
        bits[TFBS-1] = 0;
        statesArray[(*numStates)] = convertToDecimal(bits, TFBS);
        (*numStates)++;
        int i;
        for(i=0; i<TFBS; i++){
          // printf("%d", bits[i]);
       
        }
        //printf("\n");
        if(!isHindered(bindSite, bits,startPos)){
           bits[TFBS-1]=1;
           convertToDecimal(bits, TFBS);
           statesArray[(*numStates)] = convertToDecimal(bits, TFBS);
           (*numStates)++;
           for(i=0; i<TFBS; i++){
             // printf("%d", bits[i]);
           }
          // printf("\n");
        }
      } 
          //system("PAUSE");
      
  }

int main(int argc, char *argv[])
{
  FILE *fpkdis;
  char fperrors_name[80];
  int i, j, k;
  Genotype indiv;
  float initProteinConc[NGENES], kdis[NUM_K_DISASSEMBLY];

  int c, directory_success;
  int hold_genotype_constant = 0;
  int curr_seed;

  verbose = 0;

  /* change to get a different genotype */
  dummyrun = 4;

  /* parse command-line options */
  while ((c = getopt (argc, argv, "hvgd:r:p:t:c:")) != -1) {
    switch (c)
      {
      case 'd':
        output_directory = optarg;
        break;
      case 'r':
        dummyrun = atoi(optarg);
        break;
      case 'p':
        current_ploidy = atoi(optarg);
        break;
      case 't':
        tdevelopment = atof(optarg);
        break;
      case 'c':
        critical_size = atof(optarg);
        break;
      case 'g':
        hold_genotype_constant = 1;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'h':
        fprintf(stderr, "%s [-d DIRECTORY] [-r DUMMYRUN] [-h] [-g] [-p PLOIDY] [-t DEVELOPMENTTIME] [-c CRITICALSIZE]\n", argv[0]);
        exit(0);
        break;
      default:
        abort();
      }
  }

/* create output directory if needed */
#ifdef __unix__
  directory_success = mkdir(output_directory, S_IRUSR|S_IWUSR|S_IXUSR);
#else 
#ifdef __WIN32__
  directory_success = mkdir(output_directory);
#endif
#endif

  if (directory_success==-1) 
    if (errno == EEXIST) {
      fprintf(stderr, "directory '%s' already exists\n", output_directory);
    } else {
      fprintf(stderr, "directory '%s' cannot be created\n", output_directory);
      exit(-1);
    }

  sprintf(fperrors_name, "%s/netsimerrors.txt", output_directory);
  fperrors = fopen(fperrors_name, "w");

  /* get the kdis.txt values */
  fpkdis = fopen("kdis.txt","r");
  for (j = 0; j < NUM_K_DISASSEMBLY; j++) {
    fscanf(fpkdis,"%f", &kdis[j]);
  }
  fclose(fpkdis); 

  /* change random number */
  for (curr_seed=0; curr_seed<dummyrun; curr_seed++) ran1(&seed);


  /**********************************************************************
   *  modifiable part starts here
   **********************************************************************/

  /* initialize protein and mRNA concentrations */
  /* Jasmin: you can use these initial concentrations for computing your kons */
  for (i=0; i<NGENES; i++) {
    initProteinConc[i] = exp(1.25759*gasdev(&seed)+7.25669);
  }

  /* create sequences and binding site matrix */
  initialize_genotype(&indiv, kdis);
  
  /* print binding sites */
  /*print_all_binding_sites(indiv.ploidy, indiv.allBindingSites, indiv.bindSiteCount, 
			  indiv.transcriptionFactorSeq, indiv.cisRegSeq); 
  printf("tfsPerGene = %d", indiv.tfsPerGene);

  /* Jasmin:
     pseudo-code for looping (see also print_all_binding_sites in netsim.c):*/
     //system("PAUSE");
    int sitePos[10];
    int transFactor[10];
    int *startPos;
    startPos=malloc(indiv.tfsPerGene[0]*sizeof(int));
    int TFBS;
    TFBS = 20;
    
    int *bits = calloc(TFBS, sizeof(int));
    int *viableStates;
    //struct Ttype *arrayT;
    
    viableStates = malloc((100)*sizeof(int));
     //arrayT = malloc((pow(2,TFBS)+1)*sizeof(struct Ttype));
     int array = 0;
    
    
    qsort((void *) &(indiv.allBindingSites[0]), indiv.tfsPerGene[0],                                 
           sizeof(struct AllTFBindingSites),(compfn)compare );
   printf("tfsPerGene = %d", indiv.tfsPerGene[0]);
     printf("\n");
     int lem;
     for(lem =0; lem<TFBS; lem++){
             startPos[lem] = indiv.allBindingSites[lem].leftEdgePos;
             
             printf("%d\n", indiv.allBindingSites[lem].leftEdgePos);
     }
     printf("\n");
     configure(0,bits,&array,viableStates,TFBS,startPos);
     printf("%d\n", array);
     system("PAUSE");
      
  /* free dynamically allocated all binding sites list */
  free(indiv.allBindingSites);
   /*int d;
     for (d=0; d<array; d++) {
       free(arrayT[d].col);
  }   */ 
  //free(arrayT);
  free(viableStates);
  free(bits);
  free(startPos);
  
  /* close error file */
  fclose(fperrors);
}
