/* -*- Mode: C; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster, Jasmin Uribe, Kun Xiong
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <omp.h>
#include "numerical.h"
#include "lib.h"
#include "netsim.h"
#include "RngStream.h"
#include "mutation.h"

#define INITIALIZATION -1
#define DO_NOTHING -2
#define CUT_OFF_MISMATCH_ON_EFFECTOR 2
#define CUT_OFF_MISMATCH_SIGNAL_TO_TF 2
#define CUT_OFF_MISMATCH_TF_TO_TF 2
#define CUT_OFF_MISMATCH_WHOLE_NETWORK 2

int MAXELEMENTS=100; 

#if OLD_MODEL
const float TRANSLATION_TIME=2.0; 
const float TRANSCRIPTION_TIME=2.0;
const float PROB_ACTIVATING=0.62;
const float TRANSCRIPTINIT=6.75; 
const float MAX_REP_TO_INT_RATE=0.92;
const float BASAL_REP_TO_INT_RATE=0.15;
const float MAX_INT_TO_REP_RATE=4.11;
const float BASAL_INT_TO_REP_RATE=0.67;
const float MAX_INT_TO_ACT_RATE=3.3; 
const float BASAL_INT_TO_ACT_RATE=0.025;
const float MEAN_PROTEIN_DECAY_RATE=-4.58;
const float SD_PROTEIN_DECAY_RATE=0.67;
const float MEAN_ACT_TO_INT_RATE=2.92;
const float SD_ACT_TO_INT_RATE=0.52;
const float MEAN_MRNA_DECAY_RATE=-3.43;
const float SD_MRNA_DECAY_RATE=0.63;
const float MEAN_TRANSLATION_INIT_RATE=0.94;
const float SD_TRANSLATION_INIT_RATE=0.013;
#else
const float TRANSLATION_TIME=2.0; 
const float TRANSCRIPTION_TIME=2.0;
const float PROB_ACTIVATING=0.62;
const float TRANSCRIPTINIT=6.75; 
const float TRANSLATION_INITIATION_TIME=0.5; //min
const float TRANSCRIPTION_TERMINATION_TIME=1.0; //min
const float TRANSCRIPTION_ELONGATION_RATE=1800.0; //nt/min
const float TRANSLATION_ELONGATION_RATE=330.0; //codon/s
const float MAX_REP_TO_INT_RATE=0.92;
const float BASAL_REP_TO_INT_RATE=0.15;
const float MAX_INT_TO_REP_RATE=4.11;
const float BASAL_INT_TO_REP_RATE=0.67;
const float MAX_INT_TO_ACT_RATE=3.3; 
const float BASAL_INT_TO_ACT_RATE=0.025;
const float MEAN_PROTEIN_DECAY_RATE=-2.83;
const float SD_PROTEIN_DECAY_RATE=0.286;
const float MEAN_ACT_TO_INT_RATE=1.27;
const float SD_ACT_TO_INT_RATE=0.226;
const float MEAN_MRNA_DECAY_RATE=-1.49;
const float SD_MRNA_DECAY_RATE=0.267;
const float MEAN_TRANSLATION_INIT_RATE=0.408;
const float SD_TRANSLATION_INIT_RATE=0.384;
#endif
const float MEAN_PROTEIN_LENGTH=2.568; //log10(aa)
const float SD_PROTEIN_LENGTH=0.34;
const float MIN_Kd=1.0e-9;
const float MAX_Kd=1.0e-6;
const float log_MIN_Kd=-9.0;
const float log_MAX_Kd=-6.0;
const float NS_Kd=1.0e-5;
const float KD2APP_KD=1.8e10;
const float DEFAULT_UPDATE_INTERVAL=10.0; /*min*/
const float MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING=0.01;


/*Mutations*/
float SUBSTITUTION = 3.5e-10; /* susbstitution rate per site per cell division. Lynch 2008*/
/*the following mutations are enabled after burn-in*/
float DUPLICATION = 0.0;   
float SILENCING = 0.0;      

float MUTKINETIC = 3.8e-9; /* 1% subs in a gene (~1kb, including introns) will change kinetic rates and binding seq */    
float MUT_LOCUS_LENGTH_RATE=0.1e-10;
float proportion_mut_binding_seq = 1.0;
float proportion_mut_identity = 1.0;
float proportion_mut_koff=1.0;
float proportion_mut_kdis=1.0;
float proportion_mut_mRNA_decay=1.0;
float proportion_mut_protein_decay=1.0;
float proportion_mut_translation_rate=1.0;
float proportion_mut_cooperation=0.0;
float mutational_regression_rate=0.5;

#if OLD_MODEL
float sigma_ACT_TO_INT_RATE=1.79; 
float sigma_mRNA_decay=1.79; 
float sigma_protein_decay=1.79; 
float sigma_translation_init=1.79; 
float sigma_Kd=1.79;
float miu_ACT_TO_INT_RATE=5.22;
float miu_mRNA_decay=-1.13;
float miu_protein_decay=-4.58;
float miu_translation_init=-1.36;
float miu_Kd=-11.51;
#else
float sigma_ACT_TO_INT_RATE=0.773; 
float sigma_mRNA_decay=0.396; 
float sigma_protein_decay=1.109; 
float sigma_translation_init=0.695; 
float sigma_Kd=0.78;
float miu_ACT_TO_INT_RATE=1.57;
float miu_mRNA_decay=-1.189;
float miu_protein_decay=-2.83;
float miu_translation_init=0.107;
float miu_Kd=-5.0;
#endif

const float MAX_ACT_TO_INT_RATE=64.7;
const float MIN_ACT_TO_INT_RATE=0.59;
const float MAX_MRNA_DECAY=0.54;
const float MIN_MRNA_DECAY=7.5e-4;
const float MAX_PROTEIN_DECAY=0.69;
const float MIN_PROTEIN_DECAY=3.9e-6;
const float MAX_TRANSLATION_RATE=52.8;
const float MIN_TRANSLATION_RATE=2.7e-3;
const float MAX_KD=1.0e-5;
const float MIN_KD=0.0;
const float MAX_PROTEIN_LENGTH=5000.0; //aa
const float MIN_PROTEIN_LENGTH= 50.0; //aa

/*fitness*/
float Ne_saturate = 100.0;
float c_transl=2.0e-6;//2.0e-6;
float bmax=1.0; 
float duration_of_burn_in_growth_rate; /* allow cells to reach (possiblly) steady gr*/
int recalc_new_fitness; /*calculate the growth rate of the new genotype four more times to increase accuracy*/   

/*initial conditions*/
int init_TF_genes=6;
int init_N_act=3;
int init_N_rep=3;


/*Environments*/
float background_signal_strength = 0.0; /* number of TF0 protein */
float signal_off_strength=0.0; /* number of signal protein at basal level */
float signal_profile_matrix[N_THREADS][200][15];
float env1_signal_strength;
float env2_signal_strength;
float env1_t_development;
float env2_t_development;
float env1_t_signal_on;    
float env1_t_signal_off;     
float env2_t_signal_on;
float env2_t_signal_off;
char env1_effect_of_effector;
char env2_effect_of_effector;
int env1_fixed_effector_effect;    
int env2_fixed_effector_effect;
float env1_occurence;             /* one environment can occur more frequently than*/                                
float env2_occurence;             /* the other*/ 

char set_base_pair(float x) 
{
    char base;
    if (x<0.25)
        base = 'a';
    else if (x<0.5)
        base = 'c';
    else if (x<0.75)
        base = 'g';
    else 
        base = 't';  
    return base;
}

void initialize_sequence(char *Seq, 
                         int len,                         
                         int num_elements,
                         RngStream RS)
{
    float x;
    int i;
    int current_element = len/num_elements;
    int pos_n;

    for (i=0; i<len; i++) 
    {
        pos_n = (i / current_element)*current_element + i % current_element;  
        x = RngStream_RandU01(RS);     
        Seq[pos_n] = set_base_pair(x);
    }
}

/*This function initialize kinetic constants for gene expression, as well the identities of TFs*/
void initialize_genotype_fixed(Genotype *genotype, RngStream RS)
{
    int i;
    /* the first N_SIGNAL_TF genes encode the sensor TFs. The concentration of a sensor TF
     * is determined by certain environmental signal*/
    genotype->total_loci_length=0.0;
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {  
        #if RANDOM_COOPERATION_LOGIC        
            genotype->min_act_to_transc[i]=RngStream_RandInt(RS,1,2); //if one activator is sufficient to induce expression, the gene is regualted by OR gate.
        #else
            genotype->min_N_activator_to_transc[i]=1; 
            genotype->min_N_activator_to_transc[genotype->ngenes-1]=2;
        #endif

        /* tf affinity */
        genotype->Kd[i]=pow(10.0,(log_MAX_Kd-log_MIN_Kd)*RngStream_RandU01(RS)+log_MIN_Kd); 
#if OLD_MODEL
        /* mRNA decay */
        genotype->mRNA_decay_rate[i] = exp(SD_MRNA_DECAY_RATE*gasdev(RS)+MEAN_MRNA_DECAY_RATE);
        /* protein decay */
        genotype->protein_decay_rate[i] = exp(SD_PROTEIN_DECAY_RATE*gasdev(RS)+MEAN_PROTEIN_DECAY_RATE);        
        /* translation rate */
        genotype->translation_rate[i] = exp(SD_TRANSLATION_INIT_RATE*gasdev(RS)+MEAN_TRANSLATION_INIT_RATE);      
        /*PIC disassembly rate*/
        genotype->active_to_intermediate_rate[i]=exp(SD_ACT_TO_INT_RATE*gasdev(RS)+MEAN_ACT_TO_INT_RATE);             
#else
        /* mRNA decay */
        genotype->mRNA_decay_rate[i] = pow(10.0,SD_MRNA_DECAY_RATE*gasdev(RS)+MEAN_MRNA_DECAY_RATE);
        if(genotype->mRNA_decay_rate[i]>MAX_MRNA_DECAY)
            genotype->mRNA_decay_rate[i]=MAX_MRNA_DECAY;
        if(genotype->mRNA_decay_rate[i]<MIN_MRNA_DECAY)
            genotype->mRNA_decay_rate[i]=MIN_MRNA_DECAY;
//        while(genotype->mRNA_decay_rate[i]>=MAX_MRNA_DECAY || genotype->mRNA_decay_rate[i]<=MIN_MRNA_DECAY)
//            genotype->mRNA_decay_rate[i] = pow(10.0,SD_MRNA_DECAY_RATE*gasdev(RS)+MEAN_MRNA_DECAY_RATE);
        /* protein decay */
        genotype->protein_decay_rate[i] = pow(10.0,SD_PROTEIN_DECAY_RATE*gasdev(RS)+MEAN_PROTEIN_DECAY_RATE);  
        if(genotype->protein_decay_rate[i]>MAX_PROTEIN_DECAY)
            genotype->protein_decay_rate[i]=MAX_PROTEIN_DECAY;
        if(genotype->protein_decay_rate[i]<MIN_PROTEIN_DECAY)
            genotype->protein_decay_rate[i]=MIN_PROTEIN_DECAY;
//        while(genotype->protein_decay_rate[i]>=MAX_PROTEIN_DECAY || genotype->protein_decay_rate[i]<=MIN_PROTEIN_DECAY)
//            genotype->protein_decay_rate[i] = pow(10.0,SD_PROTEIN_DECAY_RATE*gasdev(RS)+MEAN_PROTEIN_DECAY_RATE);
        /* translation rate */
        genotype->translation_rate[i] = pow(10.0,SD_TRANSLATION_INIT_RATE*gasdev(RS)+MEAN_TRANSLATION_INIT_RATE);  
        if(genotype->translation_rate[i]>MAX_TRANSLATION_RATE)
            genotype->translation_rate[i]=MAX_TRANSLATION_RATE;
        if(genotype->translation_rate[i]<MIN_TRANSLATION_RATE)
            genotype->translation_rate[i]=MIN_TRANSLATION_RATE;
//        while(genotype->translation_rate[i]>=MAX_TRANSLATION_RATE || genotype->translation_rate[i]<=MIN_TRANSLATION_RATE)
//            genotype->translation_rate[i] = pow(10.0,SD_TRANSLATION_INIT_RATE*gasdev(RS)+MEAN_TRANSLATION_INIT_RATE); 
        /*PIC disassembly rate*/
        genotype->active_to_intermediate_rate[i]=pow(10.0,SD_ACT_TO_INT_RATE*gasdev(RS)+MEAN_ACT_TO_INT_RATE);  
        if(genotype->active_to_intermediate_rate[i]>MAX_ACT_TO_INT_RATE)
            genotype->active_to_intermediate_rate[i]=MAX_ACT_TO_INT_RATE;
        if(genotype->active_to_intermediate_rate[i]<MIN_ACT_TO_INT_RATE)
            genotype->active_to_intermediate_rate[i]=MIN_ACT_TO_INT_RATE;
//        while(genotype->active_to_intermediate_rate[i]>=MAX_ACT_TO_INT_RATE || genotype->active_to_intermediate_rate[i]<=MIN_ACT_TO_INT_RATE)
//            genotype->active_to_intermediate_rate[i]=pow(10.0,SD_ACT_TO_INT_RATE*gasdev(RS)+MEAN_ACT_TO_INT_RATE); 
#endif
#if MUT_LOCUS_LENGTH
        genotype->locus_length[i]=round(pow(10.0,SD_PROTEIN_LENGTH*gasdev(RS)+MEAN_PROTEIN_LENGTH));
        if(genotype->locus_length[i]>MAX_PROTEIN_LENGTH)
            genotype->locus_length[i]=MAX_PROTEIN_LENGTH;
        if(genotype->locus_length[i]<MIN_PROTEIN_LENGTH)
            genotype->locus_length[i]=MIN_PROTEIN_LENGTH;
        genotype->locus_length[i]*=3.0;
        genotype->total_loci_length+=genotype->locus_length[i];
#else
        genotype->locus_length[i]=1.0;
        genotype->total_loci_length+=genotype->locus_length[i];
#endif
    }          
    /* assign tf identity*/
    genotype->N_act=0;
    genotype->N_rep=0;    
    if(init_N_rep==-1 && init_N_act==-1) /*randomly generate activators and repressors*/
    {
        for(i=N_SIGNAL_TF;i<genotype->ntfgenes;i++)
        {   
            if (RngStream_RandU01(RS)<PROB_ACTIVATING) 
            {
                genotype->N_act++; 
                genotype->protein_identity[i] = ACTIVATOR;
            }
            else 
            {
                genotype->N_rep++;
                genotype->protein_identity[i]= REPRESSOR;
            }
        }
    }
    else
    {
        genotype->N_act=init_N_act;
        genotype->N_rep=init_N_rep;
        for(i=N_SIGNAL_TF;i<N_SIGNAL_TF+init_N_act;i++)            
            genotype->protein_identity[i]=ACTIVATOR;            
        for(i=N_SIGNAL_TF+init_N_act;i<genotype->ntfgenes;i++)
            genotype->protein_identity[i]=REPRESSOR;
    }
    /* parameterize sensor TF*/
    for(i=0;i<N_SIGNAL_TF;i++)
    {
        genotype->mRNA_decay_rate[i]=0.0; // we assume environmental signal toggles the state of sensor TF between active and inactive 
        genotype->protein_decay_rate[i]=0.0; // the concentration of sensor TF is constant.
        genotype->translation_rate[i]=0.0;
        genotype->active_to_intermediate_rate[i]=0.0; 
        genotype->protein_identity[i]=ACTIVATOR; /*make sensor TF an activator*/
        genotype->N_act++;
        genotype->Kd[i]=pow(10.0,(log_MAX_Kd-log_MIN_Kd)*RngStream_RandU01(RS)+log_MIN_Kd); 
    }
#if RANDOMIZE_SIGNAL2
    #if N_SIGNAL_TF==2
        if(RngStream_RandU01(RS)<=0.5) // we assume there is a background "on" signal, which is sensor TF 0, in the network.
            genotype->protein_identity[1]=ACTIVATOR; // Other sensor TFs can be either activators or repressors.
        else
        {
            genotype->protein_identity[1]=REPRESSOR;
            genotype->N_act--;
            genotype->N_rep++;
        }
    #endif
#endif
}

/*
 * initialize the genotype, this initializes random cis-regulatory
 * sequences for each individual, etc.  (full list below)
 */
void initialize_genotype(Genotype *genotype, RngStream RS)
{ 
    int i,k;

    genotype->ngenes=1+N_SIGNAL_TF+init_TF_genes; /*including the signal genes and 1 selection gene*/
    genotype->ntfgenes=N_SIGNAL_TF+init_TF_genes; /*including the signal genes*/
    genotype->nproteins=genotype->ngenes;  /*at initialization, each protein is encoded by one copy of gene*/   
    genotype->nTF_families=genotype->nproteins-1;
    /*at initialization, each copy of gene should have a unique cis-regulatory sequence*/
    for(i=0;i<genotype->ngenes;i++)
    {    
        genotype->which_cluster[i]=i; 
        genotype->cisreg_cluster[i][0]=i;
    }     
    /* initially, each protein has only one copy of gene*/    
    for(i=0;i<genotype->nproteins;i++)
    {
        genotype->protein_pool[i][0][0]=1;
        genotype->protein_pool[i][1][0]=i;
        genotype->which_protein[i]=i;        
    }  
    for(i=0;i<genotype->nTF_families;i++)
    {
        genotype->which_TF_family[i]=i;
        genotype->TF_family_pool[i][0][0]=1;
        genotype->TF_family_pool[i][0][0]=i;
    }
    initialize_sequence((char *)genotype->cisreg_seq, CISREG_LEN*NGENES, genotype->ngenes, RS);  // initialize cis-reg sequence
    initialize_sequence((char *)genotype->tf_seq, TF_ELEMENT_LEN*TFGENES, genotype->ntfgenes, RS);    //initialize binding sequence of TFs    
    /* We now generate the complementary sequence of BS that are on the non-template strand.
     * The complementary sequence is used to search for BS that on the non-template strand.  
     * We also assume that all the TFs can work on both strands, but cen induce expression in one direction.*/  
    for(i=0;i< genotype->ntfgenes;i++)
    {        
        for(k=0;k<TF_ELEMENT_LEN;k++)
        {
            switch (genotype->tf_seq[i][TF_ELEMENT_LEN-k-1])
            {
                case 'a': genotype->tf_seq_rc[i][k]='t'; break;
                case 't': genotype->tf_seq_rc[i][k]='a'; break;
                case 'c': genotype->tf_seq_rc[i][k]='g'; break;
                case 'g': genotype->tf_seq_rc[i][k]='c'; break;
            }
        }        
    }     
    initialize_genotype_fixed(genotype, RS);     
    calc_all_binding_sites(genotype);
}

/*
 * compute the list binding sites for specified gene and gene copy
 */
void calc_all_binding_sites_copy(Genotype *genotype, int gene_id)
{
    int i, j, k;
    int match,match_rc; // number of nucleotide that matches the binding sequence of TF, in a binding site in the coding and in the non-coding strand.    
    int N_hindered_BS=0;   
    int N_binding_sites=0;
    int start_TF;
    FILE *fperror;
    genotype->N_act_BS[gene_id]=0;
    genotype->N_rep_BS[gene_id]=0;
    genotype->max_hindered_sites[gene_id]=0;  
    //some helper pointer 
    char *tf_seq;
    char *cis_seq;
    char *tf_seq_rc; 
    cis_seq=&(genotype->cisreg_seq[gene_id][0]); 
  
    for(i=3; i < CISREG_LEN-TF_ELEMENT_LEN-3; i++) /* scan promoter */
    {  
        /*calc the number of BS within the hindrance range*/
        N_hindered_BS=0;        
        if(N_binding_sites>0)
        {
            for(j=0;j<N_binding_sites;j++)
            {
               if(genotype->all_binding_sites[gene_id][j].BS_pos> i-TF_ELEMENT_LEN-2*HIND_LENGTH)
                    N_hindered_BS++;
            }
        }  
        /* loop through TF proteins */        
#if !DIRECT_REG 
        if(genotype->which_protein[gene_id]==genotype->nproteins-1) // if the gene is an effector gene
            start_TF=N_SIGNAL_TF;// the environmental signals cannot directly regulate the selection gene
        else
            start_TF=0;
#else
        start_TF=0;
#endif        
        for (k=start_TF; k < genotype->nproteins-1; k++) 
        { 
            tf_seq=&(genotype->tf_seq[k][0]);
            tf_seq_rc=&(genotype->tf_seq_rc[k][0]);            
            /*find BS on the template strand*/
            match=0;
            for (j=i; j < i+TF_ELEMENT_LEN; j++) /*calculate the number of nucleotides that match in each [i,i+TF_ELEMENT_LEN] window. The window slides by 1 each time when scanning the promoter*/
                if (cis_seq[j] == tf_seq[j-i]) match++; 
            if (match >= NMIN)
            {  
                if (N_binding_sites + 1 >= genotype->N_allocated_elements) 
                {  
                    while(genotype->N_allocated_elements<=N_binding_sites+1)
                        genotype->N_allocated_elements+=100;
                   
                    for(j=0;j<NGENES;j++)
                    {
                        genotype->all_binding_sites[j] = realloc(genotype->all_binding_sites[j], genotype->N_allocated_elements*sizeof(AllTFBindingSites));
                        if (!genotype->all_binding_sites[j]) 
                        {
                            fperror=fopen(error_file,"a+");
                            LOG("error in calc_all_binding_sites_copy\n");
                            fclose(fperror);
                            exit(-1);                                       
                        }     
                    }                    
                }                
                genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;                      
                genotype->all_binding_sites[gene_id][N_binding_sites].Kd=KD2APP_KD*genotype->Kd[k]*pow(NS_Kd/genotype->Kd[k],(float)(TF_ELEMENT_LEN-match)/(TF_ELEMENT_LEN-NMIN+1));
                genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i ; 
                genotype->all_binding_sites[gene_id][N_binding_sites].mis_match = TF_ELEMENT_LEN-match;             
                genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                N_hindered_BS++;              
                N_binding_sites++;
                if(genotype->protein_identity[k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
            }
            else /*find BS on the non-template strand.*/
            {
                match_rc=0;
                for (j=i; j < i+TF_ELEMENT_LEN; j++)                
                    if (cis_seq[j] == tf_seq_rc[j-i]) match_rc++;
                if (match_rc >= NMIN)
                {
                    /**********************************************************************/     
                    if (N_binding_sites + 1 >= genotype->N_allocated_elements) 
                    {  
                        while(genotype->N_allocated_elements<=N_binding_sites+1)
                            genotype->N_allocated_elements+=100;

                        for(j=0;j<NGENES;j++)
                        {
                            genotype->all_binding_sites[j] = realloc(genotype->all_binding_sites[j], genotype->N_allocated_elements*sizeof(AllTFBindingSites));
                            if (!genotype->all_binding_sites[j]) 
                            {
                                fperror=fopen(error_file,"a+");
                                LOG("error in calc_all_binding_sites_copy\n");
                                fclose(fperror);
                                exit(-1);                                       
                            }     
                        }                    
                    }
                    /************************************************************************************************************/
                    genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;                                     
                    genotype->all_binding_sites[gene_id][N_binding_sites].Kd=KD2APP_KD*genotype->Kd[k]*pow(NS_Kd/genotype->Kd[k],(float)(TF_ELEMENT_LEN-match_rc)/(TF_ELEMENT_LEN-NMIN+1));
                    genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i;
                    genotype->all_binding_sites[gene_id][N_binding_sites].mis_match = TF_ELEMENT_LEN-match_rc;
                    genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                    N_hindered_BS++;                  
                    N_binding_sites++;  //two binding sites on different strands can also hinder each other                  
                    if(genotype->protein_identity[k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
                }
            } 
        }/* looping through TFs ends */
    }/*end of promoter scanning*/ 
    
    genotype->binding_sites_num[gene_id]=N_binding_sites;  
    genotype->N_rep_BS[gene_id]=N_binding_sites-(genotype->N_act_BS[gene_id]);
    /* calculate max_hindered_sites */    
    for(i=0;i<genotype->binding_sites_num[gene_id];i++)
    {
        genotype->max_hindered_sites[gene_id]=(genotype->max_hindered_sites[gene_id] > genotype->all_binding_sites[gene_id][i].N_hindered)?
                                      genotype->max_hindered_sites[gene_id] : genotype->all_binding_sites[gene_id][i].N_hindered;           
    }  
    /* calculate max_unhindered_sites */
    /* max_unhindered_sites is maximum number of TFs that can bind to a cis-reg sequence at the same time*/
    /* We use it to faciliate the calculation of Pact and Prep. See calc_TF_dist_from_all_BS for its usage.*/
    int act_BS[MAXELEMENTS][2],rep_BS[MAXELEMENTS][2];
    int N_act_BS,N_rep_BS;    
    N_act_BS=1;
    N_rep_BS=1;
    for(i=0;i<genotype->binding_sites_num[gene_id];i++) /* make lists BS by their types*/    
    {
        if(genotype->protein_identity[genotype->all_binding_sites[gene_id][i].tf_id]==ACTIVATOR)
        {
            act_BS[N_act_BS][0]=i;
            N_act_BS++;
        } 
        else
        {
            rep_BS[N_rep_BS][0]=i;
            N_rep_BS++;
        }
    }
    /*calculate the maximum number of activator binding sites that do not hinder each other*/
    /*Assuming that when site n is bound, at most x sites among site 1 to site n can be bound at the same time 
    *If site n+1 does not hinder n, then when site n+1 is bound, at most x+1 sites among site 1 to n+1 can be bound
    *If site n+1 hinders n, we check if it hinders site n-1,n-2.. until one site n-m which is not hindered. 
    *Obviously, site n-m+1 must be hindered by site n, therefore at most x-1 or x site among site 1 to n-m can be bound,
    *depending on whether n-m is hindered by n. This means at most x or x+1 sites among site 1 to n+1 can be bound. 
    *This means as n increases the maximum number of sites that can be bound at the same time won't decrease. 
    *We will see the maximum number of binding sites that won't hinder each other when n=N_act_BS.*/
    act_BS[0][0]=-1;
    act_BS[0][1]=0; 
    for(i=1;i<N_act_BS;i++) 
    {
        j=i-1;
        while(j!=0 && genotype->all_binding_sites[gene_id][act_BS[i][0]].BS_pos - genotype->all_binding_sites[gene_id][act_BS[j][0]].BS_pos<TF_ELEMENT_LEN+2*HIND_LENGTH)j--;
        act_BS[i][1]=act_BS[j][1]+1;
    }
    /*calculate the maximum number of repressor binding sites that do not hinder each other*/
    rep_BS[0][0]=-1;
    rep_BS[0][1]=0;
    for(i=1;i<N_rep_BS;i++) 
    {
        j=i-1;
        while(j!=0 && genotype->all_binding_sites[gene_id][rep_BS[i][0]].BS_pos - genotype->all_binding_sites[gene_id][rep_BS[j][0]].BS_pos<TF_ELEMENT_LEN+2*HIND_LENGTH)j--;
        rep_BS[i][1]=rep_BS[j][1]+1;
    }
    genotype->max_unhindered_sites[gene_id][1]=act_BS[N_act_BS-1][1];
    genotype->max_unhindered_sites[gene_id][2]=rep_BS[N_rep_BS-1][1];
}

/*
 * compute the list of binding sites for the specified number of gene
 * copies
 */
void calc_all_binding_sites(Genotype *genotype)
{    
    int gene_id;
    if(genotype->N_allocated_elements<MAXELEMENTS)
    {
        for(gene_id=0;gene_id<NGENES;gene_id++)
            genotype->all_binding_sites[gene_id]=realloc(genotype->all_binding_sites[gene_id], MAXELEMENTS*sizeof(AllTFBindingSites));
        genotype->N_allocated_elements=MAXELEMENTS;
    }
    for(gene_id=N_SIGNAL_TF;gene_id < genotype->ngenes;gene_id++)
    {        
        if(genotype->recalc_TFBS[gene_id]) /* do not calculate the binding sites if there's no mutation in the promoter or in TF binding seq*/
        {            
            calc_all_binding_sites_copy(genotype,gene_id);
            genotype->recalc_TFBS[gene_id]=NO;
        }
    }
}
/*calc_all_binding_sites2 is a copy of calc_all_binding_site. 
 *use this function and calc_all_binding_sites_copy2 to plot networks 
 * at different mismatch cutoff.
 */
void calc_all_binding_sites2(Genotype *genotype)
{    
    int gene_id;
    if(genotype->N_allocated_elements<MAXELEMENTS)
    {
        for(gene_id=0;gene_id<NGENES;gene_id++)
            genotype->all_binding_sites[gene_id]=realloc(genotype->all_binding_sites[gene_id], MAXELEMENTS*sizeof(AllTFBindingSites));
        genotype->N_allocated_elements=MAXELEMENTS;
    }
    for(gene_id=N_SIGNAL_TF;gene_id < genotype->ngenes;gene_id++)
    {        
        if(genotype->recalc_TFBS[gene_id]) /* do not calculate the binding sites if there's no mutation in the promoter or in TF binding seq*/
        {            
            calc_all_binding_sites_copy2(genotype,gene_id);
            genotype->recalc_TFBS[gene_id]=NO;
        }
    }
}

void calc_all_binding_sites_copy2(Genotype *genotype, int gene_id)
{
    int i, j, k;
    int match,match_rc; // number of nucleotide that matches the binding sequence of TF, in a binding site in the coding and in the non-coding strand.    
    int N_hindered_BS=0;   
    int N_binding_sites=0;
    int start_TF;
    FILE *fperror;
    genotype->N_act_BS[gene_id]=0;
    genotype->N_rep_BS[gene_id]=0;
    genotype->max_hindered_sites[gene_id]=0;  
    //some helper pointer 
    char *tf_seq;
    char *cis_seq;
    char *tf_seq_rc; 
    cis_seq=&(genotype->cisreg_seq[gene_id][0]); 
  
    for(i=3; i < CISREG_LEN-TF_ELEMENT_LEN-3; i++) /* scan promoter */
    {  
        /*calc the number of BS within the hindrance range*/
        N_hindered_BS=0;        
        if(N_binding_sites>0)
        {
            for(j=0;j<N_binding_sites;j++)
            {
               if(genotype->all_binding_sites[gene_id][j].BS_pos> i-TF_ELEMENT_LEN-2*HIND_LENGTH)
                    N_hindered_BS++;
            }
        }  
        /* loop through TF proteins */        
#if !DIRECT_REG 
        if(genotype->which_protein[gene_id]==genotype->nproteins-1) // if the gene is an effector gene
            start_TF=N_SIGNAL_TF;// the environmental signals cannot directly regulate the selection gene
        else
            start_TF=0;
#else
        start_TF=0;
#endif        
        for (k=start_TF; k < genotype->nproteins-1; k++) 
        { 
            tf_seq=&(genotype->tf_seq[k][0]);
            tf_seq_rc=&(genotype->tf_seq_rc[k][0]);            
            /*find BS on the template strand*/
            match=0;
            for (j=i; j < i+TF_ELEMENT_LEN; j++) /*calculate the number of nucleotides that match in each [i,i+TF_ELEMENT_LEN] window. The window slides by 1 each time when scanning the promoter*/
                if (cis_seq[j] == tf_seq[j-i]) match++; 
            if (match >= TF_ELEMENT_LEN-CUT_OFF_MISMATCH_WHOLE_NETWORK)
            {  
                if (N_binding_sites + 1 >= genotype->N_allocated_elements) 
                {  
                    while(genotype->N_allocated_elements<=N_binding_sites+1)
                        genotype->N_allocated_elements+=100;
                   
                    for(j=0;j<NGENES;j++)
                    {
                        genotype->all_binding_sites[j] = realloc(genotype->all_binding_sites[j], genotype->N_allocated_elements*sizeof(AllTFBindingSites));
                        if (!genotype->all_binding_sites[j]) 
                        {
                            fperror=fopen(error_file,"a+");
                            LOG("error in calc_all_binding_sites_copy\n");
                            fclose(fperror);
                            exit(-1);                                       
                        }     
                    }                    
                }                
                genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;                      
                genotype->all_binding_sites[gene_id][N_binding_sites].Kd=KD2APP_KD*genotype->Kd[k]*pow(NS_Kd/genotype->Kd[k],(float)(TF_ELEMENT_LEN-match)/(TF_ELEMENT_LEN-NMIN+1));
                genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i ; 
                genotype->all_binding_sites[gene_id][N_binding_sites].mis_match = TF_ELEMENT_LEN-match;             
                genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                N_hindered_BS++;              
                N_binding_sites++;
                if(genotype->protein_identity[k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
            }
            else /*find BS on the non-template strand.*/
            {
                match_rc=0;
                for (j=i; j < i+TF_ELEMENT_LEN; j++)                
                    if (cis_seq[j] == tf_seq_rc[j-i]) match_rc++;
                if (match_rc >= TF_ELEMENT_LEN-CUT_OFF_MISMATCH_WHOLE_NETWORK)
                {
                    /**********************************************************************/     
                    if (N_binding_sites + 1 >= genotype->N_allocated_elements) 
                    {  
                        while(genotype->N_allocated_elements<=N_binding_sites+1)
                            genotype->N_allocated_elements+=100;

                        for(j=0;j<NGENES;j++)
                        {
                            genotype->all_binding_sites[j] = realloc(genotype->all_binding_sites[j], genotype->N_allocated_elements*sizeof(AllTFBindingSites));
                            if (!genotype->all_binding_sites[j]) 
                            {
                                fperror=fopen(error_file,"a+");
                                LOG("error in calc_all_binding_sites_copy\n");
                                fclose(fperror);
                                exit(-1);                                       
                            }     
                        }                    
                    }
                    /************************************************************************************************************/
                    genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;                                     
                    genotype->all_binding_sites[gene_id][N_binding_sites].Kd=KD2APP_KD*genotype->Kd[k]*pow(NS_Kd/genotype->Kd[k],(float)(TF_ELEMENT_LEN-match_rc)/(TF_ELEMENT_LEN-NMIN+1));
                    genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i;
                    genotype->all_binding_sites[gene_id][N_binding_sites].mis_match = TF_ELEMENT_LEN-match_rc;
                    genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                    N_hindered_BS++;                  
                    N_binding_sites++;  //two binding sites on different strands can also hinder each other                  
                    if(genotype->protein_identity[k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
                }
            } 
        }/* looping through TFs ends */
    }/*end of promoter scanning*/ 
    
    genotype->binding_sites_num[gene_id]=N_binding_sites;  
    genotype->N_rep_BS[gene_id]=N_binding_sites-(genotype->N_act_BS[gene_id]);
    /* calculate max_hindered_sites */    
    for(i=0;i<genotype->binding_sites_num[gene_id];i++)
    {
        genotype->max_hindered_sites[gene_id]=(genotype->max_hindered_sites[gene_id] > genotype->all_binding_sites[gene_id][i].N_hindered)?
                                      genotype->max_hindered_sites[gene_id] : genotype->all_binding_sites[gene_id][i].N_hindered;           
    }  
    /* calculate max_unhindered_sites */
    /* max_unhindered_sites is maximum number of TFs that can bind to a cis-reg sequence at the same time*/
    /* We use it to faciliate the calculation of Pact and Prep. See calc_TF_dist_from_all_BS for its usage.*/
    int act_BS[MAXELEMENTS][2],rep_BS[MAXELEMENTS][2];
    int N_act_BS,N_rep_BS;    
    N_act_BS=1;
    N_rep_BS=1;
    for(i=0;i<genotype->binding_sites_num[gene_id];i++) /* make lists BS by their types*/    
    {
        if(genotype->protein_identity[genotype->all_binding_sites[gene_id][i].tf_id]==ACTIVATOR)
        {
            act_BS[N_act_BS][0]=i;
            N_act_BS++;
        } 
        else if(genotype->protein_identity[genotype->all_binding_sites[gene_id][i].tf_id]==REPRESSOR)
        {
            rep_BS[N_rep_BS][0]=i;
            N_rep_BS++;
        }
    }
    /*calculate the maximum number of activator binding sites that do not hinder each other*/
    /*Assuming that when site n is bound, at most x sites among site 1 to site n can be bound at the same time 
    *If site n+1 does not hinder n, then when site n+1 is bound, at most x+1 sites among site 1 to n+1 can be bound
    *If site n+1 hinders n, we check if it hinders site n-1,n-2.. until one site n-m which is not hindered. 
    *Obviously, site n-m+1 must be hindered by site n, therefore at most x-1 or x site among site 1 to n-m can be bound,
    *depending on whether n-m is hindered by n. This means at most x or x+1 sites among site 1 to n+1 can be bound. 
    *This means as n increases the maximum number of sites that can be bound at the same time won't decrease. 
    *We will see the maximum number of binding sites that won't hinder each other when n=N_act_BS.*/
    act_BS[0][0]=-1;
    act_BS[0][1]=0; 
    for(i=1;i<N_act_BS;i++) 
    {
        j=i-1;
        while(j!=0 && genotype->all_binding_sites[gene_id][act_BS[i][0]].BS_pos - genotype->all_binding_sites[gene_id][act_BS[j][0]].BS_pos<TF_ELEMENT_LEN+2*HIND_LENGTH)j--;
        act_BS[i][1]=act_BS[j][1]+1;
    }
    /*calculate the maximum number of repressor binding sites that do not hinder each other*/
    rep_BS[0][0]=-1;
    rep_BS[0][1]=0;
    for(i=1;i<N_rep_BS;i++) 
    {
        j=i-1;
        while(j!=0 && genotype->all_binding_sites[gene_id][rep_BS[i][0]].BS_pos - genotype->all_binding_sites[gene_id][rep_BS[j][0]].BS_pos<TF_ELEMENT_LEN+2*HIND_LENGTH)j--;
        rep_BS[i][1]=rep_BS[j][1]+1;
    }
    genotype->max_unhindered_sites[gene_id][1]=act_BS[N_act_BS-1][1];
    genotype->max_unhindered_sites[gene_id][2]=rep_BS[N_rep_BS-1][1];
}

/*Calculate probability of binding configurations*/
void calc_TF_dist_from_all_BS( AllTFBindingSites *BS_info,
                                int nproteins,
                                int max_N_hindered_BS,
                                int N_BS,                                                                     
                                int activating[NPROTEINS],
                                int max_unhindered_sites[3],
                                float protein_number[NGENES],
                                int min_N_activator_to_transcr,
                                float *Pact,
                                float *Prep,
								float *PaNr,
                                float *Pno_TF)
{
    int max_N_binding_act=max_unhindered_sites[1]+1; //Binding configurations can contain at most x activators, plus 1 type of configurations that don't have activators at all. 
    int max_N_binding_rep=max_unhindered_sites[2]+1; //Binding configurations can contain at most y repressors, plus 1 type of configurations that don't have repressors at all. 
    double ratio_matrices[N_BS+1][max_N_binding_rep][max_N_binding_act]; 
    double sum;    
    double product_of_freq; 
    int pos_of_last_record;  
    int pos_of_mat_nH;
    int pos_next_record;
    int i,j,k,m,n;
    double temp;
    
    /* initializing matrices to all zeros */
    for(i=0;i<max_N_hindered_BS+1;i++)
    {
        for(j=0;j<max_N_binding_rep;j++) 
        {
            for(k=0;k<max_N_binding_act;k++)
                ratio_matrices[i][j][k]=0.0;
        }
    }   
    /* body of the forward algorithm*/    
    pos_next_record=0; //where in the ratio_matrices to put the next record
    ratio_matrices[pos_next_record][0][0]=BS_info[0].Kd;   
    /*calculate distribution based on the first BS*/
    if(activating[BS_info[0].tf_id]==1) // if a activator binds to BS 0   
        ratio_matrices[pos_next_record][0][1]=protein_number[BS_info[0].tf_id];
    else    
        ratio_matrices[pos_next_record][1][0]=protein_number[BS_info[0].tf_id]; 
    /*keep calculating distribution from the remaining BS*/
    for(m=1;m<N_BS;m++)
    {
        pos_next_record++;
        /*If binding of site m blocks other binding sites*/
        product_of_freq = protein_number[BS_info[m].tf_id]; 
        if(BS_info[m].N_hindered!=0) 
        {
            for(n=m-BS_info[m].N_hindered;n<=m-1;n++)
                product_of_freq*=BS_info[n].Kd;            
        }
        /*Check whether m is a site of activator or repressor*/
        switch(activating[BS_info[m].tf_id])
        {
            case ACTIVATOR: // a BS of activators              
               if(m-BS_info[m].N_hindered!=0)//if binding of m does not block all of the BS evaluated before
                {      
                    /*find matrix(n-H)*/
                    pos_of_mat_nH=pos_next_record-BS_info[m].N_hindered-1; 
                   /*find matrix(n)*/
                    pos_of_last_record=pos_next_record-1;
                    for(i=0;i<max_N_binding_rep;i++)
                    {                
                        /*Note that it is possible for pos_of_mat_nH=pos_next_record, 
                         *which means we will be reading and writing to the same
                         *matrix. Updating the x column of matrix(n+1) uses
                         *the x-1 column of matrix(n-H). To avoiding changing the values
                         *of matrix(n-H) before using it to update matrix(n+1), we
                         *need to update matrix(n+1) from its last column.*/
                        for(j=max_N_binding_act-1;j>0;j--) 
                            ratio_matrices[pos_next_record][i][j]=BS_info[m].Kd*ratio_matrices[pos_of_last_record][i][j]+product_of_freq*ratio_matrices[pos_of_mat_nH][i][j-1];
                    }  
                    for(i=0;i<max_N_binding_rep;i++)
                        ratio_matrices[pos_next_record][i][0]=BS_info[m].Kd*ratio_matrices[pos_of_last_record][i][0];                    
                }
                else
                {
                    /*find matrix(n)*/
                    pos_of_last_record=pos_next_record-1;  //find last record     
                    for(i=0;i<max_N_binding_rep;i++)
                    {
                        for(j=0;j<max_N_binding_act;j++) 
                            ratio_matrices[pos_next_record][i][j]=BS_info[m].Kd*ratio_matrices[pos_of_last_record][i][j];
                    }
                    ratio_matrices[pos_next_record][0][1]+=product_of_freq;                    
                }
                break;

            case REPRESSOR: // a BS of repressors            
                if(m-BS_info[m].N_hindered!=0)
                {
                    /*find matrix(n-H)*/
                    pos_of_mat_nH=pos_next_record-BS_info[m].N_hindered-1; 
                    /*find matrix(n)*/
                    pos_of_last_record=pos_next_record-1; 
                    /*Similar problem when pos_of_mat_nH=pos_next_record*/
                    for(i=max_N_binding_rep-1;i>0;i--)
                    {
                        for(j=0;j<max_N_binding_act;j++)
                            ratio_matrices[pos_next_record][i][j]=BS_info[m].Kd*ratio_matrices[pos_of_last_record][i][j]+product_of_freq*ratio_matrices[pos_of_mat_nH][i-1][j];
                            
                    }  
                    for(j=0;j<max_N_binding_act;j++)
                        ratio_matrices[pos_next_record][0][j]=BS_info[m].Kd*ratio_matrices[pos_of_last_record][0][j];     
                }
                else
                {
                    /*find matrix(n)*/
                    pos_of_last_record=pos_next_record-1;  
                    for(i=0;i<max_N_binding_rep;i++)
                    {
                        for(j=0;j<max_N_binding_act;j++)
                            ratio_matrices[pos_next_record][i][j]=BS_info[m].Kd*ratio_matrices[pos_of_last_record][i][j];                                                   
                    }
                    ratio_matrices[pos_next_record][1][0]+=product_of_freq;
                } 
                break;
        }
    }

    sum=0.0;

    for(i=0;i<max_N_binding_rep;i++)
    {
        for(j=0;j<max_N_binding_act;j++)
//        for(j=0;j<((max_N_binding_act<max_unhindered_sites[0]+1-i)?max_N_binding_act:max_unhindered_sites[0]+1-i);j++)        
            sum+=ratio_matrices[pos_next_record][i][j];
    }   

    temp=0.0;
    for(i=0;i<max_N_binding_rep;i++)
    {  
        for(j=min_N_activator_to_transcr;j<max_N_binding_act;j++)
//        for(;j<((max_N_binding_act<max_unhindered_sites[0]+1-i)?max_N_binding_act:max_unhindered_sites[0]+1-i);j++)        
            temp+=ratio_matrices[pos_next_record][i][j];        
    } 
    *Pact=(float)(temp/sum);
    
    temp=0.0;
    for(i=1;i<max_N_binding_rep;i++)
    {
        for(j=0;j<max_N_binding_act;j++)
            temp+=ratio_matrices[pos_next_record][i][j];
    }    
    *Prep=(float)(temp/sum);
   
	temp=0.0;
	for(j=min_N_activator_to_transcr;j<max_N_binding_act;j++)
		temp+=ratio_matrices[pos_next_record][0][j];
	*PaNr = (float)(temp / sum);
    
    temp=0.0;
    for(j=0;j<min_N_activator_to_transcr;j++)
        temp+=ratio_matrices[pos_next_record][0][j];
    *Pno_TF=(float)(temp/sum);
    // end of the forward algorithm   
}

/*Add fixed event to queue*/
int add_fixed_event(int i,
                    float t,
                    FixedEvent **start,
                    FixedEvent **last)
{
    FixedEvent *newtime;
    int pos;
    FILE *fp;

    newtime = malloc(sizeof*newtime);
    if (!newtime) 
    {
        fp=fopen("output.txt","a+");
        fprintf(fp,"error in add_fixed_event\n");
        fclose(fp); 
        exit(1);
    }
    newtime->event_id = i;  
    newtime->time = t;
    pos = sls_store(newtime, start, last);
    return pos;
}

/**returns 0 if new fixed event won't happen concurrently with any exisiting event*/
int check_concurrence(  float t, //new fixed event is scheduled to happen at t
                        FixedEvent *translation, 
                        FixedEvent *transcription, 
                        FixedEvent *signal_on, 
                        FixedEvent *signal_off,
                        FixedEvent *burn_in,
                        float t_to_update_probability_of_binding,
                        FixedEvent *change_signal)
{   
    while(translation!=NULL)
    {
        if(t==translation->time)            
            return 1;
        translation=translation->next;
    }
    while(transcription!=NULL)
    {
        if(t==transcription->time)        
            return 1;
        transcription=transcription->next;
    }
    while(signal_on!=NULL)
    {
        if(t==signal_on->time)        
            return 1; 
        signal_on=signal_on->next;
    }
    while(signal_off!=NULL)
    {
        if(t==signal_off->time)            
            return 1;
        signal_off=signal_off->next;
    }
    while(burn_in!=NULL)
    {
        if(t==burn_in->time)        
            return 1; 
        burn_in=burn_in->next;
    }
    while(change_signal!=NULL)
    {
        if(t==change_signal->time)        
            return 1; 
        change_signal=change_signal->next;        
    }
    if(t==t_to_update_probability_of_binding)
        return 1;        
    return 0;
}

/*delete linked table from anywhere*/
/*This function is used only to pick up mRNA that is under translation
 *initiation*/
void delete_fixed_event(int gene_x,                      
                        int mRNA_y_of_gene_x,
                        FixedEvent **head,
                        FixedEvent **tail)
{
    FixedEvent *info, *lastinfo=NULL;
    int j, done;
    FILE *fp;
  
    j = -1;
    done = 0;
    info = *head;
    while (info) 
    {
        if (info->event_id==gene_x) 
        {
            j++; //we are looking for the yth mRNA of gene x.
            if (j == mRNA_y_of_gene_x) 
            {
                if (info == *head) //if we found mRNA y in the head of the queue
                {
                    *head = info->next; //the original 2nd in queue becomes the new head
                    if (info == *tail) //if there is only one event in queue, deleting the event leaves an empty queue. 
                        *tail = NULL; // Therefore the tail of the queue is NULL (the head is set to null in the upper line),                       
                } 
                else 
                {
                    lastinfo->next = info->next;
                    if (info == *tail) //if we are going to delete the tail of the queue
                        *tail = lastinfo; //the original second to last event becomes the new tail
                }                
                done = 1; //found mRNA y!
                break;
            } 
            else 
            {
                lastinfo = info;
                info = info->next;
            }
        } 
        else 
        {
            lastinfo = info;
            info = info->next;
        }
    }
    if (done == 0) //if could not find mRNA y
    {
        fp=fopen("output.txt","a+");
        fprintf(fp,"error in delete_fixed_event\n");
        fclose(fp);        
        exit(1);
    }
    free(info);
}

void delete_fixed_event_from_head(FixedEvent **head,FixedEvent **tail)
{
    FixedEvent *info;

    info = *head;
    *head = info->next;
    if (*tail == info) 
        *tail = NULL;
    free(info);
}

/*
 * initialize the cell state with the specified initial protein
 * concentration, mean mRNA number and mRNA decay
 */
void initialize_cell(   Genotype *genotype,
                        CellState *state,  
                        int init_mRNA_number[NGENES],
                        float init_protein_number[NPROTEINS],                        
                        float tdevelopment)
{
    int i, j;
    /* start fitness at 1.0 */
    state->cumulative_fitness = 0.0;     
    state->cumulative_fitness_after_burn_in = 0.0;
    /* initialize growth rate to zero (could also be based on 120 min doubling, i.e. 0.00578) */
    state->instantaneous_fitness = 0.0;    
    /* initialize linked tables*/
    state->mRNA_transcr_time_end_head = NULL;
    state->mRNA_transcr_time_end_tail = NULL;
    state->mRNA_transl_init_time_end_head = NULL;
    state->mRNA_transl_init_time_end_tail = NULL;
    state->signal_off_head = NULL;
    state->signal_off_tail = NULL;
    state->signal_on_head = NULL;
    state->signal_on_tail = NULL;
    state->burn_in_growth_rate_head =NULL;
    state->burn_in_growth_rate_tail=NULL;
    state->sampling_point_end_head=NULL;
    state->sampling_point_end_tail=NULL;
    state->last_event_t=0.0;  
    state->change_signal_strength_head=NULL;
    state->change_signal_strength_tail=NULL;
    state->t_to_update_probability_of_binding=TIME_INFINITY;
    state->cell_activated=0;
    /*initialize gene state, mRNA number*/
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {
        state->transcriptional_state[i]=REPRESSED;       
        state->mRNA_aft_transl_delay_num[i]=init_mRNA_number[i];
        state->mRNA_under_transl_delay_num[i]=0;
        state->mRNA_under_transc_num[i]=0;
        state->last_P_A[i]=0.0;
        state->last_P_R[i]=0.0;
        state->last_P_A_no_R[i]=0.0;;
        state->last_P_NotA_no_R[i]=0.0;
        state->protein_synthesis_index[i]=(float)state->mRNA_aft_transl_delay_num[i]*genotype->translation_rate[i]/genotype->protein_decay_rate[i];
    }       
    /* initiate protein concentration*/
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
        state->gene_specific_protein_number[i] = init_protein_number[i];    
    for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
    {
        state->protein_number[i]=0.0;        
        for(j=0;j<genotype->protein_pool[i][0][0];j++)
            state->protein_number[i]+=state->gene_specific_protein_number[genotype->protein_pool[i][1][j]];
    }    
    /* deal with the sensor tf*/
    for(i=0;i<N_SIGNAL_TF;i++)
    {
        state->mRNA_aft_transl_delay_num[i]=0;        
        state->mRNA_under_transc_num[i]=0;
        state->mRNA_under_transl_delay_num[i]=0;
        state->gene_specific_protein_number[i]=0.0;
    }        
    /*mark when to start calculating average growth rate*/
    if(duration_of_burn_in_growth_rate!=0.0)
        add_fixed_event(-1,duration_of_burn_in_growth_rate,&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail)); 
    else
        add_fixed_event(-1,(float)TIME_INFINITY,&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail));                
    /*plot protein concentration and fitness vs time*/
    #if JUST_PLOTTING 
        float t;
        int N_data_points; 
        t=0.0+TIME_OFFSET;
        N_data_points=(int)tdevelopment;
        for(i=0;i<N_data_points;i++)
        {
            add_fixed_event(-1,t,&(state->sampling_point_end_head),&(state->sampling_point_end_tail)); //get a timepoint each minute
            t+=1.0;            
        } 
    #endif    
}

/*
 * Set how the environmental signal should change
 */
void set_signal(CellState *state, 
                float duration_signal_on, 
                float duration_signal_off, 
                float *signal_profile,                 
                float tdevelopment,
                float signal_on_strength)
{
    float t=0.0;     
    char flag;    
    if(signal_profile==NULL)   
    {
        flag='o'; 
        state->protein_number[N_SIGNAL_TF-1]=signal_on_strength;    //always start with signal on
        #if N_SIGNAL_TF==2
            state->protein_number[0]=background_signal_strength;
        #endif      
        while(t<tdevelopment)
        {
            if(flag=='o')
            {
                add_fixed_event(-1,t+duration_signal_on,&(state->signal_off_head),&(state->signal_off_tail));
                flag='f';
                t=t+duration_signal_on;                    
            }    
            else
            {
                add_fixed_event(-1,t+duration_signal_off,&(state->signal_on_head),&(state->signal_on_tail));
                flag='o';
                t=t+duration_signal_off;
            }
        } 
    }
    else
    {
        int time_point=1;
        state->protein_number[N_SIGNAL_TF-1]=signal_profile[0];
        t=10.0;
        while(t<tdevelopment)
        {
            add_fixed_event(time_point,t,&(state->change_signal_strength_head),&(state->change_signal_strength_tail));
            time_point++;
            t+=10.0;
        } 
    }
}

/*
 * Calculate the rates of all Gillespie events
 */
void calc_all_rates(Genotype *genotype,
                    CellState *state,
                    GillespieRates *rates,
                    float developmental_time,
                    int UPDATE_WHAT,
                    int thread_id)
{
    int i,cluster_id,gene_id;
    int concurrent;
    float t_to_update_probability_of_binding,interval_to_update_probability_of_binding;  
    float diff_PA,diff_PR,diff_PnotAnoR,diff_PAnoR,diff_max;  
    /* reset rates */
    rates->total_mRNA_decay_rate=0.0;
    rates->total_active_to_intermediate_rate=0.0;
    rates->total_repressed_to_intermediate_rate=0.0;
    rates->total_intermediate_to_repressed_rate=0.0;
    rates->total_N_gene_transcript_initiated=0;
    rates->total_intermediate_to_active_rate=0.0;
    rates->total_Gillespie_rate=0.0;    
    for(i=0;i<genotype->ngenes;i++)
    {
        rates->repressed_to_intermediate_rate[i]=0.0;
        rates->intermediate_to_repressed_rate[i]=0.0;
        rates->intermediate_to_active_rate[i]=0.0;
        rates->active_to_intermediate_rate[i]=0.0;
        rates->mRNA_decay_rate[i]=0.0;
        rates->transcript_initiation_state[i]=0;
        state->P_A[i]=0.0;
        state->P_R[i]=0.0;
        state->P_A_no_R[i]=0.0;
        state->P_NotA_no_R[i]=0.0;
    } 
    /* update mRNA decay rates*/
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
        rates->mRNA_decay_rate[i] = genotype->mRNA_decay_rate[i] * (state->mRNA_aft_transl_delay_num[i] + state->mRNA_under_transl_delay_num[i]);
        rates->total_mRNA_decay_rate += rates->mRNA_decay_rate[i];        
    }
    /*update probability of binding configurations that activates expression
     * and use it to update other rates*/
    for(i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {    
        #if NO_REGULATION //if we manually turn off the expression of non-sensor TFs
            if(genotype->protein_identity[genotype->which_protein[i]]==NON_TF) // we only calculate the binding configurations for effector genes
            {
                cluster_id=genotype->which_cluster[i];        
                if(genotype->cisreg_cluster[cluster_id][0]!=i)  /*if this gene does not have a unique cis-reg sequence*/
                {                
                    state->P_A[i]=state->P_A[genotype->cisreg_cluster[cluster_id][0]]; /* copy TF distribution from elsewhere*/
                    state->P_R[i]=state->P_R[genotype->cisreg_cluster[cluster_id][0]];
                    state->P_A_no_R[i]=state->P_A_no_R[genotype->cisreg_cluster[cluster_id][0]];
                    state->P_NotA_no_R[i]=state->P_NotA_no_R[genotype->cisreg_cluster[cluster_id][0]];
                }
                else /* otherwise, we need to calc the ratio*/
                {
                    if(genotype->N_act_BS[i]!=0 || genotype->N_rep_BS[i]!=0)                
                        calc_TF_dist_from_all_BS(   genotype->all_binding_sites[i],
                                                    genotype->nproteins,
                                                    genotype->max_hindered_sites[i],
                                                    genotype->binding_sites_num[i],                                                                
                                                    genotype->protein_identity, 
                                                    genotype->max_unhindered_sites[i],
                                                    state->protein_number,
                                                    genotype->min_N_activator_to_transc[i],
                                                    &(state->P_A[i]),
                                                    &(state->P_R[i]),
													&(state->P_A_no_R[i]),
                                                    &(state->P_NotA_no_R[i]));
                    else
                    {
                        state->P_A[i]=0.0;     
                        state->P_R[i]=0.0;
						state->P_A_no_R[i] = 0.0;
                        state->P_NotA_no_R[i]=0.0;
                    }
                }
            }
            else
            {
               state->P_A[i]=0.0; 
               state->P_R[i]=0.0;
			   state->P_A_no_R[i] = 0.0;
               state->P_NotA_no_R[i]=0.0;
            }
        #else
            cluster_id=genotype->which_cluster[i];        
            if(genotype->cisreg_cluster[cluster_id][0]!=i)  /*if this gene does not have a unique cis-reg sequence*/
            {                
                state->P_A[i]=state->P_A[genotype->cisreg_cluster[cluster_id][0]]; /* copy TF distribution from elsewhere*/
                state->P_R[i]=state->P_R[genotype->cisreg_cluster[cluster_id][0]];
                state->P_A_no_R[i]=state->P_A_no_R[genotype->cisreg_cluster[cluster_id][0]];
                state->P_NotA_no_R[i]=state->P_NotA_no_R[genotype->cisreg_cluster[cluster_id][0]];
            }
            else /* otherwise, we need to calc the ratio*/
            {
                if(genotype->N_act_BS[i]!=0 || genotype->N_rep_BS[i]!=0)                
                    calc_TF_dist_from_all_BS(   genotype->all_binding_sites[i],
                                                genotype->nproteins,
                                                genotype->max_hindered_sites[i],
                                                genotype->binding_sites_num[i],                                                                
                                                genotype->protein_identity, 
                                                genotype->max_unhindered_sites[i],
                                                state->protein_number,
                                                genotype->min_N_activator_to_transc[i],
                                                &(state->P_A[i]),
                                                &(state->P_R[i]),
												&(state->P_A_no_R[i]),
                                                &(state->P_NotA_no_R[i]));
                else
                {
                    state->P_A[i]=0.0;     
                    state->P_R[i]=0.0;
					state->P_A_no_R[i] = 0.0;
                    state->P_NotA_no_R[i]=0.0;
                }
            }
        #endif
        /* calc other rates*/
        switch (state->transcriptional_state[i])
        {
            case REPRESSED:
                rates->repressed_to_intermediate_rate[i]=state->P_A[i]*(MAX_REP_TO_INT_RATE-BASAL_REP_TO_INT_RATE)+BASAL_REP_TO_INT_RATE;             
                rates->total_repressed_to_intermediate_rate+=rates->repressed_to_intermediate_rate[i];
                rates->intermediate_to_repressed_rate[i]=0.0;
                rates->intermediate_to_active_rate[i]=0.0;
                rates->active_to_intermediate_rate[i]=0.0;
                break;                
            case INTERMEDIATE:
                rates->intermediate_to_repressed_rate[i]=state->P_R[i]*(MAX_INT_TO_REP_RATE-BASAL_INT_TO_REP_RATE)+BASAL_INT_TO_REP_RATE;
                rates->total_intermediate_to_repressed_rate+=rates->intermediate_to_repressed_rate[i];
                rates->intermediate_to_active_rate[i]=MAX_INT_TO_ACT_RATE*state->P_A_no_R[i]+BASAL_INT_TO_ACT_RATE*state->P_NotA_no_R[i];
                rates->total_intermediate_to_active_rate+=rates->intermediate_to_active_rate[i]; 
                rates->active_to_intermediate_rate[i]=0.0;
                rates->repressed_to_intermediate_rate[i]=0.0;
                break;                
            case ACTIVE: 
                rates->active_to_intermediate_rate[i]=genotype->active_to_intermediate_rate[i];
                rates->total_active_to_intermediate_rate+=rates->active_to_intermediate_rate[i]; 
                rates->transcript_initiation_state[i]= 1;
                rates->total_N_gene_transcript_initiated+=1;
                rates->intermediate_to_repressed_rate[i]=0.0; 
                rates->repressed_to_intermediate_rate[i]=0.0;
                rates->intermediate_to_active_rate[i]=0.0;
                break;
        }
    }
    rates->total_Gillespie_rate+=rates->total_intermediate_to_repressed_rate;
    rates->total_Gillespie_rate+=rates->total_intermediate_to_active_rate;
    rates->total_Gillespie_rate+=rates->total_repressed_to_intermediate_rate;
    rates->total_Gillespie_rate+=rates->total_mRNA_decay_rate;
    rates->total_Gillespie_rate+=rates->total_active_to_intermediate_rate;
    rates->total_Gillespie_rate+=(float)rates->total_N_gene_transcript_initiated*TRANSCRIPTINIT;  
    
    /*Check if Pact needs to be updated more or less often*/ 
    if(UPDATE_WHAT!=INITIALIZATION && state->cell_activated==1)
    {  
        diff_max=0.0;
        cluster_id=1;    
        while(genotype->cisreg_cluster[cluster_id][0]!=NA) //check if Pact changes too much
        {
            gene_id=genotype->cisreg_cluster[cluster_id][0];
            diff_PA=fabs(state->P_A[gene_id]-state->last_P_A[gene_id]);            
            diff_PAnoR=fabs(state->P_A_no_R[gene_id]-state->last_P_A_no_R[gene_id]); 
            diff_PnotAnoR=fabs(state->P_NotA_no_R[gene_id]-state->last_P_NotA_no_R[gene_id]);  
            diff_PR=fabs(state->P_R[gene_id]-state->last_P_R[gene_id]);
            diff_max=(diff_max>diff_PA)?diff_max:diff_PA;
            diff_max=(diff_max>diff_PR)?diff_max:diff_PR;
            diff_max=(diff_max>diff_PAnoR)?diff_max:diff_PAnoR;
            diff_max=(diff_max>diff_PnotAnoR)?diff_max:diff_PnotAnoR;
            cluster_id++;
        }
        if(diff_max<EPSILON)
            interval_to_update_probability_of_binding=DEFAULT_UPDATE_INTERVAL;
        else
            interval_to_update_probability_of_binding=MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING/diff_max*(state->t-state->last_event_t);
        if(UPDATE_WHAT!=DO_NOTHING)          
            calc_leaping_interval(genotype,state,&interval_to_update_probability_of_binding,developmental_time,UPDATE_WHAT);  
    
        /*Update the next time that Pact will be updated mandatorily*/
        t_to_update_probability_of_binding=state->t+interval_to_update_probability_of_binding;
        concurrent=check_concurrence(   t_to_update_probability_of_binding,
                                        state->mRNA_transl_init_time_end_head,
                                        state->mRNA_transcr_time_end_head,
                                        state->signal_on_head,
                                        state->signal_off_head,
                                        state->burn_in_growth_rate_head,
                                        state->t_to_update_probability_of_binding,
                                        state->change_signal_strength_head);
        while(concurrent)//if the time to update overlaps with existing events, add a tiny offset
        {
            t_to_update_probability_of_binding+=TIME_OFFSET;
            concurrent=check_concurrence(   t_to_update_probability_of_binding,
                                            state->mRNA_transl_init_time_end_head,
                                            state->mRNA_transcr_time_end_head,
                                            state->signal_on_head,
                                            state->signal_off_head,
                                            state->burn_in_growth_rate_head,
                                            state->t_to_update_probability_of_binding,
                                            state->change_signal_strength_head);        
        }
        state->t_to_update_probability_of_binding=t_to_update_probability_of_binding;
    }
    /*Keep a copy of Pact and time for comparison at next time Pact is updated*/
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
        state->last_P_A[i]=state->P_A[i];     
        state->last_P_R[i]=state->P_R[i];
        state->last_P_A_no_R[i] = state->P_A_no_R[i];
        state->last_P_NotA_no_R[i]=state->P_NotA_no_R[i];
    }
    state->last_event_t=state->t;  
}

/* 
 * check to see if a fixed event ends within dt
 *
 * returns:
 *  0 if there is no fixed event occuring before time t+dt
 *  1 if a transcription event happens first before time t+dt
 *  2 if a translation initiation event  
 *  3 if the environmental signal is turned on 
 *  4 if the signal is turned off
 *  5 if burn-in growth ends
 *  6 if time to update Pact mandatorily
 */
 
int does_fixed_event_end(CellState *state, float t) 
{
    int retval=0;
    float t1;
    float t2;
    float t3;
    float t4;
    float t5;  
    float t6;    
    float t7;
    t1 = state->mRNA_transcr_time_end_head ? state->mRNA_transcr_time_end_head->time : TIME_INFINITY;
    t2 = state->mRNA_transl_init_time_end_head ? state->mRNA_transl_init_time_end_head->time : TIME_INFINITY;
    t3 = state->signal_off_head? state->signal_off_head->time : TIME_INFINITY;
    t4 = state->signal_on_head ? state->signal_on_head->time : TIME_INFINITY;
    t5 = state->burn_in_growth_rate_head ? state->burn_in_growth_rate_head->time : TIME_INFINITY;
    t6 = state->t_to_update_probability_of_binding;
    t7= state->change_signal_strength_head ? state->change_signal_strength_head->time : TIME_INFINITY;
    if((t1 <= t2) && (t1 <= t) && (t1 <= t3) && (t1 <= t4) && (t1<=t5) &&(t1<=t6) && (t1<=t7))
    {
        retval = 1;	
    }
    else if ((t2 <= t1) && (t2 <= t) && (t2 <= t3) && (t2 <= t4) && (t2<=t5)&&(t2<=t6)&&(t2<=t7))
    { 
        retval = 2;
    }  
    else if ((t3 <= t1) && (t3 <= t) && (t3 <= t2) && (t3 <= t4) && (t3<=t5)&&(t3<=t6) &&(t3<=t7)) 
    {
        retval = 3;
    }
    else if ((t4 <= t1) && (t4 <= t) && (t4 <= t2) && (t4 <= t3) && (t4<=t5)&&(t4<=t6)&&(t4<=t7)) 
    {
        retval = 4;
    }               
    else if((t5 <= t1) && (t5 <= t) && (t5 <= t2) && (t5 <= t3) && (t5<=t4)&&(t5<=t6)&&(t5<=t7))
    {
        retval = 5;
    }             
    else if((t6 <= t1) && (t6 <= t) && (t6 <= t2) && (t6 <= t3) && (t6<=t4)&&(t6<=t5)&&(t6<=t7))
    {
        retval=6;
    }
    else if((t7 <= t1) && (t7 <= t) && (t7 <= t2) && (t7 <= t3) && (t7<=t4)&&(t7<=t5)&&(t7<=t6))
    {
        retval=7;
    }
    else
    {
        retval = 0;
    }
    return retval;
}

/*
 * end transcription: update the mRNAs ready for translation initiation
 * etc. accordingly and delete the event from the queue
 */
void fixed_event_end_transcription( float *dt,                                    
                                    CellState *state,
                                    GillespieRates *rates,
                                    Genotype *genotype,
                                    int *end_state,                        
                                    int mut_step,
                                    Mutation *mut_record,
                                    char *effect_of_effector)
{
    int gene_id;
    float concurrent;
    float endtime;
    /* recompute the delta-t based on difference between now and the time of transcription end */
    *dt = state->mRNA_transcr_time_end_head->time - state->t;   
    /* update cell growth and protein concentration during dt*/
    update_protein_number_and_fitness(genotype, state, rates, *dt, *effect_of_effector, end_state, mut_step, mut_record);
    /* get the gene which is ending transcription */
    gene_id = state->mRNA_transcr_time_end_head->event_id;    
    /* increase number of mRNAs that are initializing translation*/
    (state->mRNA_under_transl_delay_num[gene_id])++;
    /* decrease the number of mRNAs undergoing transcription */
    (state->mRNA_under_transc_num[gene_id])--;
    /* delete the fixed even which has just occurred */
    delete_fixed_event_from_head(&(state->mRNA_transcr_time_end_head), &(state->mRNA_transcr_time_end_tail));   
    /*add transcription initialization event*/  
#if MUT_LOCUS_LENGTH
    endtime=state->t+*dt+genotype->locus_length[gene_id]/3.0/TRANSLATION_ELONGATION_RATE+TRANSLATION_INITIATION_TIME;
#else
    endtime=state->t+*dt+TRANSLATION_TIME;
#endif
    concurrent=check_concurrence(   endtime,
                                    state->mRNA_transl_init_time_end_head,
                                    state->mRNA_transcr_time_end_head,
                                    state->signal_on_head,
                                    state->signal_off_head,
                                    state->burn_in_growth_rate_head,
                                    state->t_to_update_probability_of_binding,
                                    state->change_signal_strength_head);
    while(concurrent)//if the time to update overlaps with existing events, add a tiny offset
    {
        endtime+=TIME_OFFSET;
        concurrent=check_concurrence(   endtime,
                                        state->mRNA_transl_init_time_end_head,
                                        state->mRNA_transcr_time_end_head,
                                        state->signal_on_head,
                                        state->signal_off_head,
                                        state->burn_in_growth_rate_head,
                                        state->t_to_update_probability_of_binding,
                                        state->change_signal_strength_head);        
    }  
    /*add to translation initiation event*/
    add_fixed_event(gene_id, endtime, &(state->mRNA_transl_init_time_end_head), &(state->mRNA_transl_init_time_end_tail));
}

/*
 * end translation initiation
 */
int fixed_event_end_translation_init(   Genotype *genotype, 
                                        CellState *state, 
                                        GillespieRates *rates, 
                                        float *dt,                                      
                                        int *end_state,                          
                                        int mut_step,
                                        Mutation *mut_record,
                                        char *effect_of_effector)
{
    int gene_id;    
    /* calc the remaining time till translation initiation ends */
    *dt = state->mRNA_transl_init_time_end_head->time - state->t;         
    /* update cell growth and protein concentration during dt*/
    update_protein_number_and_fitness(genotype, state, rates, *dt, *effect_of_effector, end_state, mut_step, mut_record);
    /* get identity of gene that has just finished translating */
    gene_id=state->mRNA_transl_init_time_end_head->event_id; 
    /* there is one less mRNA that is initializing translation */
    (state->mRNA_under_transl_delay_num[gene_id])--;  
    /* delete the event that just happened */
    delete_fixed_event_from_head(&(state->mRNA_transl_init_time_end_head), &(state->mRNA_transl_init_time_end_tail));    
    /* there is one more mRNA that produces protein */
    (state->mRNA_aft_transl_delay_num[gene_id])++;   
    /* update protein synthesis rate*/
    state->protein_synthesis_index[gene_id]= (float)state->mRNA_aft_transl_delay_num[gene_id]*genotype->translation_rate[gene_id]/genotype->protein_decay_rate[gene_id];
    
    if(genotype->which_protein[gene_id]==genotype->nproteins-1)//if the mRNA encodes a non-sensor TF, there could be a huge change in TF concentration
        return DO_NOTHING;
    else
        return gene_id; //so update Pact and reset updating interval to MIN_INTERVAL_TO_UPDATE_PACT
}

/*
 * compute t' factor used in the integration of growth rate
 * t' is the time the effector protein increases or decreases to a given amount
 */
float calc_tprime(Genotype *genotype, CellState* state, float *number_of_selection_protein_bf_dt, float dt, float given_amount, int protein_id) 
{
    int n_copies;
    int i;          
    n_copies=genotype->protein_pool[protein_id][0][0];
    float protein_synthesis_rate[n_copies],protein_decay_rate[n_copies];
    for(i=0;i<n_copies;i++)
    {
        protein_decay_rate[i]=genotype->protein_decay_rate[genotype->protein_pool[protein_id][1][i]];
        protein_synthesis_rate[i]=state->protein_synthesis_index[genotype->protein_pool[protein_id][1][i]]*protein_decay_rate[i];
    }       
    return rtsafe(  &calc_fx_dfx,
                    n_copies,
                    given_amount,
                    number_of_selection_protein_bf_dt,
                    protein_synthesis_rate,
                    protein_decay_rate,
                    0.0,
                    dt,
                    0.01); //rtsafe is in numerical.c
}

/*
 * calculate f(x)-Pp_s and f'(x),
 * f(x) is the number of effector protein molecules at time x
 */
void calc_fx_dfx(float x, int n_copies, float given_amount, float *intial_protein_number, float *protein_synthesis_rate, float *protein_decay_rate,float *fx, float *dfx)
{
    int i;    
    *fx=0;
    *dfx=0;    
    for(i=0;i<n_copies;i++)
    {
        *fx+=(intial_protein_number[i]-protein_synthesis_rate[i]/protein_decay_rate[i])*exp(-protein_decay_rate[i]*x)+protein_synthesis_rate[i]/protein_decay_rate[i];
        *dfx+=(protein_synthesis_rate[i]-protein_decay_rate[i]*intial_protein_number[i])*exp(-protein_decay_rate[i]*x);
    }
    *fx-=given_amount;    
}

/*
 * calculate F(delta_t)/Ne_sat. F(x) is the integral of f(x) over delta_t.
 * f(x) is the number of effector protein molecules at time x
 */
float calc_integral(Genotype *genotype, CellState *state, float *initial_protein_number, float dt, float saturate_protein_number)
{
    int i,n_copies,gene_ids[NPROTEINS];
    float integral=0.0,ect_minus_one;    
   
    n_copies=genotype->protein_pool[genotype->nproteins-1][0][0];
    for(i=0;i<n_copies;i++)
        gene_ids[i]=genotype->protein_pool[genotype->nproteins-1][1][i];
        
    for(i=0;i<n_copies;i++)
    {
        ect_minus_one=exp(-genotype->protein_decay_rate[gene_ids[i]]*dt)-1.0;    
        integral+=(state->protein_synthesis_index[gene_ids[i]]*ect_minus_one/genotype->protein_decay_rate[gene_ids[i]]-
                initial_protein_number[i]*ect_minus_one/genotype->protein_decay_rate[gene_ids[i]]+ 
                state->protein_synthesis_index[gene_ids[i]]*dt);
    }
    return integral/saturate_protein_number;
}

/*
 * return the instantaneous growth rate given the current cell state and environment,
 * also return the integrated growth rate as a pointer
 */
float calc_fitness(float *integrated_fitness,
                            Genotype *genotype,
                            CellState *state,
                            float* number_of_selection_protein_bf_dt,
                            float dt,
                            char effect_of_effector,
                            int *end_state,                         
                            int mut_step,
                            Mutation *mut_record)
{
    int i;
    float instantaneous_fitness=0.0;  /* this is returned from the function */
    float total_translation_rate = 0.0;    
    float dt_prime;
    float penalty=bmax/Ne_saturate;   
    float cost_of_expression;     
    float Ne_next=state->protein_number[genotype->nproteins-1];   
    float Ne=0.0;
    for(i=0;i<genotype->protein_pool[genotype->nproteins-1][0][0];i++) 
       Ne+=number_of_selection_protein_bf_dt[i];
                      
    /* compute the total cost of translation across all genes  */
#if NO_REGULATION_COST
    for(i=N_SIGNAL_TF; i < genotype->ngenes; i++)        
    {    
        if(genotype->which_protein[i]==genotype->nproteins-1)           
            total_translation_rate += genotype->translation_rate[i]*(float)state->mRNA_aft_transl_delay_num[i]+
                                    0.5*genotype->translation_rate[i]*(float)state->mRNA_under_transl_delay_num[i];
    } 
#else
    for(i=N_SIGNAL_TF; i < genotype->ngenes; i++)        
    {     
        total_translation_rate += genotype->translation_rate[i]*(float)state->mRNA_aft_transl_delay_num[i]+
                                    0.5*genotype->translation_rate[i]*(float)state->mRNA_under_transl_delay_num[i];
    } 
#endif  
    cost_of_expression=total_translation_rate*c_transl;

    switch (effect_of_effector)
    {
        case 'b': /* effector is beneficial!*/
            if(Ne>Ne_next)//decrease in effector protein
            {
                if(Ne_next>=Ne_saturate) //too many effector throughout
                {
                    *integrated_fitness =dt*(bmax-cost_of_expression);
                }
                else if(Ne<=Ne_saturate) // not enough effector throughout
                {
                    *integrated_fitness = bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)
                                            -cost_of_expression*dt;
                }
                else // bf dt_prime, the benefit saturates
                {
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate,genotype->nproteins-1); 
                    *integrated_fitness = bmax*dt_prime+bmax*(calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)-
                                                  calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate))-
                                                    cost_of_expression*dt;                    
                }                    
            }
            else // increase in effector protein
            {
                if(Ne>=Ne_saturate) //too many effector throughout
                {
                    *integrated_fitness =dt*(bmax-cost_of_expression);
                }   
                else if(Ne_next<=Ne_saturate)// not enough effector throughout
                {
                    *integrated_fitness = bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)
                                                    -cost_of_expression*dt;
                }
                else //Aft dt_prime, the benefit saturates
                {
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate,genotype->nproteins-1); 
                    *integrated_fitness = bmax*(dt-dt_prime)+bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate)-
                                                    cost_of_expression*dt;
                }                
            } 
            /* compute instantaneous growth rate at t */
            if (Ne_next < Ne_saturate)
                instantaneous_fitness = bmax*Ne_next/Ne_saturate;
            else
                instantaneous_fitness = bmax;
            break;
    
        case 'd': /* effector is deleterious! */      
#if !NO_PENALTY
            if(Ne>Ne_next)//decrease in effector protein
            {
                if(Ne_next>=Ne_saturate) //too many effector throughout
                {
                    *integrated_fitness =0.0-dt*cost_of_expression;
                }
                else if(Ne<=Ne_saturate) // not enough effector throughout
                {
                    *integrated_fitness = bmax*dt-penalty*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0)
                                                    -cost_of_expression*dt;
                }
                else // aft dt_prime, the benefit becomes positive
                {
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate,genotype->nproteins-1); 
                    *integrated_fitness = bmax*(dt-dt_prime)-penalty*(calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0)-
                                                  calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, 1.0))-
                                                    cost_of_expression*dt;                    
                }                    
            }
            else // increase in effector protein
            {
                if(Ne>=Ne_saturate) //too many effector throughout
                {
                    *integrated_fitness =0.0-dt*cost_of_expression;
                }   
                else if(Ne_next<=Ne_saturate)// not enough effector throughout
                {
                    *integrated_fitness = bmax*dt-penalty*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0)
                                                    -cost_of_expression*dt;
                }
                else //Aft dt_prime, the benefit becomes zero
                {
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate,genotype->nproteins-1); 
                    *integrated_fitness = bmax*dt_prime-penalty*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, 1.0)-
                                                    cost_of_expression*dt;
                }                
            } 
            if(Ne_next<Ne_saturate)
                instantaneous_fitness = bmax - penalty*Ne_next;
            else
                instantaneous_fitness = 0.0;            
#else
            *integrated_fitness=(bmax-cost_of_expression)*dt;
            instantaneous_fitness=bmax;            
#endif
            break; 
    } 
    /* and instantaneous integrated rate */
    instantaneous_fitness -= cost_of_expression;
    /* return the instantaneous growth rate */
    return instantaneous_fitness;
}


/* 
 * update both the protein concentration and current cell size *
 * 
 */
void update_protein_number_and_fitness( Genotype *genotype,
                                        CellState *state,                                   
                                        GillespieRates *rates,
                                        float dt,                     
                                        char effect_of_effector,
                                        int *end_state,                                   
                                        int mut_step,
                                        Mutation *mut_record)
{
    int i,j;
    float ct, ect, one_minus_ect;
    float N_effector_molecules_bf_dt[genotype->protein_pool[genotype->nproteins-1][0][0]];
    float instantaneous_fitness = 0.0;
    float integrated_fitness = 0.0;
    FILE *fperror;
  
    /* store the numbers of the effector proteins encoded by each copy of gene before updating*/   
    for(i=0;i<genotype->protein_pool[genotype->nproteins-1][0][0];i++)
        N_effector_molecules_bf_dt[i]=state->gene_specific_protein_number[genotype->protein_pool[genotype->nproteins-1][1][i]];
    /* update protein numbers*/
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {     
        ct=genotype->protein_decay_rate[i]*dt;
        ect = exp(-ct);
        if (fabs(ct)<EPSILON) one_minus_ect=ct;
        else one_minus_ect = 1.0-ect;      
        /* get the new protein concentration for this gene */
        state->gene_specific_protein_number[i]=ect*state->gene_specific_protein_number[i]+state->protein_synthesis_index[i]*one_minus_ect;        
    }    
    /* now, use protein_pool to pool gene specific protein number*/
    for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
    {
        state->protein_number[i]=0.0;        
        for(j=0;j<genotype->protein_pool[i][0][0];j++)
            state->protein_number[i]+=state->gene_specific_protein_number[genotype->protein_pool[i][1][j]];
    }   
    /* now find out the protein numbers at end of dt interval and compute growth rate */   
    instantaneous_fitness = calc_fitness(&integrated_fitness, 
                                            genotype, 
                                            state, 
                                            N_effector_molecules_bf_dt, 
                                            dt,                                                                                                                     
                                            effect_of_effector,
                                            end_state,                                                   
                                            mut_step,
                                            mut_record);  
    if(*end_state==0)/*use 0 to indicate abnormal behavior of the program.
                   *I am expecting rounding error to raise this flag.
                   *In case this flag is raised, quit the current replicate
                   *of growth and rerun a replicate.*/
    {
        fperror=fopen(error_file,"a+");
        LOG("at mut step %d",mut_step);
        fclose(fperror);
        return;
    } 
    /* use the integrated growth rate to compute the cell size in the next timestep */
    state->cumulative_fitness += integrated_fitness;    
    /* update the instantaneous growth rate for the beginning of the next timestep */
    state->instantaneous_fitness = instantaneous_fitness;      
}

/*
 * 
 * Functions that handle each possible Gillespie event 
 *
 */
int Gillespie_event_mRNA_decay(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    int gene_id;
    float x;
    int mRNA_id;
    while(1)/*in case of numerical error*/    
    {           
        x=RngStream_RandU01(RS)*rates->total_mRNA_decay_rate;
        gene_id=N_SIGNAL_TF-1;
        /* loop through mRNA products, to choose the mRNA with the
        proportionally higher decay rate */
        while (gene_id < genotype->ngenes-1 && x > 0.0) 
        {
            gene_id++;
            x-= rates->mRNA_decay_rate[gene_id];
        }
        /*rarely, numerical error picks up a gene that has no mRNA at all*/
        if((state->mRNA_aft_transl_delay_num[gene_id]+state->mRNA_under_transl_delay_num[gene_id])>=1)
            break;
    }    
    /* assume mRNAs are equally likely to be degraded */
    x = RngStream_RandInt(RS,1,state->mRNA_aft_transl_delay_num[gene_id] + state->mRNA_under_transl_delay_num[gene_id]);    
    /* decay mRNA in cytoplasm */
    if (x <= state->mRNA_aft_transl_delay_num[gene_id])
    {
        /* remove the mRNA from the cytoplasm count */
        (state->mRNA_aft_transl_delay_num[gene_id])--;  
        /*update protein synthesis rate*/
        state->protein_synthesis_index[gene_id] = (float)state->mRNA_aft_transl_delay_num[gene_id]*genotype->translation_rate[gene_id]/genotype->protein_decay_rate[gene_id];
        if(genotype->which_protein[gene_id]==genotype->nproteins-1)
        	return DO_NOTHING;
    	else // an mRNA of transcription factor is degraded, which can cause fluctuation in transcription factor concentrations.
        	return gene_id;    
    } 
    else 
    {
        /* decay mRNA in process of translation initialization */       
        mRNA_id = RngStream_RandInt(RS,0,state->mRNA_under_transl_delay_num[gene_id]-1);
        /* delete this fixed event: this mRNA will never be translated */
        delete_fixed_event(gene_id, mRNA_id, &(state->mRNA_transl_init_time_end_head), &(state->mRNA_transl_init_time_end_tail));       
        /* remove the mRNA from the count */
        (state->mRNA_under_transl_delay_num[gene_id])--; 
        return DO_NOTHING;
    }
}

void Gillespie_event_repressed_to_intermediate(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    int gene_id;
    float x;
    while(1)
    {
        x= RngStream_RandU01(RS)*rates->total_repressed_to_intermediate_rate;
        gene_id=N_SIGNAL_TF-1;
        while(gene_id<genotype->ngenes-1 && x>0.0)
        {
            gene_id++;
            x-=rates->repressed_to_intermediate_rate[gene_id];
        }
        if(rates->repressed_to_intermediate_rate[gene_id]>0.0)
            break;
    }
    /* set state */
    state->transcriptional_state[gene_id]=INTERMEDIATE;
}

void Gillespie_event_intermediate_to_repressed(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    int gene_id; 
    float x; 
    while(1)
    {    
        x= RngStream_RandU01(RS)*rates->total_intermediate_to_repressed_rate;
        gene_id=N_SIGNAL_TF-1;
        /* choose a particular gene copy to change state */
        while(gene_id<genotype->ngenes-1 && x>0.0)
        {
            gene_id++;
            x-=rates->intermediate_to_repressed_rate[gene_id];
        }
        if(rates->intermediate_to_repressed_rate[gene_id]>0.0)
            break;
    }
    /* set state */
    state->transcriptional_state[gene_id]=REPRESSED;
}

void Gillespie_event_intermediate_to_active(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    float x; 
    int gene_id;
    while(1)
    {    
        x= RngStream_RandU01(RS)*rates->total_intermediate_to_active_rate;
        gene_id=N_SIGNAL_TF-1;
        /* choose a particular gene copy to change state */
        while(gene_id<genotype->ngenes-1 && x>0.0)
        {
            gene_id++;
            x-=rates->intermediate_to_active_rate[gene_id];
        }
        if(rates->intermediate_to_active_rate[gene_id]>0.0)
            break;
    }
    /* set state */
    state->transcriptional_state[gene_id] =ACTIVE;  
}

void Gillespie_event_active_to_intermediate(Genotype *genotype, CellState *state,GillespieRates *rates, RngStream RS)
{
    int gene_id;
    float x;
    while(1)
    {
        x=RngStream_RandU01(RS)*rates->total_active_to_intermediate_rate;   
        gene_id=N_SIGNAL_TF-1;       
        while(gene_id < genotype->ngenes-1 && x>0.0) 
        {
            gene_id++;       
            x -= rates->active_to_intermediate_rate[gene_id];
        }
        if(rates->active_to_intermediate_rate[gene_id]>0.0)
            break;
    }
    state->transcriptional_state[gene_id]=INTERMEDIATE;   
}

void Gillespie_event_transcription_init(GillespieRates *rates, CellState *state, Genotype *genotype, float dt, RngStream RS)
{
    int gene_id;  
    int x;
    float candidate_t;
    int concurrent;    
    gene_id=N_SIGNAL_TF-1;    
    x=RngStream_RandInt(RS,1,rates->total_N_gene_transcript_initiated);
    while(gene_id<genotype->ngenes-1 && x>0)
    {
        gene_id++;
        x-=rates->transcript_initiation_state[gene_id];
    }     
    /* now that transcription of gene has been initiated, 
     * we add the timepoint at which the transcription ends, 
     * which is dt+time-of-transcription from now */
#if MUT_LOCUS_LENGTH
    candidate_t=state->t+dt+genotype->locus_length[gene_id]/TRANSCRIPTION_ELONGATION_RATE+TRANSCRIPTION_TERMINATION_TIME;
#else
    candidate_t=state->t+dt+TRANSCRIPTION_TIME;
#endif
    concurrent=check_concurrence(   candidate_t,
                                    state->mRNA_transl_init_time_end_head,
                                    state->mRNA_transcr_time_end_head,
                                    state->signal_on_head,
                                    state->signal_off_head,
                                    state->burn_in_growth_rate_head,
                                    state->t_to_update_probability_of_binding,
                                    state->change_signal_strength_head);
    while(concurrent)//if the time to update overlaps with existing events, add a tiny offset
    {
        candidate_t+=TIME_OFFSET;
        concurrent=check_concurrence(   candidate_t,
                                        state->mRNA_transl_init_time_end_head,
                                        state->mRNA_transcr_time_end_head,
                                        state->signal_on_head,
                                        state->signal_off_head,
                                        state->burn_in_growth_rate_head,
                                        state->t_to_update_probability_of_binding,                                        
                                        state->change_signal_strength_head);        
    }    
    add_fixed_event(gene_id, candidate_t,&(state->mRNA_transcr_time_end_head), &(state->mRNA_transcr_time_end_tail));
    /* increase the number mRNAs being transcribed */
    (state->mRNA_under_transc_num[gene_id])++;                      
}
/*
 * END
 * Functions that handle each possible Gillespie event 
 */

/*copy genotype from the acestor to offsprings*/
void clone_genotype(Genotype *genotype_templet, Genotype *genotype_clone)
{
    int i, j;           
    /*reset which_cluster for the clone*/
    for(i=0;i<NGENES;i++)
        genotype_clone->which_cluster[i]=NA;
    /*copy which_cluster and cis-reg sequence*/
    for(i=0; i< genotype_templet->ngenes;i++)
    {
        genotype_clone->which_cluster[i]=genotype_templet->which_cluster[i];            
        memcpy(&genotype_clone->cisreg_seq[i][0],&genotype_templet->cisreg_seq[i][0],CISREG_LEN*sizeof(char));                    
        genotype_clone->recalc_TFBS[i]=YES;                
    }    
    /*reset clone's cisreg_cluster*/
    i=0;
    while(genotype_clone->cisreg_cluster[i][0]!=-1)
    {
        j=0;
        while(genotype_clone->cisreg_cluster[i][j]!=-1)
        {
            genotype_clone->cisreg_cluster[i][j]=-1;
            j++;
        }
        i++;
    }        
    /*then copy from templet*/
    i=0;
    while(genotype_templet->cisreg_cluster[i][0]!=-1)
    {
        j=0;
        while(genotype_templet->cisreg_cluster[i][j]!=-1)
        {
            genotype_clone->cisreg_cluster[i][j]=genotype_templet->cisreg_cluster[i][j];
            j++;
        }
        i++;
    }
    /*reset clone's information*/
    for(i=0;i<NGENES;i++)
    {
        genotype_clone->which_protein[i]=NA;
        genotype_clone->min_N_activator_to_transc[i]=MAX_BINDING+1;   
    }    
    /*reset clone's tf_family_pool and which_tf_family*/
    for(i=0;i<NPROTEINS;i++)
    {        
        for(j=0;j<NPROTEINS;j++)
            genotype_clone->TF_family_pool[i][1][j]=NA;
        genotype_clone->TF_family_pool[i][0][0]=0;
        genotype_clone->which_TF_family[i]=NA;
    }
    /*copy from templet's tf_family_pool*/
    for(i=0;i<genotype_templet->nTF_families;i++)
    {
        genotype_clone->TF_family_pool[i][0][0]=genotype_templet->TF_family_pool[i][0][0];
        for(j=0;j<genotype_templet->TF_family_pool[i][0][0];j++)
            genotype_clone->TF_family_pool[i][1][j]=genotype_templet->TF_family_pool[i][1][j];
    }
    /*reset clone's protein_pool*/
    for(i=0;i<NPROTEINS;i++)
    {            
        for(j=0;j<NGENES;j++)
            genotype_clone->protein_pool[i][1][j]=NA;
        genotype_clone->protein_pool[i][0][0]=0;            
    }
    /*copy from templet's protein_pool and tf_family_pool*/
    for(i=0;i<genotype_templet->nproteins;i++)
    {            
        genotype_clone->which_TF_family[i]=genotype_templet->which_TF_family[i];
        genotype_clone->protein_pool[i][0][0]=genotype_templet->protein_pool[i][0][0];            
        for(j=0;j<genotype_templet->protein_pool[i][0][0];j++)
            genotype_clone->protein_pool[i][1][j]=genotype_templet->protein_pool[i][1][j];                     
    }    
    /* copy binding sites' sequences*/  
    for(i=0; i < genotype_templet->ntfgenes; i++) 
    {          
        for(j=0;j<TF_ELEMENT_LEN;j++)
        {    
            genotype_clone->tf_seq[i][j]=genotype_templet->tf_seq[i][j];
            genotype_clone->tf_seq_rc[i][j]=genotype_templet->tf_seq_rc[i][j];
        }
    }
    /*copy kinetic constants*/
    for(i=0; i < genotype_templet->ngenes; i++) 
    {            
        genotype_clone->mRNA_decay_rate[i]=genotype_templet->mRNA_decay_rate[i];
        genotype_clone->protein_decay_rate[i]=genotype_templet->protein_decay_rate[i];
        genotype_clone->translation_rate[i]=genotype_templet->translation_rate[i];            
        genotype_clone->active_to_intermediate_rate[i]=genotype_templet->active_to_intermediate_rate[i];
        genotype_clone->which_protein[i]=genotype_templet->which_protein[i];
        genotype_clone->locus_length[i]=genotype_templet->locus_length[i];
        genotype_clone->min_N_activator_to_transc[i]=genotype_templet->min_N_activator_to_transc[i];   
    } 
    /* copy TF information*/
    for(i=0;i<NPROTEINS;i++)
    {
        genotype_clone->protein_identity[i]=genotype_templet->protein_identity[i];
        genotype_clone->Kd[i]=genotype_templet->Kd[i];
    }    
    /* copy gene and protein numbers*/
    genotype_clone->ngenes=genotype_templet->ngenes;
    genotype_clone->ntfgenes=genotype_templet->ntfgenes;
    genotype_clone->nproteins=genotype_templet->nproteins;
    genotype_clone->nTF_families=genotype_templet->nTF_families;
    genotype_clone->N_act=genotype_templet->N_act;
    genotype_clone->N_rep=genotype_templet->N_rep;
    genotype_clone->total_loci_length=genotype_templet->total_loci_length;
}

/*
 * run the model for a specified cell for a single timestep:
 */
void do_single_timestep(Genotype *genotype, 
                        CellState *state,                         
                        GillespieRates *rates,
                        char *effect_of_effector, 
                        float signal_on_strength,
                        float duration_signal_on,
                        float duration_signal_off,
                        char init_effector_effect,
                        int fixed_effector_effect,
                        RngStream RS,                       
                        int mut_step,
                        Mutation *mut_record,
                        int *end_state,                    
                        int thread_id,
                        float tdevelopment,
                        float *signal_profile) 
{    
    int event, UPDATE_WHAT;     
    float fixed_time; 
    float dt;
    float x;
    FILE *fperror;

    x = expdev(RS);        /* draw random number */
    dt = x/rates->total_Gillespie_rate;
    if (dt < 0.0) 
    {	
        fperror=fopen(error_file,"a+");
        LOG("negative dt at mut_step %d\n",mut_step);
        fclose(fperror);
        *end_state=0; /*use 0 to indicate abnormal behavior of the program.
                       *I am expecting rounding error to raise this flag.
                       *In case this flag is raised, quit the current replicate
                       *of growth and rerun a replicate.*/
        return;    
    }
    /* first check to see if a fixed event occurs in current t->dt window, or in tdevelopment if running for a fixed development time */
    fixed_time = (state->t+dt<tdevelopment)?(state->t+dt):tdevelopment;
    event = does_fixed_event_end(state, fixed_time);
    while(event!=0)
    {           
        /*after doing fixed event, return a flag to indicate whether mandatorily update Pact or Prep*/
        UPDATE_WHAT=do_fixed_event( genotype, 
                                    state, 
                                    rates, 
                                    &dt,                                     
                                    event, 
                                    signal_on_strength,
                                    duration_signal_on, 
                                    duration_signal_off, 
                                    effect_of_effector, 
                                    init_effector_effect,
                                    fixed_effector_effect, 
                                    end_state, 
                                    mut_step, 
                                    mut_record,
                                    signal_profile);        
        if(*end_state==0) 
            return; 
        state->t += dt;  /* advance time by the dt */        
#if CAUTIOUS //watch for overlapping fixed events -- one of which will be skipped and cause error
        float dt_copy,x_copy,rate_total_Gillespie_rate_copy,t_copy;  
        int event2;           
        fixed_time=*t+dt;
        event2=does_fixed_event_end(state, fixed_time);
        if(event2!=0)
        {
            fperror=fopen(error_file,"a+");
            LOG("overlapping fixed events at mut_step %d", mut_step);
            fclose(fperror);
        }                       
        t_copy=*t; 
        dt_copy=dt;
        x_copy=x; 
        rate_total_Gillespie_rate_copy=rates->total_Gillespie_rate; 
#endif        
        x -= dt*rates->total_Gillespie_rate;  /* we've been running with rates->total_Gillespie_rate for dt, so substract it from x*/        
        calc_all_rates(genotype, state, rates, tdevelopment,UPDATE_WHAT,thread_id);  /* update rates->total_Gillespie_rate and re-compute a new dt */      
        dt = x/rates->total_Gillespie_rate;
        /*deal with rounding error*/
        if(dt<0.0)
        {  	
#if CAUTIOUS // this rounding error can happen very often, therefore the error_log can be huge
            fperror=fopen(error_file,"a+");
            LOG("rounding error in dt at mut_step %d\n",mut_step);            	
            fclose(fperror);
#endif 
            dt=TIME_OFFSET; 
        }
        fixed_time = (state->t+dt<tdevelopment)?(state->t+dt):tdevelopment; 
        /* check to see there aren't more fixed events to do */
        event = does_fixed_event_end(state, fixed_time);                                    
    } 
  /* no remaining fixed events to do in dt, now do stochastic events */  
  /* if we haven't already reached end of development with last
     delta-t, if there is no fixed development time, we always execute
     this  */          
    if (state->t+dt < tdevelopment)
    {
        update_protein_number_and_fitness(genotype, state, rates, dt, *effect_of_effector, end_state, mut_step, mut_record);        
        if(*end_state==0)
            return; 
        UPDATE_WHAT=do_Gillespie_event(genotype, state, rates, dt, RS, end_state, mut_step, mut_record);
        if(*end_state==0)
            return;  
        /* Gillespie step: advance time to next event at dt */
        state->t += dt;
        calc_all_rates(genotype,state,rates,tdevelopment,UPDATE_WHAT,thread_id);        
    } 
    else 
    { 
        /* do remaining dt */
        dt = tdevelopment - state->t;
        /* final update of protein concentration */
        update_protein_number_and_fitness(genotype, state, rates, dt, *effect_of_effector, end_state, mut_step, mut_record); 
        if(*end_state==0)
            return;         
        /* advance to end of development (this exits the outer while loop) */
        state->t = tdevelopment;
    }
}

/* while there are fixed events
     occuring in current t->dt window */
int do_fixed_event(Genotype *genotype, 
                    CellState *state, 
                    GillespieRates *rates, 
                    float *dt,                           
                    int event,
                    float signal_on_strength,
                    float duration_signal_on,
                    float duration_signal_off,
                    char *effect_of_effector,
                    char init_effector_effect,
                    int fixed_effector_effect,                
                    int *end_state,                  
                    int mut_step,
                    Mutation *mut_record,
                    float *signal_profile)
{     
    int return_value;
    return_value=DO_NOTHING;
    switch (event) 
    {
        case 1:     /* a transcription event ends */
            fixed_event_end_transcription(dt, state, rates, genotype,end_state,mut_step,mut_record,effect_of_effector); 
            break;
        case 2:     /* a translation initialization event ends */ 
            return_value=fixed_event_end_translation_init(genotype, state, rates, dt, end_state,mut_step,mut_record,effect_of_effector);
            state->cell_activated=1;
            break;
        case 3:     /* turn signal off*/ 
            *dt = state->signal_off_head->time - state->t;     
            update_protein_number_and_fitness(genotype, state, rates, *dt, *effect_of_effector, end_state, mut_step, mut_record); 
            delete_fixed_event_from_head(&(state->signal_off_head),&(state->signal_off_tail));
            if(fixed_effector_effect)
                *effect_of_effector=init_effector_effect;
            else
                *effect_of_effector='d';
            state->protein_number[N_SIGNAL_TF-1]=signal_off_strength;         
            break;
        case 4:     /*turn signal on*/
            *dt = state->signal_on_head->time - state->t;   
            update_protein_number_and_fitness(genotype, state, rates, *dt, *effect_of_effector, end_state, mut_step, mut_record);  
            delete_fixed_event_from_head(&(state->signal_on_head),&(state->signal_on_tail));
            state->protein_number[N_SIGNAL_TF-1]=signal_on_strength;
            if(fixed_effector_effect)                               
                *effect_of_effector=init_effector_effect;            
            else                 
                *effect_of_effector='b';            
            break;	
        case 5: /* finishing burn-in growth rate*/
            *dt=duration_of_burn_in_growth_rate-state->t;     
            update_protein_number_and_fitness(genotype, state, rates, *dt, *effect_of_effector, end_state, mut_step, mut_record);
            state->cumulative_fitness_after_burn_in=state->cumulative_fitness;           
            delete_fixed_event_from_head(&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail));
            break;
        case 6: /* mandatorily updating Pact and Prep*/
            *dt=state->t_to_update_probability_of_binding-state->t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, *effect_of_effector, end_state, mut_step, mut_record);          
            break;
        case 7: /* update signal strength */
            *dt=state->change_signal_strength_head->time-state->t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, *effect_of_effector, end_state, mut_step, mut_record);
            state->protein_number[N_SIGNAL_TF-1]=signal_profile[state->change_signal_strength_head->event_id];
            delete_fixed_event_from_head(&(state->change_signal_strength_head),&(state->change_signal_strength_tail));
            break;      
    }  
    return return_value;
}

/*
 *Gillespie events
 */
int do_Gillespie_event(Genotype *genotype,
                        CellState *state,
                        GillespieRates *rates,
                        float dt,                        
                        RngStream RS,
                        int *end_state,                        
                        int mut_step,
                        Mutation *mut_record)
{
    float x; 
    int return_value;
    return_value=DO_NOTHING;
    while(1)
    {    
        x=RngStream_RandU01(RS)*(rates->total_Gillespie_rate);          
        if (x <= rates->total_mRNA_decay_rate)  /*STOCHASTIC EVENT: an mRNA decay event */
        {    
            return_value=Gillespie_event_mRNA_decay(rates, state, genotype, RS);                
            break;
        } 
        else 
        {               
            x -= rates->total_mRNA_decay_rate;             
            if (x <= rates->total_active_to_intermediate_rate) /* STOCHASTIC EVENT: active transcriptional state to intermediate*/
            {                   
                Gillespie_event_active_to_intermediate(genotype, state, rates, RS);
                break;
            } 
            else 
            {  
                x -= rates->total_active_to_intermediate_rate;  
                if (x <= rates->total_repressed_to_intermediate_rate)  /* repressed to intermediate*/
                {    
                    Gillespie_event_repressed_to_intermediate(rates, state, genotype, RS);
                    break;
                } 
                else 
                { 
                    x-= rates->total_repressed_to_intermediate_rate;                  
                    if (x <= rates->total_intermediate_to_repressed_rate)/* intermediate to repressed */ 
                    {   
                        Gillespie_event_intermediate_to_repressed(rates, state, genotype, RS);
                        break;
                    } 
                    else 
                    {
                        x -= rates->total_intermediate_to_repressed_rate;                        
                        if (x <= rates->total_intermediate_to_active_rate)/* intermediate to active*/
                        {  
                            Gillespie_event_intermediate_to_active(rates, state, genotype, RS); 
                            break;
                        } 
                        else 
                        {
                            x -= rates->total_intermediate_to_active_rate;                            
                            if (x <= rates->total_N_gene_transcript_initiated * TRANSCRIPTINIT) /* STOCHASTIC EVENT: transcription initiation */
                            {   
                                Gillespie_event_transcription_init(rates, state, genotype, dt, RS);
                                break;
                            }   
                        }
                    }
                }
            }       
        }
    }
    return return_value;
}

/*Free linked tables*/
void free_fixedevent(CellState *state)
{
    FixedEvent *temp1, *temp2;
    /*signal_on_starts_end*/
    temp1=state->signal_on_head;
    while(temp1){		
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->signal_on_head=NULL;
    state->signal_on_tail=NULL;
    /*signal_off_starts_end*/
    temp1=state->signal_off_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->signal_off_head=NULL;
    state->signal_off_tail=NULL;
    /*mRNA_transcr_time_end*/
    temp1=state->mRNA_transcr_time_end_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;           
            free(temp2);	
    }
    state->mRNA_transcr_time_end_head=NULL;
    state->mRNA_transcr_time_end_tail=NULL;
    /*mRNA_transl_init_time_end*/
    temp1=state->mRNA_transl_init_time_end_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->mRNA_transl_init_time_end_head=NULL;
    state->mRNA_transl_init_time_end_tail=NULL;
    /*burn_in_growth_rate_end*/
    temp1=state->burn_in_growth_rate_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->burn_in_growth_rate_head=NULL;
    state->burn_in_growth_rate_tail=NULL;
    /*change_singal_strength_head*/
    temp1=state->change_signal_strength_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->change_signal_strength_head=NULL;
    state->change_signal_strength_tail=NULL;
}

/**
 *Calculate the fintess of a given genotype.
 *Essentially calling do_single_timestep until tdevelopment and calculate 
 *average growth rate over tdevelopment.
 */
void calc_avg_growth_rate(  Genotype *genotype,
                            int init_mRNA[NGENES],
                            float init_protein_number[NPROTEINS],
                            RngStream RS_parallel[N_THREADS],                           
                            int mut_step,
                            float GR1[N_REPLICATES],
                            float GR2[N_REPLICATES],                        
                            Mutation *mut_record)       
{ 
    /*Making clones of a genotype, and have the clones run in parallel*/
    #pragma omp parallel num_threads(N_THREADS) 
    {
        int thread_ID=omp_get_thread_num();
//        int thread_ID=0;
        int i,j,k;
        int N_replicates_per_thread=N_REPLICATES/N_THREADS;
        int end_state;
        char effect_of_effector;
        Genotype genotype_clone;
        CellState state_clone;
        GillespieRates rate_clone;
        int init_mRNA_clone[NGENES]; 
        float init_protein_number_clone[NGENES];
        float gr1[N_replicates_per_thread],gr2[N_replicates_per_thread];       
        int mRNA[genotype->ngenes];
        float protein[genotype->ngenes];  
        float *signal_profile;
#if CAUTIOUS
        FILE *fperror;
#endif
        /*alloc space for linked tables and set default values for parameters, in genotype*/
        initialize_cache(&genotype_clone);  
        /*initialize growth rate to 0*/
        for(i=0;i<N_replicates_per_thread;i++)
        {
            gr1[i]=0.0;
            gr2[i]=0.0;
        }
        /*clone genotype and initial mRNA and protein numbers*/
        #pragma omp critical
        { 
            genotype_clone.ngenes=genotype->ngenes;
            genotype_clone.ntfgenes=genotype->ntfgenes;
            genotype_clone.nproteins=genotype->nproteins;
            clone_genotype(genotype, &genotype_clone);  
            for(j=0; j < NGENES; j++) 
            {  
                init_mRNA_clone[j] = init_mRNA[j];
                init_protein_number_clone[j] = init_protein_number[j];
            } 
        } 
        calc_all_binding_sites(&genotype_clone); 
#if PLOT_ALTERNATIVE_FITNESS  /*modify motifs*/       
        for(i=0;i<39;i++)
            genotype_clone.N_motifs[i]=genotype->N_motifs[i];
        for(i=0;i<NGENES;i++)
        {
            genotype_clone.gene_in_core_C1ffl[i]=genotype->gene_in_core_C1ffl[i];
            for(j=0;j<NPROTEINS;j++)
                genotype_clone.TF_in_core_C1ffl[i][j]=genotype->TF_in_core_C1ffl[i][j];
        }
#if FORCE_DIAMOND
//        remove_edges_iteratively(&genotype_clone);    
        modify_topology(&genotype_clone);
#else
        modify_topology(&genotype_clone);
#endif
#endif        
        /*Set initial mRNA and protein number using given values*/
        for(j=N_SIGNAL_TF; j < genotype_clone.ngenes; j++)  //expression of the sensor TF is not considered in the model           
            mRNA[j] = init_mRNA_clone[j];                       
        for(j=N_SIGNAL_TF; j<genotype_clone.nproteins;j++)
        {
            for(k=0;k<genotype_clone.protein_pool[j][0][0];k++)
                protein[genotype_clone.protein_pool[j][1][k]]=(float)init_protein_number_clone[j]/genotype_clone.protein_pool[j][0][0]; //split the initial protein number equally to different copies
                                                                                                                                        //this is to make sure all proteins have equal initial numbers
        }
        /* now calc growth rate under two environments*/
        for(i=0;i<N_replicates_per_thread;i++) /* env 1, usually a constant signal that matches env*/
        {	 
            effect_of_effector=env1_initial_effect_of_effector; // whether the effector is beneficial or deleterious
            end_state=1; //flag to show whether do_single_time_step encounters an error. 1 means no error, 0 means error         
#if EXTERNAL_SIGNAL
            j=RngStream_RandInt(RS_parallel[thread_ID],0,99);
            signal_profile=&(signal_profile_matrix[thread_ID][j][0]);
#else
            signal_profile=NULL;             
#endif            
            /*initialize mRNA and protein numbers, and gene states, in cell*/
            initialize_cell(&genotype_clone,
                            &state_clone,
                            mRNA, 
                            protein,                     
                            env1_t_development);
            /*set how the environment signal should change during simulation*/
            set_signal(&state_clone,
                        env1_t_signal_on,
                        env1_t_signal_off,
                        signal_profile,
                        env1_t_development,
                        env1_signal_strength);          
            /*growth starts at 0 min*/
            state_clone.t = 0.0;
            /*calcualte the rates of cellular activity based on the initial cellular state*/
            calc_all_rates(&genotype_clone, &state_clone, &rate_clone, env1_t_development,INITIALIZATION,thread_ID); 
            /*run growth simulation until tdevelopment or encounter an error*/
            while(state_clone.t<env1_t_development && end_state==1) 
            {
                 do_single_timestep(&genotype_clone, 
                                    &state_clone, 
                                    &rate_clone,                                                                                                                      
                                    &effect_of_effector,
                                    env1_signal_strength,
                                    env1_t_signal_on,
                                    env1_t_signal_off,
                                    env1_initial_effect_of_effector,
                                    env1_fixed_effector_effect,
                                    RS_parallel[thread_ID],                                   
                                    mut_step,  
                                    mut_record,
                                    &end_state,                                    
                                    thread_ID,
                                    env1_t_development,
                                    signal_profile);
            } 
            if(end_state==1) // no error
            {	
                /*calculate average growth rate*/
                gr1[i]=(state_clone.cumulative_fitness-state_clone.cumulative_fitness_after_burn_in)/(env1_t_development-duration_of_burn_in_growth_rate);              
            #if CAUTIOUS                
                if(gr1[i]<0.0)
                {
                    fperror=fopen(error_file,"a+");
                    LOG("negative growth rate at replicate %d in test 1 thread %d at mut step %d\n", i, thread_ID, mut_step);
                    fclose(fperror);
                }              
            #endif
            }
            else
            {                
            	/* if do_single_timestep throws out an error. Re-run replicates*/
                i--;
#if CAUTIOUS
                fperror=fopen(error_file,"a+");
                LOG("Rerun replicates at replicate %d in test 1 thread %d at mut step %d\n", i, thread_ID, mut_step);
                fclose(fperror);
#endif
            }
            /*free linked tables*/
            free_fixedevent(&state_clone);           
        }            
        /*******************env2*********************/
        for(i=0;i<N_replicates_per_thread;i++) 
        {	 
            effect_of_effector=env2_initial_effect_of_effector;  
            end_state=1;          
            signal_profile=NULL;          
            initialize_cell(&genotype_clone,
                            &state_clone,                        
                            mRNA, 
                            protein,                         
                            env2_t_development);
            set_signal( &state_clone,
                        env2_t_signal_on,   
                        env2_t_signal_off,
                        signal_profile,
                        env2_t_development,
                        env2_signal_strength);  
            state_clone.t = 0.0;
            calc_all_rates( &genotype_clone, 
                            &state_clone, 
                            &rate_clone,                            
                            env2_t_development,
                            INITIALIZATION,
                            thread_ID);	
            while(state_clone.t<env2_t_development && end_state==1)
            {
                do_single_timestep(&genotype_clone, 
                                    &state_clone, 
                                    &rate_clone,                                                                                                                  
                                    &effect_of_effector,
                                    env2_signal_strength,
                                    env2_t_signal_on,
                                    env2_t_signal_off,
                                    env2_initial_effect_of_effector,
                                    env2_fixed_effector_effect,                                  
                                    RS_parallel[thread_ID],                                  
                                    mut_step,
                                    mut_record,
                                    &end_state,                               
                                    thread_ID,
                                    env2_t_development,
                                    signal_profile);
            } 
            if(end_state==1)
            {
                gr2[i]=(state_clone.cumulative_fitness-state_clone.cumulative_fitness_after_burn_in)/(env2_t_development-duration_of_burn_in_growth_rate);             
#if CAUTIOUS                
                if(gr2[i]<0.0)
                {
                    fperror=fopen(error_file,"a+");
                    LOG("negative growth rate at replicate %d in test 2 thread %d at mut step %d\n", i, thread_ID, mut_step);
                    fclose(fperror);
                }               
#endif
            }
            else
            {                
                i--;
#if CAUTIOUS
                fperror=fopen(error_file,"a+");
                LOG("Rerun replicates at replicate %d in test 2 thread %d at mut step %d\n", i, thread_ID, mut_step);
                fclose(fperror);
#endif
            }            
            free_fixedevent(&state_clone);            
        } 

        /*free linked tables*/
        for(j=0;j<NGENES;j++)
            free(genotype_clone.all_binding_sites[j]);
        /*pool growth rates from each thread*/
        #pragma omp critical
        {
            j=0;
            for(i=thread_ID*N_replicates_per_thread;i<(thread_ID+1)*N_replicates_per_thread;i++)
            {
                GR1[i]=gr1[j];
                GR2[i]=gr2[j];
                j++;
            }
        }        
    } 
}

/*To get protein number and fitness at given timepoint, we add a special fixed event (
 *sampling_point_end)to take a snapshot of cell state at given timepoint. The additional fixed event 
 *requires some revision to do_singe_timestep etc. The functions that require revision
 *are compiled only if necessary.
 */
#if JUST_PLOTTING 
int does_fixed_event_end_plotting(  CellState *state,
                                    float t) 
{
    int retval=0;
    float t1;
    float t2;
    float t3;
    float t4;
    float t5;    
    float t6;
    float t7;
    float t8;
    t1 = state->mRNA_transcr_time_end_head ? state->mRNA_transcr_time_end_head->time : TIME_INFINITY;
    t2 = state->mRNA_transl_init_time_end_head ? state->mRNA_transl_init_time_end_head->time : TIME_INFINITY;
    t3 = state->signal_off_head? state->signal_off_head->time : TIME_INFINITY;
    t4 = state->signal_on_head ? state->signal_on_head->time : TIME_INFINITY;
    t5 = state->burn_in_growth_rate_head ? state->burn_in_growth_rate_head->time : TIME_INFINITY;    
    t6 = state->sampling_point_end_head ? state->sampling_point_end_head->time : TIME_INFINITY;
    t7 = state->t_to_update_probability_of_binding;
    t8 = state->change_signal_strength_head ? state->change_signal_strength_head->time : TIME_INFINITY;
    if((t1 <= t2) && (t1 <= t) && (t1 <= t3) && (t1 <= t4) && (t1<=t5) && (t1<=t6) && (t1<=t7)&&(t1<=t8))
    {
        retval = 1;	
    }
    else if ((t2 <= t1) && (t2 <= t) && (t2 <= t3) && (t2 <= t4) && (t2<=t5) && (t2<=t6) && (t2<=t7)&&(t2<=t8)) 
    {
        retval = 2;
    }            
    else if ((t3 <= t1) && (t3 <= t) && (t3 <= t2) && (t3 <= t4) && (t3<=t5) && (t3<=t6) && (t3<=t7)&&(t3<=t8)) 
    {
        retval = 3;
    }
    else if ((t4 <= t1) && (t4 <= t) && (t4 <= t2) && (t4 <= t3) && (t4<=t5) && (t4<=t6) && (t4<=t7)&&(t4<=t8)) 
    {
        retval = 4;
    }                    
    else if((t5 <= t1) && (t5 <= t) && (t5 <= t2) && (t5 <= t3) && (t5<=t4) && (t5<=t6) && (t5<=t7)&&(t5<=t8))
    {
        retval = 5;
    }
    else if((t6 <= t1) && (t6 <= t) && (t6 <= t2) && (t6 <= t3) && (t6<=t4) && (t6<=t5) && (t6<=t7)&&(t6<=t8))
    {
        retval=6;
    }
    else if((t7 <= t1) && (t7 <= t) && (t7 <= t2) && (t7 <= t3) && (t7<=t4) && (t7<=t5) && (t7<=t6)&&(t7<=t8))
    {
        retval=7;
    }
    else if((t8 <= t1) && (t8 <= t) && (t8 <= t2) && (t8 <= t3) && (t8<=t4) && (t8<=t5) && (t8<=t6)&&(t8<=t7))            
    {
        retval=8;
    }
    else    
    {
        retval = 0;
    }                       
    return retval;
}

int do_fixed_event_plotting(    Genotype *genotype, 
                                CellState *state, 
                                GillespieRates *rates, 
                                float *dt,                                      
                                int event, 
                                float signal_on_strength,
                                float duration_signal_on,
                                float duration_signal_off,
                                char *effect_of_effector,
                                char init_effector_effect,
                                int fixed_effector_effect,                              
                                int *end_state,
                                float *signal_profile)
{      
    int return_value;
    return_value=DO_NOTHING;
    switch (event) 
    {        
         case 1:     /* a transcription event ends */
            fixed_event_end_transcription(dt, state, rates, genotype,end_state,0,NULL,effect_of_effector); 
            break;
        case 2:     /* a translation initialization event ends */ 
            return_value=fixed_event_end_translation_init(genotype, state, rates, dt, end_state,0,NULL,effect_of_effector);
            state->cell_activated=1;
            break;
        case 3:     /* turn off signal*/ 
            *dt = state->signal_off_head->time - state->t;     
            update_protein_number_and_fitness(genotype, state, rates, *dt, *effect_of_effector, end_state, 0,NULL); 
            delete_fixed_event_from_head(&(state->signal_off_head),&(state->signal_off_tail));            
            state->protein_number[N_SIGNAL_TF-1]=signal_off_strength;
            if(fixed_effector_effect)
                *effect_of_effector=init_effector_effect;
            else
                *effect_of_effector='d';
            state->protein_number[N_SIGNAL_TF-1]=signal_off_strength;         
            break;
        case 4:     /* turn on signal*/
            *dt = state->signal_on_head->time - state->t;   
            update_protein_number_and_fitness(genotype, state, rates, *dt, *effect_of_effector, end_state, 0,NULL);  
            delete_fixed_event_from_head(&(state->signal_on_head),&(state->signal_on_tail));
            state->protein_number[N_SIGNAL_TF-1]=signal_on_strength;
            if(fixed_effector_effect)                               
                *effect_of_effector=init_effector_effect;            
            else                 
                *effect_of_effector='b';
            break;	
        case 5: /* finishing burn-in growth rate*/
            *dt=duration_of_burn_in_growth_rate-state->t;     
            update_protein_number_and_fitness(genotype, state, rates, *dt, *effect_of_effector, end_state, 0,NULL);
            state->cumulative_fitness_after_burn_in=state->cumulative_fitness;           
            delete_fixed_event_from_head(&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail));
            break;
        case 6:
            *dt=state->sampling_point_end_head->time-state->t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, *effect_of_effector, end_state, 0, NULL);
            delete_fixed_event_from_head(&(state->sampling_point_end_head),&(state->sampling_point_end_tail));
            break;
        case 7: /* update force to update Pact and Prep*/
            *dt=state->t_to_update_probability_of_binding-state->t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, *effect_of_effector, end_state, 0,NULL);          
            break;  
        case 8: /* update signal strength */
            *dt=state->change_signal_strength_head->time-state->t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, *effect_of_effector, end_state, 0, NULL);
            state->protein_number[N_SIGNAL_TF-1]=signal_profile[state->change_signal_strength_head->event_id];
            delete_fixed_event_from_head(&(state->change_signal_strength_head),&(state->change_signal_strength_tail));
            break;
    }    
    return return_value;
}

void do_single_timestep_plotting(   Genotype *genotype, 
                                    CellState *state,                         
                                    GillespieRates *rates,                                                           
                                    char *effect_of_effector,  
                                    float signal_on_strength,
                                    float duration_signal_on,
                                    float duration_signal_off,
                                    int fixed_effector_effect,
                                    char initial_effect_of_effector,    
                                    float (*phenotype)[N_TIMP_POINTS],                                    
                                    float fitness[N_TIMP_POINTS], 
                                    int N_update[4],
                                    float max_change[4],
                                    RngStream RS,                                    
                                    int *timepoint,                                   
                                    int *end_state,
                                    int thread_id,
                                    float tdevelopment,
                                    float *signal_profile)
{    
    int event,i,UPDATE_WHAT;      
    float fixed_time; 
    float dt;
    float x;    
    x = expdev(RS);        
    dt = x/rates->total_Gillespie_rate;    
    fixed_time = (state->t+dt<tdevelopment)?(state->t+dt):tdevelopment;
    event = does_fixed_event_end_plotting(state, fixed_time);
    while(event!=0)
    {           
        UPDATE_WHAT=do_fixed_event_plotting(genotype, 
                                            state, 
                                            rates, 
                                            &dt, 
                                            event, 
                                            signal_on_strength,
                                            duration_signal_on, 
                                            duration_signal_off, 
                                            effect_of_effector,
                                            initial_effect_of_effector,
                                            fixed_effector_effect, 
                                            end_state,
                                            signal_profile);
        if(*end_state==0)
            return;        
        if(event==6)/*time to take a snapshot of fitness and protein number*/
        {    
            for(i=0;i<genotype->nproteins;i++)
                phenotype[i][*timepoint]=state->protein_number[i];            
            fitness[*timepoint]=state->instantaneous_fitness;          
            (*timepoint)++;    
        }   
        fixed_time=state->t+dt; 
        state->t += dt;                     
        x -= dt*rates->total_Gillespie_rate;       
        calc_all_rates(genotype, state, rates, tdevelopment,UPDATE_WHAT,thread_id); 
        dt = x/rates->total_Gillespie_rate;	
        if(dt<0.0)	
            dt=TIME_OFFSET;         
        /* check to see there aren't more fixed events to do */
        fixed_time = (state->t+dt<tdevelopment)?(state->t+dt):tdevelopment;         
        event = does_fixed_event_end_plotting(state, fixed_time);                                    
    } 
    /* no remaining fixed events to do in dt, now do stochastic events */  
    /* if we haven't already reached end of development with last
     delta-t, if there is no fixed development time, we always execute
     this  */          
    if (state->t+dt < tdevelopment)
    { 
        update_protein_number_and_fitness(genotype, state, rates, dt, *effect_of_effector, end_state, 0, NULL);  
        if(*end_state==0)
            return; 
        UPDATE_WHAT=do_Gillespie_event(genotype, state, rates, dt, RS, end_state, 0, NULL);
        if(*end_state==0)
            return; 
        /* Gillespie step: advance time to next event at dt */
        state->t += dt;
        calc_all_rates(genotype,state,rates,tdevelopment,UPDATE_WHAT,thread_id);    
    } 
    else 
    { 
        /* do remaining dt */
        dt = tdevelopment - state->t;
        /* final update of protein concentration */
        update_protein_number_and_fitness(genotype, state, rates, dt, *effect_of_effector, end_state, 0, NULL); 
       
        if(*end_state==0)
            return;         
        /* advance to end of development (this exits the outer while loop) */
        state->t = tdevelopment;
    }
}

void calc_avg_growth_rate_plotting( Genotype *genotype,
                                    int init_mRNA[NGENES],
                                    float init_protein_number[NPROTEINS],
                                    RngStream RS[N_THREADS])
{     
    float phenotypeA[N_REPLICATES][genotype->nproteins][N_TIMP_POINTS];
    float phenotypeB[N_REPLICATES][genotype->nproteins][N_TIMP_POINTS];
    float fitnessA[N_REPLICATES][N_TIMP_POINTS];
    float fitnessB[N_REPLICATES][N_TIMP_POINTS];
    int N_updateA[N_REPLICATES][4];
    int N_updateB[N_REPLICATES][4];
    float max_changeA[N_REPLICATES][4];
    float max_changeB[N_REPLICATES][4];
    int l,m,n;
    char filename1[32],filename2[32];
    FILE *fp1,*fp2;
    
    for(l=0;l<N_REPLICATES;l++)
    {
        for(m=0;m<genotype->nproteins;m++)
        {
            for(n=0;n<N_TIMP_POINTS;n++)
            {
                phenotypeA[l][m][n]=0.0;
                phenotypeB[l][m][n]=0.0;               
            }
        }
    }    
    for(l=0;l<N_REPLICATES;l++)
    {   
        for(n=0;n<N_TIMP_POINTS;n++)
        {
            fitnessA[l][n]=0.0;
            fitnessB[l][n]=0.0;               
        }        
    }

    for(l=0;l<N_REPLICATES;l++)
    {
        for(m=0;m<4;m++)
        {
            N_updateA[l][m]=0;
            N_updateB[l][m]=0;    
            max_changeA[l][m]=0.0;
            max_changeB[l][m]=0.0;
        }
    }   
    #pragma omp parallel num_threads(N_THREADS)
    {
        int thread_ID=omp_get_thread_num();
//        int thread_ID=0;
        int i,j,timepoint,k;    
        int end_state;  
        int N_replicates_per_thread=N_REPLICATES/N_THREADS;
//        int N_replicates_per_thread=40;
        char effect_of_effector;
        Genotype genotype_clone;
        CellState state_clone;
        GillespieRates rate_clone;
        float init_mRNA_clone[NGENES]; 
        float init_protein_number_clone[NGENES];    
        int mRNA[genotype->ngenes];
        float protein[genotype->ngenes];
        float *signal_profile;
        
        initialize_cache(&genotype_clone);
        #pragma omp critical
        {  
            clone_genotype(genotype, &genotype_clone);
            for(j=0; j < NGENES; j++) 
            {  
                init_mRNA_clone[j] = init_mRNA[j];
                init_protein_number_clone[j] = init_protein_number[j];
            }        
        }
        
        for(j=N_SIGNAL_TF; j < genotype_clone.ngenes; j++)             
            mRNA[j] = init_mRNA_clone[j];                       
        for(j=N_SIGNAL_TF; j<genotype_clone.nproteins;j++)
        {
            for(k=0;k<genotype_clone.protein_pool[j][0][0];k++)
                protein[genotype_clone.protein_pool[j][1][k]]=(float)init_protein_number_clone[j]/genotype_clone.protein_pool[j][0][0];                
        }           
        calc_all_binding_sites(&genotype_clone);  
        for(i=0;i<N_replicates_per_thread;i++)        
        {             
            effect_of_effector=env1_initial_effect_of_effector;
            end_state=1;
            signal_profile=NULL;  
            initialize_cell(&genotype_clone,
                            &state_clone,
                            mRNA, 
                            protein,                     
                            env1_t_development);
            set_signal(&state_clone,env1_t_signal_on,env1_t_signal_off,signal_profile,env1_t_development,env1_signal_strength);        
            state_clone.t = 0.0;
            calc_all_rates(&genotype_clone, &state_clone, &rate_clone, env1_t_development,INITIALIZATION, thread_ID);//,N_updateA[thread_ID*N_replicates_per_thread+i],max_changeA[thread_ID*N_replicates_per_thread+i]);
            timepoint=0;
            while(state_clone.t<env1_t_development && end_state==1)
            {    
                do_single_timestep_plotting(&genotype_clone, 
                                            &state_clone, 
                                            &rate_clone,                                                                                                                             
                                            &effect_of_effector,
                                            env1_signal_strength,
                                            env1_t_signal_on,
                                            env1_t_signal_off,
                                            env1_fixed_effector_effect,  
                                            env1_initial_effect_of_effector,
                                            phenotypeA[thread_ID*N_replicates_per_thread+i],                                           
                                            fitnessA[thread_ID*N_replicates_per_thread+i],  
                                            N_updateA[thread_ID*N_replicates_per_thread+i], 
                                            max_changeA[thread_ID*N_replicates_per_thread+i],
                                            RS[thread_ID],                                           
                                            &timepoint,                                           
                                            &end_state,
                                            thread_ID,
                                            env1_t_development,
                                            signal_profile);        
            }
            if(end_state==0)
            {                  
                for(j=0;j<4;j++)
                {
                    N_updateA[thread_ID*N_replicates_per_thread+i][j]=0;
                    max_changeA[thread_ID*N_replicates_per_thread+i][j]=0.0;
                }
                i--; 
            }        
            free_fixedevent(&state_clone);
        }
        for(i=0;i<N_replicates_per_thread;i++)        
        {              
            effect_of_effector=env2_initial_effect_of_effector; 
            end_state=1;
            signal_profile=NULL;            
            initialize_cell(&genotype_clone,
                            &state_clone,
                            mRNA, 
                            protein,                     
                            env1_t_development);
            set_signal(&state_clone,env2_t_signal_on,env2_t_signal_off,signal_profile,env2_t_development,env2_signal_strength);        
            state_clone.t = 0.0;
            calc_all_rates(&genotype_clone, &state_clone, &rate_clone, env2_t_development,INITIALIZATION,thread_ID);//N_updateB[thread_ID*N_replicates_per_thread+i],max_changeB[thread_ID*N_replicates_per_thread+i]);        
            timepoint=0;
            while(state_clone.t<env2_t_development && end_state==1)
            {
                do_single_timestep_plotting(&genotype_clone, 
                                            &state_clone, 
                                            &rate_clone,                                                                                                                                       
                                            &effect_of_effector,
                                            env2_signal_strength,
                                            env2_t_signal_on,
                                            env2_t_signal_off,
                                            env2_fixed_effector_effect,  
                                            env2_initial_effect_of_effector,
                                            phenotypeB[thread_ID*N_replicates_per_thread+i],
                                            fitnessB[thread_ID*N_replicates_per_thread+i], 
                                            N_updateB[thread_ID*N_replicates_per_thread+i],
                                            max_changeB[thread_ID*N_replicates_per_thread+i],
                                            RS[thread_ID],                                          
                                            &timepoint,                                            
                                            &end_state,
                                            thread_ID,
                                            env2_t_development,
                                            signal_profile);
            }
            if(end_state==0)
            {                  
                for(j=0;j<4;j++)
                {
                    N_updateB[thread_ID*N_replicates_per_thread+i][j]=0;
                    max_changeB[thread_ID*N_replicates_per_thread+i][j]=0.0;
                }
                i--; 
            }             
            free_fixedevent(&state_clone);
        }        
    }//end of parallel section   
    /*Output protein numbers at each timepoint*/
    /*Each protein has its own file: A is env1 and B env2*/
    /*In each file, columns are replicates, rows are timepoints*/
    for(l=0;l<genotype->nproteins;l++)
    {
        snprintf(filename1,sizeof(char)*32,"phenotypeA_%i",l);
        snprintf(filename2,sizeof(char)*32,"phenotypeB_%i",l);
        fp1=fopen(filename1,"w");
        fp2=fopen(filename2,"w");
        for(m=0;m<N_TIMP_POINTS;m++)
        {
            for(n=0;n<N_REPLICATES;n++)
            {
                fprintf(fp1,"%f ",phenotypeA[n][l][m]);
                fprintf(fp2,"%f ",phenotypeB[n][l][m]);
            }
            fprintf(fp1,"\n");
            fprintf(fp2,"\n");
        }
        fclose(fp1);
        fclose(fp2);
    }
    /*Output growth rate*/
    /*Each environment has its file: A is env1 and B env2*/
    /*In each file, columns are replicates, rows are timepoints*/
    fp1=fopen("fitnessA","w");
    fp2=fopen("fitnessB","w");    
    for(m=0;m<N_TIMP_POINTS;m++)
    {
        for(n=0;n<N_REPLICATES;n++)
        {
            fprintf(fp1,"%f ",fitnessA[n][m]);
            fprintf(fp2,"%f ",fitnessB[n][m]);
        }
        fprintf(fp1,"\n");
        fprintf(fp2,"\n");
    }  
    fclose(fp1);
    fclose(fp2);
           
//    fp1=fopen("N_update_A","w");
//    fp2=fopen("N_update_B","w");
//    for(m=0;m<4;m++)
//    {
//        for(n=0;n<N_REPLICATES;n++)
//        {
//            fprintf(fp1,"%d ",N_updateA[n][m]);
//            fprintf(fp2,"%d ",N_updateB[n][m]);
//        }
//        fprintf(fp1,"\n");
//        fprintf(fp2,"\n");
//    }
//    fclose(fp1);
//    fclose(fp2);  
//    
//    fp1=fopen("max_change_A","w");
//    fp2=fopen("max_change_B","w");
//    for(m=0;m<4;m++)
//    {
//        for(n=0;n<N_REPLICATES;n++)
//        {
//            fprintf(fp1,"%f ",max_changeA[n][m]);
//            fprintf(fp2,"%f ",max_changeB[n][m]);
//        }
//        fprintf(fp1,"\n");
//        fprintf(fp2,"\n");
//    }
//    fclose(fp1);
//    fclose(fp2);   
    
}
#endif /*End of plotting functions*/

/*
 *Set default values and allocate space for variables in Genotype
 */
void initialize_cache(Genotype *genotype)
{
    int j,k;
    FILE *fperror;  
    /*Initialize variables that applies to loci*/
    for(j=0;j<NGENES;j++)
    {
        genotype->which_protein[j]=NA;         
        genotype->recalc_TFBS[j]=YES;
        genotype->which_cluster[j]=NA; 
        genotype->min_N_activator_to_transc[j]=MAX_BINDING+1; /*by default a gene cannot be turned on. 
                                                       *MAX_BINDING is the maximum number of tf that 
                                                       *can bind to a cis-reg sequence.*/        
        genotype->Kd[j]=-1.0;
        for(k=0;k<NGENES;k++)        
            genotype->cisreg_cluster[j][k]=NA;
    }    
    for(j=0;j<NGENES;j++)
        genotype->cisreg_cluster[NGENES][j]=NA;
    /* initialize variables that applies to protein */
    for(j=0;j<NPROTEINS;j++)
    {
        genotype->which_TF_family[j]=NA;
        genotype->protein_pool[j][0][0]=0;
        genotype->TF_family_pool[j][0][0]=0;
        for(k=0;k<NGENES;k++)        
            genotype->protein_pool[j][1][k]=NA; 
        for(k=0;k<NPROTEINS;k++) 
            genotype->TF_family_pool[j][1][k]=NA;
        genotype->protein_identity[j]=NON_TF;
    }
    /* alloc space for binding sites*/
    genotype->N_allocated_elements=MAXELEMENTS;
    for(j=0;j<NGENES;j++)
    {
        genotype->all_binding_sites[j] = malloc(MAXELEMENTS*sizeof(AllTFBindingSites));
        if (!(genotype->all_binding_sites[j])) 
        {
            fperror=fopen(error_file,"a+");
            LOG("Failed to allocate space\n");
            fclose(fperror);            
            exit(-1);
        }
    }
    /*Initialize binding sites summary*/
    for(j=N_SIGNAL_TF;j<NGENES;j++)
    {
        genotype->N_act_BS[j]=0;
        genotype->N_rep_BS[j]=0;
        genotype->binding_sites_num[j]=0;
    }    
}


/**
 * Given the fitness of the resident and a mutant, decide whether the mutant can replace the resident
 */
float try_fixation(Genotype *resident, Genotype *mutant, int N_measurement_resident, int N_measurement_mutant, int *fixation, RngStream RS)
{ 
    float s;
    s=(mutant->fitness-resident->fitness)/fabs(resident->fitness);
    if(s>=minimal_selection_coefficient)
        *fixation=1;
    else          
        *fixation=0;     
    return s;
}

void replay_mutations(  Genotype *genotype_ori,
                        Genotype *genotype_ori_copy,
                        FILE *file_mutation,
                        Mutation *mut_record,
                        int replay_N_steps,
                        RngStream RS)
{
    int i;
//    remove("proportion_c1ffl.txt"); 
//    remove("summary_BS.txt");    
    for(i=0;i<replay_N_steps;i++)
    {
        calc_all_binding_sites(genotype_ori);
//        if(i%OUTPUT_INTERVAL==0)
//            summarize_binding_sites(genotype_ori,i); 
//        find_ffl(genotype_ori);
//        print_core_c1ffls(genotype_ori);   
        clone_genotype(genotype_ori,genotype_ori_copy);
        fscanf(file_mutation,"%c %d %d %s %d %a\n",&(mut_record->mut_type),
                                                    &(mut_record->which_gene),
                                                    &(mut_record->which_nucleotide), 
                                                    mut_record->nuc_diff,               
                                                    &(mut_record->kinetic_type),
                                                    &(mut_record->kinetic_diff));
        reproduce_mutate(genotype_ori_copy,mut_record,RS);        
        clone_genotype(genotype_ori_copy,genotype_ori); 
    }
    printf("Reproduce mutations successfully!\n");
//    calc_all_binding_sites(genotype_ori);
//    find_ffl(genotype_ori);
//    print_core_c1ffls(genotype_ori);
}

#if NEUTRAL
void evolve_neutrally(  Genotype *genotype_ori,
                        Genotype *genotype_ori_copy,                      
                        Mutation *mut_record,                      
                        RngStream RS_main)
{
    int i;    
    FILE *fp;
    DUPLICATION=5.25e-9;                 
    SILENCING =5.25e-9;
    N_EFFECTOR_GENES=2;
    N_TF_GENES=9; 
    miu_ACT_TO_INT_RATE=0.62;
    miu_mRNA_decay=-3.43;
    miu_protein_decay=-4.58;
    miu_translation_init=0.94;   
    
    for(i=0;i<BURN_IN_I;i++)
    {       
        clone_genotype(genotype_ori,genotype_ori_copy);
        mutate(genotype_ori_copy,RS_main,mut_record);        
           
        fp=fopen(mutation_file,"a+");
        fprintf(fp,"%c %d %d '%s' %d %a\n",
                mut_record->mut_type,    
                mut_record->which_gene,
                mut_record->which_nucleotide,
                mut_record->nuc_diff,
                mut_record->kinetic_type,
                mut_record->kinetic_diff);
        fclose(fp);        
        clone_genotype(genotype_ori_copy,genotype_ori);         
        calc_all_binding_sites(genotype_ori);
        find_ffl(genotype_ori); 
        print_core_c1ffls(genotype_ori);        
        /*output network topology every OUTPUT_INTERVAL steps*/ 
        if(i%OUTPUT_INTERVAL==0 && i!=0)
            summarize_binding_sites(genotype_ori,i);        
        /*output a summary of simulation every step*/
        output_genotype(genotype_ori, i);
    } 
    
    DUPLICATION=1.5e-7;                 
    SILENCING = 1.5e-7;
    N_EFFECTOR_GENES=EFFECTOR_GENES;
    N_TF_GENES=TFGENES;
    miu_ACT_TO_INT_RATE=5.22;
    miu_mRNA_decay=-1.13;
    miu_protein_decay=-4.58;
    miu_translation_init=-1.36;    
    
    for(;i<MAX_MUT_STEP;i++)
    {
        clone_genotype(genotype_ori,genotype_ori_copy);
        mutate(genotype_ori_copy,RS_main,mut_record);        
           
        fp=fopen(mutation_file,"a+");
        fprintf(fp,"%c %d %d '%s' %d %a\n",
                mut_record->mut_type,    
                mut_record->which_gene,
                mut_record->which_nucleotide,
                mut_record->nuc_diff,
                mut_record->kinetic_type,
                mut_record->kinetic_diff);
        fclose(fp);        
        clone_genotype(genotype_ori_copy,genotype_ori);         
        calc_all_binding_sites(genotype_ori);
        find_ffl(genotype_ori); 
        print_core_c1ffls(genotype_ori);        
        /*output network topology every OUTPUT_INTERVAL steps*/ 
        if(i%OUTPUT_INTERVAL==0 && i!=0)
            summarize_binding_sites(genotype_ori,i);        
        /*output a summary of simulation every step*/
        output_genotype(genotype_ori, i);
    }
}

#endif

#if JUST_PLOTTING
void run_plotting(  Genotype *genotype_ori,
                    Genotype *genotype_ori_copy,
                    int init_mRNA[NGENES],
                    float init_protein_number[NGENES],
                    RngStream RS_parallel[N_THREADS],
                    Mutation *mut_record,
                    FILE *fp,
                    int replay_N_steps)
{      
    replay_mutations(genotype_ori, genotype_ori_copy, fp, mut_record, replay_N_steps,RS_parallel[0]); 
    calc_all_binding_sites(genotype_ori); 
    print_mutatable_parameters(genotype_ori);
    summarize_binding_sites(genotype_ori,1);   
   // exit(0);
    /* conditions under which the phenotype and fitness is measured */    
    env1_t_development=90.1;
    env2_t_development=90.1;
    env1_signal_strength=10000.0;
    env2_signal_strength=10000.0;
    duration_of_burn_in_growth_rate = 0.0; 
    env1_t_signal_on=200.0;     
    env1_t_signal_off=0.0;
    env2_t_signal_on=10.0;
    env2_t_signal_off=130.0;
    env1_initial_effect_of_effector='b';    
    env2_initial_effect_of_effector='d'; 
    env1_fixed_effector_effect=0;
    env2_fixed_effector_effect=1;
    recalc_new_fitness=1;
    env1_occurence=0.4;
    env2_occurence=0.6;            
    calc_avg_growth_rate_plotting(  genotype_ori, 
                                    init_mRNA, 
                                    init_protein_number,
                                    RS_parallel);    
}
#endif

#if PLOT_ALTERNATIVE_FITNESS
void plot_alternative_fitness(  Genotype *genotype_ori,
                                Genotype *genotype_ori_copy,
                                int init_mRNA[NGENES],
                                float init_protein_number[NGENES],
                                RngStream RS_parallel[N_THREADS],
                                Mutation *mut_record,
                                FILE *file_mutation,
                                FILE *fitness_record,
                                int replay_N_steps)
{
    int i,j,k;
    int steps_to_be_recalc[MAX_MUT_STEP],number_of_steps;
    FILE *alternative_fitness,*original_fitness;   
    char buffer[600];
    env1_t_development=89.9;
    env2_t_development=89.9;
    duration_of_burn_in_growth_rate=0.0;
    env1_signal_strength=10000.0;
    env2_signal_strength=10000.0;
    env1_t_signal_on=200.0;    
    env1_t_signal_off=0.0;     
    env2_t_signal_on=15.0;
    env2_t_signal_off=130.0;
    env1_initial_effect_of_effector='b';
    env2_initial_effect_of_effector='d';
    env1_fixed_effector_effect=0;    
    env2_fixed_effector_effect=1; 
    recalc_new_fitness=5; // make sure its value is smaller than MAX_RECALC_FITNESS
    env1_occurence=0.33;
    env2_occurence=0.67;      
    float GR1[recalc_new_fitness][N_REPLICATES],GR2[recalc_new_fitness][N_REPLICATES];  
    fgets(buffer,600,fitness_record);
    fgets(buffer,600,fitness_record);
    original_fitness=fopen("original_fitness.txt","w");
    for(i=1;i<=replay_N_steps;i++)
    {               
        clone_genotype(genotype_ori,genotype_ori_copy);
        fscanf(file_mutation,"%c %d %d %s %d %a\n",&(mut_record->mut_type),
                                                    &(mut_record->which_gene),
                                                    &(mut_record->which_nucleotide), 
                                                    mut_record->nuc_diff,               
                                                    &(mut_record->kinetic_type),
                                                    &(mut_record->kinetic_diff));
        reproduce_mutate(genotype_ori_copy,mut_record,NULL);        
        clone_genotype(genotype_ori_copy,genotype_ori);       
        fgets(buffer,600,fitness_record);
        if(i>=32201)
        {            
            calc_all_binding_sites(genotype_ori);        
            find_ffl(genotype_ori); 
#if DIRECT_REG
            if(genotype_ori->N_motifs[5]!=0 && genotype_ori->N_motifs[5]==genotype_ori->N_motifs[0])
#else          
    #if FORCE_SINGLE_FFL
            if(genotype_ori->N_motifs[23]!=0 && genotype_ori->N_motifs[9]==0 && genotype_ori->N_motifs[23]==genotype_ori->N_motifs[18]) //force single ffl
    #elif FORCE_DIAMOND //networks contain only AND-gated double C1ffl
            if(genotype_ori->N_motifs[23]!=0 && 
                genotype_ori->N_motifs[9]==0 && 
                genotype_ori->N_motifs[23]==genotype_ori->N_motifs[18] && 
                genotype_ori->N_motifs[27]==0 &&
                genotype_ori->N_motifs[36]==0 &&
                genotype_ori->N_motifs[38]==0) 
    #else //FORCE_OR_GATE. Networks contain only AND-gated double C1ffl   
            if(genotype_ori->N_motifs[14]!=0 && genotype_ori->N_motifs[18]==0 && genotype_ori->N_motifs[14]==genotype_ori->N_motifs[9] && genotype_ori->N_motifs[27]==0) 
    #endif        
#endif             
            {
                for(j=0;j<recalc_new_fitness;j++)                    
                {    
                    calc_avg_growth_rate(   genotype_ori, 
                                            init_mRNA,
                                            init_protein_number,
                                            RS_parallel,                                        
                                            0,
                                            GR1[j],
                                            GR2[j],
                                            mut_record); 
                }

                calc_fitness_stats( genotype_ori,
                                    &(GR1[0]),
                                    &(GR2[0]),
                                    recalc_new_fitness);  

                alternative_fitness=fopen("alternative_fitness.txt","a+");
                fprintf(alternative_fitness,"%d %.10f %.10f %.10f %.10f %.10f %.10f\n",i,
                        genotype_ori->fitness,
                        sqrt(genotype_ori->sq_SE_fitness),
                        genotype_ori->avg_GR1,
                        genotype_ori->avg_GR2,
                        sqrt(genotype_ori->sq_SE_GR1),
                        sqrt(genotype_ori->sq_SE_GR2));
                fclose(alternative_fitness);
                fputs(buffer,original_fitness);                
            }    
            for(j=0;j<NGENES;j++)
            {
                genotype_ori->gene_in_core_C1ffl[j]=0;
                for(k=0;k<NPROTEINS;k++)
                    genotype_ori->TF_in_core_C1ffl[j][k]=0;
            }        
        }
    }  
    fclose(original_fitness);
}
#endif

void run_simulation(    Genotype *genotype_ori, 
                        Genotype *genotype_ori_copy, 
                        float init_protein_number[NPROTEINS],
                        int init_mRNA[NGENES],                    
                        Mutation *mut_record, 
                        int init_N_tot_mutations,
                        int init_step,
                        RngStream RS_main,
                        RngStream RS_parallel[N_THREADS])
{
    FILE *fp;
    int i;
    int max_mut_steps,run_burn_in,N_tot_trials,first_step,end_state=0; 
    first_step=init_step;
    N_tot_trials=init_N_tot_mutations; 
 
    /* first, run burn-in */
    if(BURN_IN_I)
    {
        run_burn_in=1;      
        max_mut_steps=BURN_IN_I;    
        env1_t_development=89.9;
        env2_t_development=89.9;                 
        duration_of_burn_in_growth_rate=0.0;
        env1_signal_strength=10000.0;
        env2_signal_strength=10000.0;
        env1_t_signal_on=200.0;    
        env1_t_signal_off=0.0;     
        env2_t_signal_on=10.0;
        env2_t_signal_off=200.0;
        env1_initial_effect_of_effector='b';
        env2_initial_effect_of_effector='d';
        env1_fixed_effector_effect=0;    
        env2_fixed_effector_effect=1;            
        recalc_new_fitness=5;              
        env1_occurence=0.67;                 
        env2_occurence=0.33;                
        DUPLICATION=5.25e-9;                 
        SILENCING =5.25e-9;
        N_EFFECTOR_GENES=2;
        N_TF_GENES=9; 
        miu_ACT_TO_INT_RATE=0.762; //10%
        miu_Kd=-7.5;       
        miu_translation_init=0.832; //10%
//        sigma_ACT_TO_INT_RATE=0.78; 
//        sigma_mRNA_decay=0.78; 
//        sigma_protein_decay=0.78; 
//        sigma_translation_init=0.78; 
//        sigma_Kd=0.78;
//        N_EFFECTOR_GENES=2;
//        N_TF_GENES=9; 
//        miu_ACT_TO_INT_RATE=0.62;
//        miu_mRNA_decay=-3.43;       
//        miu_translation_init=0.94;
        
        float GR1[recalc_new_fitness][N_REPLICATES],GR2[recalc_new_fitness][N_REPLICATES];
      
        fp=fopen(RuntimeSumm,"a+");
        fprintf(fp,"**********Burn-in_I conditions**********\n");
        fprintf(fp,"BURN_IN_I=%d\n",BURN_IN_I);                
        fprintf(fp,"N_REPLICATES=%d\n",N_REPLICATES);        
        fprintf(fp,"N_recalc_fitness=%d\n",recalc_new_fitness);
        fprintf(fp,"env1_t_development=%f, env2_t_development=%f\n",env1_t_development,env2_t_development);
        fprintf(fp,"Duration of burn-in growth rate=%f\n",duration_of_burn_in_growth_rate);        
        fprintf(fp,"env1: signal on duration=%f min, signal off duration=%f min, initial effector effect=%c, always_deleterious_effector:%d occurrence=%f\n",env1_t_signal_on, env1_t_signal_off, env1_initial_effect_of_effector, env1_fixed_effector_effect, env1_occurence);
        fprintf(fp,"env2: signal on duration=%f min, signal off duration=%f min, initial effector effect=%c, always_deleterious_effector:%d occurrence=%f\n",env2_t_signal_on, env2_t_signal_off, env2_initial_effect_of_effector, env2_fixed_effector_effect, env2_occurence);
        fprintf(fp,"env1 init effecto effect %c, fixed? %d. env2 init effecto effect %c, fixed? %d.\n",env1_initial_effect_of_effector,env1_fixed_effector_effect,env2_initial_effect_of_effector,env2_fixed_effector_effect);
        fclose(fp); 
        
        end_state=evolve_N_steps( genotype_ori, 
                        genotype_ori_copy,
                        &first_step, 
                        max_mut_steps, 
                        &N_tot_trials,                      
                        init_protein_number,
                        init_mRNA,                         
                        mut_record, 
                        RS_main,
                        RS_parallel,
                        run_burn_in);         
       
        fp=fopen(RuntimeSumm,"a+");
        fprintf(fp,"Burn_in_I completes after the %dth step.\n",first_step);
        fclose(fp);        
        if(end_state==-1)
            return;
        if(init_step<BURN_IN_I)
        {
            /*calcualte fitness of the current genotype under post burn_in condition*/
            env1_signal_strength=10000.0;
            env2_signal_strength=10000.0;
            env1_t_signal_on=200.0;    
            env1_t_signal_off=0.0;     
            env2_t_signal_on=10.0;
            env2_t_signal_off=200.0;  
            env1_occurence=0.33;                     
            env2_occurence=0.67;  
            for(i=0;i<recalc_new_fitness;i++)      
            {
                calc_avg_growth_rate(   genotype_ori, 
                                        init_mRNA,
                                        init_protein_number,
                                        RS_parallel,  
                                        BURN_IN_I,
                                        GR1[i],
                                        GR2[i],                                     
                                        mut_record);               
            }
            calc_fitness_stats(genotype_ori,&(GR1[0]),&(GR2[0]),recalc_new_fitness);   
            output_genotype(genotype_ori, BURN_IN_I);
        #if OUTPUT_RNG_SEEDS
            unsigned long seeds[6];      
            RngStream_GetState(RS_main,seeds);
            fp=fopen("RngSeeds.txt","a+");
            fprintf(fp,"%lu %lu %lu %lu %lu %lu ",seeds[0],seeds[1],seeds[2],seeds[3],seeds[4],seeds[5]);            
            for(i=0;i<N_THREADS;i++)
            {
                RngStream_GetState(RS_parallel[i],seeds);
                fprintf(fp,"%lu %lu %lu %lu %lu %lu ",seeds[0],seeds[1],seeds[2],seeds[3],seeds[4],seeds[5]); 
            }
            fprintf(fp,"\n");
            fclose(fp);        
        #endif
            /*output precise fitness*/        
            fp=fopen("precise_fitness.txt","a+");
            fprintf(fp,"%d %d %a %a %a %a %a %a\n",N_tot_trials, 
                                                mut_record->N_hit_bound,
                                                genotype_ori->fitness,
                                                genotype_ori->sq_SE_fitness,
                                                genotype_ori->avg_GR1,
                                                genotype_ori->avg_GR2,
                                                genotype_ori->sq_SE_GR1,
                                                genotype_ori->sq_SE_GR2);
            fclose(fp);        
            /* marks the last step at which all state of the program has been output*/
            fp=fopen("saving_point.txt","w");
            fprintf(fp,"%d %d\n",BURN_IN_I,N_tot_trials);
            fclose(fp);
        }
    }    
    
    /* post-burn-in simulations*/
    run_burn_in=0;
    max_mut_steps=MAX_MUT_STEP;    
    env1_t_development=89.9;
    env2_t_development=89.9;                    
    duration_of_burn_in_growth_rate=0.0;      
    env1_signal_strength=10000.0;
    env2_signal_strength=10000.0;
    env1_t_signal_on=200.0;    
    env1_t_signal_off=0.0;     
    env2_t_signal_on=10.0;
    env2_t_signal_off=200.0;
    env1_initial_effect_of_effector='b';
    env2_initial_effect_of_effector='d';
    env1_fixed_effector_effect=0;    
    env2_fixed_effector_effect=1;                
    recalc_new_fitness=5;                        
    env1_occurence=0.33;                     
    env2_occurence=0.67;                   
    DUPLICATION=1.5e-7;                 
    SILENCING = 1.5e-7;
    N_EFFECTOR_GENES=EFFECTOR_GENES;
    N_TF_GENES=TFGENES;
    miu_ACT_TO_INT_RATE=1.57;
    miu_Kd=-5;       
    miu_translation_init=0.107; 
//    sigma_ACT_TO_INT_RATE=0.78; 
//    sigma_mRNA_decay=0.78; 
//    sigma_protein_decay=0.78; 
//    sigma_translation_init=0.78; 
//    sigma_Kd=0.78;
    
//    miu_ACT_TO_INT_RATE=5.22;
//    miu_mRNA_decay=-1.13;       
//    miu_translation_init=-1.36; 

    fp=fopen(RuntimeSumm,"a+");
    fprintf(fp,"**********Post-burn-in conditions**********\n");
    fprintf(fp,"second phase steps=%d\n",max_mut_steps);                
    fprintf(fp,"N_replicates=%d\n",N_REPLICATES);        
    fprintf(fp,"N_recalc_fitness=%d\n",recalc_new_fitness);
    fprintf(fp,"env1_t_development=%f,env2_t_development=%f\n",env1_t_development,env2_t_development);        
    fprintf(fp,"Duration of burn-in growth rate=%f\n",duration_of_burn_in_growth_rate);         
    fprintf(fp,"env1: signal on duration=%f min, signal off duration=%f min, initial effector effect=%c, always_deleterious_effector:%d occurrence=%f\n",env1_t_signal_on, env1_t_signal_off, env1_initial_effect_of_effector, env1_fixed_effector_effect, env1_occurence);
    fprintf(fp,"env2: signal on duration=%f min, signal off duration=%f min, initial effector effect=%c, always_deleterious_effector:%d occurrence=%f\n",env2_t_signal_on, env2_t_signal_off, env2_initial_effect_of_effector, env2_fixed_effector_effect, env2_occurence);       
    fprintf(fp,"Background signal strength=%f\n",background_signal_strength);
    fprintf(fp,"Signal off strength=%f, env1 signal on strength=%f, env2 signal on strength=%f \n",signal_off_strength,env1_signal_strength,env2_signal_strength);
    fprintf(fp,"env1 init effecto effect %c, fixed? %d. env2 init effecto effect %c, fixed? %d.\n",env1_initial_effect_of_effector,env1_fixed_effector_effect,env2_initial_effect_of_effector,env2_fixed_effector_effect);
    fclose(fp);
    
    end_state=evolve_N_steps( genotype_ori, 
                    genotype_ori_copy,
                    &first_step, 
                    max_mut_steps, 
                    &N_tot_trials,                   
                    init_protein_number,
                    init_mRNA,                 
                    mut_record, 
                    RS_main,
                    RS_parallel,
                    run_burn_in); 
    calc_all_binding_sites(genotype_ori);
    summarize_binding_sites(genotype_ori,max_mut_steps); /*snapshot of the final distribution binding sites */
}

void continue_simulation(   Genotype *genotype_ori, 
                            Genotype *genotype_ori_copy,                            
                            int replay_N_steps,                       
                            float init_protein_number[NPROTEINS],
                            int init_mRNA[NGENES],                
                            Mutation *mut_record, 
                            RngStream RS_main,
                            RngStream RS_parallel[N_THREADS])
{
    int i,j,N_tot_mutations;    
    unsigned long rng_seeds[N_THREADS+1][6];
    char buffer[200]; 
    FILE *fp,*fperror;
    
    /*delete the incomplete lines in the output files*/
    tidy_output_files(output_file,mutation_file);
    
    /* set genotype based on previous steps*/
    fp=fopen(mutation_file,"r");
    if(fp!=NULL)
        replay_mutations(genotype_ori, genotype_ori_copy, fp, mut_record, replay_N_steps,RS_main);
    else
    {
        fperror=fopen(error_file,"a+");
        LOG("cannot open mutation_file\n");
        fclose(fperror);
        exit(-2);
    }
    fclose(fp);

    /* set random number seeds*/
    fp=fopen("RngSeeds.txt","r");
    if(fp!=NULL)
    {
        for(i=0;i<replay_N_steps;i++)
        {
            for(j=0;j<N_THREADS;j++)        
            {
                fscanf(fp,"%lu %lu %lu %lu %lu %lu ",
                        &(rng_seeds[j][0]),
                        &(rng_seeds[j][1]),
                        &(rng_seeds[j][2]),
                        &(rng_seeds[j][3]),
                        &(rng_seeds[j][4]),
                        &(rng_seeds[j][5]));
            }
            fscanf(fp,"%lu %lu %lu %lu %lu %lu \n",
                    &(rng_seeds[N_THREADS][0]),
                    &(rng_seeds[N_THREADS][1]),
                    &(rng_seeds[N_THREADS][2]),
                    &(rng_seeds[N_THREADS][3]),
                    &(rng_seeds[N_THREADS][4]),
                    &(rng_seeds[N_THREADS][5]));
        }
    }
    else
    {
        fperror=fopen(error_file,"a+");
        LOG("cannot open RngSeeds.txt\n");
        fclose(fperror);
        exit(-2);
    }
    fclose(fp);
    RngStream_SetSeed(RS_main,rng_seeds[0]);
    for(i=0;i<N_THREADS;i++)
        RngStream_SetSeed(RS_parallel[i],rng_seeds[i+1]);
    
    /* load fitness,N_tot_mutations,N_hit_boundary*/
    fp=fopen("precise_fitness.txt","r");
    if(fp!=NULL)
    {  
        for(i=0;i<replay_N_steps-1;i++)
            fgets(buffer,200,fp);
        fscanf(fp,"%d %d %a %a %a %a %a %a\n",&N_tot_mutations, 
                                            &(mut_record->N_hit_bound),
                                            &(genotype_ori->fitness),
                                            &(genotype_ori->sq_SE_fitness),
                                            &(genotype_ori->avg_GR1),
                                            &(genotype_ori->avg_GR2),
                                            &(genotype_ori->sq_SE_GR1),
                                            &(genotype_ori->sq_SE_GR2));
    }
    else
    {
        fperror=fopen(error_file,"a+");
        LOG("cannot open precise_fitness.txt\n");
        fclose(fperror);
        exit(-2);
    }        
    fclose(fp);
    /*continue running simulation*/
    run_simulation( genotype_ori, 
                    genotype_ori_copy,  
                    init_protein_number,
                    init_mRNA,                           
                    mut_record, 
                    N_tot_mutations,
                    replay_N_steps+1,
                    RS_main,
                    RS_parallel);    
}

void calc_fitness_stats(Genotype *genotype,
                        float (*GR1)[N_REPLICATES],
                        float (*GR2)[N_REPLICATES],                
                        int N_recalc_fitness)
{
    float avg_GR1=0.0;
    float avg_GR2=0.0;       
    float sum_sq_diff_GR1=0.0;
    float sum_sq_diff_GR2=0.0;   
    float sum_sq_diff_mean_GR=0.0;
    float diff_GR1,diff_GR2,sq_SE_GR1,sq_SE_GR2;    
    int counter=0;
    int i,j;

    for(i=0;i<N_recalc_fitness;i++)
    {
        for(j=0;j<N_REPLICATES;j++)
        {

            avg_GR1+=GR1[i][j];
            avg_GR2+=GR2[i][j];          
            genotype->fitness_measurement[counter]=env1_occurence*GR1[i][j]+env2_occurence*GR2[i][j];
            counter++;
        }
    }
    avg_GR1=avg_GR1/(N_recalc_fitness*N_REPLICATES);
    avg_GR2=avg_GR2/(N_recalc_fitness*N_REPLICATES);  
    
    for(i=0;i<N_recalc_fitness;i++)
    {
        for(j=0;j<N_REPLICATES;j++)
        {
            diff_GR1=GR1[i][j]-avg_GR1;
            diff_GR2=GR2[i][j]-avg_GR2;
            sum_sq_diff_GR1+=pow(diff_GR1,2.0);
            sum_sq_diff_GR2+=pow(diff_GR2,2.0);
            sum_sq_diff_mean_GR+=pow(diff_GR1*env1_occurence+diff_GR2*env2_occurence,2.0);
        }
    }
    sq_SE_GR1=sum_sq_diff_GR1/(N_recalc_fitness*N_REPLICATES*(N_recalc_fitness*N_REPLICATES-1));
    sq_SE_GR2=sum_sq_diff_GR2/(N_recalc_fitness*N_REPLICATES*(N_recalc_fitness*N_REPLICATES-1));     
    genotype->avg_GR1=avg_GR1;
    genotype->avg_GR2=avg_GR2;
    genotype->sq_SE_GR1=sq_SE_GR1*N_recalc_fitness*N_REPLICATES;
    genotype->sq_SE_GR2=sq_SE_GR2*N_recalc_fitness*N_REPLICATES;     
    genotype->fitness=env1_occurence*avg_GR1+env2_occurence*avg_GR2;
    genotype->sq_SE_fitness=sum_sq_diff_mean_GR/(N_recalc_fitness*N_REPLICATES-1)/(N_recalc_fitness*N_REPLICATES); 
}

int evolve_N_steps(Genotype *genotype_ori, 
                    Genotype *genotype_ori_copy,
                    int *init_step, 
                    int max_steps, 
                    int *N_tot_trials,                 
                    float init_protein_number[NPROTEINS],
                    int init_mRNA[NGENES],                 
                    Mutation *mut_record, 
                    RngStream RS_main,
                    RngStream RS_parallel[N_THREADS],
                    int run_burn_in)
{
    int i,j;
    int fixation;
    int N_trials;
    float GR1[recalc_new_fitness][N_REPLICATES],GR2[recalc_new_fitness][N_REPLICATES];
    float score;  
    FILE *fp;
 
    for(i=(*init_step);i<=max_steps;i++)
    {             
        fixation=0;      
        N_trials=0;          
        
        if(*N_tot_trials>MAX_MUTATIONS) /*This is an alternative termination condition, although we usually use MAX_MUT_STEP as termination condition*/
        {
            *init_step=i;
            return 0;
        }
        while(1) /*try mutations until one replaces the current genotype*/
        {	
            N_trials++;
            (*N_tot_trials)++;
            if(N_trials>MAX_TRIALS) /*Tried too many mutation in one step.*/
            {
                fp=fopen(output_file,"a+");                              
                fprintf(fp,"Tried %d mutations, yet none could fix\n",MAX_TRIALS);
                fclose(fp); 
                summarize_binding_sites(genotype_ori,i);
                return -1;
            }
            /*do mutation on a copy of the current genotype*/
            clone_genotype(genotype_ori,genotype_ori_copy); 
            mutate(genotype_ori_copy,RS_main,mut_record);
#if OUTPUT_MUTANT_DETAILS
            /*record every mutation*/
            fp=fopen("MUT_Detail.txt","a+");
            fprintf(fp,"%d %d %c %d %d '%s' %d %a\n",
                    i,
                    *N_tot_trials,
                    mut_record->mut_type,
                    mut_record->which_gene,
                    mut_record->which_nucleotide,
                    mut_record->nuc_diff,
                    mut_record->kinetic_type,
                    mut_record->kinetic_diff);
            fclose(fp); 
#endif
            calc_all_binding_sites(genotype_ori_copy);           
            MAXELEMENTS=genotype_ori_copy->N_allocated_elements;
            /*calculate the fitness of the mutant at low resolution*/
            calc_avg_growth_rate(   genotype_ori_copy,
                                    init_mRNA,
                                    init_protein_number,
                                    RS_parallel,                                   
                                    i,
                                    GR1[0],
                                    GR2[0],
                                    mut_record);
            calc_fitness_stats( genotype_ori_copy,
                                &(GR1[0]),
                                &(GR2[0]),
                                1); // N_recalc_fitness=1 implicitly when calculating the fitness of a mutant.
#if OUTPUT_MUTANT_DETAILS
            /*record low-resolution fitness*/
            fp=fopen("Mut_detail_fitness.txt","a+");
            fprintf(fp,"%.10f %.10f\n",
                    genotype_ori_copy->fitness,
                    genotype_ori_copy->sq_SE_fitness);
            fclose(fp);     
#endif
            /*Can the mutant replace the current genotype?*/
            score=try_fixation(genotype_ori, genotype_ori_copy, N_REPLICATES*recalc_new_fitness, N_REPLICATES, &fixation, RS_main); 
            /*If yes, record relevant info*/
            if(fixation==1)
            {                    
                fp=fopen(output_file,"a+");
                fprintf(fp,"%d %d %d %d %c %f ",i, *N_tot_trials, N_trials,mut_record->N_hit_bound,mut_record->mut_type,score);
                fclose(fp);
                fp=fopen(mutation_file,"a+");
                fprintf(fp,"%c %d %d '%s' %d %a\n",
                        mut_record->mut_type,    
                        mut_record->which_gene,
                        mut_record->which_nucleotide,
                        mut_record->nuc_diff,
                        mut_record->kinetic_type,
                        mut_record->kinetic_diff);
                fclose(fp);
                break;
            }
        }
        /*replace the current genotype by overwriting it*/
        clone_genotype(genotype_ori_copy,genotype_ori);
        calc_all_binding_sites(genotype_ori); 
        /*increase the accuracy of the fitness of the new genotype*/ 
        if(i!=BURN_IN_I) // If we are at the last step of BURN_IN, we should calculate the fitness under the post-burn-in condition, which is done in run_simulation.
        {
            for(j=1;j<recalc_new_fitness;j++)   

                calc_avg_growth_rate(   genotype_ori, 
                                        init_mRNA,
                                        init_protein_number,
                                        RS_parallel,                                    
                                        i,
                                        GR1[j],
                                        GR2[j],
                                        mut_record);  
            calc_fitness_stats( genotype_ori,
                                &(GR1[0]),
                                &(GR2[0]),
                                recalc_new_fitness); 
        }
        /*calculate the number of c1-ffls every step*/
        find_ffl(genotype_ori); 
        print_core_c1ffls(genotype_ori); 
        /*output network topology every OUTPUT_INTERVAL steps*/
        if(i%OUTPUT_INTERVAL==0 && i!=0) 
            summarize_binding_sites(genotype_ori,i);        
        /*output a summary of simulation every step*/
        if(i!=BURN_IN_I)
            output_genotype(genotype_ori, i);        
        /* output rng seeds*/
        #if OUTPUT_RNG_SEEDS
        unsigned long seeds[6];
        if(i!=BURN_IN_I)
        {
            RngStream_GetState(RS_main,seeds);
            fp=fopen("RngSeeds.txt","a+");
            fprintf(fp,"%lu %lu %lu %lu %lu %lu ",seeds[0],seeds[1],seeds[2],seeds[3],seeds[4],seeds[5]);            
            for(j=0;j<N_THREADS;j++)
            {
                RngStream_GetState(RS_parallel[j],seeds);
                fprintf(fp,"%lu %lu %lu %lu %lu %lu ",seeds[0],seeds[1],seeds[2],seeds[3],seeds[4],seeds[5]); 
            }
            fprintf(fp,"\n");
            fclose(fp);
        }
        #endif
        /*output precise fitness*/
        if(i!=BURN_IN_I)
        {
            fp=fopen("precise_fitness.txt","a+");
            fprintf(fp,"%d %d %a %a %a %a %a %a\n",*N_tot_trials, 
                                                mut_record->N_hit_bound,
                                                genotype_ori->fitness,
                                                genotype_ori->sq_SE_fitness,
                                                genotype_ori->avg_GR1,
                                                genotype_ori->avg_GR2,
                                                genotype_ori->sq_SE_GR1,
                                                genotype_ori->sq_SE_GR2);
            fclose(fp);        
            /* marks the last step at which all state of the program has been output*/
            fp=fopen("saving_point.txt","w");
            fprintf(fp,"%d %d\n",i,*N_tot_trials);
            fclose(fp);
        }
    } 
    *init_step=i;
    return 0;
}

int init_run_pop(unsigned long int seeds[6], int CONTINUE)
{  
    int i,j,k,buffer_int;
    Genotype genotype_ori;
    Genotype genotype_ori_copy;    
    int init_mRNA[NGENES]; 
    float init_protein_number[NGENES];   
    int init_step=0;
    Mutation mut_record;
    FILE *fp;    
    RngStream RS_main,RS_parallel[N_THREADS];    
//    RngStream RS_main, RS_parallel[N_THREADS*5+1];
    /*create threads*/
    omp_set_num_threads(N_THREADS);    
    /* initialize random number seeds*/
    RngStream_SetPackageSeed(seeds);    
    RS_main=RngStream_CreateStream("Main");
    for(i=0; i < N_THREADS; i++)
        RS_parallel[i]=RngStream_CreateStream("");
    
//    for(i=0; i < 5*N_THREADS+1; i++)
//        RS_parallel[i]=RngStream_CreateStream("");

    /* set initial numbers of mRNA and protein*/
    for(i=N_SIGNAL_TF; i < NGENES; i++) /* loop through tf genes*/
    {
        init_mRNA[i]=0;
        init_protein_number[i]=0.0;
    }
    /* initialize genotype */
    initialize_cache(&genotype_ori);
    initialize_cache(&genotype_ori_copy);    
    initialize_genotype(&genotype_ori, RS_main);
    genotype_ori_copy.ngenes=genotype_ori.ngenes;
    genotype_ori_copy.ntfgenes=genotype_ori.ntfgenes;
    genotype_ori_copy.nproteins=genotype_ori.nproteins; 
    /* initialize mut_record */
    mut_record.kinetic_diff=0.0;
    mut_record.kinetic_type=-1;
    mut_record.mut_type='\0';
    mut_record.nuc_diff[0]='\0';
    mut_record.nuc_diff[1]='\0';
    mut_record.nuc_diff[2]='\0';
    mut_record.which_gene=-1;
    mut_record.which_nucleotide=-1;
    mut_record.N_hit_bound=0;
    /*different modes of simulation*/    
    if(CONTINUE) /* continue a simulation from a previously saved state*/
    {
        int replay_N_steps=0;
        fp=fopen("saving_point.txt","r");
        fscanf(fp,"%d %d",&replay_N_steps,&buffer_int);
        fclose(fp);
        /*Load external signal profile*/
        fp=fopen("noisy_signal.txt","r");
        if(fp!=NULL)
        {        
            for(j=0;j<200;j++)
            {
                for(k=0;k<15;k++)
                {           
                    fscanf(fp,"%f\n",&(signal_profile_matrix[0][j][k]));
                    for(i=1;i<N_THREADS;i++)
                        signal_profile_matrix[i][j][k]=signal_profile_matrix[0][j][k];                    
                }
            }
            fclose(fp);
        } 
        if(replay_N_steps!=0)
        {                      
            continue_simulation(&genotype_ori, 
                                &genotype_ori_copy,                                
                                replay_N_steps,                                
                                init_protein_number,
                                init_mRNA,                                
                                &mut_record, 
//                                RS_parallel[N_THREADS],
//                                &(RS_parallel[N_THREADS+1])
                                RS_main,
                                RS_parallel
                                );                                                                    
            
        }  
    }
    else /* otherwise the simulation starts over from beginning*/
    {   
        /* record the initial network topology*/
        summarize_binding_sites(&genotype_ori,init_step); /*snapshot of the initial (0) distribution binding sites */   
        find_ffl(&genotype_ori); 
        print_core_c1ffls(&genotype_ori);
#if JUST_PLOTTING 
        fp=fopen("MUT.txt","r");    
        if(fp!=NULL)
        {
            printf("LOAD MUTATION RECORD SUCCESSFUL!\n"); 
            run_plotting(   &genotype_ori,
                            &genotype_ori_copy,
                            init_mRNA,
                            init_protein_number,
                            RS_parallel,
                            &mut_record,                                
                            fp,
                            MAX_MUT_STEP);
            fclose(fp);
        }
#elif PLOT_ALTERNATIVE_FITNESS
        FILE *fp2;
        fp=fopen("MUT.txt","r");    
        fp2=fopen("output.txt","r");
        if(fp!=NULL && fp2!=NULL)
        {
            printf("LOAD MUTATION RECORD SUCCESSFUL!\n");
            plot_alternative_fitness(&genotype_ori,&genotype_ori_copy,init_mRNA,init_protein_number,
                                        RS_parallel,&mut_record,fp,fp2,MAX_MUT_STEP);
            fclose(fp);
        }
#else     
        if(!SKIP_INITIAL_GENOTYPE)/* get the fitness of the initial genotype */ 
        {                     
            env1_t_development=89.9;
            env2_t_development=89.9;
            duration_of_burn_in_growth_rate=0.0;
            env1_signal_strength=10000.0;
            env2_signal_strength=10000.0;
            env1_t_signal_on=200.0;    
            env1_t_signal_off=0.0;     
            env2_t_signal_on=15.0;
            env2_t_signal_off=200.0;
            env1_initial_effect_of_effector='b';
            env2_initial_effect_of_effector='d';
            env1_fixed_effector_effect=0;    
            env2_fixed_effector_effect=1; 
            recalc_new_fitness=5; 
            env1_occurence=0.67;
            env2_occurence=0.33;    
            float GR1[recalc_new_fitness][N_REPLICATES],GR2[recalc_new_fitness][N_REPLICATES];
            /*Load external signal profile if there is one*/
            fp=fopen("noisy_signal.txt","r");
            if(fp!=NULL)
            {        
                for(j=0;j<200;j++)
                {
                    for(k=0;k<15;k++)
                    {           
                        fscanf(fp,"%f\n",&(signal_profile_matrix[0][j][k]));
                        for(i=1;i<N_THREADS;i++)
                            signal_profile_matrix[i][j][k]=signal_profile_matrix[0][j][k];                    
                    }
                }
                fclose(fp);
            } 
            for(i=0;i<recalc_new_fitness;i++)      
            {
                calc_avg_growth_rate(   &genotype_ori, 
                                        init_mRNA,
                                        init_protein_number,
                                        RS_parallel,     
//                                        &(RS_parallel[N_THREADS+1]),
                                        0,
                                        GR1[i],
                                        GR2[i],
                                        &mut_record);   

            }
                                      
            calc_fitness_stats(&genotype_ori,&(GR1[0]),&(GR2[0]),recalc_new_fitness);
            fp=fopen(output_file,"a+");
            fprintf(fp,"step N_tot_mut_tried N_mut_tried_this_step N_hit_bound fixed_mutation score fitness se_fitness avg_GR1 avg_GR2 std_GR1 std_GR2 N_genes N_proteins N_activator N_repressor\n");
            fprintf(fp,"0 0 0 0 na na %.10f %.10f %.10f %.10f %.10f %.10f %d %d %d %d \n",  
                    genotype_ori.fitness,
                    sqrt(genotype_ori.sq_SE_fitness),
                    genotype_ori.avg_GR1,
                    genotype_ori.avg_GR2,
                    sqrt(genotype_ori.sq_SE_GR1),
                    sqrt(genotype_ori.sq_SE_GR2),
                    genotype_ori.ngenes,
                    genotype_ori.nproteins,
                    genotype_ori.N_act,
                    genotype_ori.N_rep);
            fclose(fp); 
        }
#if RUN_FULL_SIMULATION    
        run_simulation( &genotype_ori, 
                        &genotype_ori_copy,                      
                        init_protein_number,
                        init_mRNA,                        
                        &mut_record, 
                        0, // this is the number of total mutations that have been tried
                        1, // this tells the program from which step the simulation begins
//                        RS_parallel[N_THREADS],
//                        &(RS_parallel[N_THREADS+1])
                        RS_main,
                        RS_parallel
                                );
    
#else // no selection at all, just randomly shuffle network topology and kinetic parameters      
        fp=fopen(output_file,"a+");
        fprintf(fp,"step N_tot_mut_tried N_mut_tried_this_step fixed_mutation score fitness se_fitness avg_GR1 avg_GR2 std_GR1 std_GR2 N_genes N_proteins N_activator N_repressor\n");
        fprintf(fp,"0 0 0 na na 0.0 0.0 0.0 0.0 0.0 0.0 0 0 0 0 \n");
        fclose(fp);      
        evolve_neutrally(   &genotype_ori,
                            &genotype_ori_copy,                          
                            &mut_record,                                               
                            RS_main);     
#endif
#endif
    }
    print_mutatable_parameters(&genotype_ori);
    /*release memory*/
    release_memory(&genotype_ori,&genotype_ori_copy,&RS_main, RS_parallel);
    return 1;	
}

void release_memory(Genotype *genotype_ori,Genotype *genotype_ori_copy,RngStream *RS_main,RngStream RS_parallel[N_THREADS])
{
    int i;
    
    for(i=0;i<NGENES;i++)
    {
        free(genotype_ori->all_binding_sites[i]);
        free(genotype_ori_copy->all_binding_sites[i]);
    }
    RngStream_DeleteStream (RS_main);
    for(i=0;i<N_THREADS;i++)
        RngStream_DeleteStream (&RS_parallel[i]);    
}

void output_genotype(Genotype *genotype, int step_i)
{
    FILE *OUTPUT;
    OUTPUT=fopen(output_file,"a+");
    fprintf(OUTPUT,"%.10f %.10f %.10f %.10f %.10f %.10f %d %d %d %d \n",  
            genotype->fitness,
            sqrt(genotype->sq_SE_fitness),
            genotype->avg_GR1,
            genotype->avg_GR2,
            sqrt(genotype->sq_SE_GR1),
            sqrt(genotype->sq_SE_GR2),
            genotype->ngenes,
            genotype->nproteins,
            genotype->N_act,
            genotype->N_rep);
    fclose(OUTPUT);   
}


void print_core_c1ffls(Genotype *genotype)
{
    FILE *fp; 
    fp=fopen("proportion_c1ffl.txt","a+");
    fprintf(fp,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
                                                                                                genotype->N_motifs[0],
                                                                                                genotype->N_motifs[1],
                                                                                                genotype->N_motifs[2],
                                                                                                genotype->N_motifs[3],
                                                                                                genotype->N_motifs[4],
                                                                                                genotype->N_motifs[5],
                                                                                                genotype->N_motifs[6],
                                                                                                genotype->N_motifs[7],
                                                                                                genotype->N_motifs[8],
                                                                                                genotype->N_motifs[9],
                                                                                                genotype->N_motifs[10],
                                                                                                genotype->N_motifs[11],
                                                                                                genotype->N_motifs[12],
                                                                                                genotype->N_motifs[13],
                                                                                                genotype->N_motifs[14],
                                                                                                genotype->N_motifs[15],
                                                                                                genotype->N_motifs[16],
                                                                                                genotype->N_motifs[17],
                                                                                                genotype->N_motifs[18],
                                                                                                genotype->N_motifs[19],
                                                                                                genotype->N_motifs[20],
                                                                                                genotype->N_motifs[21],
                                                                                                genotype->N_motifs[22],
                                                                                                genotype->N_motifs[23],
                                                                                                genotype->N_motifs[24],
                                                                                                genotype->N_motifs[25],
                                                                                                genotype->N_motifs[26],
                                                                                                genotype->N_motifs[27],
                                                                                                genotype->N_motifs[28],
                                                                                                genotype->N_motifs[29],
                                                                                                genotype->N_motifs[30],
                                                                                                genotype->N_motifs[31],
                                                                                                genotype->N_motifs[32],
                                                                                                genotype->N_motifs[33],
                                                                                                genotype->N_motifs[34],
                                                                                                genotype->N_motifs[35], 
                                                                                                genotype->N_motifs[36],
                                                                                                genotype->N_motifs[37],
                                                                                                genotype->N_motifs[38], 
                                                                                                genotype->N_act_genes,
                                                                                                genotype->N_act_genes_reg_by_env,
                                                                                                genotype->N_act_genes_not_reg_by_env,
                                                                                                genotype->protein_pool[genotype->nproteins-1][0][0],
                                                                                                genotype->ntfgenes-genotype->N_act_genes-N_SIGNAL_TF,
                                                                                                genotype->normalizer1,
                                                                                                genotype->normalizer2,
                                                                                                genotype->normalizer3,
                                                                                                genotype->normalizer4);
        fclose(fp); 
}

void summarize_binding_sites(Genotype *genotype, int step_i)
{
    FILE *OUTPUT1;
    int i,j;
    int table[NGENES][NGENES];

//    int table_reduced_sites[NGENES][NGENES];
        
    for(i=0;i<genotype->ngenes;i++)
    {
        for(j=0;j<genotype->ngenes;j++)
        {
            table[i][j]=0; 
//            table_reduced_sites[i][j]=0;
        }
    }
    
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)        
    {
        //print_binding_sites_distribution(genotype,i,step_i);
//        resolve_overlapping_sites(genotype,i,table_reduced_sites[i]);
        for(j=0;j<genotype->binding_sites_num[i];j++)
        {            
            table[i][genotype->all_binding_sites[i][j].tf_id]++; /*the numbers of the binding sites of each tf on promoter i*/
        }    
    }
    
    /*Output all binding sites*/ 
    OUTPUT1=fopen("summary_BS.txt","a+");
    fprintf(OUTPUT1,"step %d\n",step_i);
    fprintf(OUTPUT1,"Promoter ");    
    for(i=0;i<genotype->nproteins-1;i++)
    {
        if(genotype->protein_identity[i]==ACTIVATOR)
        {
            fprintf(OUTPUT1," A%d ",i);
//            N_act++;
        }
        if(genotype->protein_identity[i]==REPRESSOR)
        {
            fprintf(OUTPUT1," R%d ",i);
//            N_rep++;
        }
    }
    fprintf(OUTPUT1,"MIN_ACT_TO_TURN_ON\n");
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
        if(i<9)
            fprintf(OUTPUT1,"%d        ",i+1);
        else
            fprintf(OUTPUT1,"%d       ",i+1);
        
        for(j=0;j<genotype->nproteins-1;j++)
        {
            if(table[i][j]<10)
                fprintf(OUTPUT1," %d  ",table[i][j]);
            else
                fprintf(OUTPUT1," %d ",table[i][j]);
        }
        if(genotype->which_protein[i]==genotype->nproteins-1)
            fprintf(OUTPUT1,"S");
        else
        {
            if(genotype->protein_identity[genotype->which_protein[i]]==ACTIVATOR)
                fprintf(OUTPUT1,"A%d",genotype->which_protein[i]); 
            if(genotype->protein_identity[genotype->which_protein[i]]==REPRESSOR)
                fprintf(OUTPUT1,"R%d",genotype->which_protein[i]);
        }
        fprintf(OUTPUT1," %d \n",genotype->min_N_activator_to_transc[i]);
    }
    fprintf(OUTPUT1,"\n"); 
    fprintf(OUTPUT1,"\n"); 
    fclose(OUTPUT1);
    
    /*Output non-overlapping sites*/
//    OUTPUT1=fopen("summary_non_overlapping_BS","a+");
//    fprintf(OUTPUT1,"step %d\n",step_i);
//    fprintf(OUTPUT1,"Promoter ");
//    for(i=0;i<genotype->nproteins-1;i++)
//    {
//        if(genotype->protein_identity[i]==1)
//        {
//            fprintf(OUTPUT1," A%d ",i);
////            N_act++;
//        }
//        if(genotype->protein_identity[i]==0)
//        {
//            fprintf(OUTPUT1," R%d ",i);
////            N_rep++;
//        }
//    }
//    fprintf(OUTPUT1,"\n");
//    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
//    {
//        if(i<9)
//            fprintf(OUTPUT1,"%d        ",i+1);
//        else
//            fprintf(OUTPUT1,"%d       ",i+1);
//        
//        for(j=0;j<genotype->nproteins-1;j++)
//        {
//            if(table[i][j]<10)
//                fprintf(OUTPUT1," %d  ",table[i][j]);
//            else
//                fprintf(OUTPUT1," %d ",table[i][j]);
//        }
//        if(genotype->which_protein[i]==genotype->nproteins-1)
//            fprintf(OUTPUT1,"S");
//        else
//        {
//            if(genotype->protein_identity[genotype->which_protein[i]]==1)
//                fprintf(OUTPUT1,"A%d",genotype->which_protein[i]); 
//            if(genotype->protein_identity[genotype->which_protein[i]]==0)
//                fprintf(OUTPUT1,"R%d",genotype->which_protein[i]);
//        }
//        fprintf(OUTPUT1,"\n");
//    }
//    fprintf(OUTPUT1,"\n"); 
//    fprintf(OUTPUT1,"\n"); 
//    fclose(OUTPUT1);    
}

void summarize_binding_sites2(Genotype *genotype,int step_i)
{
    FILE *OUTPUT1;
    int i,j;
    int table[NGENES][NGENES];
    Genotype genotype_copy;
    initialize_cache(&genotype_copy);
    genotype_copy.ngenes=genotype->ngenes;
    genotype_copy.ntfgenes=genotype->ntfgenes;
    genotype_copy.nproteins=genotype->nproteins;    
    clone_genotype(genotype,&genotype_copy);
    calc_all_binding_sites2(&genotype_copy);
        
    for(i=0;i<genotype_copy.ngenes;i++)    
        for(j=0;j<genotype_copy.ngenes;j++)
            table[i][j]=0; 

    for(i=N_SIGNAL_TF;i<genotype_copy.ngenes;i++)        
        for(j=0;j<genotype_copy.binding_sites_num[i];j++)                    
            table[i][genotype_copy.all_binding_sites[i][j].tf_id]++; /*the numbers of the binding sites of each tf on promoter i*/
   
    /*Output all binding sites*/ 
    OUTPUT1=fopen("summary_BS.txt","a+");
    fprintf(OUTPUT1,"step %d\n",step_i);
    fprintf(OUTPUT1,"Promoter ");    
    for(i=0;i<genotype_copy.nproteins-1;i++)
    {
        if(genotype_copy.protein_identity[i]==ACTIVATOR)
            fprintf(OUTPUT1," A%d ",i);
        if(genotype_copy.protein_identity[i]==REPRESSOR)
            fprintf(OUTPUT1," R%d ",i);
    }
    fprintf(OUTPUT1,"MIN_ACT_TO_TURN_ON\n");
    for(i=N_SIGNAL_TF;i<genotype_copy.ngenes;i++)
    {
        if(i<9)
            fprintf(OUTPUT1,"%d        ",i+1);
        else
            fprintf(OUTPUT1,"%d       ",i+1);
        
        for(j=0;j<genotype_copy.nproteins-1;j++)
        {
            if(table[i][j]<10)
                fprintf(OUTPUT1," %d  ",table[i][j]);
            else
                fprintf(OUTPUT1," %d ",table[i][j]);
        }
        if(genotype_copy.which_protein[i]==genotype_copy.nproteins-1)
            fprintf(OUTPUT1,"S");
        else
        {
            if(genotype_copy.protein_identity[genotype_copy.which_protein[i]]==ACTIVATOR)
                fprintf(OUTPUT1,"A%d",genotype_copy.which_protein[i]); 
            if(genotype_copy.protein_identity[genotype_copy.which_protein[i]]==REPRESSOR)
                fprintf(OUTPUT1,"R%d",genotype_copy.which_protein[i]);
        }
        fprintf(OUTPUT1," %d \n",genotype_copy.min_N_activator_to_transc[i]);
    }
    fprintf(OUTPUT1,"\n"); 
    fprintf(OUTPUT1,"\n"); 
    fclose(OUTPUT1);   
    for(i=0;i<NGENES;i++)
        free(genotype_copy.all_binding_sites[i]);
}

/*find subtypes of C1-FFLs, FFL-in-diamond, and diamonds*/
void find_ffl(Genotype *genotype)
{
    int i,j,k,l,cluster_size;
    int found_bs;
    int gene_id,gene_id_copy,site_id,protein_id,N_copies,N_activators;
    int master_TF,aux_TF;
    int activators[NPROTEINS];
    int copies_reg_by_env[NGENES],copies_not_reg_by_env[NGENES],N_copies_reg_by_env,N_copies_not_reg_by_env;
    int hindrance[NPROTEINS][NPROTEINS];
    int pos_binding_sites_of_j[MAXELEMENTS],pos_binding_sites_of_k[MAXELEMENTS],N_binding_sites_of_j,N_binding_sites_of_k; 
    int N_motifs;
    int N_copies_reg_by_env_per_protein, N_copies_not_reg_by_env_per_protein;
    int N2_copies_reg_by_env,N2_copies_not_reg_by_env;
    
    /*record which genes and TFs form the motifs of interest. We use this information to perturb motifs*/
#if PLOT_ALTERNATIVE_FITNESS
    for(i=0;i<NGENES;i++)
    {
        genotype->gene_in_core_C1ffl[i]=0;
        for(j=0;j<NPROTEINS;j++)
            genotype->TF_in_core_C1ffl[i][j]=0;
    }
#endif
    /*look for motifs only if the genotype contains non-signal activators*/
    if(genotype->N_act>N_SIGNAL_TF)
    {        
        for(i=0;i<39;i++)
            genotype->N_motifs[i]=0; 
        /*loop through each cisreg_cluster. Genes in a cluster have the same binding sites*/
        i=0;   
        while(genotype->cisreg_cluster[i][0]!=NA) //not an empty cluster
        {              
            /*get a gene in the current cis-reg cluster*/
            gene_id=genotype->cisreg_cluster[i][0];
            gene_id_copy=gene_id;
            if(genotype->which_protein[gene_id]==genotype->nproteins-1) // is a effector gene
            {                 
                /*find how many genes in this cluster*/
                cluster_size=0;
                while(genotype->cisreg_cluster[i][cluster_size]!=NA)
                    cluster_size++;        
                /*reset table for recording activators that regulate effector gene*/
                for(j=0;j<NPROTEINS;j++)
                    activators[j]=0;
                /*scan binding sites for tfs that regulate gene_id*/
                for(j=0;j<genotype->binding_sites_num[gene_id];j++)
                {
                    protein_id=genotype->all_binding_sites[gene_id][j].tf_id;
                    if(genotype->protein_identity[protein_id]==ACTIVATOR && genotype->all_binding_sites[gene_id][j].mis_match<=CUT_OFF_MISMATCH_ON_EFFECTOR) // is a binding site of an activator 
                        activators[protein_id]=1;
                }
                /* move non-zeros entries in activators to the front. */
                k=0; //marks the entry to which a record is copied
                j=0;
                N_activators=0;
                while(j<genotype->nproteins)
                {
                    if(activators[j]!=0)               
                    {
                        activators[k]=j;
                        if(k!=j)
                            activators[j]=0;
                        k++; 
                        N_activators++; //total number of activators that regulate effector
                    }    
                    j++;
                } 
                
                /*build a table to show hindrance between binding sites on effector gene*/
                /*if every binding site of j can hinder all the binding site of k, denote hindrance[j][k]=1. Otherwise 0*/
                /*hindrance[j][j]=1 means at most one TF j can bind to the effector gene */
                /*A strict AND gate between j and K should have H[j][j]=H[k][k]=1, and H[j][k]=0*/
                /*set hindrance to 1 to begin with*/
                for(j=0;j<NPROTEINS;j++)
                {
                    for(l=0;l<NPROTEINS;l++)
                        hindrance[j][l]=1;
                }
                for(j=0;j<N_activators;j++)
                {
                    /*reset position table*/
                    for(l=0;l<MAXELEMENTS;l++)
                        pos_binding_sites_of_j[l]=-CISREG_LEN; 
                    /*list the positions of all the binding sites of activator j*/
                    N_binding_sites_of_j=0;
                    gene_id=genotype->cisreg_cluster[i][0];
                    for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                    {
                        if(genotype->all_binding_sites[gene_id][site_id].tf_id==activators[j] && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_ON_EFFECTOR)
                        {
                            pos_binding_sites_of_j[N_binding_sites_of_j]=genotype->all_binding_sites[gene_id][site_id].BS_pos;
                            N_binding_sites_of_j++;
                        }
                    } 
                    if(N_binding_sites_of_j==1 && genotype->min_N_activator_to_transc[gene_id_copy]==2) // there is only one binding site                     
                    {
                        hindrance[activators[j]][activators[j]]=1;
                    }
                    else
                    {
                        /*if the first bs hinders the last bs, we can be sure at most one binding site of j can be bound by TF j */
                        if(pos_binding_sites_of_j[N_binding_sites_of_j-1]-pos_binding_sites_of_j[0]<TF_ELEMENT_LEN+2*HIND_LENGTH)                     
                            hindrance[activators[j]][activators[j]]=1;
                        else
                            hindrance[activators[j]][activators[j]]=0;
                    }
                    
                    /*list the positions of all the binding sites of activator j+1...*/                  
                    N_binding_sites_of_k=0;
                    for(k=j+1;k<N_activators;k++)
                    {
                        for(l=0;l<MAXELEMENTS;l++)
                            pos_binding_sites_of_k[l]=-CISREG_LEN;
                        N_binding_sites_of_k=0;
                        gene_id=genotype->cisreg_cluster[i][0];                     
                        for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                        {
                            if(genotype->all_binding_sites[gene_id][site_id].tf_id==activators[k] && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_ON_EFFECTOR)
                            {
                                pos_binding_sites_of_k[N_binding_sites_of_k]=genotype->all_binding_sites[gene_id][site_id].BS_pos;
                                N_binding_sites_of_k++;
                            } 
                        }  
                        /*determine the relation between j and j+1, ...*/
                        if(abs(pos_binding_sites_of_j[N_binding_sites_of_j-1]-pos_binding_sites_of_k[0])>= TF_ELEMENT_LEN+2*HIND_LENGTH ||
                            abs(pos_binding_sites_of_k[N_binding_sites_of_k-1]-pos_binding_sites_of_j[0])>= TF_ELEMENT_LEN+2*HIND_LENGTH)
                        {
                            hindrance[activators[j]][activators[k]]=0; 
                            hindrance[activators[k]][activators[j]]=0;
                        }
                    }
                }
                
                /*make lists of gene copies that are regulated by environmental signal and of those that are not*/
                /*reset counters of tf genes*/
                N_copies_reg_by_env=0;
                N_copies_not_reg_by_env=0;
                for(j=0;j<NGENES;j++)
                {
                    copies_reg_by_env[j]=-1;
                    copies_not_reg_by_env[j]=-1;
                }
                for(j=0;j<N_activators;j++)
                {
                    if(activators[j]>=N_SIGNAL_TF)
                    {  
                        N_copies=genotype->protein_pool[activators[j]][0][0]; // the number of genes encoding j
                        for(k=0;k<N_copies;k++)
                        {
                            gene_id=genotype->protein_pool[activators[j]][1][k];
                            /*reset flag*/
                            found_bs=0;
                            for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                            {
                                if(genotype->all_binding_sites[gene_id][site_id].tf_id==N_SIGNAL_TF-1 && genotype->all_binding_sites[gene_id][site_id].tf_id<=CUT_OFF_MISMATCH_SIGNAL_TO_TF)
                                {
                                    found_bs=1;
                                    break;
                                }
                            } 
                            if(found_bs)
                            {
                                copies_reg_by_env[N_copies_reg_by_env]=gene_id;
                                N_copies_reg_by_env++;                                
                            }
                            else
                            {
                                copies_not_reg_by_env[N_copies_not_reg_by_env]=gene_id;
                                N_copies_not_reg_by_env++;                                 
                            }
                        }                        
                    }
                } 

                /*******************************count motifs ************************/                
#if DIRECT_REG  
                /*count c1-ffls formed by signal and one signal-regulated copy */
                if(activators[0]==N_SIGNAL_TF-1 || activators[1]==N_SIGNAL_TF-1)
                {
                    N_motifs=cluster_size*N_copies_reg_by_env;
                    genotype->N_motifs[0]+=N_motifs;           

                    for(j=0;j<N_copies_reg_by_env;j++)
                    {
                        protein_id=genotype->which_protein[copies_reg_by_env[j]];
                        if(hindrance[N_SIGNAL_TF-1][protein_id]) 
                        {
                            if(hindrance[N_SIGNAL_TF-1][N_SIGNAL_TF-1]) // signal alone cannot activate effector
                            {
                                if(hindrance[protein_id][protein_id])
                                    genotype->N_motifs[1]+=cluster_size; // effectively no regulation
                                else
                                    genotype->N_motifs[2]+=cluster_size; // I3
                            }   
                            else
                            {
                                if(hindrance[protein_id][protein_id])
                                    genotype->N_motifs[3]+=cluster_size; //I1
                                else
                                    genotype->N_motifs[4]+=cluster_size; //impossible
                            }
                        }
                        else
                        {
                            if(hindrance[N_SIGNAL_TF-1][N_SIGNAL_TF-1])
                            {
                                if(hindrance[protein_id][protein_id])
                                {
                                    genotype->N_motifs[5]+=cluster_size; // AND-gated 
#if FORCE_OR_GATE
                                    genotype->gene_in_core_C1ffl[gene_id_copy]=1;
#if !FORCE_MASTER_CONTROLLED
                                    genotype->TF_in_core_C1ffl[gene_id_copy][protein_id]=1;
#endif
#elif FORCE_DIAMOND
                                    genotype->gene_in_core_C1ffl[copies_reg_by_env[j]]=1;                                        
#endif
                                }
                                else
                                    genotype->N_motifs[6]+=cluster_size; // aux. tf controlled
                            }   
                            else
                            {
                                if(hindrance[protein_id][protein_id])
                                    genotype->N_motifs[7]+=cluster_size; // signal controll
                                else
                                    genotype->N_motifs[8]+=cluster_size; // OR-gated
                            }
                        }
                    }
                } 
#else
                /*count c1-ffl formed by one env-regulated copy and one unregulated copy*/
                for(j=0;j<N_copies_reg_by_env;j++)
                {
                    for(k=0;k<N_copies_not_reg_by_env;k++)
                    {
                        /*j and k must encode TFs of different families*/
                        if(genotype->which_TF_family[genotype->which_protein[copies_reg_by_env[j]]]!=genotype->which_TF_family[genotype->which_protein[copies_not_reg_by_env[k]]])
                        {
                            /*search bs of k on j*/
                            gene_id=copies_reg_by_env[j];
                            found_bs=0;
                            protein_id=genotype->which_protein[copies_not_reg_by_env[k]];
                            for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                            {
                                if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_TF_TO_TF)
                                {
                                    found_bs=1;
                                    break;
                                }
                            }
                            if(!found_bs) // j is not regulated by k  
                            {
                                /*search bs of j on k*/
                                gene_id=copies_not_reg_by_env[k];
                                found_bs=0;
                                protein_id=genotype->which_protein[copies_reg_by_env[j]];
                                for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                                {
                                    if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_TF_TO_TF)
                                    {
                                        found_bs=1;
                                        break;
                                    }
                                }
                                if(found_bs) // k is regulated by j, then we found an isolated C1FFL                           
                                {    
                                    genotype->N_motifs[9]+=cluster_size; //any logic
                                    master_TF=genotype->which_protein[copies_reg_by_env[j]];
                                    aux_TF=genotype->which_protein[copies_not_reg_by_env[k]];
                                    if(hindrance[master_TF][aux_TF])
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[10]+=cluster_size; //no regulation
                                            else  
                                                genotype->N_motifs[11]+=cluster_size; //emergent I3
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[12]+=cluster_size; //emergent I1
                                            else  
                                                genotype->N_motifs[13]+=cluster_size; //no possible
                                        }    
                                    }
                                    else
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                            {
                                                genotype->N_motifs[14]+=cluster_size; //AND-gated 
#if FORCE_OR_GATE                                   
                                                genotype->gene_in_core_C1ffl[gene_id_copy]=1;
#if FORCE_MASTER_CONTROLLED
                                                genotype->TF_in_core_C1ffl[gene_id_copy][master_TF]=1;
#else
                                                genotype->TF_in_core_C1ffl[gene_id_copy][aux_TF]=1;
#endif                                    
#elif FORCE_DIAMOND                                    
                                                genotype->gene_in_core_C1ffl[copies_not_reg_by_env[k]]=1;
                                                genotype->TF_in_core_C1ffl[copies_not_reg_by_env[k]][master_TF]=1;                                                          
#endif
                                            }
                                            else  
                                                genotype->N_motifs[15]+=cluster_size; //aux. TF controlled
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[16]+=cluster_size; //master TF controlled
                                            else  
                                                genotype->N_motifs[17]+=cluster_size; //OR-gated
                                        }                                     
                                    } 
                                }
                            }
                            else //j is regulated by k
                            {
                                /*search bs of j on k*/
                                gene_id=copies_not_reg_by_env[k];
                                found_bs=0;
                                protein_id=genotype->which_protein[copies_reg_by_env[j]];
                                for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                                {
                                    if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_TF_TO_TF)
                                    {
                                        found_bs=1;
                                        break;
                                    }
                                }
                                if(found_bs)
                                    genotype->N_motifs[36]+=cluster_size; // strongly-connected
                                else
                                    genotype->N_motifs[37]+=cluster_size; // ffl w/o source
                            }
                        }
                    }
                }
                /*count motifs formed by two env-regulated copies*/
                /*also count parallel structures*/
                for(j=0;j<N_copies_reg_by_env;j++)
                {
                    for(k=j+1;k<N_copies_reg_by_env;k++)
                    {
                        /*k and j must encode different TFs*/
                        if(k!=j && genotype->which_TF_family[genotype->which_protein[copies_reg_by_env[j]]]!=genotype->which_TF_family[genotype->which_protein[copies_reg_by_env[k]]])
                        {                    
                            /*search bs of k on j*/
                            gene_id=copies_reg_by_env[j];
                            found_bs=0;
                            protein_id=genotype->which_protein[copies_reg_by_env[k]];
                            for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                            {
                                if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_TF_TO_TF)
                                {
                                    found_bs=1;
                                    break;
                                }
                            }
                            if(!found_bs) // j is not regulated by k
                            {
                                /*search bs of j on k*/
                                gene_id=copies_reg_by_env[k];
                                found_bs=0;
                                protein_id=genotype->which_protein[copies_reg_by_env[j]];
                                for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                                {
                                    if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_TF_TO_TF)
                                    {
                                        found_bs=1;
                                        break;
                                    }
                                }
                                if(found_bs)  // k is regulated by j, then we've found a FFL-in-diamond                          
                                {                                    
                                    genotype->N_motifs[18]+=cluster_size; //any logic
                                    master_TF=genotype->which_protein[copies_reg_by_env[j]];
                                    aux_TF=genotype->which_protein[copies_reg_by_env[k]];
                                    if(hindrance[master_TF][aux_TF])
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[19]+=cluster_size; //no regulation
                                            else  
                                                genotype->N_motifs[20]+=cluster_size; //emergent I3
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[21]+=cluster_size; //emergent I1
                                            else  
                                                genotype->N_motifs[22]+=cluster_size; //impossible
                                        }    
                                    }
                                    else
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                            {
                                                genotype->N_motifs[23]+=cluster_size; //AND-gated
#if FORCE_OR_GATE                                        
                                                genotype->gene_in_core_C1ffl[gene_id_copy]=1;
#if FORCE_MASTER_CONTROLLED
                                                genotype->TF_in_core_C1ffl[gene_id_copy][master_TF]=1;
#else
                                                genotype->TF_in_core_C1ffl[gene_id_copy][aux_TF]=1;
#endif                                        
#elif FORCE_DIAMOND                                           
                                                genotype->gene_in_core_C1ffl[copies_reg_by_env[k]]=1;
                                                genotype->TF_in_core_C1ffl[copies_reg_by_env[k]][master_TF]=1;                                           
#elif FORCE_SINGLE_FFL                                           
                                                genotype->gene_in_core_C1ffl[copies_not_reg_by_env[k]]=1;  
#endif 
                                            }
                                            else  
                                                genotype->N_motifs[24]+=cluster_size; //aux. TF controlled
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[25]+=cluster_size; //master TF controlled
                                            else  
                                                genotype->N_motifs[26]+=cluster_size;  //OR-gated
                                        }                                     
                                    }
                                }
                                else // k is not regulated by j either, then we've found a diamond:)
                                {
                                    genotype->N_motifs[27]+=cluster_size;  //any logic
                                    master_TF=genotype->which_protein[copies_reg_by_env[j]];
                                    aux_TF=genotype->which_protein[copies_reg_by_env[k]];
                                    if(hindrance[master_TF][aux_TF])
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[28]+=cluster_size; //no regulation
                                            else  
                                                genotype->N_motifs[29]+=cluster_size; //emergent I3
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[30]+=cluster_size; //emergent I1
                                            else  
                                                genotype->N_motifs[31]+=cluster_size; //impossible
                                        }    
                                    }
                                    else
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])                                            
                                            {
                                                genotype->N_motifs[32]+=cluster_size; //AND-gated
#if FORCE_OR_GATE                                        
                                                genotype->gene_in_core_C1ffl[gene_id_copy]=1;
#if FORCE_MASTER_CONTROLLED
                                                genotype->TF_in_core_C1ffl[gene_id_copy][master_TF]=1;
#else
                                                genotype->TF_in_core_C1ffl[gene_id_copy][aux_TF]=1;
#endif                             
#endif 
                                            }
                                            else  
                                                genotype->N_motifs[33]+=cluster_size; //Aux. TF controlled
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[34]+=cluster_size; //master TF controlled
                                            else  
                                                genotype->N_motifs[35]+=cluster_size; //OR-gated
                                        }                                     
                                    } 
                                }
                            }
                            else // j is regulated by k
                            {
                                /*search bs of j on k*/
                                gene_id=copies_reg_by_env[k];
                                found_bs=0;
                                protein_id=genotype->which_protein[copies_reg_by_env[j]];
                                for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                                {
                                    if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_TF_TO_TF)
                                    {
                                        found_bs=1;
                                        break;
                                    }
                                }
                                if(!found_bs) // k is not regulated by j, we've found another FFL-in-diamond                           
                                {
                                    genotype->N_motifs[18]+=cluster_size; 
                                    master_TF=genotype->which_protein[copies_reg_by_env[k]];
                                    aux_TF=genotype->which_protein[copies_reg_by_env[j]];
                                    if(hindrance[master_TF][aux_TF])
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[19]+=cluster_size;
                                            else  
                                                genotype->N_motifs[20]+=cluster_size;
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[21]+=cluster_size;
                                            else  
                                                genotype->N_motifs[22]+=cluster_size;
                                        }    
                                    }
                                    else
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                            {
                                                genotype->N_motifs[23]+=cluster_size; //AND-gated
#if FORCE_OR_GATE                                      
                                                genotype->gene_in_core_C1ffl[gene_id_copy]=1;
#if FORCE_MASTER_CONTROLLED
                                                genotype->TF_in_core_C1ffl[gene_id_copy][master_TF]=1;
#else
                                                genotype->TF_in_core_C1ffl[gene_id_copy][aux_TF]=1;
#endif                                       
#elif FORCE_DIAMOND                                       
                                                genotype->gene_in_core_C1ffl[copies_reg_by_env[j]]=1;
                                                genotype->TF_in_core_C1ffl[copies_reg_by_env[j]][master_TF]=1;                                        
#elif FORCE_SINGLE_FFL                                       
                                                genotype->gene_in_core_C1ffl[copies_reg_by_env[j]]=1; 
#endif
                                            }
                                            else  
                                                genotype->N_motifs[24]+=cluster_size; //aux. TF controlled
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[25]+=cluster_size; //master controlled
                                            else  
                                                genotype->N_motifs[26]+=cluster_size; //OR-gated
                                        }                                     
                                    }                                   
                                }
                                else
                                    genotype->N_motifs[38]+=cluster_size; // strong-connected
                            }
                        }
                    }
                }  
#endif
            }
            i++;
        }
        
        /*count the number of activator genes, excluding the signal*/
        genotype->N_act_genes=0;
        for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
        {
            if(genotype->protein_identity[i]==ACTIVATOR)
            {
                genotype->N_act_genes+=genotype->protein_pool[i][0][0];
            }
        }
        /*count the numbers of activator genes that are regulated by the signal and of those that are not*/
        genotype->N_act_genes_reg_by_env=0;
        genotype->N_act_genes_not_reg_by_env=0;
        for(gene_id=N_SIGNAL_TF;gene_id<genotype->ngenes;gene_id++)
        {
            protein_id=genotype->which_protein[gene_id];
            if(genotype->protein_identity[protein_id]==ACTIVATOR)
            {
                found_bs=0;
                for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                {
                    if(genotype->all_binding_sites[gene_id][site_id].tf_id==N_SIGNAL_TF-1)
                    {
                        found_bs=1;
                        break;
                    }
                }
                if(found_bs)
                    genotype->N_act_genes_reg_by_env++;
                else
                    genotype->N_act_genes_not_reg_by_env++;                            
            }
        }  
        /*calculate normalizers of the numbers of motifs*/
#if DIRECT_REG
        genotype->normalizer1=genotype->N_act_genes;  
        genotype->normalizer2=genotype->N_act_genes;
        genotype->normalizer3=genotype->N_act_genes;
        genotype->normalizer4=genotype->N_act_genes;
#else
        genotype->normalizer1=0;
        genotype->normalizer2=0;
        genotype->normalizer3=0;
        genotype->normalizer4=0;
        N_copies_reg_by_env=0;
        N_copies_not_reg_by_env=0;
        N2_copies_reg_by_env=0;
        N2_copies_not_reg_by_env=0;
        for(i=N_SIGNAL_TF;i<genotype->nproteins-1;i++)
        {            
            if(genotype->protein_identity[i]==ACTIVATOR)
            {
                N_copies_reg_by_env_per_protein=0;
                N_copies_not_reg_by_env_per_protein=0;
                for(j=0;j<genotype->protein_pool[i][0][0];j++)
                {
                    gene_id=genotype->protein_pool[i][1][j];
                    found_bs=0;
                    for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                    {
                        if(genotype->all_binding_sites[gene_id][site_id].tf_id==N_SIGNAL_TF-1)
                        {
                            found_bs=1;
                            break;
                        }
                    }
                    if(found_bs)
                        N_copies_reg_by_env_per_protein++;
                    else
                        N_copies_not_reg_by_env_per_protein++;
                }
                if(N_copies_reg_by_env_per_protein>N_copies_not_reg_by_env_per_protein)
                {
                    genotype->normalizer1-=genotype->protein_pool[i][0][0]*(genotype->protein_pool[i][0][0]-1)/2;
                    N_copies_reg_by_env+=genotype->protein_pool[i][0][0];
                    genotype->normalizer2-=genotype->protein_pool[i][0][0]*(genotype->protein_pool[i][0][0]-1)/2;
                    N2_copies_reg_by_env+=genotype->protein_pool[i][0][0];
                }
                else if(N_copies_reg_by_env_per_protein==N_copies_not_reg_by_env_per_protein)
                {
                    genotype->normalizer1-=genotype->protein_pool[i][0][0]*(genotype->protein_pool[i][0][0]-1)/2;
                    N_copies_reg_by_env+=genotype->protein_pool[i][0][0];                    
                    N2_copies_not_reg_by_env+=genotype->protein_pool[i][0][0];
                }
                else
                {
                    N_copies_not_reg_by_env+=genotype->protein_pool[i][0][0];
                    N2_copies_not_reg_by_env+=genotype->protein_pool[i][0][0];
                }
            }
        }  
        genotype->normalizer3=genotype->normalizer1+N_copies_reg_by_env*(N_copies_reg_by_env-1)/2;
        genotype->normalizer1+=N_copies_reg_by_env*N_copies_not_reg_by_env+N_copies_reg_by_env*(N_copies_reg_by_env-1)/2;
        genotype->normalizer2+=N2_copies_reg_by_env*N2_copies_not_reg_by_env+N2_copies_reg_by_env*(N2_copies_reg_by_env-1)/2;
        genotype->normalizer4=N_copies_reg_by_env*N_copies_not_reg_by_env;
#endif
        genotype->normalizer1*=genotype->protein_pool[genotype->nproteins-1][0][0];
        genotype->normalizer2*=genotype->protein_pool[genotype->nproteins-1][0][0];
        genotype->normalizer3*=genotype->protein_pool[genotype->nproteins-1][0][0];
        genotype->normalizer4*=genotype->protein_pool[genotype->nproteins-1][0][0];
    }    
    else // if the genotype does not contain non-signal activator
    {
        for(i=0;i<39;i++)        
            genotype->N_motifs[i]=0;  
        genotype->N_act_genes=0; 
        genotype->N_act_genes_not_reg_by_env=0;
        genotype->N_act_genes_reg_by_env=0;
        genotype->normalizer1=0;
        genotype->normalizer2=0;
        genotype->normalizer3=0;
        genotype->normalizer4=0;
    }
}

#if PLOT_ALTERNATIVE_FITNESS
void modify_topology(Genotype *genotype)
{
    int gene_id;
    for(gene_id=N_SIGNAL_TF;gene_id < genotype->ngenes;gene_id++)
    {        
#if FORCE_OR_GATE
        add_binding_site(genotype, gene_id);
#else
        remove_binding_sites(genotype, gene_id);
#endif
    }
}

void add_binding_site(Genotype *genotype, int gene_id)
{    
    float temp;  
    int i,j;
    FILE *fperror;
    /*make sure there is enough room*/
    if (genotype->binding_sites_num[gene_id] + 1 >= genotype->N_allocated_elements) 
    {  
        while(genotype->N_allocated_elements<=genotype->binding_sites_num[gene_id]+1)
            genotype->N_allocated_elements+=100;

        for(i=0;i<NGENES;i++)
        {
            genotype->all_binding_sites[i] = realloc(genotype->all_binding_sites[i], genotype->N_allocated_elements*sizeof(AllTFBindingSites));
            if (!genotype->all_binding_sites[i]) 
            {
                fperror=fopen(error_file,"a+");
                LOG("error in calc_all_binding_sites_copy\n");
                fclose(fperror);
                exit(-1);                                       
            }     
        }                    
    }  
    
#if DIRECT_REG
    if(genotype->N_motifs[0]!=0 && genotype->gene_in_core_C1ffl[gene_id]==1) //modify Bs only if gene_id is regulated by AND gate  
#else
    if(genotype->gene_in_core_C1ffl[gene_id]==1)
#endif
    {
#if DIRECT_REG
        temp=1.0;
    #if FORCE_MASTER_CONTROLLED /* add a binding site of the signal to gene_id.*/        
        for(i=0;i<genotype->binding_sites_num[gene_id];i++)
        {
            if(genotype->all_binding_sites[gene_id][i].tf_id==N_SIGNAL_TF-1) 
                temp=(temp>genotype->all_binding_sites[gene_id][i].Kd)?genotype->all_binding_sites[gene_id][i].Kd:temp;  // find a strong binding site                      
        }
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].tf_id = N_SIGNAL_TF-1;
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].Kd = temp;                   
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].BS_pos = 2*CISREG_LEN; //no way to block any existing BS
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].mis_match = 0;
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].N_hindered = 0;
        genotype->binding_sites_num[gene_id]++;
        genotype->N_act_BS[gene_id]++;
    #endif
#endif
        for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
        {
            if(genotype->TF_in_core_C1ffl[gene_id][i]==1)
            {
                temp=1.0;
                for(j=0;j<genotype->binding_sites_num[gene_id];j++)
                {
                    if(genotype->all_binding_sites[gene_id][j].tf_id==i) 
                        temp=(temp>genotype->all_binding_sites[gene_id][j].Kd)?genotype->all_binding_sites[gene_id][j].Kd:temp;  // find a strong binding site                      
                }
                genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].tf_id = i;
                genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].Kd = temp;                   
                genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].BS_pos = (2+i)*CISREG_LEN;
                genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].mis_match = 0;
                genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].N_hindered = 0;
                genotype->binding_sites_num[gene_id]++;
                genotype->N_act_BS[gene_id]++;                
            }
        }        
    }       
    /*the maximum number of activator binding sites that do not exclude each other is increased*/
    genotype->max_unhindered_sites[gene_id][1]++;    
}

/* ignore TF x when searching binding sites on gene y*/
/* this function is almost the same as calc_all_binding_sites_copy*/
void remove_binding_sites(Genotype *genotype, int gene_id)
{
    int i, j, k;
    int match,match_rc;  
    int N_hindered_BS=0;   
    int N_binding_sites=0;
    int start_TF; 
    genotype->N_act_BS[gene_id]=0;
    genotype->N_rep_BS[gene_id]=0;
    genotype->max_hindered_sites[gene_id]=0;  
    //some helper pointer 
    char *tf_seq;
    char *cis_seq;
    char *tf_seq_rc; 
    cis_seq=&(genotype->cisreg_seq[gene_id][0]); 
  
    for(i=3; i < CISREG_LEN-TF_ELEMENT_LEN-3; i++) /* scan promoter */
    {  
        /*calc the number of BS within the hindrance range*/
        N_hindered_BS=0;        
        if(N_binding_sites>0)
        {
            for(j=0;j<N_binding_sites;j++)
            {
               if(genotype->all_binding_sites[gene_id][j].BS_pos> i-TF_ELEMENT_LEN-2*HIND_LENGTH)
                    N_hindered_BS++;
            }
        }  
        /* loop through TF proteins */        
        #if !DIRECT_REG 
            if(genotype->which_protein[gene_id]==genotype->nproteins-1) // if the gene is an effector gene
                start_TF=N_SIGNAL_TF;// the environmental signals cannot directly regulate the selection gene
            else
                start_TF=0;
        #else
            start_TF=0;
        #endif
        
        for (k=start_TF; k < genotype->nproteins-1; k++) 
        { 
#if FORCE_DIAMOND
    #if !DIRECT_REG
            if(!(genotype->gene_in_core_C1ffl[gene_id]==1 && genotype->TF_in_core_C1ffl[gene_id][k]==1))
    #else
            if(!(  genotype->N_motifs[0]!=0 
                && genotype->gene_in_core_C1ffl[gene_id]==1
                && k==start_TF))
    #endif
            {
#elif FORCE_SINGLE_FFL        
            if(!(  genotype->N_motifs[18]!=0 
                && genotype->gene_in_core_C1ffl[gene_id]==1
                && k==start_TF))       
            {    
#endif
            tf_seq=&(genotype->tf_seq[k][0]);
            tf_seq_rc=&(genotype->tf_seq_rc[k][0]);            
            /*find BS on the template strand*/
            match=0;
            for (j=i; j < i+TF_ELEMENT_LEN; j++) 
                if (cis_seq[j] == tf_seq[j-i]) match++; 
            if (match >= NMIN)
            {                              
                genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;                      
                genotype->all_binding_sites[gene_id][N_binding_sites].Kd=KD2APP_KD*genotype->Kd[k]*pow(NS_Kd/genotype->Kd[k],(float)(TF_ELEMENT_LEN-match)/(TF_ELEMENT_LEN-NMIN+1));
                genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i ; 
                genotype->all_binding_sites[gene_id][N_binding_sites].mis_match = TF_ELEMENT_LEN-match;             
                genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                N_hindered_BS++;              
                N_binding_sites++;
                if(genotype->protein_identity[k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
            }
            else /*find BS on the non-template strand.*/
            {
                match_rc=0;
                for (j=i; j < i+TF_ELEMENT_LEN; j++)                
                    if (cis_seq[j] == tf_seq_rc[j-i]) match_rc++;
                if (match_rc >= NMIN)
                {                   
                    genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;                                     
                    genotype->all_binding_sites[gene_id][N_binding_sites].Kd=KD2APP_KD*genotype->Kd[k]*pow(NS_Kd/genotype->Kd[k],(float)(TF_ELEMENT_LEN-match_rc)/(TF_ELEMENT_LEN-NMIN+1));
                    genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i;
                    genotype->all_binding_sites[gene_id][N_binding_sites].mis_match = TF_ELEMENT_LEN-match_rc;
                    genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                    N_hindered_BS++;                  
                    N_binding_sites++;                 
                    if(genotype->protein_identity[k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
                }
            } 
#if FORCE_DIAMOND
            }
#elif FORCE_SINGLE_FFL
            }
#endif
        }/* looping through TFs ends */
    }/*end of promoter scanning*/ 
    
    genotype->binding_sites_num[gene_id]=N_binding_sites;  
    genotype->N_rep_BS[gene_id]=N_binding_sites-(genotype->N_act_BS[gene_id]);
    /* calculate max_hindered_sites */    
    for(i=0;i<genotype->binding_sites_num[gene_id];i++)
    {
        genotype->max_hindered_sites[gene_id]=(genotype->max_hindered_sites[gene_id] > genotype->all_binding_sites[gene_id][i].N_hindered)?
                                      genotype->max_hindered_sites[gene_id] : genotype->all_binding_sites[gene_id][i].N_hindered;           
    }    
   
    int act_BS[MAXELEMENTS][2],rep_BS[MAXELEMENTS][2];
    int N_act_BS,N_rep_BS;    
    N_act_BS=1;
    N_rep_BS=1;
    for(i=0;i<genotype->binding_sites_num[gene_id];i++) /* make lists BS by their types*/    
    {
        if(genotype->protein_identity[genotype->all_binding_sites[gene_id][i].tf_id]==ACTIVATOR)
        {
            act_BS[N_act_BS][0]=i;
            N_act_BS++;
        } 
        else if(genotype->protein_identity[genotype->all_binding_sites[gene_id][i].tf_id]==REPRESSOR)
        {
            rep_BS[N_rep_BS][0]=i;
            N_rep_BS++;
        }
    }  
    act_BS[0][0]=-1;
    act_BS[0][1]=0; 
    for(i=1;i<N_act_BS;i++) 
    {
        j=i-1;
        while(j!=0 && genotype->all_binding_sites[gene_id][act_BS[i][0]].BS_pos - genotype->all_binding_sites[gene_id][act_BS[j][0]].BS_pos<TF_ELEMENT_LEN+2*HIND_LENGTH)j--;
        act_BS[i][1]=act_BS[j][1]+1;
    } 
    rep_BS[0][0]=-1;
    rep_BS[0][1]=0;
    for(i=1;i<N_rep_BS;i++) 
    {
        j=i-1;
        while(j!=0 && genotype->all_binding_sites[gene_id][rep_BS[i][0]].BS_pos - genotype->all_binding_sites[gene_id][rep_BS[j][0]].BS_pos<TF_ELEMENT_LEN+2*HIND_LENGTH)j--;
        rep_BS[i][1]=rep_BS[j][1]+1;
    }
    genotype->max_unhindered_sites[gene_id][1]=act_BS[N_act_BS-1][1];
    genotype->max_unhindered_sites[gene_id][2]=rep_BS[N_rep_BS-1][1];
}


void remove_edges_iteratively(Genotype *genotype)
{
    int i, j, N_FFLs;
    int copy_TF_in_C1ffl[NGENES][NPROTEINS];
    int copy_gene_in_C1ffl[NGENES];
    
    N_FFLs=0;
    
    for(i=0;i<NGENES;i++)
    {
        copy_gene_in_C1ffl[i]=genotype->gene_in_core_C1ffl[i];
        N_FFLs+=genotype->gene_in_core_C1ffl[i];
        for(j=0;j<NPROTEINS;j++)
            copy_TF_in_C1ffl[i][j]=genotype->TF_in_core_C1ffl[i][j];
    }
    
    while(N_FFLs)
    {
        modify_topology(genotype);
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            genotype->recalc_TFBS[i]=1;
        find_ffl(genotype);
        N_FFLs=0;
        for(i=0;i<NGENES;i++)
        {
            N_FFLs+=genotype->gene_in_core_C1ffl[i];
            copy_gene_in_C1ffl[i]=(copy_gene_in_C1ffl[i]>genotype->gene_in_core_C1ffl[i])?
                                    copy_gene_in_C1ffl[i]:genotype->gene_in_core_C1ffl[i];            
            genotype->gene_in_core_C1ffl[i]=copy_gene_in_C1ffl[i];
            for(j=0;j<NPROTEINS;j++)
            {
                copy_TF_in_C1ffl[i][j]=(copy_TF_in_C1ffl[i][j]>genotype->TF_in_core_C1ffl[i][j])?
                                    copy_TF_in_C1ffl[i][j]:genotype->TF_in_core_C1ffl[i][j];
                genotype->TF_in_core_C1ffl[i][j]=copy_TF_in_C1ffl[i][j];                    
            }
        }
    }
}

#endif

void tidy_output_files(char *file_genotype_summary, char *file_mutations)
{
    int i,replay_N_steps,N_tot_mutations;
    char buffer[600];    
    FILE *fp1,*fp2;
    
    fp1=fopen("saving_point.txt","r");
    fscanf(fp1,"%d %d",&replay_N_steps,&N_tot_mutations);
    fclose(fp1);
    
    /*Basically delete the last line of the file if it is not complete*/
    fp1=fopen(file_genotype_summary,"r");
    fp2=fopen("temp","w");
    for(i=0;i<replay_N_steps+2;i++)
    {
        fgets(buffer,600,fp1);
        fputs(buffer,fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove(file_genotype_summary);
    rename("temp",file_genotype_summary);
    
    fp1=fopen(file_mutations,"r");
    fp2=fopen("temp","w");
    for(i=0;i<replay_N_steps;i++)
    {
        fgets(buffer,600,fp1);
        fputs(buffer,fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove(file_mutations);
    rename("temp",file_mutations);
    
    fp1=fopen("proportion_c1ffl.txt","r");
    fp2=fopen("temp","w");
    for(i=0;i<replay_N_steps;i++)
    {
        fgets(buffer,600,fp1);
        fputs(buffer,fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove("proportion_c1ffl.txt");
    rename("temp","proportion_c1ffl.txt");
#if OUTPUT_MUTANT_DETAILS 
    fp1=fopen("MUT_Detail.txt","r");
    fp2=fopen("temp","w");
    for(i=0;i<N_tot_mutations;i++)
    {
        fgets(buffer,600,fp1);
        fputs(buffer,fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove("MUT_Detail.txt");
    rename("temp","MUT_Detail.txt");
    
    fp1=fopen("Mut_detail_fitness.txt","r");
    fp2=fopen("temp","w");
    for(i=0;i<N_tot_mutations;i++)
    {
        fgets(buffer,600,fp1);
        fputs(buffer,fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove("Mut_detail_fitness.txt");
    rename("temp","Mut_detail_fitness.txt");
#endif
    fp1=fopen("precise_fitness.txt","r");
    fp2=fopen("temp","w");
    for(i=0;i<replay_N_steps;i++)
    {
        fgets(buffer,600,fp1);
        fputs(buffer,fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove("precise_fitness.txt");
    rename("temp","precise_fitness.txt");
}

void calc_leaping_interval(Genotype *genotype, CellState *state, float *minimal_interval, float t_unreachable, int which_gene)
{
    int protein_id,j;
    float dt;
    float N_proteins_cause_change,Kd,N_at_end_of_simulation;
    float overall_rate;
    float P_binding;
    float t_remaining;
    float N_proteins_at_now[genotype->ntfgenes];
    int gene_id; 
    float ct, ect, one_minus_ect;
    
    t_remaining=t_unreachable-state->t;
    dt=t_unreachable;  
 	protein_id=genotype->which_protein[which_gene];
    Kd=KD2APP_KD*genotype->Kd[protein_id];
    P_binding=state->protein_number[protein_id]/(state->protein_number[protein_id]+Kd);        

    /*determine whether the protein tends to increase or decrease concentration*/
    overall_rate=0.0;
    for(j=0;j<genotype->protein_pool[protein_id][0][0];j++) 
    {
        gene_id=genotype->protein_pool[protein_id][1][j];
        overall_rate+=(state->protein_synthesis_index[gene_id]-state->gene_specific_protein_number[gene_id])*genotype->protein_decay_rate[gene_id];
    }

    if(overall_rate>0.0) //tend to increase
    {            
        if(P_binding>=(1.0-MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING)) //concentration already very high
            dt=t_unreachable;
        else
        {
            /* check if much change is possible within the duration of simulation*/
            N_proteins_cause_change=Kd*(P_binding+MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING)/(1.0-P_binding-MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING); 
            /* calc N_protein at the end of simulation*/
            N_at_end_of_simulation=0.0;
            for(j=0;j<genotype->protein_pool[protein_id][0][0];j++) 
            {  
                gene_id=genotype->protein_pool[protein_id][1][j];
                ct=genotype->protein_decay_rate[gene_id]*t_remaining;
                ect = exp(-ct);
                if (fabs(ct)<EPSILON) one_minus_ect=ct;
                else one_minus_ect = 1.0-ect;   
                N_at_end_of_simulation+=ect*state->gene_specific_protein_number[gene_id]+state->protein_synthesis_index[gene_id]*one_minus_ect;        
            }               
            if(N_at_end_of_simulation<N_proteins_cause_change)
               dt=t_unreachable; 
            else /*We need to solve an equation*/
            {                    
                for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)
                {
                    gene_id=genotype->protein_pool[protein_id][1][j];
                    N_proteins_at_now[j]=state->gene_specific_protein_number[gene_id];
                }
                dt=calc_tprime(genotype,state,N_proteins_at_now,t_remaining,N_proteins_cause_change,protein_id); 
            }
        }
        *minimal_interval=(*minimal_interval<dt)?*minimal_interval:dt; 
    }
    else //tends to decrease
    {
        if(P_binding<=MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING) //concentration already very low
            dt=t_unreachable;
        else
        {
            /* first check if much change is possible within the duration of simulation*/
            N_proteins_cause_change=Kd*(P_binding-MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING)/(1.0-P_binding+MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING); 
            /* calc N_protein at the end of simulation*/
            N_at_end_of_simulation=0.0;
            for(j=0;j<genotype->protein_pool[protein_id][0][0];j++) 
            {  
                gene_id=genotype->protein_pool[protein_id][1][j];
                ct=genotype->protein_decay_rate[gene_id]*t_remaining;
                ect = exp(-ct);
                if (fabs(ct)<EPSILON) one_minus_ect=ct;
                else one_minus_ect = 1.0-ect;   
                N_at_end_of_simulation+=ect*state->gene_specific_protein_number[gene_id]+state->protein_synthesis_index[gene_id]*one_minus_ect;        
            }               
            if(N_at_end_of_simulation>N_proteins_cause_change)
               dt=t_unreachable; 
            else
            {
                for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)
                {
                    gene_id=genotype->protein_pool[protein_id][1][j];
                    N_proteins_at_now[j]=state->gene_specific_protein_number[gene_id];
                }
                dt=calc_tprime(genotype,state,N_proteins_at_now,t_remaining,N_proteins_cause_change,protein_id); 
            }
        }
        *minimal_interval=(*minimal_interval<dt)?*minimal_interval:dt;  
    }
}

void print_mutatable_parameters(Genotype *genotype)
{
    int i;
    FILE *fp;    
    fp=fopen("mutatable_parameters.txt","w");
    for(i=0;i<genotype->ngenes;i++)
    {
        fprintf(fp,"%f %f %f %f ",genotype->active_to_intermediate_rate[i],
                                    genotype->mRNA_decay_rate[i],
                                    genotype->translation_rate[i],
                                    genotype->protein_decay_rate[i]);
        if(genotype->protein_identity[genotype->which_protein[i]]!=-1) //is a tf gene
            fprintf(fp,"%f\n",log10(genotype->Kd[genotype->which_protein[i]]));
        else
             fprintf(fp,"na\n");
    }
    fclose(fp);
}

void make_meta_genotype(Genotype *genotype)
{
    int i,j,which_protein;
    Gene *gene_id;
    Protein *protein_id;
    ProteinFamily *family_id;
    
    /*make protein_family*/
    family_id=genotype->functional_group;
    for(i=0;i<genotype->N_protein_family;i++)
    {
        genotype->TF_family_pool[i][0][0]=family_id->N_members;
        protein_id
        for(j=0;j<family_id->N_members;j++)
        {
            genotype->TF_family_pool[i][1][0][j]=1;
        }
    }
    
    /**/
    
    
    
    
    
    
    
    
    
    
    /*make protein_pool, which_protein */
    for(i=0;i<NPROTEINS;i++)
    {
        genotype->protein_pool[i][0][0]=0;
        for(j=0;j<NGENES;j++)
            genotype->protein_pool[i][1][j]=NA;
    }
    for(i=0;i<NGENES;i++)
    {
        genotype->which_protein[i]=NA;            
    }
    i=0;
    which_gene=genotype->genome;
    while(i<genotype->ngenes)
    {    
        which_protein=which_gene->protein_id;
        genotype->protein_pool[which_protein][1][genotype->protein_pool[which_protein]]
    }    
}