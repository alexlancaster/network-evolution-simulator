/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-  */
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster, Jasmin Uribe
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
 */
#ifndef FILE_NETSIM_SEEN
#define FILE_NETSIM_SEEN

#include <stdio.h>
#include "RngStream.h"

/*Simulation mode*/
#define JUST_PLOTTING 0
#define PLOT_ALTERNATIVE_FITNESS 0
#define NEUTRAL 0
#define RUN_FULL_SIMULATION 1
#define SKIP_INITIAL_GENOTYPE 0
#define SET_BS_MANUALLY 0
#define QUICK_BURN_IN 0
#define EXTERNAL_SIGNAL 0

/*Runtime control*/  
#define MAX_MUT_STEP 55000
#define BURN_IN_I 5000
#define BURN_IN_II 0
#define MAX_MUTATIONS 800000
#define MAX_TRIALS 2000
#define N_THREADS 10
#define N_REPLICATES 200
#define OUTPUT_INTERVAL 10

/*Miscellaneous settings*/
#define MAXIT 100          /* maximum number of iterations for Newtown-Raphson */
#define EPSILON 1.0e-6       /* original code used EPSILON 10^-6 */
#define RT_SAFE_EPSILON 1e-6
#define OUTPUT_RNG_SEEDS 1
#define TIME_INFINITY 9.99e10
#define TIME_OFFSET 0.01
#define CAUTIOUS 0

/*Biology and evolution settings*/
#define ROUND_UP_NEGATIVE_FITNESS 0
#define DIRECT_REG 0
#define NO_PENALTY 0
#define ADD_2_PATHWAYS 0
#define FORCE_OR_GATE 0
#define RANDOM_COOPERATION_LOGIC 0
#define N_SIGNAL_TF 1 // the 1st TF enables basal activity in TFN. The 2nd is the actual signal TF. 
#define NO_REGULATION_COST 0
#define NO_REGULATION 0 // this locks the state of transcription factors to NUC_NO_PIC
//#define IGNORE_SOURCE_SIGNAL_REG_ON_EFFECTOR 1 // ignores env tf when calculating the Pact of effector WHEN PLOTTING
#define ADJUST_FITNESS 0 // allows manually adjust the fitness of a phenotype
#if ADJUST_FITNESS
#define ADJUST_FITNESS_CONDITION 3 // more or less than 3 copies of regular transcription factor genes receives penaly in fitness
#define ADJUSTMENT 1.0e-3 //selection coefficient
#endif
#define REGRESSIVE_MUTATION 1
#define IGNORE_BS_OVERLAPPING 0
#define SIMPLE_SUBSTITUTION 1
#define RANDOMIZE_SIGNAL2 0
#define ALPHA 0.2
#define MAX_RECALC_FITNESS 10
#define RULE_OF_REPLACEMENT 4 /* 0 for z-score, 1 for Wilcoxon, 2 for larger-fitness-fixes, 3 for larger-than-epsilon-fixes, 4 for s>minimal_selection_coefficient */
#if RULE_OF_REPLACEMENT==4
#define minimal_selection_coefficient 1.0e-8
#endif
#ifndef MAX_COPIES
#define MAX_COPIES 2       /* each gene can have at most two copies*/
#endif
/* Because mutation can change the number of genes, the numbers defined here are used to allocate storage space only.
 * Set the numbers to be 8 folds of the initial ngenes and ntfgenes, so that we can have two whole genome duplications*/
#ifndef TFGENES             /* number of genes encoding TFs */
#define TFGENES 20         /* the initial value is set in initiate_genotype*/
#endif
#ifndef NGENES
#define NGENES 25  /* total number of genes: add the (non-TF) selection gene to the total (default case) */
#endif
#ifndef NPROTEINS           
#define NPROTEINS 25
#endif
#ifndef EFFECTOR_GENES
#define EFFECTOR_GENES 5  /* this is the upper limit of effector gene copies*/
#endif
#define CISREG_LEN 150        /* length of cis-regulatory region in base-pairs */
#define TF_ELEMENT_LEN 8      /* length of binding element on TF */
#define NMIN 6                /* minimal number of nucleotide that matches the binding sequence of a TF in its binding site*/
                              /* DO NOT MAKE NMIN<TF_ELEMENT_LEN/2, OTHERWISE calc_all_binding_sites_copy will make mistake*/  
//#define NUM_K_DISASSEMBLY 131 /* number of differents for PIC disassembly from data file  */
#ifndef HIND_LENGTH
#define HIND_LENGTH 3         /* default length of hindrance on each side of the binding site,i.e. a tf occupies TF_ELEMENT_LEN+2*HIND_LENGTH */
                              /* the binding of Lac repressor blockes 12 bp. Record MT 1981*/
#endif
#define MAX_BINDING 10  /* MAX_MODE is the max number of tf that can bind to a promoter plus 1*/
//#define MAX_BS_IN_CLUSTER 100


/* 
 * define macros for logging warning/errors 
 */
#ifndef LOGGING_OFF
  #define LOG(...) { fprintf(fperrors, "%s: ", __func__); fprintf (fperrors, __VA_ARGS__) ; fflush(fperrors); } 
#else
  #define LOG 
#endif

/*
 * primary data structures for model
 */
enum { KON_PROTEIN_DECAY_INDEX = 0, KON_SALPHAC_INDEX = 1 };

/*
 * enum for 'CellState'->active indices
 */

enum { NUC_NO_PIC = 0,
       NO_NUC_NO_PIC =1,
       PIC_NO_NUC = 3,};

/*
 * Rates for Gillespie algorithm
 *
 * Events with exponentially-distributed waiting times are:
 * - transcription initiation
 * - mRNA transport and decay
 * - acetylation and deacetylation, 
 * - and PIC assembly and disassembly.
 *
 */
typedef struct GillespieRates GillespieRates;
struct GillespieRates {
  float mRNAdecay;         /* rates[2] */
  float mRNAdecay_rate[NGENES];
  float pic_disassembly;    /* rates[3] */
  float pic_disassembly_rate[NGENES]; 
  float nuc_disassembly;
  float nuc_disassembly_rate[NGENES]; 
  float nuc_assembly;
  float nuc_assembly_rate[NGENES];  
  float pic_assembly;
  float pic_assembly_rate[NGENES];  
  int transcript_init; 
  int transcript_init_rate[NGENES];   
  float subtotal;
};

typedef struct AllTFBindingSites AllTFBindingSites;
struct AllTFBindingSites {
  int tf_id;         /* transcription factor */
//  float Koff;        /* replacing hamming_dist */
  float Kd;
  int mis_match;     /* number of mismatched nuc */
  int BS_pos;        /* start position of BS on DNA, always with reference to forward strand */                     
  int N_hindered;    /* the number of BS hindered by this TF when it binds to the current BS */  
};

typedef struct Genotype Genotype;
struct Genotype {

    int ngenes;                                             /* the number of loci */
    int ntfgenes;                                           /* the number of tf loci */
    int nproteins;                                          /* because of gene duplication, the number of proteins and mRNAs can be different 
                                                               from the number of loci. nprotein is the number of elements in protein_pool*/
    int which_protein[NGENES];                              /* in case of gene duplication, this array tells the protein corresponding to a given gene id */   
    char cisreg_seq[NGENES][CISREG_LEN];    
    
    /*these apply to protein, not loci*/
    int N_act;                                              /* number of activators*/
    int N_rep;                                              /* number of repressors*/    
    int activating[NPROTEINS];                              /* 1 for activator, 0 for repressor, -1 for non-tf */ 
    char tf_seq[NPROTEINS][TF_ELEMENT_LEN];
    char tf_seq_rc[NPROTEINS][TF_ELEMENT_LEN];                /* reversed complementary sequence of BS. Used to identify BS on the non-template strand*/
    int protein_pool[NPROTEINS][2][NGENES];                 /* element 1 record how many genes/mRNAs producing this protein,ele 2 stores which genes/mRNAs*/
    float Kd[NPROTEINS];
    
    /*these apply to loci*/
    float mRNAdecay[NGENES];                                /* kinetic rates*/
    float proteindecay[NGENES];                             /* kinetic rates*/
    float translation[NGENES];                              /* kinetic rates*/   
    float pic_disassembly[NGENES];                          /* kinetic rates*/    
    int min_act_to_transc[NGENES];                          /* 1 for OR GATE, at leat 2 FOR AND GATE */ 
 
    /* binding sites related data, applying to loci*/   
    int cisreg_cluster[NGENES+1][NGENES];                     /* For genes having the same cis-reg, tf distribution can be shared.
                                                             * Genes having the same cis-reg are clustered.
                                                               1st dim stores cluster ids, 2nd dim stores gene_ids in a cluster.
                                                               cisreg_cluster works with which_cluster*/
    int which_cluster[NGENES];                              /* which_cluster stores the cluster id of a gene*/                                                           
    int recalc_TFBS[NGENES];                                /* whether to recalc the TFBS*/
    
    int binding_sites_num[NGENES];                          /* total number of binding sites */
    int max_unhindered_sites[NGENES][3];                    /* maximal number of binding sites that do not hinder each other. 0 for activator BS, 1 for repressor BS*/  
    int max_hindered_sites[NGENES];                        /* maximal number of BSs a BS can hinder*/ 

    int N_act_BS[NGENES];                                   /* total number of binding sites of activating TF */
    int N_rep_BS[NGENES];                                   /* total number of binding sites of repressing TF */
    int avg_N_BS_in_cluster[NGENES][NGENES][2];

    AllTFBindingSites *all_binding_sites[NGENES];   


    float avg_GR1;
    float avg_GR2;
    float sq_SE_GR1;
    float sq_SE_GR2;
    float fitness;
    float sq_SE_fitness;
    float fitness_measurement[MAX_RECALC_FITNESS*N_REPLICATES];
    int N_motifs[27];    
//     float proportion_motifs[18]; 
    int TF_in_core_C1ffl[NGENES][NPROTEINS];
    int gene_in_core_C1ffl[NGENES];
    int N_act_genes; 
    int N_act_genes_reg_by_env;
    int N_act_genes_not_reg_by_env;
};

/* 
 * transcription/translation delays are sorted linked lists.  Deleting
 * the head each time, and tack new stuff on the end.  Linked lists
 * are easy to create pre-sorted.
 */
typedef struct FixedEvent FixedEvent;
struct FixedEvent {
  int event_id; 
  float time;
  FixedEvent *next;
};

typedef struct CellState CellState;
struct CellState {   
    float cumulative_fitness;                    /* size of cell */
    float cumulative_fitness_after_burn_in;          
    float instantaneous_fitness;                  /* total growth rate in the previous deltat */   
    int mRNA_aft_transl_delay_num[NGENES];          /* mRNAs that have finished the translational delay */

    int mRNA_under_transl_delay_num[NGENES];   /* mRNAs that are still under the translational delay (they do not contribute to protein 
                                                * turnover, but contribute to the cost of translation)  */
    int mRNA_under_transc_num[NGENES];       /* mRNAs which haven't finished transcription */

    FixedEvent *mRNA_transl_init_time_end_head;   /* times when mRNAs become fully loaded with ribosomes and start producing protein */
    FixedEvent *mRNA_transl_init_time_end_tail;  
    FixedEvent *mRNA_transcr_time_end_head;  /* times when transcription is complete and an mRNA is available to move to cytoplasm */
    FixedEvent *mRNA_transcr_time_end_tail;
    FixedEvent *signal_off_head;          /* times when env=A ends. Note, this event is not gene- or copy-specific. I just use the structure of FixedEvent for convenience.*/
    FixedEvent *signal_off_tail;   
    FixedEvent *signal_on_head;          /* times when env=A ends. Note, this event is not gene- or copy-specific. I just use the structure of FixedEvent for convenience.*/
    FixedEvent *signal_on_tail; 
    FixedEvent *burn_in_growth_rate_head;
    FixedEvent *burn_in_growth_rate_tail;  
    FixedEvent *sampling_point_end_head;
    FixedEvent *sampling_point_end_tail;
    FixedEvent *change_signal_strength_head;
    FixedEvent *change_signal_strength_tail;

    float t_to_update_Pact_or_Prep;

    float Pact[NGENES];
    float Prep[NGENES];
    float Pact_No_rep[NGENES];
    float Pno_TF[NGENES];
    float last_Prep[NGENES];
    float last_Pact[NGENES];
    float last_Pact_No_rep[NGENES];
    float last_Pno_TF[NGENES];
    float last_event_t;

    float protein_number[NPROTEINS];     /* pooled protein number from gene_specific_protein_conc */
    float gene_specific_protein_number[NGENES]; /* stores the "protein" number for each gene.
                                               * can be considered temporary data. Make muation easier to
                                               * deal with. */
    float konvalues[NGENES][2];        /* moved from KonState*/  
    int active[NGENES];       /* gives the state of each of the genes, according to figure
                                          *  NUC_NO_PIC = 0,
                                          *  NO_NUC_NO_PIC =1,
                                          *  PIC_NO_NUC = 3,
                                          */ 
};

typedef struct TimeCourse TimeCourse;
struct TimeCourse
{
  float concentration;
  float time;
  TimeCourse *next;
};    

typedef struct Mutation Mutation;
struct Mutation
{
    char mut_type;
    int which_gene;
    int which_nucleotide;
    char nuc_diff[3];    /*the first three elements store the nuc after mutation (max_indel=3)*/
    int kinetic_type;    /*0 for pic_disassembly, 1 for mRNA_decay, 2 for translation, 3 for protein_decay*/
    float kinetic_diff;
};

/*
 * global variables
 */

/* see netsim.c for documentation for these global constant variables */
int MAXELEMENTS;
const float TRANSLATION_TIME;
const float TRANSCRIPTION_TIME;
const float TRANSCRIPTINIT; 
const float NUC_DISASSEMBLY;
const float NUC_ASSEMBLY;
const float BASE_NUC_ASSEMLY;
const float BASE_NUC_DISASSEMBLY;
const float PIC_ASSEMBLY;
const float BASE_PIC_ASSEMBLY;

/*output files*/
char error_file[32];
char mutation_file[32];
char RuntimeSumm[32];
char output_file[32];

/* see netsim.c for documentation for these global variables */
float tdevelopment; 
float growth_rate_scaling; 
float Ne_saturate;
float c_transl;
float bmax;
int init_TF_genes;
float signal_off_strength;
float background_signal_strength;
float env1_signal_strength;
float env2_signal_strength;
float env1_t_signal_on;    
float env1_t_signal_off;     
float env2_t_signal_on;
float env2_t_signal_off;
char env1_initial_effect_of_effector;
char env2_initial_effect_of_effector;
int env1_fixed_effector_effect;    
int env2_fixed_effector_effect;
int env1_signalA_mismatches; 
int env2_signalA_mismatches;
float env1_occurence;
float env2_occurence;
int init_N_act;
int init_N_rep;
int recalc_new_fitness;
char init_env1;
char init_env2;
int min_act_to_transcr_selection_protein;
float cost_term;
float penalty;
int N_replicates;
float duration_of_burn_in_growth_rate;
/* file output parameters */
char *output_directory ;
int verbose ;
FILE *fperrors;
FILE *fp_cellsize[2];
#if 0 
FILE *fp_koff[2];
FILE *fp_growthrate[2];
FILE *fp_tfsbound[2];
FILE *fp_rounding[2];
#endif

/* function prototypes */

extern char set_base_pair(float);

extern void initialize_growth_rate_parameters();

extern void initialize_sequence(char [], int, int, RngStream);

extern void print_genotype(Genotype *, int);

extern void print_all_binding_sites(int [NGENES],
                                    AllTFBindingSites *, 
                                    int ,
                                    char [TFGENES][TF_ELEMENT_LEN],
                                    char [NGENES][CISREG_LEN]
//                                    int [NGENES][2]
									);

extern void print_tf_occupancy(CellState *,
                               AllTFBindingSites *,
                               float);

extern void initialize_genotype(Genotype *, RngStream) ;

extern void initialize_genotype_fixed(Genotype *, RngStream);

extern void calc_all_binding_sites_copy(Genotype *, int);

extern void calc_all_binding_sites(Genotype *);

extern void set_binding_sites(Genotype *,int *, int (*)[5], int (*)[5], float (*)[5]);

extern void cluster_BS_cluster(Genotype *, int);

extern double calc_flux(AllTFBindingSites *,int,int,double [MAX_BINDING]);

extern float calc_ratio_act_to_rep(AllTFBindingSites *,
                                    int ,
                                    int ,
                                    int ,
                                    int ,
                                    int , 
                                    int [NGENES],                                    
                                    int ,
                                    int ,
                                    int *,
                                    float [NGENES]);

extern void calc_TF_dist_from_all_BS(  AllTFBindingSites *,
                                        int ,
                                        int ,
                                        int ,                                         
                                        int [NPROTEINS],                                    
                                        int [3],
                                        float [NGENES],
                                        int,
                                        float *,
                                        float *,
					float *,
                                        float *);

extern int add_fixed_event(int,                           
                           float,
                           FixedEvent **,
                           FixedEvent **);

extern void add_time_point(float,
                           float,
                           TimeCourse **,
                           TimeCourse **);

extern void add_fixed_event_end(int,                                
                                float,
                                FixedEvent **,
                                FixedEvent **);

extern void delete_fixed_event(int,                               
                               int,
                               FixedEvent **,
                               FixedEvent **);

extern void delete_fixed_event_start(FixedEvent **,
                                     FixedEvent **);

extern void initialize_cell(CellState *,
                            int,
                            int,  
                            int [NPROTEINS][2][NGENES],
                            float [NGENES],
                            int [NGENES],
                            float [NPROTEINS],
                            RngStream,float);

extern int does_fixed_event_end(CellState*,
                                float);

extern int does_fixed_event_end_plotting(CellState *,float);

extern void update_protein_number_and_fitness(Genotype *,
                                            CellState *,
                                            GillespieRates *,
                                            float,                                        
                                            float,
                                            char,
                                            int *,                                           
                                            int,
                                            Mutation *);

extern void fixed_event_end_transcription(  float *,
                                            float,
                                            CellState *,
                                            GillespieRates *,
                                            Genotype *,
                                            int *,                              
                                            int,
                                            Mutation *,
                                            char *);
extern float calc_tprime(Genotype*, CellState*, float*, float, float);

extern float calc_integral(Genotype *, CellState *, float *, float, float);

extern float calc_fitness(  float *,
                            Genotype *,
                            CellState *,
                            float*,                                       
                            float,                                       
                            char,
                            int *,                                        
                            int,
                            Mutation *);

extern void transport_event(Genotype *,
                            CellState *,
                            GillespieRates *,
                            float,
                            float,
                            RngStream);

extern int Gillespie_event_mRNA_decay(GillespieRates *, CellState *, Genotype *, float, RngStream);

extern void Gillespie_event_histone_acteylation(GillespieRates *, CellState *, Genotype *, RngStream);

extern void Gillespie_event_histone_deacteylation(GillespieRates *, CellState *, Genotype *, RngStream); 

extern void Gillespie_event_assemble_PIC(GillespieRates *, CellState *, Genotype *, RngStream);

extern void Gillespie_event_disassemble_PIC(Genotype *, CellState *, GillespieRates *, RngStream);

extern void Gillespie_event_transcription_init(GillespieRates *, CellState *, Genotype *, float, float, RngStream);

extern void do_single_timestep( Genotype *, 
                                CellState *,
                                GillespieRates *,                              
                                float *,                                                  
                                char *,
                                float,
                                float,
                                float, 
                                char,
                                int,
                                RngStream,                             
                                int,
                                Mutation *,
                                int *,
                                int *,
                                int,
                                float,
                                float *) ;

extern void do_single_timestep_plotting(    Genotype *, 
                                            CellState *,
                                            GillespieRates *,                              
                                            float *,                                                  
                                            char *,                                            
                                            float,
                                            float,
                                            int,   
                                            char,
                                            float (*)[90],
                                            float [90],
                                            RngStream,                              
                                            int *,
                                            int *,
                                            int,
                                            float,
                                            float *) ;
							   
extern void free_fixedevent(CellState *);
 
extern void calc_avg_growth_rate(   Genotype *,                                    
                                    int [NGENES],
                                    float [NPROTEINS],                                  
                                    RngStream [N_THREADS],                                                                     
                                    int,
                                    float *,
                                    float *,
                                    Mutation *); 
  
extern int init_run_pop(unsigned long int [6], int);

extern void print_time_course(TimeCourse *, FILE *);

extern void print_all_protein_time_courses(int, TimeCourse **, TimeCourse **, FILE *);

extern void log_snapshot(Genotype *,
                         CellState *,
                         GillespieRates *,
                         float ,
                         float );

extern int mod(int, int); 

extern void calc_all_rates(Genotype *,
                            CellState *,
                            GillespieRates *, 
                            float,
                            int,int);

extern int fixed_event_end_translation_init(   Genotype *, 
                                                CellState *,    
                                                GillespieRates *, 
                                                float *, 
                                                float,
                                                int *,                                    
                                                int,
                                                Mutation *,
                                                char *);

extern int do_fixed_event(  Genotype *, 
                            CellState *, 
                            GillespieRates *, 
                            float *,
                            float ,
                            int , 
                            float,
                            float ,
                            float,  
                            char *,
                            char,
                            int,                         
                            int *,                          
                            int,
                            Mutation *,
                            float *);

extern int do_fixed_event_plotting( Genotype *, 
                                    CellState *, 
                                    GillespieRates *,                                     
                                    float *,
                                    float ,
                                    int , 
                                    float ,
                                    float,        
                                    char *,  
                                    char ,
                                    int ,                                    
                                    int *,
                                    float *);

extern int do_Gillespie_event(Genotype*, CellState *, GillespieRates *, float, float, RngStream, int *, int, Mutation *);

extern void calc_configurations(Genotype *, int);

extern void mutate(Genotype *, RngStream, Mutation *);

extern void mut_substitution(Genotype *, Mutation *, RngStream);

extern void mut_insertion(Genotype *,Mutation *, RngStream);

extern void mut_partial_deletion(Genotype *,Mutation *, RngStream);

extern void mut_whole_gene_deletion(Genotype *,Mutation *, RngStream);

extern void mut_duplication(Genotype *,Mutation *, RngStream);

extern void mut_binding_sequence(Genotype *,Mutation *, RngStream);

extern void mut_kinetic_constant(Genotype *, Mutation *, RngStream);

extern float mut_make_new_value(float, float, RngStream);

extern void mut_identity(Genotype *, Mutation *, RngStream);

extern void mut_koff(Genotype *, Mutation *, RngStream);

extern void reproduce_mutate(Genotype *, Mutation *,RngStream);

extern void reproduce_substitution(Genotype *, Mutation *);

extern void reproduce_insertion(Genotype *,Mutation *);

extern void reproduce_partial_deletion(Genotype *,Mutation *);

extern void reproduce_whole_gene_deletion(Genotype *,Mutation *);

extern void reproduce_gene_duplication(Genotype *,Mutation *);

extern void reproduce_mut_binding_sequence(Genotype *,Mutation *);

extern void reproduce_mut_kinetic_constant(Genotype *, Mutation *);

extern void reproduce_mut_identity(Genotype *, Mutation *);

extern void reproduce_mut_koff(Genotype *, Mutation *);

extern void draw_mutation(Genotype *, char *, RngStream);

extern void initialize_cache(Genotype *);

extern void update_protein_pool(Genotype *, int, int, char);

extern void update_cisreg_cluster(Genotype *, int, char, int [NGENES][NGENES], int, int);

extern void clone_genotype(Genotype *, Genotype *);

extern float try_fixation(Genotype *, Genotype *, int, int, int *, RngStream);

extern void calc_avg_growth_rate_plotting(  Genotype *,                               
                                            int [NGENES],
                                            float [NPROTEINS], 
                                            RngStream [N_THREADS]); 

extern void summarize_binding_sites(Genotype *,int);

extern void print_binding_sites_distribution(Genotype *,int, int);

extern int check_concurrence(   float , 
                                FixedEvent *, 
                                FixedEvent *, 
                                FixedEvent *, 
                                FixedEvent *,
                                FixedEvent *,
                                float,
                                FixedEvent *);

extern void set_signal(CellState *, float, float, float *, float,float);

extern void output_genotype(Genotype *, int);

extern void release_memory(Genotype *,Genotype *, RngStream *, RngStream [N_THREADS]);

extern void calc_fx_dfx(float, int, float, float*, float*, float*, float*, float*);

extern void resolve_overlapping_sites(Genotype *, int, int [NGENES]);

extern int evolve_N_steps(  Genotype *, 
                            Genotype *,
                            int *, 
                            int, 
                            int *,   
                            float [NPROTEINS],
                            int [NGENES],
                            Mutation *, 
                            RngStream,
                            RngStream [N_THREADS],
                            int );

extern void run_simulation( Genotype *, 
                            Genotype *,                     
                            float [NPROTEINS],
                            int [NGENES],
                            Mutation *,
                            int,
                            int,
                            RngStream,
                            RngStream [N_THREADS]);
extern void continue_simulation(Genotype *, 
                                Genotype *,                                
                                int, 
                                float [NPROTEINS],
                                int [NGENES],
                                Mutation *, 
                                RngStream,
                                RngStream [N_THREADS]);
extern void run_plotting(   Genotype *,
                            Genotype *,
                            int [NGENES],
                            float [NGENES],
                            RngStream [N_THREADS],
                            Mutation *,
                            FILE *,
                            int);

extern void calc_fitness_stats( Genotype *,
                                float (*)[N_REPLICATES],
                                float (*)[N_REPLICATES],
                                int);


extern void evolve_neutrally(   Genotype *,
                                Genotype *,                             
                                Mutation *,
                                int,
//                                float [NUM_K_DISASSEMBLY],
                                RngStream);


extern void replay_mutations(   Genotype *,
                                Genotype *,
                                FILE *,
                                Mutation *,
                                int,
                                RngStream);

extern void plot_alternative_fitness(   Genotype *,
                                        Genotype *,
                                        int [NGENES],
                                        float [NGENES],
                                        RngStream [N_THREADS],
                                        Mutation *,
                                        FILE *,
                                        int);

extern void find_ffl(Genotype *);

extern void tidy_output_files(char*, char*);

extern void print_core_c1ffls(Genotype *);

extern void initialization_add_regulation(Genotype *, RngStream);
#endif /* !FILE_NETSIM_SEEN */