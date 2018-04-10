/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RngStream.h"
#include "netsim.h"
#include "numerical.h"
#include "mutation.h"



/*
 ************ begin of mutation functions **************
 */
/*single nucleic acid substitution in cis-reg*/
void mut_substitution(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int which_nucleotide, which_gene,i,N_BS_bf_mutation,flag_difference;
    char *Genome,nucleotide; 
    AllTFBindingSites *container; // used to store binding sites before substitution
                                  // we compare whether substituion changes any binding sites,
                                  // in order to determine whether mutation creates a unique cis-reg     
                                  // sequence whose binding configuration always needs computation     
    /*points to the current cis-reg*/
    Genome= &genotype->cisreg_seq[0][0]; 
    /*whether to simulate bias in the frequency of different substitutions*/
    #if !SIMPLE_SUBSTITUTION 
        float random;
        while(1)
        {
            /*this is the nucleic acid that's going to be mutated*/ 
            which_nucleotide=RngStream_RandInt(RS,N_SIGNAL_TF*CISREG_LEN,genotype->ngenes*CISREG_LEN-1);
            random=RngStream_RandU01(RS);
            if (Genome[which_nucleotide]=='a'||Genome[which_nucleotide]=='t')
            {
                 if(random<=0.2)break;  // 20% of mutations hit A or T
            }	            
            else
            { 
                 if(random>0.2)break; //80% hit C or G
            }
        }  
        /*this is gene hit by the mutation*/
        which_gene=which_nucleotide/CISREG_LEN;
        /*calculate and store the distribution of binding sites before mutation*/
        calc_all_binding_sites_copy(genotype,which_gene);
        container = malloc(genotype->N_allocated_elements*sizeof(AllTFBindingSites));
        N_BS_bf_mutation=genotype->binding_sites_num[which_gene];
        for(i=0;i<N_BS_bf_mutation;i++)
        {
            container[i].BS_pos=genotype->all_binding_sites[which_gene][i].BS_pos;
            container[i].tf_id=genotype->all_binding_sites[which_gene][i].tf_id;
            container[i].mis_match=genotype->all_binding_sites[which_gene][i].mis_match;
        }
        /*generate new nucleic acid*/
        random=RngStream_RandU01(RS);
        switch (Genome[which_nucleotide])
        {
            case 'a':
                if(random<=0.2)
                    Genome[which_nucleotide]='c'; // 20% of mutation to A changes it to C
                else
                {
                    if((random-0.2)<=0.3)
                        Genome[which_nucleotide]='t'; // 30% changes it to T
                    else
                        Genome[which_nucleotide]='g'; // 50% changes it to G
                }
                break;
            case 't':
                if(random<=0.2)
                    Genome[which_nucleotide]='g';
                else
                {
                    if((random-0.2)<=0.3)
                        Genome[which_nucleotide]='a';
                    else
                        Genome[which_nucleotide]='c';
                }
                break;        
            case 'c':
                if(random<=0.2)
                    Genome[which_nucleotide]='g';
                else
                {
                    if((random-0.2)<=0.3)
                        Genome[which_nucleotide]='t';
                    else
                        Genome[which_nucleotide]='a';
                }
                break;           
            case 'g':  
                if(random<=0.2)
                    Genome[which_nucleotide]='c';
                else
                {
                    if((random-0.2)<=0.3)
                        Genome[which_nucleotide]='a';
                    else
                        Genome[which_nucleotide]='t';
                }
                break;
        }   
    #else // if do not simulate bias in substitution
        /*this is the nucleic acid that's going to be mutated*/        
        which_nucleotide=RngStream_RandInt(RS,N_SIGNAL_TF*CISREG_LEN,genotype->ngenes*CISREG_LEN-1);
        /*this is gene hit by the mutation*/
        which_gene=which_nucleotide/CISREG_LEN;
        /*calculate and store the distribution of binding sites before mutation*/
        calc_all_binding_sites_copy(genotype,which_gene);
        container = malloc(genotype->N_allocated_elements*sizeof(AllTFBindingSites));
        N_BS_bf_mutation=genotype->binding_sites_num[which_gene];
        for(i=0;i<N_BS_bf_mutation;i++)
        {
            container[i].BS_pos=genotype->all_binding_sites[which_gene][i].BS_pos;
            container[i].tf_id=genotype->all_binding_sites[which_gene][i].tf_id;
            container[i].mis_match=genotype->all_binding_sites[which_gene][i].mis_match;
        }
        /*generate new nucleic acid*/
        nucleotide=set_base_pair(RngStream_RandU01(RS));
        while(nucleotide==Genome[which_nucleotide]) // make sure we get a different nucleic acid
            nucleotide=set_base_pair(RngStream_RandU01(RS));
        Genome[which_nucleotide]=nucleotide;
    #endif    
    /*compare binding site bf and aft substitution to decide whether to update cisreg_cluster*/
    calc_all_binding_sites_copy(genotype,which_gene);    
    if(N_BS_bf_mutation!=genotype->binding_sites_num[which_gene])
    {
        update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1);  
    }
    else
    {
        /*assuming no change in binding sites*/
        flag_difference=0;
        /*comparing binding sites pair by pair*/
        for(i=0;i<N_BS_bf_mutation;i++)
        {
            if(container[i].BS_pos!=genotype->all_binding_sites[which_gene][i].BS_pos ||
                container[i].tf_id!=genotype->all_binding_sites[which_gene][i].tf_id ||
                container[i].mis_match!=genotype->all_binding_sites[which_gene][i].mis_match)
            {    
                flag_difference=1;
                break;
            }
        }
        if(flag_difference==1)
           update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1); 
    }
    free(container);
    /*record mutation info*/
    mut_record->which_nucleotide=which_nucleotide;
    mut_record->which_gene=which_gene;
    mut_record->nuc_diff[0]=Genome[which_nucleotide];
    genotype->recalc_TFBS[which_gene]=YES;  /*the binding sites on this cis-reg needs to be updated*/    
}

/*
 * Replay substitution. Used when replay the evolution for analysis
 */
void reproduce_substitution(Genotype *genotype, Mutation *mut_record)
{    
    char *Genome;
    int i, N_BS_bf_mutation,which_gene,flag_difference;    
    AllTFBindingSites *container;   
    /*get the mutated gene from record*/
    which_gene=mut_record->which_gene;
    /*calculate and store the distribution of binding sites before mutation*/
    calc_all_binding_sites_copy(genotype,which_gene);
    container = malloc(genotype->N_allocated_elements*sizeof(AllTFBindingSites));
    N_BS_bf_mutation=genotype->binding_sites_num[which_gene];
    for(i=0;i<N_BS_bf_mutation;i++)
    {
        container[i].BS_pos=genotype->all_binding_sites[which_gene][i].BS_pos;
        container[i].tf_id=genotype->all_binding_sites[which_gene][i].tf_id;
        container[i].mis_match=genotype->all_binding_sites[which_gene][i].mis_match;
    }
    /*apply the mutation from record*/
    Genome= &genotype->cisreg_seq[0][0];    
    Genome[mut_record->which_nucleotide]=mut_record->nuc_diff[1]; 
    /*compare binding site bf and aft substitution to decide whether to update cisreg_cluster*/
    calc_all_binding_sites_copy(genotype,which_gene);    
    if(N_BS_bf_mutation!=genotype->binding_sites_num[which_gene])
    {
        update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1);  
    }
    else
    {
        flag_difference=0;
        for(i=0;i<N_BS_bf_mutation;i++)
        {
            if(container[i].BS_pos!=genotype->all_binding_sites[which_gene][i].BS_pos ||
                container[i].tf_id!=genotype->all_binding_sites[which_gene][i].tf_id ||
                container[i].mis_match!=genotype->all_binding_sites[which_gene][i].mis_match)
            {    
                flag_difference=1;
                break;
            }
        }
        if(flag_difference==1)
           update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1); 
    }
    free(container);    
    genotype->recalc_TFBS[which_gene]=YES;  /*recalc TFBS*/
}

/**
 *Deletion whole cis-reg sequence
 */
void mut_whole_gene_deletion(Genotype *genotype, Mutation *mut_record, RngStream RS) // any gene can be deleted
{
    int i,j;
    int which_gene, offset, cisreg_seq_cluster_id,cisreg_seq_cluster_id_copy;
    char *temp;
    int protein_id, N_gene_of_selection_protein, puppet_gene,position_of_puppet_gene;    
    /* check which genes can be deleted*/   
    N_gene_of_selection_protein=genotype->protein_pool[genotype->nproteins-1][0][0]; //does the effector gene have multipe copies   
    if(genotype->ntfgenes-N_SIGNAL_TF>1)//if there are more than one copies of genes of non-sensor TF  
    {
        if(N_gene_of_selection_protein==1)//if the effector have only one copies
        {
            protein_id=genotype->nproteins-1;
            while(protein_id==genotype->nproteins-1) //the last copy of the effector cannot be deleted. choose a different gene
            {   
                which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
                protein_id=genotype->which_protein[which_gene];
            }
        }
        else //otherwise any gene can be deleted
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
    }
    else //if there's only one copy of non-sensor TF gene
    {        
        while(1)
        {
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
            protein_id=genotype->which_protein[which_gene];
            if(protein_id==genotype->nproteins-1)//obviously we can only delete a copy of the effector gene
                break;                           //it's impossible for deletion when the non-sensor TF and the effector each has only one copy
                                                 //function draw_mutation would have set deletion rate to 0 in such a case
        }
    }
    /*record mutation info*/
    mut_record->which_gene=which_gene;    
    /*To preserve data structure, the last gene has to be an effector gene. If this gene is selected for deletion, 
     * we overwrite it with one of its copy. Then we deleted the original copy. 
     * For example, assuming ngenes=10, gene 9 is an effector gene and is marked for deletion. 
     * We know gene 5 was duplicated from gene 9, so we overwrite gene 9 with all the info of gene 5, 
     * then we delete gene 5 (the puppet gene).*/
    if(which_gene==genotype->ngenes-1)
    {
        cisreg_seq_cluster_id=genotype->which_cluster[which_gene];
        /*check if the effector gene has a copy that has the same cis-reg sequence*/
        if(genotype->cisreg_cluster[cisreg_seq_cluster_id][1]!=NA)//there is at least one copy that meets the requirement
        {            
            puppet_gene=genotype->cisreg_cluster[cisreg_seq_cluster_id][0];//genes in cisreg_cluster is ordered acendingly
        }
        else// otherwise, this effector gene has a unique cis-reg sequence. We have to use a copy that has a different cis-reg sequence.
            // We also need to remove the cis-reg cluster of the effector gene to be deleted. 
        {     
            /*find such a copy -- call it puppet_gene*/
            protein_id=genotype->which_protein[which_gene];
            position_of_puppet_gene=0;
            if(genotype->protein_pool[protein_id][1][position_of_puppet_gene]==which_gene)position_of_puppet_gene++;
            puppet_gene=genotype->protein_pool[protein_id][1][position_of_puppet_gene];
            /*Add the effector gene to the cis_reg cluster of the puppet_gene, so that we can edit the cluster with standard procedure later*/
            cisreg_seq_cluster_id_copy=genotype->which_cluster[puppet_gene];            
            j=0;
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]!=NA)j++;//find an empty slot in cisreg_cluster
            genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]=which_gene;            
            genotype->which_cluster[which_gene]=cisreg_seq_cluster_id_copy;
            /*then delete which_gene's original cisreg cluster by shifting cisreg_cluster*/
            cisreg_seq_cluster_id_copy=cisreg_seq_cluster_id;            
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][0]!=NA)
            {
                /*reset cluster=cisreg_seq_cluster_id*/
                j=0;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]!=NA)
                {
                    genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]=NA;
                    j++;
                }
                /* copy cisreg_seq_cluster_id+1 to cisreg_seq_cluster_id, and update which_cluster*/                    
                j=0;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][j]!=NA)
                {                        
                    genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]=genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][j];
                    genotype->which_cluster[genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][j]]--;
                    j++;
                }
                cisreg_seq_cluster_id_copy++;
            }
        }       
        /*overwrite which_gene with the info of puppet_gene*/
        for(j=0;j<CISREG_LEN;j++)        
            genotype->cisreg_seq[which_gene][j]=genotype->cisreg_seq[puppet_gene][j];
        genotype->min_N_activator_to_transc[which_gene]=genotype->min_N_activator_to_transc[puppet_gene]; 
        genotype->active_to_intermediate_rate[which_gene]=genotype->active_to_intermediate_rate[puppet_gene];           
        genotype->mRNA_decay_rate[which_gene]=genotype->mRNA_decay_rate[puppet_gene];
        genotype->protein_decay_rate[which_gene]=genotype->protein_decay_rate[puppet_gene];
        genotype->translation_rate[which_gene]=genotype->translation_rate[puppet_gene]; 
        genotype->locus_length[which_gene]=genotype->locus_length[puppet_gene];
        genotype->recalc_TFBS[which_gene]=YES;             
        /*now the effector gene is same as the puppet_gene, so delete the puppet_gene*/
        which_gene=puppet_gene;        
    }
    /*points the first nucleotide of the cis-reg sequence to be deleted*/
    temp = &genotype->cisreg_seq[which_gene][0];			
    offset=CISREG_LEN;
    /* shift the cisreg array to overwrite the cis_seq to be deleted */
    for(i=0;i<CISREG_LEN*(genotype->ngenes-which_gene-1);i++) 
    {				
        *temp=*(temp+offset);         
        temp++;				
    }    
    /* if the to-be-deleted gene is a tf gene*/
    protein_id=genotype->which_protein[which_gene];   
    if(protein_id<genotype->nproteins-1)
    {  
        /*if this tf has only one copy of gene, then we'll delete the binding seq of this tf*/ 
        if(genotype->protein_pool[protein_id][0][0]==1) 
        {            
            /* shift the tf_reg array to overwrite the binding sequence to be deleted */
            temp=&genotype->tf_seq[protein_id][0];
            offset=TF_ELEMENT_LEN;
            for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
            {
                *temp=*(temp+offset);
                temp++;
            }
            /* shift the tf_reg_rc array to overwrite the binding sequence to be deleted */    
            temp=&genotype->tf_seq_rc[protein_id][0];
            for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
            {
                *temp=*(temp+offset);
                temp++;
            }   
        }
    } 
    /* remove it from PIC_assembly, mRNAdecay, proteinDecay, translation and re_calc*/
    for(i=which_gene;i<genotype->ngenes;i++)
    {
        genotype->min_N_activator_to_transc[i]=genotype->min_N_activator_to_transc[i+1];
        genotype->active_to_intermediate_rate[i]=genotype->active_to_intermediate_rate[i+1];             
        genotype->mRNA_decay_rate[i]=genotype->mRNA_decay_rate[i+1];
        genotype->protein_decay_rate[i]=genotype->protein_decay_rate[i+1];
        genotype->translation_rate[i]=genotype->translation_rate[i+1];  
        genotype->locus_length[i]=genotype->locus_length[i+1];
        genotype->recalc_TFBS[i]=YES;        
    }
    /* if the to-be-deleted gene is a tf gene, change ntfgenes as well*/
    if(protein_id<genotype->nproteins-1)
        genotype->ntfgenes--;    
    /* now change protein_pool and cisreg_cluster*/   
    update_protein_pool(genotype,protein_id,which_gene,'w'); 
    update_cisreg_cluster(genotype,which_gene,'w',NULL,-1,-1);  
    genotype->ngenes--;
}

void mut_whole_gene_deletion2(Genotype *genotype, Mutation *mut_record, RngStream RS)
{    
    int which_gene;
    int protein_id, N_gene_of_selection_protein, puppet_gene,position_of_puppet_gene;   
    Gene *target_gene,*gene;
    Protein *target_protein,*protein;
    ProteinFamily *target_protein_family,*protein_family;
    CisRegFamily *target_cisreg_family, *cisreg_family;
    /* check which genes can be deleted*/   
    N_gene_of_selection_protein=genotype->protein_pool[genotype->nproteins-1][0][0]; //does the effector gene have multipe copies   
    if(genotype->ntfgenes-N_SIGNAL_TF>1)//if there are more than one copies of genes of non-sensor TF  
    {
        if(N_gene_of_selection_protein==1)//if the effector have only one copies
        {
            protein_id=genotype->nproteins-1;
            while(protein_id==genotype->nproteins-1) //the last copy of the effector cannot be deleted. choose a different gene
            {   
                which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
                protein_id=genotype->which_protein[which_gene];
            }
        }
        else //otherwise any gene can be deleted
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
    }
    else //if there's only one copy of non-sensor TF gene
    {        
        while(1)
        {
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
            protein_id=genotype->which_protein[which_gene];
            if(protein_id==genotype->nproteins-1)//obviously we can only delete a copy of the effector gene
                break;                           //it's impossible for deletion when the non-sensor TF and the effector each has only one copy
                                                 //function draw_mutation would have set deletion rate to 0 in such a case
        }
    }
    /*record mutation info*/
    mut_record->which_gene=which_gene; 
    
    /**********************************/
    gene=genotype->genome;
    while(gene->next_gene_in_genome->id=which_gene)
        gene=gene->next_gene_in_genome;
    target_gene=gene->next_gene_in_genome;
    gene->next_gene_in_genome=target_gene->next_gene_in_genome;
    gene=gene->next_gene_in_genome;
    while(gene!=NULL)
        gene->id--;
    
    /*add target to waste*/
    if(genotype->last_empty_gene!=NULL)
        genotype->last_empty_gene->next_empty_gene=target_gene;
    else    
        genotype->last_empty_gene=target_gene;
    
    /*deal with proteins*/
    target_protein=target_gene->which_protein;
    if(target_protein->N_members==1)//protein has only gene
    {
        /*we need to delete targe protein as well*/
        target_protein_family=target_protein->which_protein_family;
        if(target_protein_family->N_members==1) // family has only one protein
        {
            /*we need to delete target protein family*/
            if(genotype->last_empty_protein_family!=NULL)
                genotype->last_empty_protein_family->next_empty_family=target_protein_family;
            else
                genotype->last_empty_protein_family=target_protein_family;
            protein_family=genotype->protein_family;
            while(protein_family->next_family!=target_protein_family)
                protein_family=protein_family->next_family;
            protein_family=target_protein_family->next_family;
            protein_family=protein_family->next_family;
            while(protein_family!=NULL)
                protein_family->id--;            
        }
        else
        {
            if(target_protein_family->first_protein==target_protein)
                target_protein_family=target_protein->next_protein_in_family;
            else                
            {
                protein=target_protein_family->first_protein;
                while(protein->next_protein_in_family!=target_protein)
                    protein=protein->next_protein_in_family;
                protein->next_protein_in_family=target_protein->next_protein_in_proteome;
            }
        }
        
        if(genotype->last_empty_protein==NULL)
            genotype->last_empty_protein=target_protein;
        else
            genotype->last_empty_protein=
        
        protein=genotype->proteome;
        while(protein->next_protein_in_proteome!=target_protein)/*by default the signal, i.e the first protein, cannot be deleted*/
            protein=protein->next_protein_in_proteome;
        protein->next_protein_in_proteome=target_protein->next_protein_in_proteome;        
        protein=protein->next_protein_in_family;
        while(protein!=NULL)
            protein->id--;  
    }
    else
    {
        gene=target_protein->first_gene;
        while(gene->next_gene_making_Protein!=target_gene)
            gene=gene->next_gene_making_Protein;
        gene->next_gene_making_Protein=target_gene->next_gene_making_Protein;
        gene=genotype->genome;
        while(gene->next_gene_in_genome!=target_gene)
            gene=gene->next_gene_making_Protein;
        gene->next_gene_making_Protein=target_gene->next_gene_making_Protein;      
    }  
    
    /*deal with cisreg family*/
    target_cisreg_family=target_gene->which_cisreg_family;
    if(target_cisreg_family->N_members==1)
    {
        cisreg_family=genotype->cisreg_family;        
        while(cisreg_family!=NULL)
            cisreg_family=cisreg_family->next_family;
        cisreg_family->next_family=target_cisreg_family->next_family;
        cisreg_family=cisreg_family->next_family;
        while(cisreg_family!=NULL)
            cisreg_family->id--;
    }
    else
    {
        
    }    
}


void reproduce_whole_gene_deletion(Genotype *genotype, Mutation *mut_record) // any gene can be deleted
{    
    int which_gene, offset, i,j, cisreg_seq_cluster_id, cisreg_seq_cluster_id_copy,puppet_gene,position_of_puppet_gene;
    char *temp;    
    int protein_id;       
    which_gene=mut_record->which_gene;    
    if(which_gene==genotype->ngenes-1)
    {
        cisreg_seq_cluster_id=genotype->which_cluster[which_gene];
        if(genotype->cisreg_cluster[cisreg_seq_cluster_id][1]!=NA)
        {            
            puppet_gene=genotype->cisreg_cluster[cisreg_seq_cluster_id][0];
        }
        else 
        {
            protein_id=genotype->which_protein[which_gene];
            position_of_puppet_gene=0;
            if(genotype->protein_pool[protein_id][1][position_of_puppet_gene]==which_gene)position_of_puppet_gene++;
            puppet_gene=genotype->protein_pool[protein_id][1][position_of_puppet_gene];            
            cisreg_seq_cluster_id_copy=genotype->which_cluster[puppet_gene];            
            j=0;
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]!=NA)j++;
            genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]=which_gene;            
            genotype->which_cluster[which_gene]=cisreg_seq_cluster_id_copy;            
            cisreg_seq_cluster_id_copy=cisreg_seq_cluster_id;            
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][0]!=NA)
            {                
                j=0;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]!=NA)
                {
                    genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]=NA;
                    j++;
                }                                    
                j=0;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][j]!=NA)
                {                        
                    genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]=genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][j];
                    genotype->which_cluster[genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][j]]--;
                    j++;
                }
                cisreg_seq_cluster_id_copy++;
            }
        }
        for(j=0;j<CISREG_LEN;j++)        
            genotype->cisreg_seq[which_gene][j]=genotype->cisreg_seq[puppet_gene][j];
        genotype->min_N_activator_to_transc[which_gene]=genotype->min_N_activator_to_transc[puppet_gene];  
        genotype->active_to_intermediate_rate[which_gene]=genotype->active_to_intermediate_rate[puppet_gene];           
        genotype->mRNA_decay_rate[which_gene]=genotype->mRNA_decay_rate[puppet_gene];
        genotype->protein_decay_rate[which_gene]=genotype->protein_decay_rate[puppet_gene];
        genotype->translation_rate[which_gene]=genotype->translation_rate[puppet_gene]; 
        genotype->locus_length[which_gene]=genotype->locus_length[puppet_gene];
        genotype->recalc_TFBS[which_gene]=YES;        
        which_gene=puppet_gene;        
    }    
    temp = &genotype->cisreg_seq[which_gene][0];			
    offset=CISREG_LEN;    
    for(i=0;i<CISREG_LEN*(genotype->ngenes-which_gene-1);i++) 
    {				
        *temp=*(temp+offset);         
        temp++;				
    }
    protein_id=genotype->which_protein[which_gene];   
    if(protein_id<genotype->nproteins-1)
    {      
        if(genotype->protein_pool[protein_id][0][0]==1) 
        {  
            temp=&genotype->tf_seq[protein_id][0];
            offset=TF_ELEMENT_LEN;
            for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
            {
                *temp=*(temp+offset);
                temp++;
            }            
            temp=&genotype->tf_seq_rc[protein_id][0];
            for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
            {
                *temp=*(temp+offset);
                temp++;
            }   
        }
    }
    for(i=which_gene;i<genotype->ngenes;i++)
    {
        genotype->min_N_activator_to_transc[i]=genotype->min_N_activator_to_transc[i+1];
        genotype->active_to_intermediate_rate[i]=genotype->active_to_intermediate_rate[i+1];            
        genotype->mRNA_decay_rate[i]=genotype->mRNA_decay_rate[i+1];
        genotype->protein_decay_rate[i]=genotype->protein_decay_rate[i+1];
        genotype->translation_rate[i]=genotype->translation_rate[i+1];  
        genotype->locus_length[i]=genotype->locus_length[i+1];
        genotype->recalc_TFBS[i]=YES;
    }
    if(protein_id<genotype->nproteins-1)
        genotype->ntfgenes--;
    update_protein_pool(genotype,protein_id,which_gene,'w'); 
    update_cisreg_cluster(genotype,which_gene,'w',NULL,-1,-1);
    genotype->ngenes--;   
}

/*
 * Duplicate a whole cis-reg sequence
 */
void mut_duplication(Genotype *genotype, Mutation *mut_record, RngStream RS) 
{
    int which_gene, which_gene_copy, i, protein_id;
    char *temp1, *temp2;    
    if(genotype->protein_pool[genotype->nproteins-1][0][0]>=N_EFFECTOR_GENES) // too many effector gene copies
    {     
        //note that it's not possible to have too many effector gene copies and too many tf gene copies at the same time
        //because that will make duplication rate 0.
        which_gene=genotype->ngenes-1;
        while(genotype->which_protein[which_gene]==genotype->nproteins-1)
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
    }
    else
    {
        if(genotype->ntfgenes>=N_TF_GENES) // too many tf gene copies
        {
            while(1)
            {
                which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
                if(genotype->which_protein[which_gene]==genotype->nproteins-1)
                    break;
            }    
        }            
        else //any gene can be duplicated
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1); 
    }   
    /*record mutation info*/
    mut_record->which_gene=which_gene;  
    /*To preserve data structure, the last gene has to be an effector gene.*/
    /*We want to shift the effector gene, and put the duplicated gene to the second last place.
     *But if it is the effector gene that gets duplicated, things are a little different */
    if(which_gene==genotype->ngenes-1)
        which_gene_copy=which_gene+1; 
    else
        which_gene_copy=which_gene;
    /* copy the promoter*/
    temp1=&genotype->cisreg_seq[which_gene_copy][0]; /* points to the gene to be duplicated*/
    temp2=&genotype->cisreg_seq[genotype->ngenes-1][CISREG_LEN-1]; /* points to the end of the effector gene */
    /* shift the sequences of the effector gene CISREG_LEN bp */
    for(i=0;i<CISREG_LEN;i++) 
    {
        *(temp2+CISREG_LEN)=*temp2;
        temp2--;
    }  
    /* point temp2 to the start of the original position of the selection gene */
    temp2=&genotype->cisreg_seq[genotype->ngenes-1][0];
    /* put the duplicated gene at the original place of the selection gene */
    for(i=0;i<CISREG_LEN;i++) 
        *temp2++=*temp1++;
    /* shift the effector genes to make a slot*/
    genotype->min_N_activator_to_transc[genotype->ngenes]=genotype->min_N_activator_to_transc[genotype->ngenes-1];
    genotype->active_to_intermediate_rate[genotype->ngenes]=genotype->active_to_intermediate_rate[genotype->ngenes-1]; 
    genotype->mRNA_decay_rate[genotype->ngenes]=genotype->mRNA_decay_rate[genotype->ngenes-1];
    genotype->protein_decay_rate[genotype->ngenes]=genotype->protein_decay_rate[genotype->ngenes-1];
    genotype->translation_rate[genotype->ngenes]=genotype->translation_rate[genotype->ngenes-1];
    genotype->locus_length[genotype->ngenes]=genotype->locus_length[genotype->ngenes-1];
    genotype->recalc_TFBS[genotype->ngenes]=YES;        
    /* copy and paste info to the slot*/
    genotype->min_N_activator_to_transc[genotype->ngenes-1]=genotype->min_N_activator_to_transc[which_gene_copy];
    genotype->active_to_intermediate_rate[genotype->ngenes-1]=genotype->active_to_intermediate_rate[which_gene_copy];
    genotype->mRNA_decay_rate[genotype->ngenes-1]=genotype->mRNA_decay_rate[which_gene_copy];
    genotype->protein_decay_rate[genotype->ngenes-1]=genotype->protein_decay_rate[which_gene_copy];
    genotype->translation_rate[genotype->ngenes-1]=genotype->translation_rate[which_gene_copy];
    genotype->locus_length[genotype->ngenes-1]=genotype->locus_length[which_gene_copy];
    genotype->recalc_TFBS[genotype->ngenes-1]=YES;
    /* update protein_pool*/
    protein_id=genotype->which_protein[which_gene];    
    update_protein_pool(genotype,protein_id,which_gene,'d'); 
    /* update cisreg_cluster*/    
    update_cisreg_cluster(genotype,which_gene,'d',NULL,-1,-1);
    /* update gene numbers*/  
    if(protein_id<genotype->nproteins-1)//note duplication do not change nproteins           
        genotype->ntfgenes++;     
    genotype->ngenes++;
}

void reproduce_gene_duplication(Genotype *genotype, Mutation *mut_record) //any gene can be duplicated
{   
    int which_gene,which_gene_copy, i, protein_id;
    char *temp1, *temp2; 
    /*get the gene to be duplicated from record*/
    which_gene=mut_record->which_gene;    
    if(which_gene==genotype->ngenes-1)
        which_gene_copy=which_gene+1; 
    else
        which_gene_copy=which_gene; 
    temp1=&genotype->cisreg_seq[which_gene_copy][0]; 
    temp2=&genotype->cisreg_seq[genotype->ngenes-1][CISREG_LEN-1];    
    for(i=0;i<CISREG_LEN;i++) 
    {
        *(temp2+CISREG_LEN)=*temp2; 
        temp2--;
    }    
    temp2=&genotype->cisreg_seq[genotype->ngenes-1][0];  
    for(i=0;i<CISREG_LEN;i++) 
        *temp2++=*temp1++;
    genotype->min_N_activator_to_transc[genotype->ngenes]=genotype->min_N_activator_to_transc[genotype->ngenes-1];
    genotype->active_to_intermediate_rate[genotype->ngenes]=genotype->active_to_intermediate_rate[genotype->ngenes-1]; 
    genotype->mRNA_decay_rate[genotype->ngenes]=genotype->mRNA_decay_rate[genotype->ngenes-1];
    genotype->protein_decay_rate[genotype->ngenes]=genotype->protein_decay_rate[genotype->ngenes-1];
    genotype->translation_rate[genotype->ngenes]=genotype->translation_rate[genotype->ngenes-1];
    genotype->locus_length[genotype->ngenes]=genotype->locus_length[genotype->ngenes-1];
    genotype->recalc_TFBS[genotype->ngenes]=YES;
    genotype->min_N_activator_to_transc[genotype->ngenes-1]=genotype->min_N_activator_to_transc[which_gene_copy];
    genotype->active_to_intermediate_rate[genotype->ngenes-1]=genotype->active_to_intermediate_rate[which_gene_copy];
    genotype->mRNA_decay_rate[genotype->ngenes-1]=genotype->mRNA_decay_rate[which_gene_copy];
    genotype->protein_decay_rate[genotype->ngenes-1]=genotype->protein_decay_rate[which_gene_copy];
    genotype->translation_rate[genotype->ngenes-1]=genotype->translation_rate[which_gene_copy]; 
    genotype->locus_length[genotype->ngenes-1]=genotype->locus_length[which_gene_copy];
    genotype->recalc_TFBS[genotype->ngenes-1]=YES;
    protein_id=genotype->which_protein[which_gene];
    update_protein_pool(genotype,protein_id,which_gene,'d');     
    update_cisreg_cluster(genotype,which_gene,'d',NULL,-1,-1);
    if(protein_id<genotype->nproteins-1)
        genotype->ntfgenes++;      
    genotype->ngenes++;     
}

/*
 * Mutation to the binding sequence of a TF gene
 */
void mut_binding_sequence(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int which_gene, which_nucleotide, protein_id, i;
    char nucleotide;      
    char *tf_seq, *tf_seq_rc,*temp1,*temp2;
    /*get which gene to mutate*/
    which_gene=RngStream_RandInt(RS,0,genotype->ngenes-2); // ngenes-1 is the selection gene, so tf genes must be 0 to ngenes-2
    while(genotype->which_protein[which_gene]==genotype->nproteins-1)
        which_gene=RngStream_RandInt(RS,0,genotype->ngenes-2);
    /*if this TF has more than one copies of gene, then the mutation adds a new tf
     *which requires a new slot in tf_seq and tf_seq_rc to store the new binding seq*/
    protein_id=genotype->which_protein[which_gene];    
    if(genotype->protein_pool[protein_id][0][0]>1)
    {
        /*points to the first nucleotide of the binding sequence*/
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
        /*points to an empty slot*/
        temp1=&genotype->tf_seq[genotype->nproteins-1][0];
        temp2=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        /*copy the binding sequences to empty slots*/
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_seq++;
            *temp2++=*tf_seq_rc++;
        }
        /*point tf_seq and tf_seq_rc to the new slot so that we can apply mutation later*/
        tf_seq=&genotype->tf_seq[genotype->nproteins-1][0];
        tf_seq_rc=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        /*update protein pool*/
        update_protein_pool(genotype,protein_id,which_gene,'c');  
    }    
    else /*if this tf has only one copy of gene, no new slot is required*/
    {    
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
    }
    /*mutate the binding sequence*/
    /*mutation only changes one nucleotide in the binding sequence*/
    which_nucleotide=RngStream_RandInt(RS,0,TF_ELEMENT_LEN-1);
    nucleotide=set_base_pair(RngStream_RandU01(RS));        
    while (nucleotide == tf_seq[which_nucleotide])
        nucleotide=set_base_pair(RngStream_RandU01(RS));
    tf_seq[which_nucleotide]=nucleotide;    
    /*record mutation info*/
    mut_record->which_gene=which_gene;
    mut_record->which_nucleotide=which_nucleotide;
    mut_record->nuc_diff[0]=nucleotide;    
    /* update the reverse complement sequence*/
    switch (nucleotide)
    {
        case 'g':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='c'; break;
        case 'c':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='g'; break;
        case 'a':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='t'; break;
        case 't':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='a'; break;
    }  
    /* The binding sites on every promoter needs recalculation */
    for(i=0;i<genotype->ngenes;i++)    
        genotype->recalc_TFBS[i]=YES;
    /*decide whether to update cisreg clusters. Mutation to binding seq may differ bs distributions among genes in a cluster*/       
    int new_clusters[NGENES][NGENES]; // 
    int genes_in_cluster[NGENES];
    int N_genes_in_cluster,no_difference,reference_gene,gene_to_be_sorted;
    int N_new_clusters,N_genes_in_new_cluster,j,k;
    calc_all_binding_sites(genotype); 
    i=N_SIGNAL_TF;
    while(genotype->cisreg_cluster[i][0]!=NA)/*check each cisreg cluster*/
    {        
        N_new_clusters=0;
        N_genes_in_cluster=0;
        for(j=0;j<NGENES;j++)
        {
            for(k=0;k<NGENES;k++)
                new_clusters[j][k]=NA;
        }    
        while(genotype->cisreg_cluster[i][N_genes_in_cluster]!=NA)
        {
            genes_in_cluster[N_genes_in_cluster]=genotype->cisreg_cluster[i][N_genes_in_cluster];
            N_genes_in_cluster++;
        }
        /*the while loop below sort genes in a cluster into groups based on whether they have the same BS distributions*/ 
        /*We use one gene in the original cluster as a reference, and sort all genes 
         *that have the same binding sites as the reference gene into a new cluster. 
         *After all genes are sorted, we check if the original cluster turns into multiple
         *new cluster.*/
        while(N_genes_in_cluster>0)//run until every gene in the original cluster has been sorted into new clusters
        {
            reference_gene=genes_in_cluster[0];//compare each gene to the first gene in cluster            
            N_genes_in_new_cluster=0;
            which_gene=0;//start comparison from the first gene in the cluster
            while(which_gene<N_genes_in_cluster)
            {
                no_difference=1;
                gene_to_be_sorted=genes_in_cluster[which_gene];                
                if(genotype->binding_sites_num[gene_to_be_sorted]==genotype->binding_sites_num[reference_gene])
                {
                    for(j=0;j<genotype->binding_sites_num[reference_gene];j++)
                    {
                        if(genotype->all_binding_sites[reference_gene][j].BS_pos!=genotype->all_binding_sites[gene_to_be_sorted][j].BS_pos ||
                            genotype->all_binding_sites[reference_gene][j].tf_id!=genotype->all_binding_sites[gene_to_be_sorted][j].tf_id ||
                            genotype->all_binding_sites[reference_gene][j].mis_match!=genotype->all_binding_sites[gene_to_be_sorted][j].mis_match)
                        {
                            no_difference=0;
                            break;
                        }
                    }
                }
                else
                    no_difference=0;                
                if(no_difference)//if the gene under comparison has the same binding sites as the reference gene
                {
                    /*put the gene into the new cluster*/
                    new_clusters[N_new_clusters][N_genes_in_new_cluster]=genes_in_cluster[which_gene];
                    N_genes_in_new_cluster++;
                    /*shift to remove the gene from the original cluster*/
                    for(j=which_gene;j<N_genes_in_cluster-1;j++)
                        genes_in_cluster[j]=genes_in_cluster[j+1];
                    N_genes_in_cluster--;                
                }
                else //otherwise, compare next gene with the reference gene                
                    which_gene++;                
            }
            N_new_clusters++;
        }
        if(N_new_clusters!=1)//if the original cluster turns into multiple clusters
            update_cisreg_cluster(genotype,-1,'c',new_clusters,N_new_clusters,i);
        i++;
    }
    /* Calling calc_all_binding_sites reset recalc_TFBS to 0, so we need to turn them on again */
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=YES;  
}

void reproduce_mut_binding_sequence(Genotype *genotype, Mutation *mut_record)
{    
    int which_gene, which_nucleotide, protein_id, i;    
    char *tf_seq, *tf_seq_rc, *temp1,*temp2;     
    which_gene=mut_record->which_gene;
    protein_id=genotype->which_protein[which_gene];  
    if(genotype->protein_pool[protein_id][0][0]>1)
    {        
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
        temp1=&genotype->tf_seq[genotype->nproteins-1][0];
        temp2=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_seq++;
            *temp2++=*tf_seq_rc++;
        }      
        tf_seq=&genotype->tf_seq[genotype->nproteins-1][0];
        tf_seq_rc=&genotype->tf_seq_rc[genotype->nproteins-1][0];        
        update_protein_pool(genotype,protein_id,which_gene,'c');  
    }    
    else 
    {    
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
    }        			
    which_nucleotide=mut_record->which_nucleotide;        
    tf_seq[which_nucleotide]=mut_record->nuc_diff[1];
    switch (tf_seq[which_nucleotide])
    {
//        case 'g':
//            tf_seq_rc[which_nucleotide]='c'; break;
//        case 'c':
//            tf_seq_rc[which_nucleotide]='g'; break;
//        case 'a':
//            tf_seq_rc[which_nucleotide]='t'; break;
//        case 't':
//            tf_seq_rc[which_nucleotide]='a'; break;
        case 'g':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='c'; break;
        case 'c':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='g'; break;
        case 'a':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='t'; break;
        case 't':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='a'; break;
    }     
    for(i=0;i<genotype->ngenes;i++)
        genotype->recalc_TFBS[i]=YES;
    calc_all_binding_sites(genotype);
    int new_clusters[NGENES][NGENES],genes_in_cluster[NGENES];
    int N_genes_in_cluster,no_difference,reference_gene,gene_to_be_sorted;
    int N_new_clusters,N_genes_in_new_cluster,j,k;
    i=N_SIGNAL_TF;
    while(genotype->cisreg_cluster[i][0]!=NA)
    {        
        N_new_clusters=0;
        N_genes_in_cluster=0;
        for(j=0;j<NGENES;j++)
        {
            for(k=0;k<NGENES;k++)
                new_clusters[j][k]=NA;
        }  
        while(genotype->cisreg_cluster[i][N_genes_in_cluster]!=NA)
        {            
            genes_in_cluster[N_genes_in_cluster]=genotype->cisreg_cluster[i][N_genes_in_cluster];
            N_genes_in_cluster++;
        }
        while(N_genes_in_cluster>0)
        {
            reference_gene=genes_in_cluster[0];
            N_genes_in_new_cluster=0;
            which_gene=0;
            while(which_gene<N_genes_in_cluster)
            {
                no_difference=1;
                gene_to_be_sorted=genes_in_cluster[which_gene];                
                if(genotype->binding_sites_num[gene_to_be_sorted]==genotype->binding_sites_num[reference_gene])
                {
                    for(j=0;j<genotype->binding_sites_num[reference_gene];j++)
                    {
                        if(genotype->all_binding_sites[reference_gene][j].BS_pos!=genotype->all_binding_sites[gene_to_be_sorted][j].BS_pos ||
                            genotype->all_binding_sites[reference_gene][j].tf_id!=genotype->all_binding_sites[gene_to_be_sorted][j].tf_id ||
                            genotype->all_binding_sites[reference_gene][j].mis_match!=genotype->all_binding_sites[gene_to_be_sorted][j].mis_match)
                        {
                            no_difference=0;
                            break;
                        }
                    }
                }
                else
                    no_difference=0;                
                if(no_difference)
                {
                    new_clusters[N_new_clusters][N_genes_in_new_cluster]=genes_in_cluster[which_gene];
                    N_genes_in_new_cluster++;
                    for(j=which_gene;j<N_genes_in_cluster-1;j++)
                        genes_in_cluster[j]=genes_in_cluster[j+1];
                    N_genes_in_cluster--;                
                }
                else 
                    which_gene++;
                
            }
            N_new_clusters++;
        }
        if(N_new_clusters!=1)
            update_cisreg_cluster(genotype,-1,'c',new_clusters,N_new_clusters,i);
        i++;
    }    
    for(i=0;i<genotype->ngenes;i++)
        genotype->recalc_TFBS[i]=YES;   
}

/* Mutations to the rate of mRNA_decay, translation, protein_decay, and pic_disassembly 
 * will be mutated. 
 */
void mut_kinetic_constant(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    float random;    
    int which_gene;
    float total_mut_flux=   proportion_mut_kdis+
                            proportion_mut_mRNA_decay+
                            proportion_mut_protein_decay+
                            proportion_mut_translation_rate+
                            proportion_mut_cooperation;                          
    
    random=RngStream_RandU01(RS);
    which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);    
    /*record mutation info*/
    mut_record->which_gene=which_gene;
    while(1)
    {
        /********************************** mut kdis *******************************/
        if(random<=proportion_mut_kdis/total_mut_flux) 
        {  
            genotype->active_to_intermediate_rate[which_gene]=mut_make_new_value(genotype->active_to_intermediate_rate[which_gene],
                                                                                miu_ACT_TO_INT_RATE,
                                                                                sigma_ACT_TO_INT_RATE,
                                                                                MAX_ACT_TO_INT_RATE,
                                                                                MIN_ACT_TO_INT_RATE,
                                                                                RS,
                                                                                mut_record);
            /*record mutation info*/
            mut_record->kinetic_type=0;
            mut_record->kinetic_diff=genotype->active_to_intermediate_rate[which_gene];
            break;
        }
        else 
        {
            random-=proportion_mut_kdis/total_mut_flux;
            /******************************** mut mRNAdecay **********************/
            if(random<=proportion_mut_mRNA_decay/total_mut_flux) 
            {       
                genotype->mRNA_decay_rate[which_gene]=mut_make_new_value(genotype->mRNA_decay_rate[which_gene],miu_mRNA_decay,sigma_mRNA_decay, MAX_MRNA_DECAY, MIN_MRNA_DECAY,RS,mut_record);
                /*record mutation info*/
                mut_record->kinetic_type=1;
                mut_record->kinetic_diff=genotype->mRNA_decay_rate[which_gene];
                break;
            }
            else 
            {
                random-=proportion_mut_mRNA_decay/total_mut_flux;
                /*************************** mut translation *********************/
                if(random<=proportion_mut_translation_rate/total_mut_flux) 
                {       
                    genotype->translation_rate[which_gene]=mut_make_new_value(genotype->translation_rate[which_gene],miu_translation_init,sigma_translation_init, MAX_TRANSLATION_RATE, MIN_TRANSLATION_RATE,RS,mut_record);        
                    mut_record->kinetic_type=2;
                    mut_record->kinetic_diff=genotype->translation_rate[which_gene];
                    break;
                }
                else 
                {  
                    random-=proportion_mut_translation_rate/total_mut_flux;
                    /********************* mut protein decay **********************/
                    if(random<=proportion_mut_protein_decay/total_mut_flux)
                    {                   
                        genotype->protein_decay_rate[which_gene]=mut_make_new_value(genotype->protein_decay_rate[which_gene],miu_protein_decay,sigma_protein_decay,MAX_PROTEIN_DECAY,MIN_PROTEIN_DECAY,RS,mut_record);                  
                        mut_record->kinetic_type=3;
                        mut_record->kinetic_diff=genotype->protein_decay_rate[which_gene];
                        break;
                    }
                    else 
                    {
                        random-=proportion_mut_protein_decay/total_mut_flux;
                        /****************** mut cooperation ***********************/
                        if(random<=proportion_mut_cooperation/total_mut_flux)
                        {                            
                            while(genotype->which_protein[which_gene]!=genotype->nproteins-1)
                                which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
                            mut_record->which_gene=which_gene;
                            if(genotype->min_N_activator_to_transc[which_gene]==1)
                                genotype->min_N_activator_to_transc[which_gene]=2;
                            else
                                genotype->min_N_activator_to_transc[which_gene]=1;
                            mut_record->kinetic_type=4;
                            mut_record->kinetic_diff=(float)genotype->min_N_activator_to_transc[which_gene];
                            update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1); 
                            break;
                        }
                    }
                } 
            }
        }
    }
}

float mut_make_new_value(float old_val, float miu, float sigma, float upper_bound, float lower_bound, RngStream RS, Mutation *mut_record)
{
    float new_val;
    new_val=old_val;
    while(new_val==old_val) 
    {
#if OLD_MODEL
        new_val=old_val*exp(sigma*gasdev(RS)+mutational_regression_rate*(miu-log(old_val)));
#else
        new_val=old_val*pow(10.0,sigma*gasdev(RS)+mutational_regression_rate*(miu-log10(old_val)));
        if(new_val>upper_bound)
        {
            new_val=upper_bound;
            mut_record->N_hit_bound++;
        }
        if(new_val<lower_bound)
        {
            new_val=lower_bound;
            mut_record->N_hit_bound++;
        }
#endif
//        while(new_val>=upper_bound || new_val<=lower_bound)
//        {
//            new_val=old_val*pow(10.0,sigma*gasdev(RS)+mutational_regression_rate*(miu-log10(old_val)));
//            mut_record->N_hit_bound++;
//        }            
    }
    return new_val;
}

void reproduce_mut_kinetic_constant(Genotype *genotype, Mutation *mut_record)
{    
    int which_gene;            
    which_gene=mut_record->which_gene; /* which gene */ 
    switch (mut_record->kinetic_type)
    {
        case 0: /* mut r_act_to_int */
            genotype->active_to_intermediate_rate[which_gene]=mut_record->kinetic_diff;
            break;
        case 1: /* mut mRNAdecay */ 
            genotype->mRNA_decay_rate[which_gene]=mut_record->kinetic_diff;
            break;
        case 2: /* mut translation */                      
            genotype->translation_rate[which_gene]=mut_record->kinetic_diff;           
            break;        
        case 3: /* mut protein decay */
            genotype->protein_decay_rate[which_gene]=mut_record->kinetic_diff;            
            break;
        case 4: /* mut cooperation */
            genotype->min_N_activator_to_transc[which_gene]=(int)mut_record->kinetic_diff;
            update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1); 
    }
}

void mut_locus_length(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int which_gene;
    float random;    
    while(1)
    {
        random=RngStream_RandU01(RS)*genotype->total_loci_length;
        which_gene=N_SIGNAL_TF-1;
        while(which_gene<genotype->ngenes-1 && random>0.0)
        {
            which_gene++;
            random-=genotype->locus_length[which_gene];
        }
        if(genotype->locus_length[which_gene]>0.0)
            break;
    }
    genotype->total_loci_length-=genotype->locus_length[which_gene];
    if(genotype->locus_length[which_gene]>=3.0*MAX_PROTEIN_LENGTH)
        genotype->locus_length[which_gene]-=3.0;
    else if(genotype->locus_length[which_gene]<=3.0*MAX_PROTEIN_LENGTH)
        genotype->locus_length[which_gene]+=3.0;
    else                            
        genotype->locus_length[which_gene]+=(RngStream_RandU01(RS)<0.5)?3.0:-3.0;   
    genotype->total_loci_length+=genotype->locus_length[which_gene];
    mut_record->mut_type='l';   
    mut_record->which_gene=which_gene;
    mut_record->kinetic_diff=genotype->locus_length[which_gene];
}

void reproduce_mut_locus_length(Genotype *genotype, Mutation *mut_record)
{
    genotype->total_loci_length-=genotype->locus_length[mut_record->which_gene];
    genotype->locus_length[mut_record->which_gene]=mut_record->kinetic_diff;
    genotype->total_loci_length+=genotype->locus_length[mut_record->which_gene];
}

void mut_identity(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int tf_id,protein_id,i;   
    char *tf_seq,*tf_seq_rc,*temp1,*temp2; 
    /*which tf gene to mutate*/
    tf_id = RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-2); // the first sensor tf must be an activator, therefore is not subject to mutation
    protein_id=genotype->which_protein[tf_id];    
    while(protein_id==genotype->nproteins-1) // copies of effector gene are not subjected to this mutation
    {
        tf_id = RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-2); 
        protein_id=genotype->which_protein[tf_id];
    }    
    /*save record*/
    mut_record->which_gene=tf_id;       
    /* if this tf gene has more than one copies, the mutation adds a new protein*/
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {
        /*give the new protein its own binding sequence*/
        tf_seq=&genotype->tf_seq[protein_id][0]; 
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
        temp1=&genotype->tf_seq[genotype->nproteins-1][0];
        temp2=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_seq++;
            *temp2++=*tf_seq_rc++;
        } 
        update_protein_pool(genotype,protein_id,tf_id,'e');
        /* update_protein_pool put the new protein at nproteins-1 and then increases nproteins by 1, 
        * so the new protein is at nproteins-2 now. Note that N_act and N_rep is updated in update_protein_pool*/ 
        genotype->protein_identity[genotype->nproteins-2]=(genotype->protein_identity[genotype->nproteins-2]==ACTIVATOR)?REPRESSOR:ACTIVATOR;         
    }
    else //otherwise we just flip the property of an exisiting TF
    {
        genotype->protein_identity[protein_id]=(genotype->protein_identity[protein_id]==ACTIVATOR)?REPRESSOR:ACTIVATOR; 
        if(genotype->protein_identity[protein_id]==ACTIVATOR)//mutate to activator
        {
            genotype->N_act++;
            genotype->N_rep--;
        }
        else if(genotype->protein_identity[protein_id]==REPRESSOR) //otherwise to repressor
        { 
            genotype->N_rep++;
            genotype->N_act--;
        }
    }    
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=YES; /* recalculate binding sites on every promoter */        
}

void reproduce_mut_identity(Genotype *genotype, Mutation *mut_record)
{
    int tf_id, protein_id,i;
    char *tf_seq,*tf_seq_rc,*temp1,*temp2;    
    tf_id = mut_record->which_gene;
    protein_id=genotype->which_protein[tf_id];
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
        temp1=&genotype->tf_seq[genotype->nproteins-1][0];
        temp2=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_seq++;
            *temp2++=*tf_seq_rc++;
        }     
        update_protein_pool(genotype,protein_id,tf_id,'e');  
        genotype->protein_identity[genotype->nproteins-2]=(genotype->protein_identity[genotype->nproteins-2]==ACTIVATOR)?REPRESSOR:ACTIVATOR;     
    }
    else
    {
        genotype->protein_identity[protein_id]=(genotype->protein_identity[protein_id]==ACTIVATOR)?REPRESSOR:ACTIVATOR; 
        if(genotype->protein_identity[protein_id]==ACTIVATOR)
        {
            genotype->N_act++;
            genotype->N_rep--;
        }
        else if(genotype->protein_identity[protein_id]==REPRESSOR)
        { 
            genotype->N_rep++;
            genotype->N_act--;
        }
    }    
    for(i=0;i<genotype->ngenes;i++)    
        genotype->recalc_TFBS[i]=YES;   
}

/*
 * mutate affinity of TF
 */
void mut_Kd(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int i;   
    float new_Kd; 
    int tf_id,protein_id;
    char *tf_seq,*tf_seq_rc,*temp1,*temp2;
    /*which TF to mutate*/
    tf_id=RngStream_RandInt(RS,0,genotype->ngenes-2);
    protein_id=genotype->which_protein[tf_id];
    while(protein_id==genotype->nproteins-1)
    {
        tf_id=RngStream_RandInt(RS,0,genotype->ngenes-2);
        protein_id=genotype->which_protein[tf_id];
    }
    /*generate a new koff */    
    new_Kd=mut_make_new_value(genotype->Kd[protein_id],miu_Kd,sigma_Kd,MAX_KD,MIN_KD,RS,mut_record);  
#if OLD_MODEL  
    while(new_Kd>NS_Kd)
        new_Kd=mut_make_new_value(genotype->Kd[protein_id],miu_Kd,sigma_Kd,MAX_KD,MIN_KD,RS,mut_record);  
#endif
    /* if this tf gene has more than one copies, the mutation adds a new protein*/    
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {    
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
        temp1=&genotype->tf_seq[genotype->nproteins-1][0];
        temp2=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_seq++;
            *temp2++=*tf_seq_rc++;
        }
        update_protein_pool(genotype,protein_id,tf_id,'f');
        /* update_protein_pool put the new protein at nproteins-1 and then increases nproteins by 1, 
         * so the new protein is at nproteins-2 now.*/
        genotype->Kd[genotype->nproteins-2]=new_Kd;  
    }                                                                                                           
    else
        genotype->Kd[protein_id]=new_Kd;     
    /*record mutation*/
    mut_record->which_gene=tf_id;
    mut_record->kinetic_diff=new_Kd; 
    /* recalculate binding sites on every promoter */
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=YES;
}

void reproduce_mut_Kd(Genotype *genotype, Mutation *mut_record)
{   
    int tf_id,protein_id,i;
    char *tf_seq, *tf_seq_rc, *temp1, *temp2;    
    tf_id=mut_record->which_gene;
    protein_id=genotype->which_protein[tf_id];    
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {    
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
        temp1=&genotype->tf_seq[genotype->nproteins-1][0];
        temp2=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_seq++;
            *temp2++=*tf_seq_rc++;
        }    
        update_protein_pool(genotype,protein_id,tf_id,'f');  
        genotype->Kd[genotype->nproteins-2]=mut_record->kinetic_diff;  
    }    
    else
        genotype->Kd[protein_id]=mut_record->kinetic_diff;
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=YES;   
}

/*
 *Wrapper of all mutation functions
 */
void mutate(Genotype *genotype, RngStream RS, Mutation *mut_record)
{
    int i;   
    draw_mutation(genotype, &(mut_record->mut_type),RS);
    /*reset records of nucleotides*/
    for(i=0;i<3;i++)
        mut_record->nuc_diff[i]='\0';
    switch (mut_record->mut_type)
    {
        case 's': //substitution in cis-reg       		
            mut_substitution(genotype,mut_record,RS);
            break;      
        case 'w': // whole gene deletion.          
            mut_whole_gene_deletion(genotype,mut_record,RS);            
            break;
        case 'd': // Whole gene duplication                  
            mut_duplication(genotype,mut_record,RS);            
            break;
        case 'c': //binding sequence 
            mut_binding_sequence(genotype,mut_record,RS);            
            break;
        case 'k': //mutations to kinetic constants        
            mut_kinetic_constant(genotype, mut_record,RS);           
            break;
        case 'e': //activator to repressor or the reverse
            mut_identity(genotype, mut_record, RS);            
            break;
        case 'f': //mutations to the Kd of a tf
            mut_Kd(genotype,mut_record,RS);            
            break;
        case 'l': //mutate locus length
            mut_locus_length(genotype,mut_record, RS);
            break;
    }
}

/* this function perform mutation indicated by input. Used to reproduce the genotype following a serial of mutations*/
void reproduce_mutate(Genotype *genotype, Mutation *mut_record,RngStream RS)
{ 
    switch (mut_record->mut_type)
    {
        case 's': //substitution        		
            reproduce_substitution(genotype,mut_record);
            break;        
        case 'w': // whole gene deletion.         
            reproduce_whole_gene_deletion(genotype,mut_record);
            break;
        case 'd': // Whole gene duplication                 
            reproduce_gene_duplication(genotype,mut_record);
            break;        
        case 'c': //binding sequence        
            reproduce_mut_binding_sequence(genotype,mut_record);
            break;
        case 'k': //mutations in kinetic constants        
            reproduce_mut_kinetic_constant(genotype, mut_record);            
            break;
        case 'e': //changing the identity of a TF
            reproduce_mut_identity(genotype, mut_record);
            break;
        case 'f':
            reproduce_mut_Kd(genotype,mut_record);
            break;  
        case 'l':
            reproduce_mut_locus_length(genotype, mut_record);
            break;
    }   
}

/* this function calculates the probability of different mutations based on
 * the current genotype, and draws a mutation accordingly.
 */
void draw_mutation(Genotype *genotype, char *mut_type, RngStream RS)
{
    float random;
    float tot_mut_rate=0.0;
    float tot_subs_rate, tot_dup_rate, tot_sil_rate, tot_mut_kin_rate, tot_mut_identity_rate, tot_mut_binding_seq_rate, tot_mut_koff_rate; 
    float tot_mut_locus_length_rate;
    
    int N_genes_to_be_duplicated; 
    /* duplication rate*/    
    tot_dup_rate=0.0;
    N_genes_to_be_duplicated=genotype->ngenes-N_SIGNAL_TF;//NA to sensor TF gene
    if(genotype->ntfgenes>=N_TF_GENES)//too many non-sensor TF genes.
        N_genes_to_be_duplicated-=genotype->ntfgenes-N_SIGNAL_TF; //do not duplicate non-sensor TF gene anymore
    if(genotype->protein_pool[genotype->nproteins-1][0][0]>=N_EFFECTOR_GENES)//too many effector gene
        N_genes_to_be_duplicated-=genotype->protein_pool[genotype->nproteins-1][0][0];//do not duplicate effector gene anymore
    tot_dup_rate=N_genes_to_be_duplicated*DUPLICATION;
    tot_mut_rate+=tot_dup_rate;   
    /* silencing rate*/ 
    if(genotype->ntfgenes-N_SIGNAL_TF>1)//if there are more than one copy of non-sensor TF gene
        tot_sil_rate=(genotype->ntfgenes-N_SIGNAL_TF)*SILENCING;    
    else
        tot_sil_rate=0.0; //the last non-sensor TF cannot be deleted
    if(genotype->protein_pool[genotype->nproteins-1][0][0]>1)//if there are more than one effector gene
        tot_sil_rate+=genotype->protein_pool[genotype->nproteins-1][0][0]*SILENCING;
    else
        tot_sil_rate+=0.0; // otherwise last copy effector gene cannot be deleted either*/
    tot_mut_rate+=tot_sil_rate;    
    /* calc total susbtitution rate*/
    tot_subs_rate=(genotype->ngenes-N_SIGNAL_TF)*CISREG_LEN*SUBSTITUTION; //NA to the sensor TF gene
    tot_mut_rate+=tot_subs_rate;
    /* mut in kinetic constants */    
    tot_mut_kin_rate=(genotype->ngenes-N_SIGNAL_TF)*MUTKINETIC*(proportion_mut_kdis+proportion_mut_mRNA_decay+proportion_mut_protein_decay+proportion_mut_translation_rate); // NA to the sensor TF gene
    tot_mut_kin_rate+=genotype->protein_pool[genotype->nproteins-1][0][0]*MUTKINETIC*proportion_mut_cooperation;
    tot_mut_rate+=tot_mut_kin_rate;
    /* mut in binding seq*/
    tot_mut_binding_seq_rate=genotype->ntfgenes*MUTKINETIC*proportion_mut_binding_seq; // NA to the effector genes
    tot_mut_rate+=tot_mut_binding_seq_rate;    
    /* mut in identity*/
    tot_mut_identity_rate=(genotype->ntfgenes-N_SIGNAL_TF)*MUTKINETIC*proportion_mut_identity; // NA to the sensor TF gene
    tot_mut_rate+=tot_mut_identity_rate;
    /* mut in koff*/
    tot_mut_koff_rate=genotype->ntfgenes*MUTKINETIC*proportion_mut_koff;  
    tot_mut_rate+=tot_mut_koff_rate;    
    /* mut locus length*/
#if MUT_LOCUS_LENGTH
    tot_mut_locus_length_rate=MUT_LOCUS_LENGTH_RATE*genotype->total_loci_length;
#else
    tot_mut_locus_length_rate=0.0;
#endif
    tot_mut_rate+=tot_mut_locus_length_rate;
    /*Draw a mutation based on the above rates*/
    while(1)
    {
        random=RngStream_RandU01(RS);    
        if(random<=tot_subs_rate/tot_mut_rate) 
        {
            *mut_type='s';                      /* substitution*/
            break;
        }
        else
        {   
            random-=tot_subs_rate/tot_mut_rate;
            if(random<= tot_dup_rate/tot_mut_rate)
            {
                *mut_type='d';                          /* gene duplication */
                break;
            }
            else
            {
                random-=tot_dup_rate/tot_mut_rate;
                if(random<=tot_sil_rate/tot_mut_rate)
                {
                    *mut_type='w';                              /* gene deletion*/  
                    break;
                }
                else
                {                  
                    random-=tot_sil_rate/tot_mut_rate; 
                    if(random<=tot_mut_locus_length_rate/tot_mut_rate)
                    {
                        *mut_type='l';      /* mut locus length */
                        break;
                    }                    
                    else
                    {                        
                        random-=tot_mut_locus_length_rate/tot_mut_rate;
                        if(random<=tot_mut_kin_rate/tot_mut_rate)
                        {
                            *mut_type='k';                  /* mut kinetic const*/
                            break;
                        }                        
                        else
                        {
                            random-=tot_mut_kin_rate/tot_mut_rate;
                            if(random<=tot_mut_binding_seq_rate/tot_mut_rate)
                            {
                                *mut_type='c';              /* mut binding seq*/
                                break;
                            }                            
                            else
                            {
                                random-=tot_mut_binding_seq_rate/tot_mut_rate;
                                if(random<=tot_mut_koff_rate/tot_mut_rate)
                                {                                
                                    *mut_type='f';          /* mut koff*/
                                    break;
                                }                               
                                else
                                {                                    
                                    random-=tot_mut_koff_rate/tot_mut_rate;
                                    if(random<=tot_mut_identity_rate/tot_mut_rate)
                                    {
                                        *mut_type='e';           /* mut identity of a TF */
                                        break;
                                    }                                    
                                }
                            }
                        }
                    }                
                }
            }
        }  
    }
}

/*
 *We use which_protein to look for the protein encoded by a given gene copy,
 *and protein_pool to look for genes encoding a given protein. These two tables
 *are updated upon mutations.
 */
void update_protein_pool(Genotype *genotype, int which_protein, int which_gene, char mut_type)
{
    int i, j, protein_id, gene_id, tf_family_id; 
    /*Protein_pool stores the numbers of gene copies that encode a given protein, and the ids of these gene copies.
     *One important thing is that the genes encoding a given protein are not stored by the order of their ids in protein_pool.
     *To delete a gene, which might be the only gene, encoding a given protein, we shift protein_pool to overwrite the to-be-deleted gene
     *We need to update the ids of the remaining genes and proteins
     *For gene duplication, the new gene is always add to the end of the list of genes encoding a given protein.
     *A new protein is also add to the end of protein_pool
     *which_protein can be updated easily. Changing which protein a gene encodes is always easy. For deletion, 
     *we just shift the array to overwrite the to-be-deleted gene and update the ids of the remaining genes.*/ 
    switch (mut_type)
    {
        case 'w':/*a whole gene deletion*/        
            if(genotype->protein_pool[which_protein][0][0]==1) /* if this is the only gene copy,we also need to delete a protein*/
            {   
                /*
                 * UPDATE protein_pool for protein>=which_protein
                 */
                protein_id=which_protein;
                /*shift protein>which_protein to overwrite the to-be-deleted protein*/
                for(i=0;i<genotype->nproteins-which_protein;i++)  
                {   
                    /*reset the portion of protein_pool to be overwritten*/
                    for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)
                        genotype->protein_pool[protein_id][1][j]=NA;
                    /*overwrite*/
                    genotype->protein_pool[protein_id][0][0]=genotype->protein_pool[protein_id+1][0][0]; //this is number of gene copies encoding a protein
                    for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)
                    {
                        gene_id=genotype->protein_pool[protein_id+1][1][j];//these are the gene copies encoding a protein
                        /*note that deletion changes the ids of the remaining genes!!! Any gene that is greater than which_gene is reduced by one*/
                        genotype->protein_pool[protein_id][1][j]=(gene_id>which_gene)?gene_id-1:gene_id;
                    }            
                    protein_id++;
                }
                /*
                 * SPECIAL CASE: if a tf PROTEIN is deleted, we need to change the number of TF proteins and update TF_family_pool
                 */                
                if(which_protein<genotype->nproteins-1) 
                {   
                    /* reduce the number of activators or that of repressors */ 
                    if(genotype->protein_identity[which_protein]==ACTIVATOR) 
                        genotype->N_act--;
                    else if(genotype->protein_identity[which_protein]==REPRESSOR)
                        genotype->N_rep--;
                    /* also remove it from activating and kd */
                    protein_id=which_protein;
                    for(i=0;i<genotype->nproteins-which_protein-1;i++)
                    {
                        genotype->protein_identity[protein_id]=genotype->protein_identity[protein_id+1];
                        genotype->Kd[protein_id]=genotype->Kd[protein_id+1];
                        protein_id++;
                    }  
                    /* update TF_family_pool*/
                    tf_family_id=genotype->which_TF_family[protein_id];
                    if(genotype->TF_family_pool[tf_family_id][0][0]==1) /*there is one member in the family, so delete the family*/
                    {
                        /*shift tf_family_id to overwrite*/
                        for(i=tf_family_id;i<genotype->nTF_families-1;i++)
                        {
                            /*reset entries of family i*/
                            for(j=0;j<genotype->TF_family_pool[i][0][0];j++)
                                genotype->TF_family_pool[i][1][j]=NA;
                            /*copy entries of family i+1 to family i. Proteins with id > protein_id needs a new id */
                            for(j=0;j<genotype->TF_family_pool[i+1][0][0];j++)
                                genotype->TF_family_pool[i][1][j]=(genotype->TF_family_pool[i+1][1][j]>protein_id)?genotype->TF_family_pool[i+1][1][j]-1:genotype->TF_family_pool[i+1][1][j];
                            genotype->TF_family_pool[i][0][0]=genotype->TF_family_pool[i+1][1][j];
                        }
                        /*reset the last family*/
                        for(i=0;i<genotype->TF_family_pool[genotype->nTF_families-1][0][0];i++)
                            genotype->TF_family_pool[genotype->nTF_families-1][1][i]=NA;
                        genotype->TF_family_pool[genotype->nTF_families-1][0][0]=0; 
                        /*update nTF_family*/
                        genotype->nTF_families--;
                        /*shift which_tf_family and update the family id of protein with family id > tf_family_id*/
                        for(i=protein_id;i<genotype->nproteins-1;i++)
                            genotype->which_TF_family[i]=(genotype->which_TF_family[i+1]>tf_family_id)?genotype->which_TF_family[i+1]-1:genotype->which_TF_family[i+1];
                    }
                    else /* the tf family has multiple proteins, no need to delete the family*/
                    {
                        /*find where is protein_id in its family*/
                        for(i=0;i<genotype->TF_family_pool[tf_family_id][0][0];i++)
                        {
                            if(genotype->TF_family_pool[tf_family_id][1][i]==protein_id)
                                break;
                        }
                        /*shift entries in its family*/
                        for(;i<genotype->TF_family_pool[tf_family_id][0][0];i++)
                            genotype->TF_family_pool[tf_family_id][1][i]=genotype->TF_family_pool[tf_family_id][1][i+1];
                        /*update protein ids in all families*/
                        for(i=0;i<genotype->nTF_families;i++)
                        {
                            for(j=0;j<genotype->TF_family_pool[i][0][0];j++)
                                genotype->TF_family_pool[i][1][j]=(genotype->TF_family_pool[i][1][j]>protein_id)?genotype->TF_family_pool[i][1][j]-1:genotype->TF_family_pool[i][1][j];
                        } 
                        /*shift which_tf_family*/
                        for(i=protein_id;i<genotype->nproteins-1;i++)
                            genotype->which_TF_family[i]=genotype->which_TF_family[i+1];
                    }
                    /* in the case, all genes need to recalc binding sites*/
                    for(i=N_SIGNAL_TF;i<which_gene;i++)                    
                        genotype->recalc_TFBS[i]=YES; /* recalc BS */                                       
                }  
                /*
                 * UPDATE which_protein
                 */
                /* update which_protein for gene<which_gene in which_protein*/
                for(i=N_SIGNAL_TF;i<which_gene;i++)
                    genotype->which_protein[i]=(genotype->which_protein[i]<which_protein)?genotype->which_protein[i]:genotype->which_protein[i]-1;//the deletion also changes the ids of proteins
                /* shift and update which_protein for gene>=which_gene in which_protein*/                
                for(i=which_gene;i<genotype->ngenes;i++)
                    genotype->which_protein[i]=(genotype->which_protein[i+1]>which_protein)?genotype->which_protein[i+1]-1:genotype->which_protein[i+1];                   
                /*one less protein*/
                genotype->nproteins--;
            }  
            else /*if the protein has more than one genes*/
            {
                /*
                 * UPDATE protein_pool for protein>=which_protein
                 */
                /* find where is this which_gene*/
                i=0;
                while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++; 
                /*shift protein_pool to overwrite which_gene, and update ids of genes*/
                j=i;
                for(;i<genotype->protein_pool[which_protein][0][0];i++)
                    genotype->protein_pool[which_protein][1][i]=(genotype->protein_pool[which_protein][1][i+1]>which_gene)?genotype->protein_pool[which_protein][1][i+1]-1:genotype->protein_pool[which_protein][1][i+1];// deletion changes the ids of genes!!!
                /*also update the ids for genes before j*/
                for(i=0;i<j;i++)
                    genotype->protein_pool[which_protein][1][i]=(genotype->protein_pool[which_protein][1][i]>which_gene)?genotype->protein_pool[which_protein][1][i]-1:genotype->protein_pool[which_protein][1][i];                            
                /*one less copy encoding which_protein*/
                genotype->protein_pool[which_protein][0][0]--;                
                /*update the ids of genes in protein>which_protein*/                
                for(i=which_protein+1;i<genotype->nproteins;i++)  
                {   
                    for(j=0;j<genotype->protein_pool[i][0][0];j++)
                        genotype->protein_pool[i][1][j]=(genotype->protein_pool[i][1][j]>which_gene)?genotype->protein_pool[i][1][j]-1:genotype->protein_pool[i][1][j];//note that deletion changes the ids of genes!!!
                }
                /*
                 * UPDATE which_protein
                 */
                /*shift which_protein to delete which_gene*/
                for(i=which_gene;i<genotype->ngenes;i++)
                    genotype->which_protein[i]=genotype->which_protein[i+1];                
            }
            /*
             * UPDATE protein_pool for protein<which_protein
             */           
            for(i=N_SIGNAL_TF;i<which_protein;i++)
            {
                for(j=0;j<genotype->protein_pool[i][0][0];j++)
                    genotype->protein_pool[i][1][j]=(genotype->protein_pool[i][1][j]<which_gene)?genotype->protein_pool[i][1][j]:genotype->protein_pool[i][1][j]-1;
            }            
            break;
        case 'd': /*a gene duplication*/
            /* add it to protein_pool, but do not change nproteins*/    
            genotype->protein_pool[which_protein][1][genotype->protein_pool[which_protein][0][0]]=genotype->ngenes-1; //the newly duplicated gene takes the original place of the effector gene
            genotype->protein_pool[which_protein][0][0]++; 
            /*update the id of the original effector gene stored in protein_pool*/           
            j=0;
            while(genotype->protein_pool[genotype->nproteins-1][1][j]!=genotype->ngenes-1)j++;
            genotype->protein_pool[genotype->nproteins-1][1][j]++; 
            /*update which_protein*/          
            genotype->which_protein[genotype->ngenes]=genotype->nproteins-1;//points to the effector protein       
            genotype->which_protein[genotype->ngenes-1]=which_protein;            
            break;
        case 'c': /*mutation in tf binding seq, creating a new tf and a new TF family*/
            /* remove this copy of gene from the original protein_pool*/
            i=0;
            while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++;
            /*shift to delete which_gene*/
            for(;i<genotype->protein_pool[which_protein][0][0];i++) 
                genotype->protein_pool[which_protein][1][i]= genotype->protein_pool[which_protein][1][i+1]; 
            /*one less gene copy to encoding which_protein*/
            genotype->protein_pool[which_protein][0][0]--;
            /* increase the protein id of the effector protein to make room for the new tf*/ 
            genotype->protein_pool[genotype->nproteins][0][0]=genotype->protein_pool[genotype->nproteins-1][0][0];
            for(j=0;j<genotype->protein_pool[genotype->nproteins-1][0][0];j++)                
            {
                genotype->protein_pool[genotype->nproteins][1][j]=genotype->protein_pool[genotype->nproteins-1][1][j]; //shift the protein pool of effector gene
                genotype->which_protein[genotype->protein_pool[genotype->nproteins-1][1][j]]++;//update which_protein
                genotype->protein_pool[genotype->nproteins-1][1][j]=-1; //reset the original protein pool of the effector 
            }                       
            /* create a new protein and link it to which_gene*/
            genotype->which_protein[which_gene]=genotype->nproteins-1; //put the new protein to the original place of the effector protein
            genotype->protein_pool[genotype->nproteins-1][0][0]=1;
            genotype->protein_pool[genotype->nproteins-1][1][0]=which_gene;
            /* update activator or repressor numbers, and activating*/
            if(genotype->protein_identity[which_protein]==ACTIVATOR) //mutation to binding seq does not change the identity of a tf
                genotype->N_act++;
            else if(genotype->protein_identity[which_protein]==REPRESSOR)
                genotype->N_rep++;
            genotype->protein_identity[genotype->nproteins-1]=genotype->protein_identity[which_protein];
            /* update Kd*/
            genotype->Kd[genotype->nproteins-1]=genotype->Kd[which_protein];
            /* make a new tf family and add in the new tf*/
            genotype->TF_family_pool[genotype->nTF_families][0][0]=1;
            genotype->TF_family_pool[genotype->nTF_families][1][0]=genotype->nproteins;
            genotype->which_TF_family[genotype->nproteins]=genotype->nTF_families;
            /* update tf family number*/
            genotype->nTF_families++;
            /* finally, update protein numbers*/
            genotype->nproteins++;            
            break;            
        case 'e': /*mutation in the identity of a TF, creating a new tf and a new tf family*/
            /* remove this copy of gene from the original protein*/
            i=0;
            while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++;
            for(;i<genotype->protein_pool[which_protein][0][0];i++) 
            {
                genotype->protein_pool[which_protein][1][i]= genotype->protein_pool[which_protein][1][i+1]; 
            }
            genotype->protein_pool[which_protein][0][0]--;
            /* increase the protein id of the effector protein*/               
            genotype->protein_pool[genotype->nproteins][0][0]=genotype->protein_pool[genotype->nproteins-1][0][0];
            for(j=0;j<genotype->protein_pool[genotype->nproteins-1][0][0];j++)                
            {
                genotype->protein_pool[genotype->nproteins][1][j]=genotype->protein_pool[genotype->nproteins-1][1][j];
                genotype->which_protein[genotype->protein_pool[genotype->nproteins-1][1][j]]++;
                genotype->protein_pool[genotype->nproteins-1][1][j]=-1;
            }            
            /* create a new protein and link it to which_gene*/
            genotype->which_protein[which_gene]=genotype->nproteins-1; 
            genotype->protein_pool[genotype->nproteins-1][0][0]=1;
            genotype->protein_pool[genotype->nproteins-1][1][0]=which_gene;
            /* update activating*/
            if(genotype->protein_identity[which_protein]==ACTIVATOR) 
                genotype->N_rep++;  /* an activator turns into a repressor */
            else if(genotype->protein_identity[which_protein]==REPRESSOR)
                genotype->N_act++;
            genotype->protein_identity[genotype->nproteins-1]=genotype->protein_identity[which_protein];
            /* update Kd*/
            genotype->Kd[genotype->nproteins-1]=genotype->Kd[which_protein];
            /* make a new tf family and add in the new tf*/
            genotype->TF_family_pool[genotype->nTF_families][0][0]=1;
            genotype->TF_family_pool[genotype->nTF_families][1][0]=genotype->nproteins;
            genotype->which_TF_family[genotype->nproteins]=genotype->nTF_families;
            /* update tf family number*/
            genotype->nTF_families++;
            /* finally, update protein numbers*/
            genotype->nproteins++;           
            break;            
        case 'f': /*mutation in tf kd, creating a new tf, BUT NOT CREATING NEW TF FAMILY*/
            /* remove this copy of gene from the original protein pool*/
            i=0;
            while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++;
            for(;i<genotype->protein_pool[which_protein][0][0];i++) 
            {
                genotype->protein_pool[which_protein][1][i]= genotype->protein_pool[which_protein][1][i+1]; /* rearrange data array */
            }
            genotype->protein_pool[which_protein][0][0]--;
            /* increase the protein id of the effector protein*/               
            genotype->protein_pool[genotype->nproteins][0][0]=genotype->protein_pool[genotype->nproteins-1][0][0];
            for(j=0;j<genotype->protein_pool[genotype->nproteins-1][0][0];j++)                
            {
                genotype->protein_pool[genotype->nproteins][1][j]=genotype->protein_pool[genotype->nproteins-1][1][j];
                genotype->which_protein[genotype->protein_pool[genotype->nproteins-1][1][j]]++;
                genotype->protein_pool[genotype->nproteins-1][1][j]=-1;
            }            
            /* create a new protein and link it to this gene*/
            genotype->which_protein[which_gene]=genotype->nproteins-1; 
            genotype->protein_pool[genotype->nproteins-1][0][0]=1;
            genotype->protein_pool[genotype->nproteins-1][1][0]=which_gene;            
            /* update activating*/
            if(genotype->protein_identity[which_protein]==ACTIVATOR) 
                genotype->N_act++;  
            else if(genotype->protein_identity[which_protein]==REPRESSOR) 
                genotype->N_rep++;
            genotype->protein_identity[genotype->nproteins-1]=genotype->protein_identity[which_protein];  
            /* finally, update protein numbers*/
            genotype->nproteins++;
            /* NOTE: this mutation does not change the number of genes*/
            break;
    }
}


/*To reduce the amount of calculation on the probability of binding distributions, we group gene copies
 *that are created by whole gene duplication. We call such a group a cis-reg cluster because gene copies 
 *in the group should have the same cis-reg sequence. For each cis-reg cluster we only need to calculate
 *the probability of binding distributions once. However, substitutions is cis-reg sequence can create/remove
 *binding sites, therefore we need to check whether a gene copy is still in the original cis-reg cluster 
 *after mutation.We use cisreg_cluster and which_cluster to track the bi-way relation between a gene and 
 *a cis-reg cluster.*/
void update_cisreg_cluster(Genotype *genotype, int which_gene, char mut_type, int new_clusters[NGENES][NGENES], int N_new_clusters, int original_cluster_id)
{
    /*In a cis-reg cluster, gene copies are ordered ascendingly by their ids. There are no empty slots in the list 
     *of gene copies. Empty slots after the list are marked by -1. We do not track the number of gene copies in a
     *cluster. In order to count the number of gene copies correctly, marking the empty slots accurately become important.
     *Therefore, we always set all the slots in a cluster to -1, when deleting or overwrting the cluster. 
     */
    int cisreg_seq_cluster_id, cisreg_seq_cluster_id_copy, i, j, last_cluster; 
    
    if(mut_type=='c')//a mutation to the binding sequence of a TF happened. this can create new clusters.
    {    
        /*reset the original cisreg_cluster*/
        i=0;
        while(genotype->cisreg_cluster[original_cluster_id][i]!=NA)
        {
            genotype->cisreg_cluster[original_cluster_id][i]=NA;
            i++;
        }
        /*assign the first new cluster to the original cluster*/
        i=0;
        while(new_clusters[0][i]!=NA)
        {
            genotype->cisreg_cluster[original_cluster_id][i]=new_clusters[0][i];
            genotype->which_cluster[new_clusters[0][i]]=original_cluster_id;
            i++;
        }
        /*find a empty slot in cisreg_cluster for the rest of new clusters*/    
        last_cluster=N_SIGNAL_TF;
        while(genotype->cisreg_cluster[last_cluster][0]!=NA)last_cluster++;//an empty cluster has all of its slots marked by -1
        /*assign the rest of new clusters into empty slots*/       
        for(i=1;i<N_new_clusters;i++)
        {
            /*reset the empty cluster, just in case*/
            j=0;
            while(genotype->cisreg_cluster[last_cluster][j]!=NA)j++;
            /*assign new clusters*/
            j=0;
            while(new_clusters[i][j]!=NA)
            {
                genotype->cisreg_cluster[last_cluster][j]=new_clusters[i][j];
                genotype->which_cluster[new_clusters[i][j]]=last_cluster;
                j++;
            }
            last_cluster++;
        }        
        return;
    } 
    /*find the cis-reg cluster of which gene*/
    cisreg_seq_cluster_id=genotype->which_cluster[which_gene];
    if(mut_type!='d')/*not a duplication*/
    {   
        if(mut_type!='w') /*not a gene deletion. Substitution (or small indel) just kicked one gene copy out of a cluster*/
        {
            if(genotype->cisreg_cluster[cisreg_seq_cluster_id][1]!=NA) //if a cluster contains more than 1 gene
            {   
                /*then remove which_gene from the original cluster*/        
                i=0;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=which_gene) i++;
                /*shift the slots to overwrite which_gene*/
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=NA)
                {
                    genotype->cisreg_cluster[cisreg_seq_cluster_id][i]=genotype->cisreg_cluster[cisreg_seq_cluster_id][i+1];
                    i++;
                }                     
                /* and create a new cluster*/
                last_cluster=N_SIGNAL_TF; // start from the first non-signal-tf gene
                while(genotype->cisreg_cluster[last_cluster][0]!=NA) last_cluster++; 
                genotype->cisreg_cluster[last_cluster][0]=which_gene;
                genotype->which_cluster[which_gene]=last_cluster;
            }
        }
        else /*is gene deletion*/
        {               
            if(genotype->cisreg_cluster[cisreg_seq_cluster_id][1]!=NA) //if a cluster contains more than 1 gene
            {   
                /*then remove which_gene from the original cluster*/        
                i=0;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=which_gene) i++;
                /*shift to overwrite which_gene, and update ids of the remaining genes*/
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=NA)
                {      
                    /*note that gene copies are ordered ascendingly*/
                    genotype->cisreg_cluster[cisreg_seq_cluster_id][i]=(genotype->cisreg_cluster[cisreg_seq_cluster_id][i+1]==NA)?NA:(genotype->cisreg_cluster[cisreg_seq_cluster_id][i+1]-1); 
                    i++;
                }               
                /*take care of clusters>cisreg_seq_cluster_id*/
                cisreg_seq_cluster_id_copy=cisreg_seq_cluster_id+1;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][0]!=NA)
                {
                    i=0;
                    while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]!=NA)
                    {
                        genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]=(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]<which_gene)?genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]:genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]-1;
                        i++;
                    }
                    cisreg_seq_cluster_id_copy++;
                }              
                /*update which_cluster*/
                for(i=which_gene;i<genotype->ngenes;i++)                              
                    genotype->which_cluster[i]=genotype->which_cluster[i+1];               
            }
            else//if a cluster contains only 1 gene
            {   
                /*need to shift cisreg_cluster*/ 
                cisreg_seq_cluster_id_copy=cisreg_seq_cluster_id;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][0]!=NA)
                {
                    /*reset cluster=cisreg_seq_cluster_id_copy*/
                    i=0;
                    while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]!=NA)
                    {
                        genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]=NA;
                        i++;
                    }
                    /*then copy from cluster=cisreg_seq_cluster_id_copy+1*/                    
                    i=0;
                    while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][i]!=NA)
                    {                        
                        genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]=(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][i]<which_gene)?genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][i]:genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][i]-1; /*note deletion changes gene thread_ID*/
                        i++;
                    }
                    cisreg_seq_cluster_id_copy++;
                }
                /*shift which_cluster and update cluster ids for gene>=which_gene*/                
                for(i=which_gene;i<genotype->ngenes;i++)                                
                    genotype->which_cluster[i]=(genotype->which_cluster[i+1]<cisreg_seq_cluster_id)?genotype->which_cluster[i+1]:genotype->which_cluster[i+1]-1;
                /*update cluster ids for gene<which_gene*/
                for(i=N_SIGNAL_TF;i<which_gene;i++)
                    genotype->which_cluster[i]=(genotype->which_cluster[i]>cisreg_seq_cluster_id)?genotype->which_cluster[i]-1:genotype->which_cluster[i];
            }
            /*update gene ids in clusters<cisreg_seq_cluster_id*/
            for(i=N_SIGNAL_TF;i<cisreg_seq_cluster_id;i++)
            {
                j=0;
                while(genotype->cisreg_cluster[i][j]!=NA)
                {
                    genotype->cisreg_cluster[i][j]=(genotype->cisreg_cluster[i][j]<which_gene)?genotype->cisreg_cluster[i][j]:genotype->cisreg_cluster[i][j]-1;
                    j++;
                }
            }                                            
        }
    }
    else /*gene duplication*/
    {      
        /*Assuming the last gene is in cis-reg cluster X, check whether the duplicated
         *gene is also in cluster X*/
        if(genotype->which_cluster[which_gene]!=genotype->which_cluster[genotype->ngenes-1]) //No
        {    
            /*find an empty slot in cisreg_seq_cluster_id*/
            i=0;
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=NA)i++;
            /*the duplicated gene is always add to ngene-1. Note that ngenes has not been 
             *updated when update_cisreg_cluster is called.*/
            genotype->cisreg_cluster[cisreg_seq_cluster_id][i]=genotype->ngenes-1;            
            /* now find the cis-reg cluster of the effector gene at ngenes-1*/           
            cisreg_seq_cluster_id_copy=genotype->which_cluster[genotype->ngenes-1];
            /*find this effector gene in the cluster*/
            i=0;
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]!=genotype->ngenes-1)i++;
            /*increase the gene id by 1*/
            genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]=genotype->ngenes;
        }
        else //the last gene is duplicated
        {   
            /*simply add ngenes to the cis-reg cluster of the last gene*/
            i=0;
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=NA)i++; //find an empty slot
            genotype->cisreg_cluster[cisreg_seq_cluster_id][i]=genotype->ngenes; //add ngenes to the slot
            genotype->cisreg_cluster[cisreg_seq_cluster_id][i-1]=genotype->ngenes-1; //this is the duplicated gene                     
        }
        /*update which_cluster*/
        genotype->which_cluster[genotype->ngenes]=genotype->which_cluster[genotype->ngenes-1]; 
        genotype->which_cluster[genotype->ngenes-1]=cisreg_seq_cluster_id;
    }   
}

/****************** end of mutation functions *********************************/

