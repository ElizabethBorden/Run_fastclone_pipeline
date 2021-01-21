#Process somatic mutation
import sys
import pandas as pd
import numpy as np

vep_file=sys.argv[1] ## Output file from the variant effects predictor
peptide_file=sys.argv[2] ## Output file from pvac-seq with peptides generated
sample_name=sys.argv[3] ## Name of the sample
out_dir=sys.argv[4] ## Directory where you want the output file
scores_file=sys.argv[5] ## Output file with fastclone scores

print("Matching mutation location with the variant from VEP")
mutation_list=[]
position=[]
variant=[]
f_mutation=open(vep_file,'r')
for ele in f_mutation:
    if ele.startswith('#'):
        continue
    else:
        line=ele.strip().split('\t')
        if line[0]!="chrMT":
            gene_info=line[7].split('missense_variant')
            for unit in gene_info:
                #print(unit)
                if len(unit.split('|')) > 3:
                    chro_pos=line[0]+':'+line[1]
                    mut_id=sample_name+":"+chro_pos
                    position.append(unit.split('|')[13])
                    variant.append(unit.split('|')[14])
                    #print(variant)
                else:
                    continue
            #print(transcript)
            #for ele in transcripts:
                mutation_list.append(mut_id)
        else:
            continue

import pandas as pd 
data_snp=pd.DataFrame()
data_snp["mutation_id"]=mutation_list
data_snp['variant']=[i+j for i,j in zip(position, variant)]

print("Matching the scores from Fastclone with the mutation position")
scores = pd.read_csv(scores_file,header=0,sep=',')
scores_count = scores[['position', 'score1', 'score2', 'score3']]
range_dic={}
for i in range(len(scores_count.position)):
    if scores_count.position[i] not in range_dic.keys():
        range_dic[scores_count.position[i]]=[]
        range_dic[scores_count.position[i]].append([scores_count['score1'][i], scores_count['score2'][i], scores_count['score3'][i]])
    else:
        range_dic[scores_count.position[i]].append([scores_count['score1'][i], scores_count['score2'][i], scores_count['score3'][i]])
score1_list=[]
score2_list=[]
score3_list=[]
for i in range(len(data_snp.mutation_id)):
    mutation_pos=data_snp.mutation_id[i]
    try:
        for ele in range_dic[mutation_pos]:
            score1=ele[0]
            score2=ele[1]
            score3=ele[2]
            score1_list.append(score1)
            score2_list.append(score2)
            score3_list.append(score3)
    except KeyError:
        score1_list.append("N/A")
        score2_list.append("N/A")
        score3_list.append("N/A")    

data_snp["score1"]=score1_list
data_snp["score2"]=score2_list
data_snp["score3"]=score3_list

print("Matching the peptide with the variant from VEP")
peptides = pd.read_csv(peptide_file,header=0,sep=',')
peptide_count = peptides[['Mutation', 'Peptide']]
range_dic={}
for i in range(len(peptide_count.Mutation)):
    if peptide_count.Mutation[i] not in range_dic.keys():
        range_dic[peptide_count.Mutation[i]]=[]
        range_dic[peptide_count.Mutation[i]].append([peptide_count['Peptide'][i]])
    else:
        continue
        #range_dic[peptide_count.Mutation[i]].append([peptide_count['Peptide'][i]])

peptide_list=[]
for i in range(len(data_snp.variant)):
    variant=data_snp.variant[i]
    try:
        for ele in range_dic[variant]:
            #print(ele[0])
            peptide=ele[0]
            peptide_list.append(peptide)
    except KeyError:
            #continue
        peptide_list.append("N/A")

data_snp["peptide"]=peptide_list
data_snp.replace('', np.nan,inplace=True)
data_snp.dropna(axis=0, how='any',inplace=True)
print("Outputing final data")
data_snp.to_csv(out_dir+"/"+sample_name+"_matched_transcript_mutation.txt",header=1,sep="\t",index=0)
