# Setting up filenames here:
from os.path import join
#configfile: "consortium.config.json"

# Files
sample = ["2"]

# Path to files
mutation_path = "/scratch/eknodel/Cancer_Genomics/01_somatic_mutation_calling/gatk_mutect2/"
sequenza_path = "/scratch/eknodel/sequenza/"
vep_path = "/scratch/eknodel/Cancer_Genomics/02_variant_annotation/"

rule all:
    input:
        expand(sequenza_path + "{sample}_sequenza2pyclone.txt", sample=sample),
        expand(sequenza_path+"Patient{sample}_fastclone/scores.csv", sample=sample),
        expand(vep_path+"Patient{sample}/{sample}_matched_transcript_mutation.txt", sample=sample)

rule convert_to_pyclone_format:
    input:
        mut = os.path.join(mutation_path, "{sample}.somatic.filtered.pass.vcf.forpyclone"),
        seg = os.path.join(sequenza_path, "{sample}_segments.txt")
    output:
        fc_input = os.path.join(sequenza_path, "{sample}_sequenza2pyclone.txt")
    params:
        out_dir = sequenza_path,
        sample = sample
    shell:
        """
        python sequenza2pyclone.py {input.mut} {input.seg} {params.sample} {params.out_dir}
        """

rule run_fastclone:
    input:
        fc_input = os.path.join(sequenza_path, "{sample}_sequenza2pyclone.txt")
    output:
        fc_output = os.path.join(sequenza_path, "Patient{sample}_fastclone/scores.csv")
    params:
        out_dir = os.path.join(sequenza_path, "Patient{sample}_fastclone")
    shell:
        """
        rm -r {params.out_dir};
        fastclone load-pyclone prop {input.fc_input} None solve {params.out_dir};
        sed -i '1 s/^.*$/position,score1,score2,score3/' {output.fc_output}
        """

rule match_results:
    input:
        mut = os.path.join(vep_path, "Patient{sample}/GATK_pt{sample}.vep"),
        pep = os.path.join(vep_path, "Patient{sample}/GATK_pt{sample}.peptides.formatted"),
        scores = os.path.join(sequenza_path, "Patient{sample}_fastclone/scores.csv")
    output:
        out = os.path.join(vep_path, "Patient{sample}/{sample}_matched_transcript_mutation.txt")
    params:
        sample = sample,
        out_dir = os.path.join(vep_path, "Patient{sample}")
    shell:
        """
        python match_genename_peptide.py {input.mut} {input.pep} {params.sample} {params.out_dir} {input.scores}
        """
        
