import os
import sys
import pandas as pd
import subprocess


annovar = "/home/jis215/tools/annovar/annotate_variation.pl"
GH_db_dir = sys.argv[3]
GH_bed = "GeneHancer_hg19_v5.12.bed"

#load all GH databases (GH_elements, Gene_association, GH_tissue, GH_TFBS) as pandas dfs
GH_elements = pd.read_csv(os.path.join(GH_db_dir,"GeneHancer_AnnotSV_elements_v5.12.txt"),sep="\t",low_memory=False)
GH_gene_interactions = pd.read_csv(os.path.join(GH_db_dir,"GeneHancer_AnnotSV_gene_association_scores_v5.12.txt"),sep="\t",low_memory=False)
GH_tissue = pd.read_csv(os.path.join(GH_db_dir,"GeneHancer_AnnotSV_tissues_v5.12.txt"),sep="\t",low_memory=False)
GH_TFBS = pd.read_csv(os.path.join(GH_db_dir,"GeneHancer_TFBSs_v5.12.txt"),sep="\t",low_memory=False)

def read_GH_avoutput(filename, colnames):
    '''
    Process the avoutput from annotating avinput
    with the GH_ID bed file and save it as pandas df
    '''

    p1 = subprocess.Popen(["cut", "--complement", "-f1", filename],stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["sed \'s/Name=//\'"], stdin=p1.stdout, stdout=subprocess.PIPE,shell=True)
    p3 = subprocess.Popen(["awk \'{print $2,$3,$4,$5,$6,$1}\'"],stdin=p2.stdout, stdout=subprocess.PIPE,shell=True)

    out,err=p3.communicate()
    row_list=out.decode().split("\n")
    row_list_2=[]
    for item in row_list:
        row_list_2.append(item.split(" "))

    #remove empty string at the end of the list
    if row_list[-1] == "":
        row_list_2.pop()
    col_names=colnames
    df = pd.DataFrame.from_records(row_list_2,columns=col_names)
    #convert all entries in df to string for merging
    df = df.astype(str)
    return df


def GH_annotation(sample_id, avinput, avoutput_dir, outfile):
    '''
    Annotate the avinput file with the GH databases (v5.12)
    and save the annotated output as the vcf outfile
    '''

    suffix = avoutput_dir + "/" + sample_id + "_GH_ID"
    #create avoutput_dir if not present
    if not os.path.exists(avoutput_dir):
        #need to use os.makedirs() to create all intermediate dirs
        os.makedirs(avoutput_dir)
    GH_anno_cmd = [annovar, "-regionanno", "-build","hg19", avinput, GH_db_dir, "-dbtype", "bed", "-bedfile", GH_bed, "-colsWanted", "4","-out",suffix]
    avinput_colnames=["Chrom","Start","End","Ref","Alt"]
    avinput_df = pd.read_csv(avinput, sep="\t",names=avinput_colnames)

    #run the annovar region-based annotation to annotate the avinput with GH_id bedfile (hg19)
    subprocess.run(GH_anno_cmd)

    #then read in the annotated bedfile as pandas df
    input_filename = suffix + ".hg19_bed"
    input_colnames = avinput_colnames + ["GH_ID"]
    GH_ID_df = read_GH_avoutput(input_filename, input_colnames)

    #add all additional info columns to the GH_ID_df per GH_ID
    element_elite_list = []
    reg_type_list = []
    gene_association_list = []
    tissue_info_list = []
    TFBS_info_list = []

    for item in GH_ID_df["GH_ID"]:

        #is_elite and reg_type info
        is_element_elite = GH_elements.loc[GH_elements["GHid"] == item, "is_elite"].values[0]
        reg_type = GH_elements.loc[GH_elements["GHid"] == item, "regulatory_element_type"].values[0]
        element_elite_list.append(is_element_elite)
        reg_type_list.append(reg_type)

        #Gene_association info
        gene_interactions_rows = GH_gene_interactions[GH_gene_interactions["GHid"] == item]
        genes=gene_interactions_rows[['symbol','combined_score','is_elite']].astype(str).apply("/".join,axis=1).tolist()
        gene_interactions_summary = ";".join(genes)
        gene_association_list.append(gene_interactions_summary)

        #Tissue info
        tissue_rows = GH_tissue[GH_tissue["GHid"] == item]
        tissue_rows= tissue_rows.fillna(".")
        tissues = tissue_rows[['source','tissue','category']].astype(str).agg("/".join,axis=1).tolist()
        tissues_summary = ";".join(tissues)
        tissue_info_list.append(tissues_summary)

        #TFBS info
        TFBS_rows = GH_TFBS[GH_TFBS["GHid"] == item]
        if TFBS_rows.empty:
            TFBS_summary = "."
        else:
            TFBS = TFBS_rows[['TF','tissues']].astype(str).agg("/".join,axis=1).tolist()
        TFBS_summary = ";".join(TFBS)
        TFBS_info_list.append(TFBS_summary)

    GH_ID_df["GH_is_element_elite"] = element_elite_list
    GH_ID_df["GH_reg_type"] = reg_type_list
    GH_ID_df["GH_gene_associations"] = gene_association_list
    GH_ID_df["GH_tissues"] = tissue_info_list
    GH_ID_df["GH_TFBS"] = TFBS_info_list

    #finally left join the GH-related columns to the avinput df and save it as the vcf outfile
    joint_df = avinput_df.astype(str)
    joint_df = joint_df.merge(GH_ID_df,on=avinput_colnames,how="left")
    joint_df = joint_df.fillna(".")
    joint_df.to_csv(outfile,sep="\t",index=False,header=True)


def main(argv):
    if len(argv) != 6:
        sys.stderr.write("usage :" + argv[0] + " <sample_id> <avinput_path> <GH_db_dir> <avoutput_dir> <out_vcf_path>")
        sys.exit(2)
    sample_id = argv[1]
    avinput = argv[2]
    avoutput_dir = argv[4]
    out_vcf_path = argv[5]
    GH_annotation(sample_id, avinput, avoutput_dir, out_vcf_path)


main(sys.argv)
