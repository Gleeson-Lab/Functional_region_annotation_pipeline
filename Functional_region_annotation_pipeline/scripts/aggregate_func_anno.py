import os
import sys
import pandas as pd
import subprocess


def read_spliceai_output(filename,colnames):
    '''
    read the output vcf file from running spliceai
    and save it as pandas df
    '''
    p1 = subprocess.Popen(["grep", "-v", "#", filename],stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["awk \'{print $1,$2,$2+length($4)-1,$4,$5,$8}\'"],stdin=p1.stdout,stdout=subprocess.PIPE,shell=True)

    out,err=p2.communicate()
    row_list=out.decode().split("\n")
    row_list_2=[]
    for item in row_list:
        row_list_2.append(item.split(" "))

    #remove empty string at the end of the list
    if row_list[-1] == "":
        row_list_2.pop()
    col_names=colnames
    row_list_2
    df = pd.DataFrame.from_records(row_list_2,columns=col_names)
    #convert all entries in df to string for merging
    df = df.astype(str)
    return df


def aggregate_vcfs(spliceai_vcf, func_score_vcf, region_anno_vcf, GH_anno_vcf, outfile):
    '''
    combine all results from the functional and region annotations
    and save it as the designated vcf outfile
    '''

    common_colnames = ["Chrom","Start","End","Ref","Alt"]
    spliceai_colnames = common_colnames + ["spliceAI"]

    ###convert all input vcfs to pandas dfs
    spliceai_df = read_spliceai_output(spliceai_vcf, spliceai_colnames)

    func_score_df = pd.read_csv(func_score_vcf, sep="\t")
    func_score_df = func_score_df.astype(str)
    #need to change the colname of chrom column so left join works
    func_score_df = func_score_df.rename(columns={"Chr": "Chrom"})

    region_anno_df = pd.read_csv(region_anno_vcf, sep="\t")
    region_anno_df = region_anno_df.astype(str)

    GH_anno_df = pd.read_csv(GH_anno_vcf, sep="\t")
    GH_anno_df = GH_anno_df.astype(str)

    ###join all dfs together and save it as the vcf outfile
    joint_df = spliceai_df.merge(func_score_df, on=common_colnames,how="left")
    joint_df = joint_df.merge(region_anno_df, on=common_colnames,how="left")
    joint_df = joint_df.merge(GH_anno_df, on=common_colnames,how="left")
    joint_df = joint_df.fillna(".")
    joint_df.to_csv(outfile,sep="\t",index=False,header=True)


def main(argv):
    if len(argv) != 6:
        sys.stderr.write("usage :" + argv[0] + " <spliceai_vcf> <func_score_vcf> <region_anno_vcf> <GH_anno_vcf> <vcf_out_path>") 
        sys.exit(2)
    spliceai_vcf = argv[1]
    func_score_vcf = argv[2]
    region_anno_vcf = argv[3]
    GH_anno_vcf = argv[4]
    vcf_out_path = argv[5]
    aggregate_vcfs(spliceai_vcf, func_score_vcf, region_anno_vcf, GH_anno_vcf, vcf_out_path)


main(sys.argv)
