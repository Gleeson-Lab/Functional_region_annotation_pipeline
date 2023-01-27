import os
import sys
import pandas as pd
import subprocess

annovar = "/home/jis215/tools/annovar/annotate_variation.pl"
annovar_db_dir = sys.argv[3]

def read_avoutput(filename,colnames):
    '''
    process the output from annovar region-based anntations
    and save them as pd dataframes for joining
    '''
   
    p1 = subprocess.Popen(["cut", "--complement", "-f1", filename],stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["sed \'s/;/\\t/\'"], stdin=p1.stdout, stdout=subprocess.PIPE,shell=True)
    p3 = subprocess.Popen(["sed \'s/Score=//\'"], stdin=p2.stdout, stdout=subprocess.PIPE,shell=True)
    p4 = subprocess.Popen(["sed \'s/Name=//\'"], stdin=p3.stdout, stdout=subprocess.PIPE,shell=True)
    p5 = subprocess.Popen(["awk \'{print $3,$4,$5,$6,$7,$1,$2}\'"],stdin=p4.stdout, stdout=subprocess.PIPE,shell=True)

    out,err=p5.communicate()
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
    #remove the name column if histone annotation
    if "Histone" in filename: 
        df = df.drop(columns=colnames[-1])
    #convert all entries in df to string for merging
    df = df.astype(str)
    return df

def region_annotation(sample_id, avinput, avoutput_dir, outfile, *cell_line_args):
    '''
    run annovar region-based annotation 
    to annotate the given avinput with the
    desired databases, then read the avoutput
    into pd dataframes and join them together 
    and save it as the outfile
    default dbs: phastCons100, ENCODE_DNase, ENCODE_TFBS
    arg dbs: ENCODE_histone_modification, specify cell lines in the arg
    '''

    annovar_cmd = [annovar, "-regionanno", "-build","hg19", avinput,
                   annovar_db_dir]
    avinput_colnames=["Chrom","Start","End","Ref","Alt"]
    histone_types = ["H3k27ac","H3k27me3","H3k4me1","H3k4me3"]
    avinput_df = pd.read_csv(avinput, sep="\t",names=avinput_colnames)
    #make avoutput_dir if it does not exist
    if not os.path.exists(avoutput_dir):
        #need to use os.makedirs() to create all intermediate dirs
        os.makedirs(avoutput_dir)
    suffix = avoutput_dir + "/" + sample_id

    #run region_annotation for pCE 100, Dnase and TFBS regions first
    subprocess.run(annovar_cmd + ["-out",suffix,"-dbtype","wgEncodeRegDnaseClusteredV3","-scorecolumn","5"])
    subprocess.run(annovar_cmd + ["-out",suffix,"-dbtype","wgEncodeRegTfbsClusteredV3","-scorecolumn","5"])
    subprocess.run(annovar_cmd + ["-out",suffix,"-dbtype","phastConsElements100way"])

    #then iterate through the *cell_line_args and run the histone annotations
    db_file_name_list = []
    for cell_line in cell_line_args:
        for histone_type in histone_types:
            db_file_name = "wgEncodeBroadHistone" + cell_line + histone_type + "StdPk"
            db_file_name_list.append(db_file_name)
            subprocess.run(annovar_cmd + ["-out",suffix,"-dbtype",db_file_name,"-scorecolumn","5"])

    #then process the output files from running annovar and read them into pd dataframes
    db_file_name_list = db_file_name_list+ ["wgEncodeRegDnaseClusteredV3","wgEncodeRegTfbsClusteredV3","phastConsElements100way"]
    df_list = []
    for db_file_name in db_file_name_list:
        input_filename = suffix + '.hg19_' + db_file_name
        input_colnames = avinput_colnames + [db_file_name+"_Score", db_file_name+"_Name"]
        df = read_avoutput(input_filename,input_colnames)
        df_list.append(df)

    #finally join them to the avinput_df via left join to obtain the aggregated vcf outfile
    joint_df = avinput_df.astype(str)
    for df in df_list:
        joint_df = joint_df.merge(df,on=avinput_colnames,how="left")
    joint_df = joint_df.fillna(".")
    joint_df.to_csv(outfile,sep="\t",index=False,header=True)



def main(argv):
    if len(argv[1:6]) != 5:
        sys.stderr.write("usage :" + argv[0] + "<sample_id> <avinput_path> <annovar_db_dir> <avoutput_dir> <out_vcf_path> <cell_lines>")
        sys.exit(2)
    sample_id = argv[1]
    avinput = argv[2]
    avoutput_dir = argv[4]
    out_vcf_path = argv[5]
    cell_lines = argv[6:]
    region_annotation(sample_id, avinput, avoutput_dir, out_vcf_path, *cell_lines)


main(sys.argv)
