import subprocess
from Bio import SeqIO

def trim_fastq(input_fastq, output_fastq, trim_start=None, trim_end=None):
    with open(input_fastq, "r") as in_handle, open(output_fastq, "w") as out_handle:
        for record in SeqIO.parse(in_handle, "fastq"):
            # Trim the first 10 nucleotides of the sequence
            if trim_start:
                trimmed_seq = record[trim_start:]
            if trim_end:
                trimmed_seq = record[:trim_end]

            # Write the trimmed record to the output file
            SeqIO.write(trimmed_seq, out_handle, "fastq")


if __name__ == "__main__":
#     # nhej
#     sra_list = ["SRR3192757"]
#     for sra in sra_list:
#         merge_reads = ["/homes8/jingyuan/bioinformatics_tools/FLASH-1.2.11/flash",
#         "./data/nhej/%s/%s_1.fastq"%(sra,sra), "./data/nhej/%s/%s_2.fastq"%(sra,sra),
#         "-o", sra, "-d", "./data/nhej/%s"%sra, "-M", "250", "-O", "-m", "20", "-x", "0.1"]
#         subprocess.call(merge_reads)

#         input_fastq = "./data/nhej/%s/%s.extendedFrags.fastq"%(sra,sra)
#         output_fastq = "./data/nhej/%s/%s_trim.fastq"%(sra,sra)
#         trim_start = 24  # Number of nucleotides to trim
#         trim_fastq(input_fastq, output_fastq, trim_start=trim_start)

#         trim_adapters = ["/homes8/jingyuan/bioinformatics_tools/pTrimmer/pTrimmer-1.3.4",
#         "--seqtype", "single", "--ampfile", "./data/nhej/%s/amplicon.txt"%sra, "--read1",
#         "./data/nhej/%s/%s_trim.fastq"%(sra,sra), "--trim1",
#         "./data/nhej/%s/%s_trim2.fastq"%(sra,sra)]
#         subprocess.call(trim_adapters)
        
 # hgRNA "SRR4842571",
    sra_list = ["SRR4842590"]
    for sra in sra_list:
        input_fastq = "./data/nhej/hgRNA/%s/%s.fastq"%(sra,sra)
        output_fastq = "./data/nhej/hgRNA/%s/%s_trim.fastq"%(sra,sra)
        trim_start = 10  # Number of nucleotides to trim
        trim_fastq(input_fastq, output_fastq, trim_start=trim_start)

        trim_adapters = ["/homes8/jingyuan/bioinformatics_tools/pTrimmer/pTrimmer-1.3.4",
        "--seqtype", "single", "--ampfile", "./data/nhej/hgRNA/%s/amplicon.txt"%sra, "--read1",
        "./data/nhej/hgRNA/%s/%s_trim.fastq"%(sra,sra), "--trim1",
        "./data/nhej/hgRNA/%s/%s_trim2.fastq"%(sra,sra)]
        print(trim_adapters)
        subprocess.call(trim_adapters)
        
        
    # # carlin
    # sra_list = ["SRR11311821","SRR11311824"]
    # for sra in sra_list:
    #     merge_reads = ["/homes8/jingyuan/bioinformatics_tools/FLASH-1.2.11/flash",
    #     "./data/nhej/Carlin/%s_1.fastq"%sra, "./data/nhej/Carlin/%s_2.fastq"%sra,
    #     "-o", sra, "-d", "./data/nhej/%s"%sra, "-M", "250", "-O", "-m", "20", "-x", "0.1"]
    #     subprocess.call(merge_reads)
    # 
    #     input_fastq = "./data/nhej/%s/%s.extendedFrags.fastq"%(sra,sra)
    #     output_fastq = "./data/nhej/%s/%s_trim.fastq"%(sra,sra)
    #     trim_start = 10  # Number of nucleotides to trim
    #     trim_fastq(input_fastq, output_fastq, trim_start=trim_start)

#         trim_adapters = ["/homes8/jingyuan/bioinformatics_tools/pTrimmer/pTrimmer-1.3.4",
#         "--seqtype", "single", "--ampfile", "./data/nhej/%s/amplicon.txt"%sra, "--read1",
#         "./data/nhej/%s/%s_trim.fastq"%(sra,sra), "--trim1",
#         "./data/nhej/%s/%s_trim2.fastq"%(sra,sra)]
#         subprocess.call(trim_adapters)      
        
      
    #hdr      
    # sra_list = ["SRR3192758","SRR3192759"]
    # for sra in sra_list:
    #     merge_reads = ["/homes8/jingyuan/bioinformatics_tools/FLASH-1.2.11/flash",
    #     "./data/hdr/%s/%s_1.fastq"%(sra,sra), "./data/hdr/%s/%s_2.fastq"%(sra,sra),
    #     "-o", sra, "-d", "./data/hdr/%s"%sra, "-M", "250", "-O", "-m", "20", "-x", "0.1"]
    #     subprocess.call(merge_reads)
    # 
    #     input_fastq = "./data/hdr/%s/%s.extendedFrags.fastq"%(sra,sra)
    #     output_fastq = "./data/hdr/%s/%s_trim.fastq"%(sra,sra)
    #     
    # 
    #     trim_start = 20
    #     trim_fastq(input_fastq, output_fastq, trim_start=trim_start)
    # 
    #     trim_adapters = ["/homes8/jingyuan/bioinformatics_tools/pTrimmer/pTrimmer-1.3.4",
    #     "--seqtype", "single", "--ampfile", "./data/hdr/%s/amplicon.txt"%sra, "--read1",
    #     "./data/hdr/%s/%s_trim.fastq"%(sra,sra), "--trim1",
    #     "./data/hdr/%s/%s_trim2.fastq"%(sra,sra)]
    #     subprocess.call(trim_adapters)
    
#     # base editor
#     sra_list = ["SRR3305543","SRR3305544","SRR3305545"]
#     for sra in sra_list:
#         input_fastq = "./data/base_editor/%s.fastq"%sra
#         output_fastq = "./data/base_editor/%s_trim.fastq"%sra
#         trim_start = 4  # Number of nucleotides to trim
#         trim_fastq(input_fastq, output_fastq, trim_start=trim_start)

#         trim_adapters = ["/homes8/jingyuan/bioinformatics_tools/pTrimmer/pTrimmer-1.3.4",
#         "--seqtype", "single", "--ampfile", "./data/base_editor/amplicon.txt", "--read1",
#         "./data/base_editor/%s_trim.fastq"%sra, "--trim1",
#         "./data/base_editor/%s_trim2.fastq"%sra]
#         subprocess.call(trim_adapters)
  


