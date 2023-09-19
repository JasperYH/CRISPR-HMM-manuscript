import subprocess

if __name__ == "__main__":
  fasterq_dump_path = "/homes8/jingyuan/bioinformatics_tools/sratoolkit.3.0.2-centos_linux64/bin/fasterq-dump"
  # sra_list = ["SRR3305543","SRR3305544","SRR3305545","SRR3192757","SRR3192758","SRR3192759"]
  # sra_list = ["SRR11311821","SRR11311824"]
  sra_list = ["SRR7042471"]
  for sra in sra_list:
    subprocess.call([fasterq_dump_path,sra,"-O","/homes8/jingyuan/Projects/CRISPR-HMM_manuscript/data"])
    
