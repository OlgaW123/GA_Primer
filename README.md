What needs to be installed?
# BLAT - (on linux)
How to install BLAT?
Download BLAT
##
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat
Make BLAT Executable
##
    chmod +x blat
Move BLAT to a directory in your PATH eg.
##
    sudo mv blat /usr/local/bin/
Test BLAT
##
    blat

# twoBitToFa
How to install twoBitToFa?
Download twoBitToFa
##
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
Make twoBitToFa Executable
##
    chmod +x twoBitToFa
Add twoBitToFa to PATH eg.
##
    echo 'export PATH="$path/to/twoBitToFa:$PATH"' >> ~/.bashrc
    source ~/.bashrc
Test twoBitToFa
##
    twoBitToFa

Necessary data
# hg38.2bit
Can be found here
##
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
    
# BLASTN - (on linux)
How to install BLASTN?
Download BLASTN
##
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz
Extract file
##
    tar -xzvf ncbi-blast-2.12.0+-x64-linux.tar.gz
Add BLASTN to your PATH eg.
##
    export PATH=$PATH:/path/to/ncbi-blast-2.12.0+/bin
To edit .bashrc
##
    nano ~/.bashrc
    #paste the path at the end
    source ~/.bashrc
Test BLASTN
##
    blastn

# BLAST Database
How to make a database?
Download human genome 
##
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
Extract data
##
    gunzip GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
Make BLAST database
##
    makeblastdb -in GCA_000001405.15_GRCh38_full_analysis_set.fna -dbtype nucl -out human_genome_db


