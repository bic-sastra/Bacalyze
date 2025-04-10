# --- Licensing and Dependency Information ---
# This software is licensed under the GNU Affero General Public License v3 (AGPL v3).
# Dependencies:
# - Python: Licensed under the Python Software Foundation License.
# - BWA: Licensed under the GPL License.
# - Fastp: Licensed under the MIT License.
# - SAMtools: Licensed under the MIT License.
# - BCFtools: Licensed under the MIT License.
# - VCFtools: Licensed under the GPL License.
# - Nextflow: Licensed under the Apache License 2.0.
# - pandas: Licensed under the BSD 3-Clause License.
# - scikit-learn: Licensed under the BSD 3-Clause License.
# 
# Users are responsible for complying with the licenses of the dependencies.
# For more details, see the relevant license files for each tool.

FROM debian:bookworm-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install -y --no-install-recommends \
    python3 python3-pip bwa fastp samtools bcftools vcftools curl openjdk-17-jdk-headless git
RUN git clone https://bitbucket.org/genomicepidemiology/resfinder_db/ /RESFINDER_DB/
RUN echo 'export CGE_RESFINDER_RESGENE_DB="/RESFINDER_DB"' >> ~/.bashrc

RUN apt-get autoremove -y
RUN apt-get autoclean -y
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*
    
# Download and install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/

# Use pip to install Python packages
RUN pip3 install --break-system-packages pandas scikit-learn mlxtend biopython resfinder MobileElementFinder && \
    ln -s /usr/bin/python3 /usr/bin/python
RUN apt-get update && apt-get install -y --no-install-recommends ncbi-blast+ && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
RUN pip3 install --break-system-packages seaborn
# Install SPAdes genome assembler
RUN apt-get update && apt-get install -y --no-install-recommends spades && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
CMD ["bash"]
