nextflow.enable.dsl=2
params.txt_in = ''
params.ref = ''
params.read1 = ''
params.read2 = ''
params.reads = ''
params.qual = '15'
params.pca = ''
params.model = ''
params.bed = ''
params.outdir = './'
params.mode = ''
params.labels = ''
params.k = '5'
params.annotation = 'no'
params.species = ''
params.fasta_reads = ''
params.assembled = ''
params.association_threshold = '10000'
params.min_support = '0.1'
params.ram_limit = '4' // in GB

process indexReference {
    input:
    path reference
    output:
    path "${reference}*"
    script:
    """
    bwa index -p $reference $reference
    echo "Indexing reference"
    """
}

//process trimming {
//    input:

//    path raw1
//    path raw2
//    output:
//    tuple path("${raw1}_trimmed.fastq"), path("${raw2}_trimmed.fastq")
//    script:
//   """
//   echo "Trimming raw data. Inputs: ${raw1} and ${raw2} Quality: ${params.qual}"
//   fastp -i $raw1 -I $raw2 -o ${raw1}_trimmed.fastq -O ${raw2}_trimmed.fastq -q ${params.qual}
//   """
//

process trimming {
// publishDir params.outdir, mode: 'copy'
    input:
    tuple path(raw1), path(raw2)
    output:
    tuple path("${raw1}_trimmed.fastq"), path("${raw2}_trimmed.fastq")
    script:
    """
    echo "Trimming raw data. Inputs: ${raw1} and ${raw2} Quality: ${params.qual}"
    fastp -i $raw1 -I $raw2 -o ${raw1}_trimmed.fastq -O ${raw2}_trimmed.fastq -q ${params.qual}
    """
}

process trimming_single_end {
// publishDir params.outdir, mode: 'copy'
    input:
    path raw1
    output:
    path "${raw1}_trimmed.fastq"
    script:
    """
    fastp -i $raw1 -o ${raw1}_trimmed.fastq -q ${params.qual}
    """
}

process alignReads {
    input:
    path ref_index
	tuple path(forward), path(reverse)
    path reference
	
    output:
    path "*.sam"

    script:
        """
	echo "Aligning Reads. Reference: ${reference}"
	bwa mem $reference $forward $reverse > ${forward}.sam
        """
}

process alignReads_single_end {
// publishDir params.outdir, mode: 'copy'
    input:
    path ref_index
    path read1
    path reference
    
    output:
    path "*.sam"
    script:
    """
    echo "Aligning Reads. Reference: ${reference}"
    bwa mem -p $reference $read1 > ${read1}.sam
    """
}

process samToBam {
// publishDir params.outdir, mode: 'copy'
    input:
    path sam
    output:
    path "*.bam"
    script:
    """
    samtools view -Sb $sam > ${sam}.bam
    """
}

process sorting {
// publishDir params.outdir, mode: 'copy'
    input:
    path bam
    output:
    path "*.bam"
    script:
    """
    samtools sort $bam -o ${bam}_sorted.bam
    """
}

process variantCalling {
// publishDir params.outdir, mode: 'copy'
    input:
    path sorted
    path reference
    output:
    path "${sorted}*"
    script:
    """
    bcftools mpileup -f $reference $sorted | bcftools call -vm -Oz > ${sorted}.vcf.gz
    """
}

process mutationAnnotation {
    // publishDir params.outdir, mode: 'copy'
	input:
    path vcf
    path bed
    output:
    path "*.vcf"
    script:
    """
	bcftools annotate -a $bed -h <(echo -e "##INFO=<ID=ANN,Number=1,Type=String,Description=\"Annotation\">") -c CHROM,FROM,TO,INFO/ANN $vcf -o ${vcf}_annotated.vcf
    """
}

process SNP {
    input:
    path variant
    output:
    path "${variant}*"
    script:
    """
    echo "${variant}"
    vcftools --gzvcf $variant --minDP 4 --max-missing 0.2 --minQ 30 --recode --recode-INFO-all --out ${variant}.filter
    """
}

process SNPFiltering {
    input:
    path snp
    output:
    path "${snp}*"
    script:
    """
    vcftools --gzvcf $snp --remove-indels --recode --recode-INFO-all --out ${snp}.filter.snps
    """
}

process SNPOutput {
    input:
    path snps
    
    output:
    path "${snps}*"
    
    script:
    """
    bcftools query -f '%CHROM %POS %REF %ALT [%TGT]\n' $snps -o ${snps}.filter.snps.extract.txt
    """
}

process posRefAlt {
    input:
    path snpout
    output:
    path "${snpout}*"
    script:
    """
    cut -d " " -f 2,4 $snpout > ${snpout}.snp
    """
}

process encode {
	// publishDir params.outdir, mode: 'copy'
    input:
    path snp
    path pca
    output:
    path "*.csv"
    script:
    """
    python3 -c "
    import pandas as pd
    import pickle
    import random
    import numpy as np
    from sklearn.decomposition import PCA
    #pd.set_option('future.no_silent_downcasting', True)
    df = pd.read_csv('${snp}', sep=' ', header=None)
    df_transposed = df.transpose()
    new_header = df_transposed.iloc[0]
    df_transposed = df_transposed[1:]
    df_transposed.columns = new_header
    replacement_dict = {'A': 1, 'G': 2, 'C': 3, 'T': 4, 'N': 0}
    df_transposed.replace(replacement_dict, inplace=True)
    df_numeric = df_transposed.apply(pd.to_numeric, errors='coerce')
    df_numeric.fillna(0, inplace=True)
    with open ('${pca}','rb') as f:
        pca = pickle.load(f)
    def decompose(new_data, expected_features):
        if isinstance(new_data, pd.DataFrame):
            new_data = new_data.values
        input_features = new_data.shape[1]
        if input_features > expected_features:
            new_data = new_data[:, :expected_features]
        elif input_features < expected_features:
            padding = np.zeros((new_data.shape[0], expected_features - input_features))
            new_data = np.hstack((new_data, padding))
        return new_data
    n_features = pca.n_features_in_
    processed_df = decompose(df_numeric,n_features)
    df_pca = pca.transform(processed_df)
    df_pca_df = pd.DataFrame(df_pca)
    df_pca_df.to_csv('${snp}.csv',index=False)
    "
    """
}

process predict {
    publishDir params.outdir, mode: 'copy'
    input:
    path new_data
    path model
	path labels
    output:
    path '*.txt'
    script:
    """
    python3 -c "
	import pandas as pd
	import pickle
	from sklearn.svm import SVC
	import numpy as np
	import warnings
	warnings.filterwarnings('ignore')
	new_data = pd.read_csv('${new_data}')
	with open('${model}','rb') as f:
		ml_model = pickle.load(f)
	input_data = new_data.iloc[:]
	pred = ml_model.predict(input_data)
	labels_list = []
	with open('${labels}','r') as f:
		for i in f:
			label = i.strip()
			labels_list.append(label)
	r = 'RESISTANT'
	s = 'SUSCEPTIBLE'
	pred_list = pred[0]
	with open('${new_data}_output.txt','w') as f:
		f.write('The predicted trait(s) for ${new_data} is:\\n')
		for i in range(0,len(pred_list)):
			out_trait = r + '\\n' if pred_list[i] == '0' else s + '\\n'
			antibiotic = labels_list[i]
			out_line = f'{antibiotic} {out_trait}'
			f.write(out_line)
    "
    """
}

process kmer {
    publishDir params.outdir, mode: 'copy'
    input:
    path fastq_file
    output:
    path '*.csv'
    script:
    """
    python3 -c "
    from Bio import SeqIO
    import gzip
    from collections import Counter
    import pandas as pd

    fastq_file_path = '${fastq_file}'
    if fastq_file_path.endswith('.gz'):
        handle = gzip.open(fastq_file_path, 'rt')
    else:
        handle = open(fastq_file_path, 'r')

    fastq_file = SeqIO.parse(handle,'fastq')

    k = int('${params.k}')

    kmer_counts = Counter()
    
    for record in fastq_file:
        sequence = str(record.seq)
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            kmer_counts[kmer] +=1

    df_kmer = pd.DataFrame(kmer_counts.items(), columns=['k-mer', 'Frequency'])
    df_kmer = df_kmer.transpose()
    df_kmer.to_csv(f'{fastq_file_path}_csv.csv',index=False, header=False)
    "
    """
}

process assemble {
    publishDir params.outdir, mode: 'copy'
    input:
    tuple path(read1), path(read2)
    output:
    path "${read1}_contigs.fasta"
    script:
    """
    spades -1 $read1 -2 $read2 -o ./
    mv contigs.fasta ${read1}_contigs.fasta
    """
}

process assemble_single_end {
	 publishDir params.outdir, mode: 'copy'
    input:
    path reads
    output:
    path "${reads}_contigs.fasta"
    script:
    """
    spades -s $reads -o ./
    mv contigs.fasta ${reads}_contigs.fasta
    """
}

process kmer_fasta {
    publishDir params.outdir, mode: 'copy'
    input:
    path fasta_file
    output:
    path '*.csv'
    script:
    """
    python3 -c "
    from Bio import SeqIO
    import gzip
    from collections import Counter
    import pandas as pd

    fasta_file_path = '${fasta_file}'
    if fasta_file_path.endswith('.gz'):
        handle = gzip.open(fasta_file_path, 'rt')
    else:
        handle = open(fasta_file_path, 'r')

    fasta_file = SeqIO.parse(handle,'fasta')

    k = int('${params.k}')

    kmer_counts = Counter()
    
    for record in fasta_file:
        sequence = str(record.seq)
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            kmer_counts[kmer] +=1

    df_kmer = pd.DataFrame(kmer_counts.items(), columns=['k-mer', 'Frequency'])
    df_kmer = df_kmer.transpose()
    df_kmer.to_csv(f'{fasta_file_path}_csv.csv',index=False, header=False)
    "
    """
}

process genomeAnnotation {
       publishDir params.outdir, mode: 'copy'
	input:
	path fasta
	output:
	path "${fasta}*"
	script:
	"""
	source ~/.bashrc
	prokka --outdir ./ --prefix $fasta $fasta --force --cpus 0
	sed -i '/##FASTA/,\$d; /^#/d' ${fasta}.gff
	"""
}

process mge_association {
    publishDir params.outdir, mode: 'copy'
    input:
    path mge
    path gff
    output:
    path "${mge}_mge_association.csv"
    script:
    """
    python3 -c "
    import pandas as pd
    mges = pd.read_csv('${mge}', comment='#')
    genes = pd.read_csv('${gff}', comment='#', sep='\\t', header=None,names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    associations = []
    mges['start'] = mges['start'].astype(int)
    mges['end'] = mges['end'].astype(int)
    genes['start'] = genes['start'].astype(int)
    genes['end'] = genes['end'].astype(int)
    for _, mge in mges.iterrows():
        for _, gene in genes.iterrows():
            if gene['seqid'] == mge['contig']:  # Check if they are on the same contig
                if abs(gene['start'] - mge['start']) <= 31000 or abs(gene['end'] - mge['end']) <= 31000:
                    associations.append({
                        'MGE_id': mge['name'],
                            'MGE_type': mge['type'],
                        'Gene': gene['attributes'],
                        'MGE_START': mge['start'],
                        'MGE_END': mge['end'],
                        'GENE_START': gene['start'],
                        'GENE_END': gene['end'],
                        'Distance': min([abs(gene['start'] - mge['start']), abs(gene['end'] - mge['end'])])
                    })
    associations_df = pd.DataFrame(associations)
    associations_df.to_csv('${mge}_mge_association.csv', index=False)
    "
    """
}

process mgefind {
    publishDir "${params.outdir}/mge", mode: 'copy', mkdirs: true
    errorStrategy 'ignore'
    //debug true
    input:
    path fasta
    output:
    path "${fasta}_mge.csv", optional: true
    
    script:
    """
    #!/bin/bash
    rm -rf /tmp/mge_finder/$fasta
    mefinder find -c $fasta ${fasta}_mge
    if [\$? -ne 0]; then
        echo "Error: MGEFinder failed to run for ${fasta}"
    fi
    """
}

process argFind {

     publishDir "${params.outdir}/arg", mode: 'copy', mkdirs: true
	errorStrategy 'ignore'
//     maxRetries 3
	//debug true
    input:
    path fasta
    output:
    path "${fasta}_results_tab.tsv"
    script:
    """
    #!/bin/bash -ux
    
    if [ -f ~/.bashrc ]; then
        source ~/.bashrc
    fi
    echo "PATH inside container: \$PATH"
    blastn_path=\$(which blastn)
    echo \$blastn_path
#    while true; do sleep 1000; done
    python3 -m resfinder -ifa $fasta -s \"${params.species}\" -acq -o ./ -b "\$blastn_path"
    mv *results_tab.txt ${fasta}_results_tab.tsv
    python3 -c"
import os
os.makedirs('${params.outdir}/arg/', exist_ok=True)
    "
	cp pheno_table.txt ${params.outdir}/arg/${fasta}_pheno.tsv
    """
}

process argAssociation {
    publishDir "${params.outdir}/associations", mode: 'copy', mkdirs: true
	errorStrategy 'ignore'
    debug true
    //cache true
    input:
    path arg
    path mge
    output:
    path "${mge}_arg_association.csv"
    script:
    """
    python3 -c "
    import pandas as pd
    import os
    import glob
    import re
    mge_basename = os.path.basename('${mge}')
    match = re.match(r'(.+)_mge.csv', mge_basename)

    if not match:
        print(f'Error: Could not extract prefix from {mge_basename}')
        exit(1)

    mge_prefix = match.group(1)  # Extracted prefix

    # Find matching ARG file
    matching_arg_files = glob.glob(f'${params.outdir}/arg/*{mge_prefix}*_results_tab.tsv')

    if not matching_arg_files:
        print(f'No matching ARG file found for {mge_prefix}')
        exit(0)

    arg_file = matching_arg_files[0]  # Use the first match
    if os.stat('${mge}').st_size == 0:
        exit(0)
    args = pd.read_csv(arg_file, sep='\\t')
    mges = pd.read_csv('${mge}', comment='#')
    args.columns = args.columns.str.strip()
    args['contig'] = args['Contig'].str.split(' ').str[0]
    args.drop(columns=['Contig'], inplace=True)
    mges['contig'] = mges['contig'].str.split(' ').str[0]
    #print(f'Processing ARG and MGE data for {mge_prefix}')
    if args.empty or mges.empty:
        exit(0)
    associations = []
    args = args[['Resistance gene', 'Position in contig','contig']]
    args[['start', 'end']] = args['Position in contig'].str.split('\\.\\.', expand=True, n=1)
    args['start'] = args['start'].astype(int)
    args['end'] = args['end'].astype(int)
    args.drop(columns=['Position in contig'], inplace=True)
    mges['start'] = mges['start'].astype(int)
    mges['end'] = mges['end'].astype(int)
    for _, arg in args.iterrows():
        for _, mge in mges.iterrows():
            if (abs(arg['start'] - mge['start']) <= ${params.association_threshold} or abs(arg['end'] - mge['end']) <= ${params.association_threshold}) and arg['contig'] == mge['contig']:

                associations.append({
                    'ARG': arg['Resistance gene'],
                    'MGE': mge['name'],
                    'MGE_type': mge['type'],
                    'MGE_START': mge['start'],
                    'MGE_END': mge['end'],
                    'ARG_START': arg['start'],
                    'ARG_END': arg['end'],
                    'CONTIG': arg['contig'],
                    'Distance': min([abs(arg['start'] - mge['start']), abs(arg['end'] - mge['end'])])
                })
    associations_df = pd.DataFrame(associations)
    associations_df.to_csv('${mge}_arg_association.csv', index=False)
    "
    """
}

process csvCombine {
    cache false
    debug true
    input:
    path mge_files
    path arg_association_files

    output:
    path "*"

    script:
    """
    python3 -c "
import os
print('Imported os')
import pandas as pd
print('Imported pandas')
from mlxtend.frequent_patterns import apriori, association_rules
from mlxtend.preprocessing import TransactionEncoder
print('Imported mlxtend')

import resource
print('Imported resource')
# Limit RAM usage to 4 GB
soft, hard = resource.getrlimit(resource.RLIMIT_AS)
resource.setrlimit(resource.RLIMIT_AS, (${params.ram_limit} * 1024 * 1024 * 1024, hard))
print(f'Set resource limit to ${params.ram_limit} GB')
csv_dir = '${params.outdir}/mge/'
mge_data = {}

for filename in os.listdir(csv_dir):
    if filename.endswith('mge.csv'):
        file_path = os.path.join(csv_dir, filename)
        
        with open(file_path, 'r') as file:
            genome_name = file.readlines()[1].strip().replace('#sample:', '').strip()
        
    
        df = pd.read_csv(file_path, comment='#')
        
        for mge_no, mge_type in zip(df['name'].unique(), df['type']):
            if mge_no not in mge_data:
                mge_data[mge_no] = {'NO_OF_GENOMES': 0, 'GENOMES': set(), 'TYPE': mge_type}
            mge_data[mge_no]['NO_OF_GENOMES'] += 1
            mge_data[mge_no]['GENOMES'].add(genome_name)

output_df = pd.DataFrame([
    {'MGE': mge_no, 'NO_OF_GENOMES': data['NO_OF_GENOMES'], 'TYPE': data['TYPE'], 'GENOMES': ','.join(data['GENOMES'])}
    for mge_no, data in mge_data.items()
])

output_df.to_csv('combined_mge.csv', index=False)
import matplotlib.pyplot as plt
# Sort the DataFrame by 'NO_OF_GENOMES' in descending order
output_df = output_df.sort_values(by='NO_OF_GENOMES', ascending=False)

# Generate bar graph with hues of red, with higher values having darker hues
plt.figure(figsize=(10, 6))
norm = plt.Normalize(output_df['NO_OF_GENOMES'].min(), output_df['NO_OF_GENOMES'].max())
colors = plt.cm.Reds(norm(output_df['NO_OF_GENOMES']))
bars = plt.bar(output_df['MGE'], output_df['NO_OF_GENOMES'], color=colors)
plt.xlabel('MGE')
plt.ylabel('Number of Genomes')
plt.title('Number of Genomes per MGE')
plt.xticks(rotation=90)
plt.tight_layout()

    # Add legend
legend_labels = [f'{mge}: {count}' for mge, count in zip(output_df['MGE'], output_df['NO_OF_GENOMES'])]
plt.legend(bars, legend_labels, title='MGE Counts', bbox_to_anchor=(1.05, 1), loc='upper left')

# Save the bar graph
plt.savefig('${params.outdir}/mge_counts_bargraph_with_legend.png', bbox_inches='tight')
print('Saved bar graph with legend for MGE counts')

#plt.savefig('mge_counts_bargraph_with_legend.png', bbox_inches='tight')

genome_mge_data = {}

for filename in os.listdir(csv_dir):
    if filename.endswith('mge.csv'):
        file_path = os.path.join(csv_dir, filename)
        
        with open(file_path, 'r') as file:
            genome_name = file.readlines()[1].strip().replace('#sample:', '').strip()
        
        df = pd.read_csv(file_path, comment='#')
        
        if genome_name not in genome_mge_data:
            genome_mge_data[genome_name] = set()
        
        for mge_no in df['name'].unique():
            genome_mge_data[genome_name].add(mge_no)

genome_mge_df = pd.DataFrame([
    {'GENOME': genome, 'MGES': ','.join(mges)}
    for genome, mges in genome_mge_data.items()
])

#print(genome_mge_df)

genome_mge_df.to_csv('genome_mge_output.csv', index=False)
#print('Written genome IDs with their MGEs to CSV file')
# Prepare the data for the apriori algorithm
genome_mge_list = genome_mge_df['MGES'].apply(lambda x: x.split(',')).tolist()

if genome_mge_df.shape[0] != 1:

    te = TransactionEncoder()
    te_ary = te.fit(genome_mge_list).transform(genome_mge_list)
    df = pd.DataFrame(te_ary, columns=te.columns_)

    df.to_csv('${params.outdir}/transaction_data.csv', index=False)
# Run the apriori algorithm
    print('Running Apriori algorithm')
    frequent_itemsets = apriori(df, min_support=${params.min_support}, use_colnames=True, low_memory=True)
    print('Apriori algorithm completed')
    frequent_itemsets.to_csv('${params.outdir}/frequent_itemsets.csv', index=False)
    if frequent_itemsets.shape[0] <= 13087:
        print('Generating association rules')
        rules = association_rules(frequent_itemsets, metric='lift', min_threshold=1)
        print('Association rules generated')

    else:
        frequent_itemsets = frequent_itemsets.sort_values(by='support', ascending=False).head(13087)
        rules = association_rules(frequent_itemsets, metric='lift', min_threshold=1,support_only=True)
        with open('${params.outdir}/warning_note.txt', 'w') as f:
            f.write('The number of frequent itemsets is too large to generate association rules. Only the top 13087 frequent itemsets were used to generate association rules.')
        print('The number of frequent itemsets is too large to generate association rules. Only the top 13087 frequent itemsets were used to generate association rules.')
        print('Association rules generated')

# Print the frequent itemsets and the association rules

# Save the frequent itemsets and the association rules to CSV files
#frequent_itemsets.to_csv('${params.outdir}/frequent_itemsets.csv', index=False)
    rules.to_csv('${params.outdir}/association_rules.csv', index=False)
    rules['antecedents'] = rules['antecedents'].apply(lambda x: ', '.join(list(x)))
    rules['consequents'] = rules['consequents'].apply(lambda x: ', '.join(list(x)))

# Save the readable association rules to a CSV file
    rules.to_csv('${params.outdir}/readable_association_rules.csv', index=False)
    print('Written association rules to CSV file')

else:
    with open('${params.outdir}/warning_note.txt', 'w') as f:
        f.write('Only one genome was found in the dataset. Association rules cannot be generated.')
    print('Only one genome was found in the dataset. Association rules cannot be generated.')
    frequent_itemsets = None
    rules = None

genome_mge_df['NUMBER_OF_MGES'] = genome_mge_df['MGES'].apply(lambda x: len(x.split(',')))
#print(genome_mge_df)
genome_mge_df.to_csv('${params.outdir}/genome_mge_output.csv', index=False)
print('Written genome IDs with their MGEs to CSV file')
# Calculate the number of MGEs per genome
num_mges = genome_mge_df['NUMBER_OF_MGES']

plt.figure(figsize=(10, 6))
box = plt.boxplot(num_mges, vert=True, patch_artist=True, notch=True, showfliers=True)

plt.xlabel('Genomes')
plt.ylabel('Number of MGEs')
plt.title('Box Plot of Number of MGEs per Genome')

# Highlight the median
for median in box['medians']:
    median.set(color='red', linewidth=2)

# Highlight the whiskers
for whisker in box['whiskers']:
    whisker.set(color='blue', linewidth=1.5)

# Highlight the boxes
for box_patch in box['boxes']:
    box_patch.set(facecolor='lightblue')

# Add labels for median and whiskers
med = box['medians'][0].get_ydata()[0]
plt.text(1.1, med, f'Median: {med:.2f}', verticalalignment='center', color='red')

whisker_low = box['whiskers'][0].get_ydata()[1]
plt.text(1.1, whisker_low, f'Lower Whisker: {whisker_low:.2f}', verticalalignment='center', color='blue')

whisker_high = box['whiskers'][1].get_ydata()[1]
plt.text(1.1, whisker_high, f'Upper Whisker: {whisker_high:.2f}', verticalalignment='center', color='blue')

plt.tight_layout()
plt.savefig('${params.outdir}/mge_boxplot.png', bbox_inches='tight')
print('Saved box plot of MGEs per genome')

# Calculate the genome(s) with the lowest and highest number of MGEs
min_mges = genome_mge_df['NUMBER_OF_MGES'].min()
max_mges = genome_mge_df['NUMBER_OF_MGES'].max()
avg_mges = genome_mge_df['NUMBER_OF_MGES'].mean()

genomes_with_min_mges = genome_mge_df[genome_mge_df['NUMBER_OF_MGES'] == min_mges]['GENOME'].tolist()
genomes_with_max_mges = genome_mge_df[genome_mge_df['NUMBER_OF_MGES'] == max_mges]['GENOME'].tolist()

# Create a new DataFrame with the results
summary_df = pd.DataFrame({
    'Metric': ['Lowest Number of MGEs', 'Highest Number of MGEs', 'Average Number of MGEs'],
    'Value': [min_mges, max_mges, avg_mges],
    'Genomes': [', '.join(genomes_with_min_mges), ', '.join(genomes_with_max_mges), 'N/A']
})

# Write the summary DataFrame to a CSV file
summary_df.to_csv('${params.outdir}/genome_mge_summary.csv', index=False)
print('Written summary of genome and MGEs to CSV file')

# Create a new DataFrame with the number of MGEs vs number of genomes with that many MGEs
mge_counts = genome_mge_df['NUMBER_OF_MGES'].value_counts().reset_index()
mge_counts.columns = ['NUMBER_OF_MGES', 'NUMBER_OF_GENOMES']

# Write the new DataFrame to a CSV file
mge_counts.to_csv('${params.outdir}/mge_counts_vs_genomes.csv', index=False)
print('Written MGE counts vs genomes to CSV file')
# Visualize the number of MGEs vs number of genomes
plt.figure(figsize=(10, 6))
norm = plt.Normalize(mge_counts['NUMBER_OF_GENOMES'].min(), mge_counts['NUMBER_OF_GENOMES'].max())
colors = plt.cm.Greens(norm(mge_counts['NUMBER_OF_GENOMES']))
bars = plt.bar(mge_counts['NUMBER_OF_MGES'], mge_counts['NUMBER_OF_GENOMES'], color=colors)
plt.xlabel('Number of MGEs')
plt.ylabel('Number of Genomes')
plt.title('Number of MGEs vs Number of Genomes')
plt.xticks(rotation=90)
plt.tight_layout()

# Create a legend with color patches
import matplotlib.patches as mpatches
legend_patches = [mpatches.Patch(color=plt.cm.Greens(norm(count)), label=f'{mge}: {count}') 
                  for mge, count in zip(mge_counts['NUMBER_OF_MGES'], mge_counts['NUMBER_OF_GENOMES'])]
plt.legend(handles=legend_patches, title='MGE Counts', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.savefig('${params.outdir}/mge_counts_vs_genomes.png', bbox_inches='tight')
print('Saved bar graph of MGE counts vs genomes')
type_counts = output_df['TYPE'].value_counts().reset_index()
type_counts.columns = ['TYPE', 'COUNT']
type_counts['TYPE'] = type_counts['TYPE'].apply(lambda x: x.upper() if x == 'mite' else x.capitalize())
# Write the new DataFrame to a CSV file
type_counts.to_csv('${params.outdir}/mge_type_counts.csv', index=False)
print('Written MGE type counts to CSV file')

# Visualize the counts of each type of MGE as a bar graph
plt.figure(figsize=(10, 6))
plt.bar(type_counts['TYPE'], type_counts['COUNT'], color='skyblue')
plt.xlabel('Type of MGE')
plt.ylabel('Count')
plt.title('Count of Each Type of MGE')
plt.xticks(rotation=90)
plt.tight_layout()
norm = plt.Normalize(type_counts['COUNT'].min(), type_counts['COUNT'].max())
colors = plt.cm.Blues(norm(type_counts['COUNT']))
bars = plt.bar(type_counts['TYPE'], type_counts['COUNT'], color=colors)

# Create a legend with color patches
legend_patches = [mpatches.Patch(color=plt.cm.Blues(norm(count)), label=f'{type_}: {count}') 
                  for type_, count in zip(type_counts['TYPE'], type_counts['COUNT'])]
plt.legend(handles=legend_patches, title='Type Counts', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig('${params.outdir}/type_counts_bargraph.png', bbox_inches='tight')
print('Saved bar graph of MGE type counts')
print('Completed processing of MGE data')

    "
    python3 -c "
print('Starting ARG data processing')
import os
print('Imported os')
import pandas as pd
print('Imported pandas')
import seaborn as sns
print('Imported seaborn')
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
print('Imported matplotlib')
# Directory containing the CSV files
csv_dir = '${params.outdir}/associations/'

# Initialize a dictionary to store ARG data
arg_data = {}
print('initialized arg_data')
# Iterate over each CSV file in the directory
try:
    for filename in os.listdir(csv_dir):
        if filename.endswith('arg_association.csv'):
            # Extract genome name from the filename
            genome_name = filename.split('.fna')[0]
        
            file_path = os.path.join(csv_dir, filename)

            try:        
            
                df = pd.read_csv(file_path)
            except Exception as e:
                print(f'Error: {e}')
                continue
        
            for index, row in df.iterrows():
                if row['Distance'] <= ${params.association_threshold}:
                    arg = row['ARG']
                    mge = row['MGE']
                    distance = row['Distance']
                    key = (arg, mge)
                
                    if key not in arg_data:
                        arg_data[key] = {'DISTANCE': [], 'GENOMES': set()}
                    arg_data[key]['DISTANCE'].append(str(distance))
                    arg_data[key]['GENOMES'].add(genome_name)

    # Convert the dictionary to a DataFrame
    output_data = [
        {
            'ARG': key[0], 
            'MGE': key[1], 
            'DISTANCE': ','.join(data['DISTANCE']), 
            'NO_OF_GENOMES': len(data['GENOMES']), 
            'GENOMES': ','.join(data['GENOMES'])
        }
        for key, data in arg_data.items()
    ]
    output_df = pd.DataFrame(output_data)

    # Sort the DataFrame in descending order of NO_OF_GENOMES
    output_df = output_df.sort_values(by='NO_OF_GENOMES', ascending=False)
    #output_df.to_csv('${params.outdir}/combined_association.csv', index=False)
    mge_csv_dir = '${params.outdir}/mge'
    arg_csv_dir = '${params.outdir}/arg'
    # Load all TSV files from the ../arg directory into a single DataFrame by concatenation
except Exception as e:
    output_df = pd.DataFrame()
    print(e)
    pass
mge_csv_dir = '${params.outdir}/mge'
arg_csv_dir = '${params.outdir}/arg'
arg_df_list = []
for filename in os.listdir(arg_csv_dir):
    if filename.endswith('tab.tsv'):
        file_path = os.path.join(arg_csv_dir, filename)
        df = pd.read_csv(file_path, sep='\t', usecols=['Resistance gene', 'Phenotype'])
        arg_df_list.append(df)

# Concatenate all DataFrames
merged_arg_df = pd.concat(arg_df_list, ignore_index=True)
# Drop all columns except 'Resistance gene' and 'Phenotype'
merged_arg_df = merged_arg_df[['Resistance gene', 'Phenotype']]

# Group by 'Resistance gene' and concatenate 'Phenotype' values
merged_arg_df = merged_arg_df.groupby('Resistance gene')['Phenotype'].apply(lambda x: ','.join(x)).reset_index()
# Remove duplicate values from within each row separated by comma
merged_arg_df['Phenotype'] = merged_arg_df['Phenotype'].apply(lambda x: ','.join(sorted(set(x.split(',')))))
#Rename 'Resistance gene' to 'ARG'
merged_arg_df = merged_arg_df.rename(columns={'Resistance gene': 'ARG'})
#print(merged_arg_df)
# Merge the two DataFrames on 'ARG'
if not output_df.empty:
    final_df = pd.merge(output_df, merged_arg_df, on='ARG', how='left')
#print(final_df)
    final_df.to_csv('${params.outdir}/combined_association.csv', index=False)
# Initialize a new dictionary to store the additional data
new_data = {}

# Iterate over each CSV file in the directory again
for filename in os.listdir(csv_dir):
    if filename.endswith('arg_association.csv'):
        # Extract genome name from the filename
        genome_name = filename.split('.fna')[0]
        file_path = os.path.join(csv_dir, filename)
        
    # Read the CSV file into a DataFrame
        
        try:
            df = pd.read_csv(file_path)
        except Exception as e:
            print(f'Error: {e}')
            continue

        

        for index, row in df.iterrows():
            if row['Distance'] <= ${params.association_threshold}:
                arg = row['ARG']
                mge = row['MGE']
                arg_start = row['ARG_START']
                arg_end = row['ARG_END']
                mge_start = row['MGE_START']
                mge_end = row['MGE_END']
                key = (arg, mge, arg_start, arg_end, mge_start, mge_end)

                if key not in new_data:
                    new_data[key] = {'GENOMES': set()}
                new_data[key]['GENOMES'].add(genome_name)

# Convert the new dictionary to a DataFrame
new_output_data = [
    {
        'ARG': key[0], 
        'MGE': key[1], 
        'ARG_START': key[2], 
        'ARG_END': key[3], 
        'MGE_START': key[4], 
        'MGE_END': key[5], 
        'GENOMES': ','.join(data['GENOMES'])
    }
    for key, data in new_data.items()
]
#print('Converting')
new_output_df = pd.DataFrame(new_output_data)
#print('Converted')
    # Merge the new DataFrame with the final_df to get the 'Phenotype' column
new_final_df = pd.merge(new_output_df, final_df[['ARG', 'Phenotype']], on='ARG', how='left')

# Save the new DataFrame to a CSV file
#new_final_df.to_csv('new_combined_output.csv', index=False)
# Function to check if ARG and MGE overlap
def check_overlap(row):
    return (row['ARG_START'] >= row['MGE_START'] and row['ARG_START'] <= row['MGE_END']) or \
           (row['ARG_END'] >= row['MGE_START'] and row['ARG_END'] <= row['MGE_END'])

# Apply the function to create the 'ON_MGE' column
new_final_df['OVERLAP'] = new_final_df.apply(check_overlap, axis=1)
# Strip white spaces from 'Phenotype' column and explode with comma delimiter
new_final_df['Phenotype'] = new_final_df['Phenotype'].str.replace(' ', '')
new_final_df = new_final_df.assign(Phenotype=new_final_df['Phenotype'].str.split(',')).explode('Phenotype')
# Reset the index of the DataFrame
new_final_df.reset_index(drop=True, inplace=True)
# Save the updated DataFrame to a CSV file
new_final_df.to_csv('${params.outdir}/combined_association_with_overlap.csv', index=False)

# Visualize a single bar plot with a part of the bar highlighted according to the number of overlaps

# Group by 'Phenotype' and count the total occurrences and overlaps
plot_data = new_final_df.groupby('Phenotype').agg(
    Total=('Phenotype', 'size'),
    Overlaps=('OVERLAP', 'sum')
).reset_index()

# Create a bar plot
plt.figure(figsize=(12, 8))
bar_plot = sns.barplot(x='Phenotype', y='Total', data=plot_data, color='blue', label='Total')

# Highlight the part of the bar that corresponds to overlaps
for index, row in plot_data.iterrows():
    bar_plot.patches[index].set_color('blue')
#bar_plot.patches[index].set_edgecolor('black')
    #bar_plot.patches[index].set_hatch('//')

# Add the overlaps as a part of the bar
bar_plot = sns.barplot(x='Phenotype', y='Overlaps', data=plot_data, color='green', label='Overlaps')

# Set plot labels and title
plt.xlabel('Phenotype')
plt.ylabel('Count')
plt.title('Total and Overlaps Grouped by Phenotype')
plt.xticks(rotation=90)
# Add legend with colors

legend_patches = [
    Patch(color='blue', label='Near MGEs'),
    Patch(color='green', label='Overlaps')
]
plt.legend(handles=legend_patches, title='Legend')

# Update the labels for the plot
bar_plot.set_xlabel('Antibiotics')
bar_plot.set_ylabel('Number of ARGs')
bar_plot.set_title('MGEs Grouped by Antibiotics')

# Show the plot
plt.tight_layout()
plt.savefig('${params.outdir}/ARG_grouped_by_antibiotics.png')
# Save the plot data to a CSV file
plot_data.to_csv('${params.outdir}/ARG_grouped_by_antibiotics.csv', index=False)
print('Saved bar graph of ARGs grouped by antibiotics')



tsv_dir = '${params.outdir}/arg'

# Initialize a list to store DataFrames
df_list = []
genome_list = []

# Iterate over each TSV file in the directory
for filename in os.listdir(tsv_dir):
    if filename.endswith('pheno.tsv'):
        file_path = os.path.join(tsv_dir, filename)
        
        # Read the TSV file into a DataFrame, skipping the first 16 lines
        df = pd.read_csv(file_path, sep='\t', skiprows=16)
        
        # Remove '# ' from the first column name
        df.columns = [col.lstrip('# ') for col in df.columns]
        #print(filename)

        genome_id = filename.split('.fna')[0]
        #print(genome_id)
        genome_list.append(genome_id)
        #print(genome_list)
        
        # Convert 'Genetic background' column to string
        df.dropna(subset=['Genetic background'], inplace=True)
        df['Genetic background'] = df['Genetic background'].astype(str)
        #df.to_csv(f'{genome_id}_test.csv', index=False)

        # Append the DataFrame to the list
        df_list.append(df)


#    df.to_csv(f'{genome_id}_counted.csv', index=False)



#print('hi')
# Iterate over each DataFrame in the list
#print(genome_list)
class_counts_list = []
import os

# Create the directory if it does not exist
output_dir = '${params.outdir}/genes_by_class'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
for genome_id, df in zip(genome_list, df_list):
    df['Genetic background'] = df.groupby('Class')['Genetic background'].transform(lambda x: ','.join(x))
    df = df.drop_duplicates(subset='Class')
    df['Genetic background'] = df['Genetic background'].apply(lambda x: ','.join(set(x.split(','))))
#    df.to_csv(f'{genome_id}_grouped.csv', index=False)
    df['Num_Genes'] = df['Genetic background'].apply(lambda x: len(x.replace(' ', '').split(',')))
    # Group by 'Class' and count the number of genes
    class_counts = df[['Class', 'Num_Genes']]
    class_counts = class_counts.sort_values(by='Num_Genes',ascending=False).reset_index(drop=True)
    #print(class_counts)
    # Drop rows where 'Class' column has the value 'under_development'
    class_counts = class_counts[class_counts['Class'] != 'under_development']
    class_counts_list.append(class_counts)
    # Create a bar plot
    plt.figure(figsize=(10, 6))
    # Normalize the number of genes for color mapping
    norm = plt.Normalize(class_counts['Num_Genes'].min(), class_counts['Num_Genes'].max())
    colors = plt.cm.Oranges(norm(class_counts['Num_Genes']))
    
    # Create a bar plot with the color mapping
    bars = plt.bar(class_counts['Class'], class_counts['Num_Genes'], color=colors)
    
    # Set plot labels and title
    plt.xlabel('Class')
    plt.ylabel('Number of Genes')
    plt.title('Number of Genes Grouped by Class')
    plt.xticks(rotation=90)
    plt.tight_layout()
    
    # Create legend labels
    legend_labels = [f'{cls}: {count}' for cls, count in zip(class_counts['Class'], class_counts['Num_Genes'])]
    bars_list = list(bars)
    plt.legend(bars_list, legend_labels, title='ARGs by Class', loc='best',ncol=2)
    # Save the plot
    plt.savefig(f'${params.outdir}/genes_by_class/genes_grouped_by_class_{genome_id}.png')
    #plt.show()

final_df = pd.concat(class_counts_list, ignore_index=True)
total_counts = final_df.groupby('Class')['Num_Genes'].sum().reset_index()
total_counts = total_counts.sort_values(by='Num_Genes',ascending=False)
total_counts = total_counts[total_counts['Class'] != 'under_development']
total_counts['Class'] = total_counts['Class'].str.capitalize()
total_counts.to_csv('${params.outdir}/total_genes_grouped_by_class.csv', index=False)
plt.figure(figsize=(10, 6))
# Normalize the number of genes for color mapping
norm = plt.Normalize(total_counts['Num_Genes'].min(), total_counts['Num_Genes'].max())
colors = plt.cm.Blues(norm(total_counts['Num_Genes']))
bars = plt.bar(total_counts['Class'], total_counts['Num_Genes'], color=colors)
plt.xlabel('Class')
plt.ylabel('Number of Genes')
plt.title('Total Number of Genes Grouped by Class')
plt.xticks(rotation=90)
plt.tight_layout()
legend_labels = [f'{cls}: {count}' for cls, count in zip(total_counts['Class'], total_counts['Num_Genes'])]
plt.legend(bars, legend_labels, title='ARGs by Class', loc='best',ncol=2)
plt.savefig('${params.outdir}/total_genes_grouped_by_class.png')
#plt.show()

# Create a DataFrame to store the number of genes for each class across all genomes
boxplot_data = pd.DataFrame()

# Iterate over each DataFrame in class_counts_list
for class_counts in class_counts_list:
    for index, row in class_counts.iterrows():
        class_name = row['Class']
        num_genes = row['Num_Genes']
        if class_name not in boxplot_data:
            boxplot_data[class_name] = []
        boxplot_data[class_name] = boxplot_data[class_name] + [num_genes]

# Convert the dictionary to a DataFrame
boxplot_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in boxplot_data.items()]))

# Create a box plot for each class
plt.figure(figsize=(12, 8))
sns.boxplot(data=boxplot_df, orient='h')
plt.xlabel('Number of Genes')
plt.ylabel('Class')
plt.title('Distribution of Number of Genes per Class')
plt.tight_layout()
plt.savefig('${params.outdir}/boxplot_genes_per_class.png')
    "
    """
}

workflow arg_mge_association {
    println "Workflow mode: To identify ARGs and MGEs associated with genes and ARGs"
    if (params.read1 != "" && params.read2 != "") {
            def paired_reads_channel = tuple(params.read1, params.read2)
        	trimming(paired_reads_channel)
    }
    else {
        paired_reads_channel = Channel.fromFilePairs("${params.reads}")
        paired_reads_channel = paired_reads_channel.map { pair -> [pair[1][0], pair[1][1]] }
        trimming(paired_reads_channel)
    }
    assemble(trimming.out)
    mgefind(assemble.out)
    argFind(assemble.out)
    argAssociation(argFind.out, mgefind.out)
    csvCombine(mgefind.out.collect(),argAssociation.out.collect())
}

workflow arg_mge_association_single_end {
    println "Workflow mode: To identify ARGs and associated MGEs (Single End reads)"
    input_channel = Channel.from(params.read1.split(/\s+/))
    trimming_single_end(input_channel)
    assemble_single_end(trimming_single_end.out)
    mgefind(assemble_single_end.out)
    argFind(assemble_single_end.out)
    argAssociation(argFind.out, mgefind.out)
    csvCombine(mgefind.out.collect(),argAssociation.out.collect())
}

workflow arg {
    println "Workflow mode: To identify ARGs"
    if (params.read1 != "" && params.read2 != "") {
            def paired_reads_channel = tuple(params.read1, params.read2)
        	trimming(paired_reads_channel)
    }
    else {
        paired_reads_channel = Channel.fromFilePairs("${params.reads}")
        paired_reads_channel = paired_reads_channel.map { pair -> [pair[1][0], pair[1][1]] }
        trimming(paired_reads_channel)
    }
    assemble(trimming.out)
    argFind(assemble.out)
}

workflow mge{
    println "Workflow mode: To generate MGEs"
    if (params.read1 != "" && params.read2 != "") {
            def paired_reads_channel = tuple(params.read1, params.read2)
        	trimming(paired_reads_channel)
    }
    else {
        paired_reads_channel = Channel.fromFilePairs("${params.reads}")
        paired_reads_channel = paired_reads_channel.map { pair -> [pair[1][0], pair[1][1]] }
        trimming(paired_reads_channel)
    }
    assemble(trimming.out)
    mgefind(assemble.out)
    genomeAnnotation(assemble.out)
    mge_association(mgefind.out, genomeAnnotation.out)
}

workflow mge_single_end{
    println "Workflow mode: To generate MGEs (Single End reads)"
    input_channel = Channel.from(params.read1.split(/\s+/))
    trimming_single_end(input_channel)
    assemble_single_end(trimming_single_end.out)
    genomeAnnotation(assemble_single_end.out)
    mgefind(assemble_single_end.out)
    mge_association(mgefind.out, genomeAnnotation.out)
}

workflow genome_annotation_single_end {
    println "Workflow mode: To generate genome annotation (Single End reads)"
	input_channel = Channel.from(params.read1.split(/\s+/))
	trimming_single_end(input_channel)
	assemble_single_end(trimming_single_end.out)
	genomeAnnotation(assemble_single_end.out)
}

workflow genome_annotation_paired_end {
    println "Workflow mode: To generate genome annotation (Paired end reads)"
    if (params.read1 != "" && params.read2 != "") {
            def paired_reads_channel = tuple(params.read1, params.read2)
        	trimming(paired_reads_channel)
    }
    else {
        paired_reads_channel = Channel.fromFilePairs("${params.reads}")
        paired_reads_channel = paired_reads_channel.map { pair -> [pair[1][0], pair[1][1]] }
        trimming(paired_reads_channel)
    }
    assemble(trimming.out)
    genomeAnnotation(assemble.out)
}

workflow kmer_single_end {
    println "Workflow mode: To generate k-mer counts (Single End reads)"
    input_channel = Channel.from(params.read1.split(/\s+/))
    trimming_single_end(input_channel)
    kmer(trimming_single_end.out)
}

workflow kmer_paired_end {
    println "Workflow mode: To generate k-mer counts (Paired end reads)"
    if (params.read1 != "" && params.read2 != "") {
            def paired_reads_channel = tuple(params.read1, params.read2)
        	trimming(paired_reads_channel)
    }
    else {
        paired_reads_channel = Channel.fromFilePairs("${params.reads}")
        paired_reads_channel = paired_reads_channel.map { pair -> [pair[1][0], pair[1][1]] }
        trimming(paired_reads_channel)
    }
    assemble(trimming.out)
    kmer_fasta(assemble.out)
}

workflow assembly_single_end {
    println "Workflow mode: To assemble single end reads"
    input_channel = Channel.from(params.read1.split(/\s+/))
    trimming_single_end(input_channel)
    assemble(trimming_single_end.out)
}

workflow assembly {
    println "Workflow mode: To assemble paired end reads"
    if (params.read1 != "" && params.read2 != "") {
            def paired_reads_channel = tuple(params.read1, params.read2)
        	trimming(paired_reads_channel)
    }
    else {
        paired_reads_channel = Channel.fromFilePairs("${params.reads}")
        paired_reads_channel = paired_reads_channel.map { pair -> [pair[1][0], pair[1][1]] }
        trimming(paired_reads_channel)
    }
    //paired_reads_channel = paired_reads_channel.map { pair ->
    //(params.read1 != "" && params.read2 != "") ? [params.read1, params.read2] : [pair[1][0], pair[1][1]]}
    assemble(trimming.out)
}

workflow annotation_singleend{
	println "Workflow mode: To annotate the variants"
	println "Warning! Only one reads file was included. It will be assumed to be a single end file"
    indexReference(params.ref)
    input_channel = Channel.from(params.read1.split(/\s+/)).view()
    trimming_single_end(input_channel)
    alignReads_single_end(indexReference.out, trimming_single_end.out, params.ref)
    samToBam(alignReads_single_end.out)
    sorting(samToBam.out)
    variantCalling(sorting.out, params.ref)
    mutationAnnotation(variantCalling.out, params.bed)
    SNP(mutationAnnotation.out)
    SNPFiltering(SNP.out)
    SNPOutput(SNPFiltering.out)
    posRefAlt(SNPOutput.out)
    encode(posRefAlt.out, params.pca)
    predict(encode.out, params.model, params.labels)
}
	
workflow no_annotation_singleend{
	println "Workflow mode: To not annotate the variants (Default)"
	println "Warning! Only one reads file was included. It will be assumed to be a single end file"
    indexReference(params.ref)
    input_channel = Channel.from(params.read1.split(/\s+/)).view()
    trimming_single_end(input_channel)
    alignReads_single_end(indexReference.out, trimming_single_end.out,params.ref)
    samToBam(alignReads_single_end.out)
    sorting(samToBam.out)
    variantCalling(sorting.out,params.ref)
    SNP(variantCalling.out)
    SNPFiltering(SNP.out)
    SNPOutput(SNPFiltering.out)
    posRefAlt(SNPOutput.out)
    encode(posRefAlt.out, params.pca)
	predict(encode.out, params.model, params.labels)
	}

workflow no_annotation{
	println "Workflow mode: To not annotate the variants and bulk processing"
	println "Using wildcard rule: ${params.reads}"
    indexReference(params.ref)
    paired_reads_channel = Channel.fromFilePairs("${params.reads}")
    //paired_reads_channel = paired_reads_channel.map { pair -> [pair[1][0], pair[1][1]] } //def reads_tuple = tuple(params.read1, params.read2)
    paired_reads_channel = paired_reads_channel.map { pair ->
    (params.read1 != "" && params.read2 != "") ? [params.read1, params.read2] : [pair[1][0], pair[1][1]]}
    trimming(paired_reads_channel)
    alignReads(indexReference.out, trimming.out,params.ref)
    samToBam(alignReads.out)
    sorting(samToBam.out)
    variantCalling(sorting.out,params.ref)
    SNP(variantCalling.out)
    SNPFiltering(SNP.out)
    SNPOutput(SNPFiltering.out)
    posRefAlt(SNPOutput.out)
    encode(posRefAlt.out, params.pca)
	predict(encode.out, params.model, params.labels)
	}

workflow annotation{
	println "Workflow mode: To annotate the variants and bulk processing"
    indexReference(params.ref)
    if (params.read1 != "" && params.read2 != "") {
            def paired_reads_channel = tuple(params.read1, params.read2)
        	trimming(paired_reads_channel)
    }
    else {
        paired_reads_channel = Channel.fromFilePairs("${params.reads}")
        paired_reads_channel = paired_reads_channel.map { pair -> [pair[1][0], pair[1][1]] }
        trimming(paired_reads_channel)
    }
    alignReads(indexReference.out, trimming.out,params.ref)
    samToBam(alignReads.out)
    sorting(samToBam.out)
    variantCalling(sorting.out,params.ref)
    mutationAnnotation(variantCalling.out, params.bed)
    SNP(mutationAnnotation.out)
    SNPFiltering(SNP.out)
    SNPOutput(SNPFiltering.out)
    posRefAlt(SNPOutput.out)
    encode(posRefAlt.out, params.pca)
	predict(encode.out, params.model, params.labels)
	}

workflow process_all {
    println "Workflow mode: To process all"
    if (params.read1 != "" && params.read2 != "") {
            def paired_reads_channel = tuple(params.read1, params.read2)
        	trimming(paired_reads_channel)
    }
    else {
        paired_reads_channel = Channel.fromFilePairs("${params.reads}")
        paired_reads_channel = paired_reads_channel.map { pair -> [pair[1][0], pair[1][1]] }
        trimming(paired_reads_channel)
    }
    alignReads(indexReference.out, trimming.out,params.ref)
    samToBam(alignReads.out)
    sorting(samToBam.out)
    variantCalling(sorting.out,params.ref)
    mutationAnnotation(variantCalling.out, params.bed)
    SNP(mutationAnnotation.out)
    SNPFiltering(SNP.out)
    SNPOutput(SNPFiltering.out)
    posRefAlt(SNPOutput.out)
    encode(posRefAlt.out, params.pca)
    predict(encode.out, params.model, params.labels)
    assemble(trimming.out)
    kmer(assemble.out)
    genomeAnnotation(assemble.out)
    mgefind(assemble.out)
    mge_association(mgefind.out, genomeAnnotation.out)
    argFind(assemble.out)
    argAssociation(argFind.out, mgefind.out)
}

workflow process_all_single_end {
    println "Workflow mode: To process all (Single End reads)"
    input_channel = Channel.from(params.read1.split(/\s+/))
    trimming_single_end(input_channel)
    alignReads_single_end(indexReference.out, trimming_single_end.out,params.ref)
    samToBam(alignReads_single_end.out)
    sorting(samToBam.out)
    variantCalling(sorting.out,params.ref)
    mutationAnnotation(variantCalling.out, params.bed)
    SNP(mutationAnnotation.out)
    SNPFiltering(SNP.out)
    SNPOutput(SNPFiltering.out)
    posRefAlt(SNPOutput.out)
    encode(posRefAlt.out, params.pca)
    predict(encode.out, params.model, params.labels)
    assemble_single_end(trimming_single_end.out)
    kmer(trimming_single_end.out)
    genomeAnnotation(assemble_single_end.out)
    mgefind(assemble_single_end.out)
    argFind(assemble_single_end.out)
    argAssociation(argFind.out, mgefind.out)
    mge_association(mgefind.out, genomeAnnotation.out)
}

workflow arg_mge_association_assembled {
    println "Workflow mode: To identify ARGs and MGEs associated with genes and ARGs (Assembled reads)"
    def isTxtInput = params.txt_in == "yes"
    if (isTxtInput){
        println"Reading input paths from the given text file. Make sure that you have given full path of the files"
        input_channel = Channel
	    .fromPath(params.read1)
    	.flatMap { file ->
    	    file.text.readLines() // Read each line from the file
    	}
    } else {
    println "Reading input paths direcly from ls command"
    input_channel = Channel.from(params.read1.split(/\s+/))
    }
    mgefind(input_channel)
    argFind(input_channel)
    argAssociation(argFind.out.collect(), mgefind.out)
    csvCombine(mgefind.out.collect(),argAssociation.out.collect())
}

workflow arg_assembled {
    println "Workflow mode: To identify ARGs (Assembled reads)"
    input_channel = Channel.from(params.read1.split(/\s+/))
    argFind(input_channel)
}

workflow mge_assembled {
    println "Workflow mode: To generate MGEs (Assembled reads)"
    input_channel = Channel.from(params.read1.split(/\s+/))
    mgefind(input_channel)
    genomeAnnotation(input_channel)
    mge_association(mgefind.out, genomeAnnotation.out)
}

workflow genome_annotation_assembled {
    println "Workflow mode: To generate genome annotation (Assembled reads)"
    input_channel = Channel.from(params.read1.split(/\s+/))
    genomeAnnotation(input_channel)
}

workflow kmer_assembled {
    println "Workflow mode: To generate k-mer counts (Assembled reads)"
    input_channel = Channel.from(params.read1.split(/\s+/))
    kmer_fasta(input_channel)
}

workflow annotation_assembled {
    println "Workflow mode: To annotate the variants (Assembled reads) and predict phenotype"
    input_channel = Channel.from(params.read1.split(/\s+/))
    indexReference(params.ref)
    alignReads(indexReference.out, input_channel,params.ref)
    samToBam(alignReads.out)
    sorting(samToBam.out)
    variantCalling(sorting.out,params.ref)
    mutationAnnotation(variantCalling.out, params.bed)
    SNP(mutationAnnotation.out)
    SNPFiltering(SNP.out)
    SNPOutput(SNPFiltering.out)
    posRefAlt(SNPOutput.out)
    encode(posRefAlt.out, params.pca)
    predict(encode.out, params.model, params.labels)
}

workflow no_annotation_assembled {
    println "Workflow mode: To not annotate the variants (Assembled reads)"
    input_channel = Channel.from(params.read1.split(/\s+/))
    indexReference(params.ref)
    alignReads(indexReference.out, input_channel,params.ref)
    samToBam(alignReads.out)
    sorting(samToBam.out)
    variantCalling(sorting.out,params.ref)
    SNP(variantCalling.out)
    SNPFiltering(SNP.out)
    SNPOutput(SNPFiltering.out)
    posRefAlt(SNPOutput.out)
    encode(posRefAlt.out, params.pca)
    predict(encode.out, params.model, params.labels)
}

workflow snp_no_predict {
    println "Workflow mode: To identify SNPs and not predict phenotype"
    indexReference(params.ref)
    if (params.read1 != "" && params.read2 != "") {
            def paired_reads_channel = tuple(params.read1, params.read2)
        	trimming(paired_reads_channel)
    }
    else {
        paired_reads_channel = Channel.fromFilePairs("${params.reads}")
        paired_reads_channel = paired_reads_channel.map { pair -> [pair[1][0], pair[1][1]] }
        trimming(paired_reads_channel)
    }
    alignReads(indexReference.out, trimming.out,params.ref)
    samToBam(alignReads.out)
    sorting(samToBam.out)
    variantCalling(sorting.out,params.ref)
    mutationAnnotation(variantCalling.out, params.bed)
    SNP(mutationAnnotation.out)
    SNPFiltering(SNP.out)
    SNPOutput(SNPFiltering.out)
}

workflow snp_no_predict_single_end {
    println "Workflow mode: To identify SNPs and not predict phenotype (Single End reads)"
    indexReference(params.ref)
    input_channel = Channel.from(params.read1.split(/\s+/))
    trimming_single_end(input_channel)
    alignReads_single_end(indexReference.out, trimming_single_end.out,params.ref)
    samToBam(alignReads_single_end.out)
    sorting(samToBam.out)
    variantCalling(sorting.out,params.ref)
    mutationAnnotation(variantCalling.out, params.bed)
    SNP(mutationAnnotation.out)
    SNPFiltering(SNP.out)
    SNPOutput(SNPFiltering.out)
}

workflow snp_no_predict_assembled {
    println "Workflow mode: To identify SNPs and not predict phenotype (Assembled reads)"
    println "Warning: The SNPs will have a dept of only one"
    indexReference(params.ref)
    input_channel = Channel.from(params.read1.split(/\s+/))
    mutationAnnotation(input_channel, params.bed)
    SNP(mutationAnnotation.out)
    SNPFiltering(SNP.out)
    SNPOutput(SNPFiltering.out)
}

workflow process_all_assembled {
    println "Workflow mode: To process all (Assembled reads)"
    input_channel = Channel.from(params.read1.split(/\s+/))
    alignReads(indexReference.out, input_channel,params.ref)
    samToBam(alignReads.out)
    sorting(samToBam.out)
    variantCalling(sorting.out,params.ref)
    mutationAnnotation(variantCalling.out, params.bed)
    SNP(mutationAnnotation.out)
    SNPFiltering(SNP.out)
    SNPOutput(SNPFiltering.out)
    posRefAlt(SNPOutput.out)
    encode(posRefAlt.out, params.pca)
    predict(encode.out, params.model, params.labels)
    assemble(input_channel)
    kmer(input_channel)
    genomeAnnotation(input_channel)
    mgefind(input_channel)
    argFind(input_channel)
    argAssociation(argFind.out, mgefind.out)
    mge_association(mgefind.out, genomeAnnotation.out)
}

workflow {
    def isPairedEnd = params.read2 != "" || params.reads != "" && params.assembled != "yes"
    def isSingleEnd = params.read2 == "" && params.reads == "" && params.assembled != "yes"
    def isAnnotation = params.annotation == "yes"
    def isSnpMode = params.mode == "snp"
    def isKmerMode = params.mode == "kmer"
    def isAssembleMode = params.mode == "assemble"
    def isGenomeAnnotationMode = params.mode == "genome_annotation"
    def isMgeMode = params.mode == "mge"
    def isSupermode = params.mode == "super"
    def isArgMgeMode = params.mode == "arg_mge" || params.mode == "mge_arg"
    def isArgMode = params.mode == "arg"
    def isAssembled = params.assembled == "yes"
    def isNoPredict = params.mode == "snp-no-predict"
    def isTxtInput = params.txt_in == "yes"

    if (isSnpMode) {
        if (isSingleEnd) {
            if (isAnnotation) {
                annotation_singleend()
            } else {
                no_annotation_singleend()
            }
        } else {
            if (isAnnotation) {
                if (params.reads != "") {
                    annotation()
                }
            } else {
                if (params.reads != "") {
                    no_annotation()
                }
            }
        }
    } else if (isKmerMode) {
        if (isPairedEnd) {
            kmer_paired_end()
        } else {
            kmer_single_end()
        }
    } else if (isAssembleMode) {
        if (isSingleEnd) {
            assembly_single_end()
        } else {
            assembly()
        }
    } else if (isGenomeAnnotationMode) {
    	if (isAssembled) {
        genome_annotation_assembled()
        }
        else if (isSingleEnd) {
            genome_annotation_single_end()
        } else if (isPairedEnd) {
            genome_annotation_paired_end()
        } 
    } else if (isMgeMode) {
        if (isSingleEnd) {
            mge_single_end()
        } else {
            mge()
        }
    } else if (isSupermode) {
        if (isSingleEnd) {
            process_all_single_end()
        } else {
            process_all()
        }
    } else if (isArgMgeMode) {
        if (isAssembled) {
            arg_mge_association_assembled()
        } else if (isSingleEnd) {
            arg_mge_association_single_end()
        } else {
            arg_mge_association()
        }
    } else if (isArgMode) {
        if (isAssembled) {
        arg_assembled()
        } else {
        arg()
        }
    } else if (isAssembled) {
        if (isSnpMode) {
            if (isAnnotation) {
                annotation_assembled()
            } else {
                no_annotation_assembled()
            }
        } else if (isKmerMode) {
            kmer_assembled()
        } else if (isAssembleMode) {
            println "Assemble mode is not applicable for assembled reads"
        } else if (isMgeMode) {
            mge_assembled()
        } else if (isArgMgeMode) {
            arg_mge_association_assembled()
        } else if (isArgMode) {
            arg_assembled()
        } else if (isNoPredict) {
            snp_no_predict_assembled()
        } else if (isSupermode) {
            process_all_assembled()
        }
        else {
            println "Invalid mode for assembled reads"
        }
    } else if (isNoPredict) {
        if (isSingleEnd) {
            snp_no_predict_single_end()
        } else {
            snp_no_predict()
        }
    } else {
        println "Invalid mode"
    }
}
