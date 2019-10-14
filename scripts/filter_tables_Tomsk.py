from lib.blister import *
import numpy as np

def splitDataFrameList(df,target_column,separator):
    ''' df = dataframe to split,
    target_column = the column containing the values to split
    separator = the symbol used to perform the split
    returns: a dataframe with each entry for the target column separated, with each element moved into a new row. 
    The values in the other columns are duplicated across the newly divided rows.
    '''
    def splitListToRows(row,row_accumulator,target_column,separator):
        split_row = row[target_column].split(separator)
        for s in split_row:
            new_row = row.to_dict()
            new_row[target_column] = s
            row_accumulator.append(new_row)
    new_rows = []
    df.apply(splitListToRows,axis=1,args = (new_rows,target_column,separator))
    new_df = pd.DataFrame(new_rows)
    return new_df

genemap = pd.read_csv("/dev/datasets/FairWind/_db/genemap2.txt", sep='\t', low_memory=False)
genemap.drop(columns=['Chromosome',  'Genomic Position Start',  'Genomic Position End', 'Cyto Location', 'Computed Cyto Location', 'Mim Number', 'Mouse Gene Symbol/ID', 'Entrez Gene ID', 'Approved Symbol',  'Ensembl Gene ID'], inplace=True)
genemap = genemap[genemap['Gene Symbols'] == genemap['Gene Symbols']]

genemap = splitDataFrameList(genemap, 'Gene Symbols', ', ')
genemap.rename(columns={'Gene Symbols' : 'SYMBOL'}, inplace=True)

def the_thread(block, output_dir):
	
	columns_list = ['AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF']
	
	order = ['HGMD', 'CHROM', 'POS', 'REF', 'ALT', 'GT', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene Name', 'Phenotypes', 'Comments', 'Gene', 'BIOTYPE', 'TOTAL_AF', 'EXON', 'INTRON', 'SIFT', 'PolyPhen','Feature_type', 'Feature', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'DISTANCE', 'STRAND', 'FLAGS', 'SYMBOL_SOURCE', 'HGNC_ID', 'CANONICAL', 'TSL', 'APPRIS', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'REFSEQ_MATCH', 'SOURCE', 'GIVEN_REF', 'USED_REF', 'BAM_EDIT', 'DOMAINS', 'HGVS_OFFSET', 'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF', 'CLIN_SIG', 'SOMATIC', 'PHENO', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'PHENOTYPES', 'Condel', 'MPC', 'AA', 'BLOSUM62', 'miRNA', '1000Gp3_AC', '1000Gp3_AF', '1000Gp3_AFR_AC', '1000Gp3_AFR_AF', '1000Gp3_AMR_AC', '1000Gp3_AMR_AF', '1000Gp3_EAS_AC', '1000Gp3_EAS_AF', '1000Gp3_EUR_AC', '1000Gp3_EUR_AF', '1000Gp3_SAS_AC', '1000Gp3_SAS_AF', 'ALSPAC_AC', 'ALSPAC_AF', 'AltaiNeandertal', 'Ancestral_allele', 'CADD_phred', 'CADD_raw', 'CADD_raw_rankscore', 'DANN_rankscore', 'DANN_score', 'Denisova', 'ESP6500_AA_AC', 'ESP6500_AA_AF', 'ESP6500_EA_AC', 'ESP6500_EA_AF', 'Eigen-PC-phred', 'Eigen-PC-raw', 'Eigen-PC-raw_rankscore', 'Eigen-phred', 'Eigen-raw', 'Eigen_coding_or_noncoding', 'Ensembl_geneid', 'Ensembl_proteinid', 'Ensembl_transcriptid', 'ExAC_AC', 'ExAC_AF', 'ExAC_AFR_AC', 'ExAC_AFR_AF', 'ExAC_AMR_AC', 'ExAC_AMR_AF', 'ExAC_Adj_AC', 'ExAC_Adj_AF', 'ExAC_EAS_AC', 'ExAC_EAS_AF', 'ExAC_FIN_AC', 'ExAC_FIN_AF', 'ExAC_NFE_AC', 'ExAC_NFE_AF', 'ExAC_SAS_AC', 'ExAC_SAS_AF', 'ExAC_nonTCGA_AC', 'ExAC_nonTCGA_AF', 'ExAC_nonTCGA_AFR_AC', 'ExAC_nonTCGA_AFR_AF', 'ExAC_nonTCGA_AMR_AC', 'ExAC_nonTCGA_AMR_AF', 'ExAC_nonTCGA_Adj_AC', 'ExAC_nonTCGA_Adj_AF', 'ExAC_nonTCGA_EAS_AC', 'ExAC_nonTCGA_EAS_AF', 'ExAC_nonTCGA_FIN_AC', 'ExAC_nonTCGA_FIN_AF', 'ExAC_nonTCGA_NFE_AC', 'ExAC_nonTCGA_NFE_AF', 'ExAC_nonTCGA_SAS_AC', 'ExAC_nonTCGA_SAS_AF', 'ExAC_nonpsych_AC', 'ExAC_nonpsych_AF', 'ExAC_nonpsych_AFR_AC', 'ExAC_nonpsych_AFR_AF', 'ExAC_nonpsych_AMR_AC', 'ExAC_nonpsych_AMR_AF', 'ExAC_nonpsych_Adj_AC', 'ExAC_nonpsych_Adj_AF', 'ExAC_nonpsych_EAS_AC', 'ExAC_nonpsych_EAS_AF', 'ExAC_nonpsych_FIN_AC', 'ExAC_nonpsych_FIN_AF', 'ExAC_nonpsych_NFE_AC', 'ExAC_nonpsych_NFE_AF', 'ExAC_nonpsych_SAS_AC', 'ExAC_nonpsych_SAS_AF', 'FATHMM_converted_rankscore', 'FATHMM_pred', 'FATHMM_score', 'GERP++_NR', 'GERP++_RS', 'GERP++_RS_rankscore', 'GM12878_confidence_value', 'GM12878_fitCons_score', 'GM12878_fitCons_score_rankscore', 'GTEx_V6p_gene', 'GTEx_V6p_tissue', 'GenoCanyon_score', 'GenoCanyon_score_rankscore', 'H1-hESC_confidence_value', 'H1-hESC_fitCons_score', 'H1-hESC_fitCons_score_rankscore', 'HUVEC_confidence_value', 'HUVEC_fitCons_score', 'HUVEC_fitCons_score_rankscore', 'Interpro_domain', 'LRT_Omega', 'LRT_converted_rankscore', 'LRT_pred', 'LRT_score', 'M-CAP_pred', 'M-CAP_rankscore', 'M-CAP_score', 'MetaLR_pred', 'MetaLR_rankscore', 'MetaLR_score', 'MetaSVM_pred', 'MetaSVM_rankscore', 'MetaSVM_score', 'MutPred_AAchange', 'MutPred_Top5features', 'MutPred_protID', 'MutPred_rankscore', 'MutPred_score', 'MutationAssessor_UniprotID', 'MutationAssessor_pred', 'MutationAssessor_score', 'MutationAssessor_score_rankscore', 'MutationAssessor_variant', 'MutationTaster_AAE', 'MutationTaster_converted_rankscore', 'MutationTaster_model', 'MutationTaster_pred', 'MutationTaster_score', 'PROVEAN_converted_rankscore', 'PROVEAN_pred', 'PROVEAN_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HDIV_rankscore', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_pred', 'Polyphen2_HVAR_rankscore', 'Polyphen2_HVAR_score', 'REVEL_rankscore', 'REVEL_score', 'Reliability_index', 'SIFT_converted_rankscore', 'SIFT_pred', 'SIFT_score', 'SiPhy_29way_logOdds', 'SiPhy_29way_logOdds_rankscore', 'SiPhy_29way_pi', 'TWINSUK_AC', 'TWINSUK_AF', 'Transcript_id_VEST3', 'Transcript_var_VEST3', 'Uniprot_aapos_Polyphen2', 'Uniprot_acc_Polyphen2', 'Uniprot_id_Polyphen2', 'VEST3_rankscore', 'VEST3_score', 'aapos', 'clinvar_clnsig', 'clinvar_golden_stars', 'clinvar_rs', 'clinvar_trait', 'codon_degeneracy', 'fathmm-MKL_coding_group', 'fathmm-MKL_coding_pred', 'fathmm-MKL_coding_rankscore', 'fathmm-MKL_coding_score', 'gnomAD_exomes_AC', 'gnomAD_exomes_AF', 'gnomAD_exomes_AFR_AC', 'gnomAD_exomes_AFR_AF', 'gnomAD_exomes_AFR_AN', 'gnomAD_exomes_AMR_AC', 'gnomAD_exomes_AMR_AF', 'gnomAD_exomes_AMR_AN', 'gnomAD_exomes_AN', 'gnomAD_exomes_ASJ_AC', 'gnomAD_exomes_ASJ_AF', 'gnomAD_exomes_ASJ_AN', 'gnomAD_exomes_EAS_AC', 'gnomAD_exomes_EAS_AF', 'gnomAD_exomes_EAS_AN', 'gnomAD_exomes_FIN_AC', 'gnomAD_exomes_FIN_AF', 'gnomAD_exomes_FIN_AN', 'gnomAD_exomes_NFE_AC', 'gnomAD_exomes_NFE_AF', 'gnomAD_exomes_NFE_AN', 'gnomAD_exomes_OTH_AC', 'gnomAD_exomes_OTH_AF', 'gnomAD_exomes_OTH_AN', 'gnomAD_exomes_SAS_AC', 'gnomAD_exomes_SAS_AF', 'gnomAD_exomes_SAS_AN', 'gnomAD_genomes_AC', 'gnomAD_genomes_AF', 'gnomAD_genomes_AFR_AC', 'gnomAD_genomes_AFR_AF', 'gnomAD_genomes_AFR_AN', 'gnomAD_genomes_AMR_AC', 'gnomAD_genomes_AMR_AF', 'gnomAD_genomes_AMR_AN', 'gnomAD_genomes_AN', 'gnomAD_genomes_ASJ_AC', 'gnomAD_genomes_ASJ_AF', 'gnomAD_genomes_ASJ_AN', 'gnomAD_genomes_EAS_AC', 'gnomAD_genomes_EAS_AF', 'gnomAD_genomes_EAS_AN', 'gnomAD_genomes_FIN_AC', 'gnomAD_genomes_FIN_AF', 'gnomAD_genomes_FIN_AN', 'gnomAD_genomes_NFE_AC', 'gnomAD_genomes_NFE_AF', 'gnomAD_genomes_NFE_AN', 'gnomAD_genomes_OTH_AC', 'gnomAD_genomes_OTH_AF', 'gnomAD_genomes_OTH_AN', 'integrated_confidence_value', 'integrated_fitCons_score', 'integrated_fitCons_score_rankscore', 'phastCons100way_vertebrate', 'phastCons100way_vertebrate_rankscore', 'phastCons20way_mammalian', 'phastCons20way_mammalian_rankscore', 'phyloP100way_vertebrate', 'phyloP100way_vertebrate_rankscore', 'phyloP20way_mammalian', 'phyloP20way_mammalian_rankscore', 'LoFtool', 'CADD_PHRED', 'CADD_RAW', 'MaxEntScan_alt', 'MaxEntScan_diff', 'MaxEntScan_ref', 'ada_score', 'rf_score', 'QUAL', 'FILTER', 'PL', 'Allele']
	
	index, input_filename = block
	output_filename = Blister.Output(input_filename, output_dir, "squeezed_named", "csv", rewrite=True, index=index)
	if not output_filename: return

	Blister.Sleep(max_level=70.0, interval=1, index=index)

	with Blister.Timestamp("Processing file", input_filename, output_filename, index=index) as start_time:
		table = pd.read_csv(input_filename, sep='\t', low_memory=False)
		for col in columns_list:
			table[col] = pd.to_numeric(table[col], errors='coerce', downcast='float')
		table['TOTAL_AF'] = table[columns_list].apply(np.nanmax, axis=1)
		table = table[((table['IMPACT'] == 'HIGH') | (table['IMPACT'] == 'MODERATE')) & ((table['TOTAL_AF'] < 0.01) | table['TOTAL_AF'].isna())]
		
		table = pd.merge(genemap, table, how='right', on=['SYMBOL'])
		
		table = table[order]
	
		table['classifier'] = table['CHROM'].str.cat(table['POS'].apply(str), sep=':')
		table = table.groupby(table['classifier'])
	
		new_table = pd.DataFrame()
		for col in order:
			new_table[col] = table[col].apply(set).apply(list).apply(lambda x: [t for t in x if str(t) != 'nan']).apply(lambda x: (x[0] if ((len(x) == 1) and (type(x) == type(list()))) else (float('nan') if ((not x) and (type(x) == type(list()))) else ', '.join([str(t) for t in x]))))
		del table
		new_table.sort_values(by=['CHROM', 'POS'], axis=0, ascending=True, inplace=True)
		new_table.to_csv(output_filename, sep='\t', index=False)
		del new_table

results = Blister.EachFile("Tomsk Filter Tables", ["/dev/datasets/FairWind/_results/bowtie/dupless_tables/tables/*.csv"], "/dev/datasets/FairWind/_results/bowtie/dupless_tables/squeezed_named/")(the_thread)()
