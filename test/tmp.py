from bioutils.io.gtf import GtfRecord, parse_gtf_attrs

if __name__ == "__main__":
    print(parse_gtf_attrs(b'gene_id "STRG.23"; transcript_id "STRG.23.1"; reference_id "XM_017000355.2"; ref_gene_id "MIB2"; ref_gene_name "MIB2"; cov "1.073021"; FPKM "196.106232"; TPM "586.075928";'))
