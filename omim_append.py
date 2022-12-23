import pandas as pd
import optparse

# helper functions
def import_omim_db(table):
    print(f"---- Reading Database: {table}\n")
    cols = ['MIM Number', 'Gene Symbols', 'Approved Gene Symbol', 'Entrez Gene ID', 'Phenotypes']
    omim_type = {'MIM Number':int, 'Gene Symbols':str, 'Approved Gene Symbol':str,'Phenotypes':str, 'Entrez Gene ID':int}
    DB = pd.read_excel(table, usecols = cols, converters = omim_type)
    return(DB)

def import_qci_table(table):
    qci_table_type = {'OMIM ID':int, 'OMIM Gene':str, 'Gene Symbol':str,'OMIM Phenotype':str}
    if table.endswith(".gz"):
        print(f"---- Reading Gzipped Table: {table}\n")
        df = pd.read_table(table, sep = "\t", converters = qci_table_type, compression = 'gzip')
        
    elif table.endswith(".tsv"):
        print(f"---- Reading Table: {table}\n")
        df = pd.read_table(table, sep = "\t", converters = qci_table_type)
        
    df_cols = list(df.columns)
    # insert specific columns for easy interpretation
    df_cols.insert(50, "HPO Link")
    df_cols.insert(50, "OMIM Link")
    df_cols.insert(50, "OMIM ID")
    df_cols.insert(9, "Gene ID")
    return(df.reindex(columns = df_cols))

def merger(qci_table, omim_db):
    cols = ['MIM Number', 'Gene Symbols', 'Approved Gene Symbol', 'Entrez Gene ID', 'Phenotypes']
    merged_df = pd.merge(qci_table, omim_db, left_on = "Gene Symbol", right_on = "Approved Gene Symbol")
    merged_df["OMIM ID"] = merged_df["MIM Number"] 
    merged_df["OMIM Gene"] = merged_df["Approved Gene Symbol"]
    merged_df["OMIM Phenotype"] = merged_df["Phenotypes"]
    merged_df["Gene ID"] = merged_df["Entrez Gene ID"]
    merged_df["OMIM Link"] = merged_df['OMIM ID'].apply(lambda x: f'https://omim.org/entry/{x}')
    merged_df["HPO Link"] = merged_df['Gene ID'].apply(lambda x: f'https://hpo.jax.org/app/browse/gene/{x}')
    return(merged_df.loc[:, ~merged_df.columns.isin(cols)])

def make_hyperlink(value, base):
        # add links to omim and hpo
        # '<a href="https://omim.org/entry/{}">https://omim.org/entry/{}</a>'
        # '=HYPERLINK("https://omim.org/entry/{}";"{}")'
        # '<a href="https://hpo.jax.org/app/browse/gene/{}">https://hpo.jax.org/app/browse/gene/{}</a>'
        # '=HYPERLINK("https://hpo.jax.org/app/browse/gene/{}";"{}")'
    url = base + "{}"
    return '=HYPERLINK("%s", "%s")' % (url.format(value), value)

def exelize(merged_df):
    merged_df["OMIM Link"] = merged_df["OMIM Link"].apply(make_hyperlink, base = "https://omim.org/entry/")
    merged_df["HPO Link"] = merged_df["HPO Link"].apply(make_hyperlink, base = "https://hpo.jax.org/app/browse/gene/")
    return(merged_df)

def Main():
    parser = optparse.OptionParser()

    # add options
    # TODO: add 'action = "append", ' to --file option
    parser.add_option('-f', '--file', dest = 'file', type = 'str', default = None, action = "append",
        help = 'TSV table exported from QCI')
    
    parser.add_option('-d', '--database', dest = 'db', type = 'str', default = None, 
        help = 'OMIM database')
    
    parser.add_option('-e', '--export', dest = 'export', type = 'str', default = "excel", 
        help = 'Output file format (Optional). Default is excel.')
    
    parser.add_option('-o', '--output', dest = 'out', type = 'str', default = None, 
        help = 'Output file name (Optional). Do not use with export flag.')

    parser.add_option('-o', '--output', dest = 'out', type = 'str', default = None, 
        help = 'Output file name (Optional). Do not use with export flag.')
    
    
    (options, args) = parser.parse_args()

    if options.db and options.file is None:
        raise Exception(f"---- Filename flag (Example: '-f sample.tsv.gz') should be filled.\n",
                        f"---- Database flag (Example: '-d OMIM.xlsx') should be filled.")

    # read database
    # col_types_omim = {'Chromosome':str, 'Genomic Position Start':int, 'Genomic Position End':int, 'Cyto Location':str,'Computed Cyto Location':str,'MIM Number':int, 'Gene Symbols':str, 'Gene Name':str,'Approved Gene Symbol':str, 'Entrez Gene ID':int, 'Ensembl Gene ID':str, 'Comments':str,'Phenotypes':str, 'Mouse Gene Symbol':str,'Mouse Gene ID':str}
    if options.db is not None:
        try:
            OMIM_DB = import_omim_db(options.db)
        
        except IOError as err:
            print(f"---- Error: File type of {options.db} does not supported.")
            print(err)
            exit(1)

    # read table
    if options.file is not None:
        try:
            df_list = [import_qci_table(file) for file in options.file]
        
        except IOError as err:
            print(f"---- Error: File type of {options.file} does not supported.")
            print(err)
            exit(1)

    # left join two tables
    try:
        final_df_list = [merger(df, OMIM_DB) for df in df_list]
    
    except IOError as err:
        print(err)
        exit(1)

    # default export behaviour
    default = options.file.split(".")[0]
    if options.export in ["excel", ".xlsx", "xlsx", ".xls", "xls"]:
        print(f"---- Exporting {default}.xlsx\n")
        final_excel = exelize(final_df)
        final_excel.to_excel(f"{default}.xlsx", index = False)

    elif options.export in ["table", ".tsv", "tab"]:
        print(f"---- Exporting {default}.tsv\n", f"---- Exporting {default}.tsv.gz\n")
        final_df.to_csv(f"{default}.tsv", index = False)
        final_df.to_csv(f"{default}.tsv.gz", index = False, compression="gzip")
        
    elif options.out is not None and options.out.endswith(".xlsx"):
        print(f"---- Exporting {options.out}\n")
        final_excel.to_excel(f"{options.out}", index = False)
    
    elif options.out is not None and options.out.endswith(".tsv"):
        print(f"---- Exporting {options.out}\n")
        final_df.to_csv(f"{options.out}.tsv", index = False)
    
    elif options.out is not None and options.out.endswith(".tsv.gz"):
        print(f"---- Exporting {options.out}\n")
        final_df.to_csv(f"{default}.tsv.gz", index = False, compression="gzip")

    elif options.export in ["excel", ".xlsx", "xlsx", ".xls", "xls"] and options.out is not None and options.out.endswith(".xlsx"):
        print(f"---- Exporting {options.out}\n")
        final_excel = exelize(final_df)
        final_excel.to_excel(f"{options.out}.xlsx", index = False)
    else:
        raise Exception("---- Something wrong with either filename or file extension.\n---- Please check help by typing:\n---- python omim-append.py --help.")
    
if __name__ == '__main__':
    Main()