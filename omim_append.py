import pandas as pd
from halo import Halo
import argparse, warnings

warnings.filterwarnings('ignore', message="^Columns.*")

# helper functions
def import_omim_db(table):
    with Halo(text=f'Reading Database: {table}', spinner='dots'):
        # col_types_omim = {'Chromosome':str, 'Genomic Position Start':int, 'Genomic Position End':int, 'Cyto Location':str,'Computed Cyto Location':str,'MIM Number':int, 'Gene Symbols':str, 'Gene Name':str,'Approved Gene Symbol':str, 'Entrez Gene ID':int, 'Ensembl Gene ID':str, 'Comments':str,'Phenotypes':str, 'Mouse Gene Symbol':str,'Mouse Gene ID':str}
        cols = ['MIM Number', 'Gene Symbols', 'Approved Gene Symbol', 'Entrez Gene ID', 'Phenotypes']
        omim_type = {'MIM Number':int, 'Gene Symbols':str, 'Approved Gene Symbol':str,'Phenotypes':str, 'Entrez Gene ID':int}
        DB = pd.read_excel(table, usecols = cols, converters = omim_type)
    return(DB)

def import_qci_table(table):
    qci_table_type = {'OMIM ID':int, 'OMIM Gene':str, 'Gene Symbol':str,'OMIM Phenotype':str}
    with Halo(text=f'Reading Table: {table}', spinner='dots'):
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

def make_hyperlink(url):
        # '<a href="https://omim.org/entry/{}">https://omim.org/entry/{}</a>'
        # '<a href="https://hpo.jax.org/app/browse/gene/{}">https://hpo.jax.org/app/browse/gene/{}</a>'
    return '=HYPERLINK("%s", "%s")' % (url, url)

def exelize(merged_df):
    merged_df["OMIM Link"] = merged_df["OMIM Link"].apply(make_hyperlink)
    merged_df["HPO Link"] = merged_df["HPO Link"].apply(make_hyperlink)
    return(merged_df)

# currently only excel or tab seperated value
def export_table(x, output, extension):
    if extension in ["excel", ".xlsx", "xlsx", ".xls", "xls"]:
        with Halo(text=f'Exporting {output}.xlsx', spinner='dots'):
            x.to_excel(f"{output}.xlsx", index = False)    

    elif extension in ["table", ".tsv", "tsv", "tab"]:
        with Halo(text=f'Exporting {output}.tsv.gz', spinner='dots'):
            x.to_csv(f"{output}.tsv.gz", index = False, compression = "gzip")
    
    else:
        raise Exception("Extension not supported.")

# MAIN FUNCTION (Delete this comment later.) -----------------------------------------

def Main():
    parser = argparse.ArgumentParser()

    # add options
    parser.add_argument('-f', '--file', dest = 'file', default = None, nargs = '+', required=True,
        help = 'TSV table exported from QCI')
    
    parser.add_argument('-d', '--database', dest = 'db', default = None, required=True,
        help = 'OMIM database')
    
    parser.add_argument('-e', '--export', dest = 'export', default = "xlsx",
        help = 'Output file format (Optional). Default is excel.')
    
    parser.add_argument('-o', '--output', dest = 'out', default = None, nargs = '+',
        help = 'Output file name (Optional). Do not use with export flag.')
    
    parser.add_argument('--duo', dest = 'duo', default = False, action = "store_true",
        help = 'DUO')
    parser.add_argument('--trio', dest = 'trio', default = False, action = "store_true",
        help = 'TRIO')
    parser.add_argument('--quadro', dest = 'quadro', default = False, action = "store_true",
        help = 'QUADRO')
    args = parser.parse_args()

    try: # read database
        OMIM_DB = import_omim_db(args.db)
    except IOError as err:
        print(f"---- Error: File type of {args.db} does not supported.")
        print(err)
        exit(1)

    try: # read table
        df_list = [import_qci_table(file) for file in args.file]
    except IOError as err:
        print(f"---- Error: File type of {args.file} does not supported.")
        print(err)
        exit(1)

    try: # left join two tables
        final_df_list = [merger(df, OMIM_DB) for df in df_list]
    except IOError as err:
        print(err)
        exit(1)
    
    try: # turn into excel file
        final_excel_list = [exelize(final_df) for final_df in final_df_list]
    except IOError as err:
        print(err)
        exit(1)

    if args.out is not None and len(args.out) == len(args.file):
        for index, output in enumerate(args.out):
            export_table(final_df_list[index], output.split(".")[0], output.split(".")[1])

    elif args.export:
        for index, excel in enumerate(final_excel_list):
            export_table(excel, args.file[index].split(".")[0], args.export)

    else:
        raise Exception("---- Something wrong with either filename, file extension or provided options.\n",
                        "---- Please check help by typing:\n",
                        "---- python omim-append.py --help.")
    
if __name__ == '__main__':
    Main()