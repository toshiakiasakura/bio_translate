import os 
import numpy as np
import pandas as pd 
import argparse
import bio_trans 


def create_default_directory():
    for d in ["xml", "json", "markdown", "html"]:
        if not os.path.exists(d):
            os.mkdir(d)

def check_env_variable():
    email= os.getenv("BIO_EMAIL")
    if email == None:
        raise Exception(""" Environmental variable "BIO_EMAIL is not specified. """)

def trans_csv(args):
    df_doi = pd.read_csv(args.csv)
    cond = df_doi["status"].isnull()
    df_run = df_doi.loc[cond]
    for i, r in df_run.iterrows():
        doi = r["doi"]
        name = r["name"]
        print(f"- doi  : {doi}\n- name : {name}")
        try:
            if args.json2html:
                const_md = bio_trans.MDConstructor(name)
                const_md.convert_json2html(args.open_html)
            else:
                print(name)
                bio_trans.doi2html(doi, name, args.open_html)
        except Exception as e:
            print(e)
            if args.access_doi:
                bio_trans.access2doi(doi)
        print("-"*50 + "\n")

def one_trans(args):
    if args.name == None:
        raise Exception(""" "--name/-n" should also be selected.""")
    if args.json2html:
        const_md = bio_trans.MDConstructor(args.name)
        const_md.convert_json2md(args.open_html)
    else:
        bio_trans.doi2html(args.doi, args.name, args.open_html)

def parse_args():
    parser = argparse.ArgumentParser(description="Tranlate article from doi into markdown/html.")
    parser.add_argument("--doi", "-d",
        help="""One doi to markdon. Should be set "--name/-n" option if this is selected.""") 
    parser.add_argument("--name","-n", 
        help="""Name for each file. Name should be without extension. Should be set with "--name/-n" option.""")
    parser.add_argument("--csv",
        help="""Process doi from csv. csv should include "doi" and "name" header.""")
    parser.add_argument("--open-html","-o",  action="store_true",
        help="""Open tab when html are produced.""")
    parser.add_argument("--access-doi","-a",  action="store_true",
        help="""Open target doi.org page when doi is not correctly converted.""")
    parser.add_argument("--json2html","-j",  action="store_true",
        help="""Convert process starts from json file. If html file does not exist, raise error.""")
    parser.add_argument("--check","-c", 
        help="""List up journals used in pubmed central database. Journal data are in jlist.csv.""") 

    args = parser.parse_args()
    return(args)

def main():
    check_env_variable()
    create_default_directory()
    args = parse_args()
    if args.check:
        bio_trans.check_journal_in_PMC(args.check)
    elif args.csv: 
        trans_csv(args)
    elif args.doi:
        one_trans(args)
    else:
        print("""Error : Use"--doi/-d" and "--name/-n", or "--csv" option.\n""")

    print("Program finished. ")

if __name__ == "__main__":
    main()

