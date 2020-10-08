import os 
import numpy as np
import pandas as pd 
import argparse
import bio_trans 


def create_default_directory():
    for d in ["xml", "json", "markdown"]:
        if not os.path.exists(d):
            os.mkdir(d)

def trans_csv(args):
    df_doi = pd.read_csv(args.csv)
    cond = df_doi["status"].isnull()
    df_run = df_doi.loc[cond]
    for i, r in df_run.iterrows():
        doi = r["doi"]
        name = r["name"]
        print(f"- doi  : {doi}\n- name : {name}")
        try:
            if args.json2markdown:
                const_md = bio_trans.MDConstructor(name)
                const_md.convert_json2md(args.open)
            else:
                bio_trans.doi2markdown(doi, name, args.open)
        except Exception as e:
            print(e)
            if args.access_doi:
                bio_trans.access2doi(doi)
        print("-"*50 + "\n")

def one_trans(args):
    if args.name == None:
        raise Exception(""" "--name/-n" should also be selected.""")
    if args.json2markdown:
        const_md = bio_trans.MDConstructor(args.name)
        const_md.convert_json2md(args.open)
    else:
        bio_trans.doi2markdown(args.doi, args.name, args.open)

def parse_args():
    parser = argparse.ArgumentParser(description="Tranlate article from doi into markdown.")

    parser.add_argument("--doi", "-d",
        help="""One doi to markdon. Should be set "--name/-n" option if this is selected.""") 
    parser.add_argument("--name","-n", 
        help="""Name for each file. Name should be without extension. Should be set with "--name/-n" option.""")
    parser.add_argument("--csv",
        help="""Process doi from csv. csv should include "doi" and "name" header.""")
    parser.add_argument("--open-markdown","-o",  action="store_true",
        help="""Open tab when markdonws are produced.""")
    parser.add_argument("--access-doi","-a",  action="store_false",
        help="""Open target doi.org page when doi is not correctly converted.""")
    parser.add_argument("--json2markdown","-j",  action="store_false",
        help="""Convert process starts from json file. If json file does not exist, raise error.""")
    parser.add_argument("--check","-c", 
        help="""List up journals used in pubmed central database. Journal data are in jlist.csv.""") 

    args = parser.parse_args()
    return(args)

def main():
    create_default_directory()
    args = parse_args()
    if args.check:
        bio_trans.check_journal_in_PMC(args.check)
    elif args.csv: 
        trans_csv(args)
    elif args.doi:
        one_trans(args)

if __name__ == "__main__":
    main()

