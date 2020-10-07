import numpy as np
import pandas as pd 
import argparse

import bio_trans 


def main():
    args = parse_args()
    if args.check:
        bio_trans.check_journal_in_PMC(args.check)
        return 


    df_doi = pd.read_csv("./doi_management.csv")
    cond = df_doi["status"].isnull()
    df_run = df_doi.loc[cond]
    for i, r in df_run.iterrows():
        doi = r["doi"]
        name = r["name"]
        print(f"- doi  : {doi}\n- name : {name}")
        try:
            bio_trans.doi2markdown(doi, name, args.open)
        except Exception as e:
            print(e)
            bio_trans.access2doi(doi)
        print("-"*50 + "\n")

def parse_args():
    parser = argparse.ArgumentParser(description="Tranlate article from doi into markdown.")
    parser.add_argument("--check","-c", help="""If this option is selected, list up journals 
used in pubmed central""", default = "")
    parser.add_argument("--open", "-o", help="""Open tab when markdonws are produced.""",
        action="store_true")

    args = parser.parse_args()
    return(args)

if __name__ == "__main__":
    main()

