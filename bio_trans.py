import os 
import requests
import re
import json
from tqdm import tqdm
import xml.etree.ElementTree as ET
import urllib.parse
import pandas as pd
import numpy as np 
import webbrowser

import pubmed_parser as pp
from Bio import Entrez
from googletrans import Translator 

Entrez.email = os.getenv("BIO_EMAIL")


def fetch_pmids(term, db):
    """Fetch pmids given the word. 

    Args: 
        term (str) : search word. 
        db (str) : database. 

    Returns:
        list : pm ids.

    References: 
        - Biopythonを使ってPMCから論文取得
        https://roy29fuku.com/natural-language-processing/paper-analysis/retrieve-pmc-full-text-with-biopython/
        - EFetch, official.  
        https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
    """
    pmids = []
    retmax = 10000 
    handle = Entrez.esearch(db=db,term=term)
    record = Entrez.read(handle)
    count = int(record['Count'])
    for retstart in tqdm(range(0, count, retmax)):
        handle = Entrez.esearch(db=db, term=term, retmax=retmax, retstart=retstart)
        record = Entrez.read(handle)
        pmids.extend(record['IdList'])
    return(pmids)

def fetch_one_abstract(pmid, db="pubmed", retmode="xml"):
    """Fetch one abstract data. 

    Args: 
        pmid (str) : pubmed id. 
        db (str) : database name.
        retmode (str) : return mode. 

    Args: 
        str : abstract data formatted in xml text format. 
    """
    handle = Entrez.efetch(db=db, id=pmid, retmode=retmode)
    text = handle.read().decode()
    return(text)

def fetch_one_full_text(pmcid, db="pmc", rettype="full", retmode="xml"):
    """Fetch one full text data. Almost similar to "fetch_one_abstract".

    Args: 
        pmcid (str) : pmcid. 
        db (str) : database name.
        rettype (str) : return type
        retmode (str) : return mode. 

    Args: 
        str : full text data formatted in xml text format. 
    """
    handle = Entrez.efetch(db=db, id=pmcid, rettype=rettype, retmode=retmode)
    text = handle.read().decode()
    return(text)


def convert_ID(id_,idtype="doi", get_id_tp="pmcid"):
    """Convert one id into another one.

    Args:
        id_ (str) : this id is converted from another one. 
            Input can be "doi", "pmid" or "pmcid".
        idtype (str) : type of input id.
        get_id_tp (str) : type for extraction id. 

    Return:
        str : id of "get_id_tp".

    References:
        - ID Converter API
        https://www.ncbi.nlm.nih.gov/pmc/pmctopmid/
    """
    #id_ = urllib.parse.quote(id_)
    tool = "service-root"
    email = os.getenv("BIO_EMAIL")
    #url = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool={}&amp;email={}&amp;email&amp;ids={}&amp;format=json'.format(tool, email, ','.join([id_]))
    url = f'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={id_}&idtype={idtype}&tool={tool}&email={email}&format=json'

    r = requests.get(url)
    j = json.loads(r.text)

    j_rec= j.get("records")
    if j_rec == None:
        raise Exception(f"""Request does not succeed.
            json data is {j}""")
    
    get_id = j_rec[0].get(get_id_tp) 
    if get_id == None:
        errmsg = j_rec[0].get("errmsg", "errormsg does not exists")
        raise Exception(f"""{get_id_tp} can not be extracted.
            Error message in response is {errmsg}
            json data is {j}""")
    return(get_id)


def print_side_by_side(a, b, size=30, space=4):
    """Print two long texts side by side. 
    
    Args: 
        a (str) : text
        b (str) : text 
        size (int) : size for print for each line.
        space (int) : space between two texts. 
    """
    while a or b:
        print(a[:size].ljust(size) + " " * space + b[:size])
        a = a[size:]
        b = b[size:] 

def remove_specific_tags(text, tag_names=None):
    """Remove specific tags from raw text
    remaining texts inside tag. 

    Args: 
        text (str) : xml formtted text. 
        tag_names (list) : list of tag names. 

    Return:
        str : tag removed text.
    """
    tag_names = tag_names if not tag_names == None else ["xref","ext-link","italic", "sub"]

    remove_list = [f"<{tag}.*?>"  for tag in tag_names] + [f"</{tag}>" for tag in tag_names]
    for k in remove_list:
        regex = re.compile(k)
        text = re.sub(regex,"",text)
    return(text)

def check_journal_in_PMC(keyword, full=True):
    """Read journal list of PMC and check the journal is in PMC Database. 
    Print out journal which is in PMC filtered by "keyword".

    Args:
        keyword (str) : used for filter.
    """
    df = pd.read_csv("jlist.csv")
    j_title = "Journal title"
    cond = df[j_title].str.contains(keyword)
    dfM = df.loc[cond]

    for i in dfM[j_title].unique():
        print(i)
    if full:
        display(dfM)

def get_xml_path(name):
    """Convert name into path to xml data.
    """
    path = f"./data/{name}.xml"
    return(path)

def get_markdown_path(name):
    """Convert name into path to markdown data.
    """
    path = f"./markdown/{name}.md"
    return(path)

def save_xml(doi, name):
    """Save xml data from doi to xml formatted texts. 

    Args: 
        doi (str) : doi 
        path (str) : path to save. 
    """
    pmc = convert_ID(doi)
    text = fetch_one_full_text(pmc)
    path = get_xml_path(name)
    with open(path, "w") as f:
        f.write(text)


class WrapText():
    fontsize = "14px"
    
    @classmethod
    def wrap_fontsize(cls,text):
        text = f"""<span style="font-size:{cls.fontsize}">{text}</span>"""
        return(text)

def wrap_two_columns(text1,text2):
    """Wrap two texts as two columns.

    Args:
        text1, text2 (str) : texts. 

    References:
        - Have two columns in Markdown
        https://stackoverflow.com/questions/30514408/have-two-columns-in-markdown
    """
    text =f"""
<div style="-webkit-column-count: 2; -moz-column-count: 2; column-count: 2; -webkit-column-rule: 1px dotted #e0e0e0; -moz-column-rule: 1px dotted #e0e0e0; column-rule: 1px dotted #e0e0e0;">
    <div style="display: inline-block;">
        {text1}
    </div>
- - - - 
    <div style="display: inline-block;">
        {text2}
    </div>
</div>
"""
    return(text)


def remove_specific_element(lis_, word): 
    """Remove specific element of list. 

    Args:
        lis_ (list) : target list of remove. 
        word (str)  : regular expression form. 
    """
    regex = re.compile(word) 
    lis_ = [ v for v in lis_ if not regex.match(v)]
    return(lis_)

class MDConstructor():
    counter = 0

    def __init__(self,name):
        """Parse xml data using pubmed_parser (OSS). 

        Args:
            path (str) : path to xml data.
        """
        self.path2xml = get_xml_path(name)
        self.path2md  = get_markdown_path(name)

        path = self.path2xml
        self.meta  = pp.parse_pubmed_xml(path)
        self.ref   = pp.parse_pubmed_references(path)
        self.paras = pp.parse_pubmed_paragraph(path, all_paragraph=True)
        self.captions = pp.parse_pubmed_caption(path)
        self.tables = pp.parse_pubmed_table(path)


    def en2jp(self,text):
        """Translate text from English to Japanese.
        Using google translation API.

        Args: 
            text (str) : text to be translated. Should be English.

        Return:
            str : translated japanese. 
        """
        self.counter += 1 
        print("google translator counter: ", self.counter)
        translator = Translator()
        trans_jp = translator.translate(text, src="en", dest="ja")
        return(trans_jp.text)

    def list_en_and_jp(self,text_en):
        """list English text and Japanese text in markdown format.
        Also fontsize is changed.

        Args:
            en_text (str) : English text. 
        """
        text_jp = self.en2jp(text_en)
        text_en = WrapText.wrap_fontsize(text_en)
        text_jp = WrapText.wrap_fontsize(text_jp)
        both_text = wrap_two_columns(text_en, text_jp)
        return(both_text)

    def paragraph_constrct(self):
        """Paragraph constrction for "convert2md" 
        """
        text = ""
        section = ""
        for p in self.paras:
            if p["section"] == "":
                continue
            if p["section"] != section:
                section = p["section"]
                text += f"\n### {section}\n"
            text += self.list_en_and_jp(p["text"])
        return(text)


    def convert2md(self, open_=True):
        """Convert parsed dictionary into markdown format. 

        Args:
            d (dict) : generated from "parse_xml". 
            path (str)  : target path for markdown file.
            open_ (bool): open markdown in webbrowser or not.

        """
        doi = self.meta["doi"]
        authors = [ " ".join(remove_specific_element(author, "Aff\d+")) for author in self.meta["author_list"] ] 
        authors_str = ", ".join(authors)

        
        text = f"""# {self.meta["full_title"]} \n
- doi : [{doi}](https://doi.org/{doi}) 
- pmc : {self.meta["pmc"]}
- Authors : {authors_str}
- Publication_date: {self.meta["publication_date"]}
"""
        abst_both = self.list_en_and_jp(self.meta["abstract"])
        text += f"""
## Abstract 
{abst_both}
"""
        text += "\nThesis Body Part\n" 
        text += self.paragraph_constrct()

        with open(self.path2md, "w") as f:
            f.write(text)
        if open_:
            webbrowser.open(self.path2md, new=2)

class XmlParser():
    def __init__(self,text):
        """Set root from xml formatted text.

        Args: 
            text (str) : xml formatted text.

        Explanation:
            self.root : xml.etree.elementtree.element object.
        """
        self.root = ET.fromstring(text)
    
    def print_all_child(self,obj=None, intend=0):
        """Print all the tags from root with modified indents. 

        Args: 
            obj (xml.etree.elementtree.element) : xml element object.
            intend (int) : number of intends, increases as recusively used.
        """
        obj = obj if obj != None else self.root

        for child in obj:
            print(" "*4*intend + child.tag)
            self.print_all_child(child, intend+1)

    def xml_print(self, iter_text, obj=None):
        """Print tab, attrib, text at once from self.root. 

        Args: 
            iter_text (str) : tag name
        """
        obj = obj if obj else self.root
        for neighbor in obj.iter(iter_text):
            print("#####", neighbor.tag, neighbor.attrib)
            print(neighbor.text)
            print()

    def xml_child_print(self,  iter_text, obj=None,):
        """Print out xml children objects.

        Args: 
            iter_text (str) : tag name
            obj (xml.etree.elementtree.element) : xml element object.
        """
        obj = obj if obj != None else self.root
        for neighbor in obj.iter(iter_text):
            for child in neighbor:
                self.xml_print(child.tag, child)
    
    def print_all_under_obj(self, obj=None):
        """Print all the texts from objects. 

        Args: 
            obj (xml.etree.elementtree.element) : xml element object.
        """
        obj = obj if obj != None else self.root
        print(obj.text)
        children = obj.getchildren()
        if not len(children):
            return()
        for child in children:
            self.print_all_under_obj(child)
            print(child.tail)
            
    def constract_list(self,obj=None, list_=None):
        obj = obj if obj != None else self.root
        list_ = list_ if list_ != None else []

        list_.append(obj.text)
        children = obj.getchildren()
        for child in children:
            self.constract_list(child, list_)
            list_.append(child.tail)
        return(list_)






