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

import tag_wrapper 

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
    path = f"./xml/{name}.xml"
    return(path)

def get_markdown_path(name):
    """Convert name into path to markdown data.
    """
    path = f"./markdown/{name}.md"
    return(path)

def get_json_path(name):
    """Convert name into path to json data.
    """
    path = f"./json/{name}.json" 
    return(path)

def save_xml_from_doi(doi, name):
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

def read_json(name):
    """Read json file produced from pubmed full text xml. 

    Args:   
        name (str) : name of target file without extensions.
    """
    path = get_json_path(name)
    with open(path, "r") as f: 
        js_ = json.load(f) 
    return(js_)

def parse_xml_into_json(path, trans=True):
    """Parse pubmed Central full text xml data and 
    convert it into json format.
    
    Args:
        path (str)   : path to xml data. 
        trans (bool) : apply translation if True. 
    """
    meta  = pp.parse_pubmed_xml(path)
    ref   = pp.parse_pubmed_references(path)
    paras = pp.parse_pubmed_paragraph(path, all_paragraph=True)
    captions = pp.parse_pubmed_caption(path)
    tables = pp.parse_pubmed_table(path)

    if trans: 
        meta["abstract_jp"] = en2jp(meta["abstract"])
        for p in tqdm(paras, desc="google trans"):
            if p["section"] == "":
                continue
            p["text"] = p["text"].rstrip().rstrip("\n")
            p["text_jp"] = en2jp(p["text"])
    js_ = {"meta" : meta, 
            "ref"  : ref, 
            "paragraphs": paras,
            "captions" : captions,
            "tables"   : tables
    }

    return(js_)

def parse_and_save_xml_into_json(name) :
    """Parse xml and save json.

    Args:   
        name (str) : name of xml. Do not need extensiono of the file.
    """
    path2xml = get_xml_path(name)
    path2json = get_json_path(name) 

    js_= parse_xml_into_json(path2xml)
    with open(path2json, "w") as f:
        json.dump(js_, f)


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
<br> <br>
    <div style="display: inline-block;">
        {text2}
    </div>
</div>
<hr>
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

def en2jp(text):
    """Translate text from English to Japanese.
    Using google translation API.

    Args: 
        text (str) : text to be translated. Should be English.

    Return:
        str : translated japanese. 
    """
    translator = Translator()
    trans_jp = translator.translate(text, src="en", dest="ja")
    return(trans_jp.text)

class MDConstructor():

    def __init__(self,name, keywords=None, model_keywords=None):
        """Parse xml data using pubmed_parser (OSS). 

        Args:
            name (str) : name of target file without extensions.
        """
        self.path2json = get_json_path(name)
        self.path2md  = get_markdown_path(name)

        with open(self.path2json, "r") as f:
            js_ = json.load(f)

        self.meta  = js_["meta"]
        self.ref   = js_["ref"]
        self.paras = js_["paragraphs"]
        self.captions = js_["captions"]
        self.tables = js_["tables"]

        self.keywords = keywords if keywords != None else ["O3", "ozone", "O 3"] 
        self.model_keywords = model_keywords if model_keywords != None \
            else ["model", "regression", "モデル", "回帰"]

    def meta_info_md(self):
        """Construct meta information for markdown. 
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
        return(text)

    def abstract_md(self):
        abst_en = self.keywords_handler(self.meta["abstract"])
        abst_jp = self.keywords_handler(self.meta["abstract_jp"], sp="。")
        abst_both = wrap_two_columns(abst_en, abst_jp)
        text = f"""
## Abstract 
{abst_both}
"""
        return(text)

    def keywords_handler(self, text, sp=". ", skip=False): 
        """Add bold tag "**" for markdown format. 

        Args:
            text (str) : text. 
            sp (str)   : separation by this. 
            skip (str) : contain sentense which do not include specific keywords.
        """
        t_lis= text.split(sp)
        t_lis_tag = []
        for t in t_lis:
            t_tag = t 
            flag = False
            for k in self.keywords:
                if k in t: 
                    t_tag = tag_wrapper.bold(t_tag)
                    t_tag = tag_wrapper.purple(t_tag)
                    flag = True
                    break 
            for k in self.model_keywords:
                if k in t:
                    k_tag = k
                    k_tag = tag_wrapper.italic(k_tag)
                    k_tag = tag_wrapper.red(k_tag)
                    t_tag = t_tag.replace(k, k_tag)

            if not skip:
                t_lis_tag.append(t_tag)
            elif flag:
                t_lis_tag.append(t_tag)

        if len(t_lis_tag) == 0:
            return("")

        text_tagged = sp.join(t_lis_tag)  
        #if sp == "。": 
        #    text_tagged += sp
        text_tagged = tag_wrapper.fontsize(text_tagged)
        
        return(text_tagged)

    def paragraph_constrct(self, skip):
        """Paragraph constrction for "convert2md" 
        """
        text = ""
        section = ""

        for i,p in enumerate(self.paras):
            if p["section"] == "":
                continue
            if p["section"] != section:
                section = p["section"]
                text += f"\n## {section}\n"

            text_en = self.keywords_handler(p["text"], skip=skip) 
            if text_en == "":
                continue
            text_jp = self.keywords_handler(p["text_jp"], sp="。", skip=skip)

            text += wrap_two_columns(text_en, text_jp)

        return(text)

    def convert2md(self, open_=True):
        """Convert parsed dictionary into markdown format. 

        Args:
            d (dict) : generated from "parse_xml". 
            path (str)  : target path for markdown file.
            open_ (bool): open markdown in webbrowser or not.

        """
        text = ""
        text += self.meta_info_md()
        text += self.abstract_md()

        #text += "\n\n# Thesis Body Part\n" 
        text += self.paragraph_constrct(skip=False)

        #text += "\n\n# Summarization Part\n" 
        #text += self.paragraph_constrct(skip=True)

        with open(self.path2md, "w") as f:
            f.write(text)
        if open_:
            webbrowser.open(self.path2md, new=2)

def doi2markdown(doi,name, open_=False):
    """Given doi, convert it into pmc id, fetch full text of the artcile, 
    parse the full text, translate english into japanse, 
    sum up into markdown.

    Args:
        doi (str) : doi.
        name (str) : name of target file without extensions.
    """
    save_xml_from_doi(doi, name)
    parse_and_save_xml_into_json(name)
    const_md = MDConstructor(name)
    const_md.convert2md(open_)
    

class XmlParser():
    """This class should be depricated.
    """
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






