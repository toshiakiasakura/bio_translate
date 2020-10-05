import os 
import requests
import re
import json
from tqdm import tqdm
import xml.etree.ElementTree as ET

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


def convert_ID(id_, get_id_tp="pmcid"):
    """Convert one id into another one.

    Args:
        id_ (str) : this id is converted from another one. 
            Input can be "doi", "pmid" or "pmcid".
        get_id_tp (str) : type for extraction id. 

    Return:
        str : id of "get_id_tp".
    """
    tool = "service-root"
    email = os.getenv("BIO_EMAIL")
    url = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool={}&amp;email={}&amp;email&amp;ids={}&amp;format=json'.format(tool, email, ','.join([id_]))

    r = requests.get(url)
    j = json.loads(r.text)

    j_rec= j.get("records")
    if j_rec == None:
        raise Exception("""Request does not succeed.
            json data is {j}""")
    
    get_id = j_rec[0].get(get_id_tp) 
    if get_id == None:
        errmsg = j_rec[0].get("errmsg", "errormsg does not exists")
        raise Exception(f"""{get_id_tp} can not be extracted.
            Error message in response is {errmsg}
            json data is {j}""")
    return(get_id)

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

def remove_specific_tags(text, tag_names):
    """Remove specific tags from raw text
    remaining texts inside tag. 

    Args: 
        text (str) : xml formtted text. 
        tag_names (list) : list of tag names. 

    Return:
        str : tag removed text.
    """
    remove_list = [f"<{tag}.*?>"  for tag in tag_names] + [f"</{tag}>" for tag in tag_names]
    for k in remove_list:
        regex = re.compile(k)
        text = re.sub(regex,"",text)
    return(text)

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






