# https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/
# example file: https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/pubmed22n1003.xml.gz
import os
import requests
from selectolax.parser import HTMLParser

base_url = 'https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/'
base_dir = '/media/yjzhang/easystore-5gb-1/research_big_data/pubmed_abstracts/'


updates_url = 'https://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/'

def download_all_files(base_url, base_dir):
    r = requests.get(base_url)
    tree = HTMLParser(r.content)
    for node in tree.css('a'):
        if '.xml.gz' in node.text():
            url = base_url + node.attributes['href']
            print(url)
            if not os.path.exists(base_dir + node.text()):
                req = requests.get(url)
                with open(os.path.join(base_dir, node.text()), 'wb') as f:
                    f.write(req.content)
            else:
                print('file exists')


download_all_files(base_url, base_dir)
download_all_files(updates_url, base_dir)
