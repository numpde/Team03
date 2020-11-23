#    HK, 2020-11-12
import requests
from bs4 import BeautifulSoup
import warnings

CLINVAR_SELECTORS = {
    'clin_significance': '#kis-variation > div.shadowbox > div > div > div > div:nth-child(1) > ul > li:nth-child(1)',
    'ref': '#kis-variation > div.shadowbox > div > div > div > div:nth-child(1) > ul > li:nth-child(2)',
    'var': '#kis-variation > div.shadowbox > div > div > div > div:nth-child(1) > ul > li:nth-child(3)',
    'var_type': '#kis-variation > div.shadowbox > div > div > div > div:nth-child(1) > ul > li:nth-child(4)',
    'organism': '#kis-variation > div.shadowbox > div > div > div > div:nth-child(2) > ul > li:nth-child(1)',
    'gene': '#kis-variation > div.shadowbox > div > div > div > div:nth-child(2) > ul > li:nth-child(2)',
    'GRCh38_p12': '#kis-variation > div.shadowbox > div > div > div > div:nth-child(2) > ul > li:nth-child(3)',
    'GRCh37_p13': '#kis-variation > div.shadowbox > div > div > div > div:nth-child(2) > ul > li:nth-child(4)'
}

SNP_SELECTORS = {
    'functional_consequence': '#maincontent > div > div:nth-child(5) > div > div.rslt > div.supp > dl > dd:nth-child(12)'
}


def get_content_from_soup(content):
    if len(content) == 1:
        for elem in content:
            kv = elem.text.replace('\n', '').split(':')
        return kv[0], ':'.join(kv[1:])
    else:
        warnings.warn(f'soup has not length 1: {content}')
        return None, None


def get_ncbi_info(rs_id: str) -> dict:
    """
    Sends a request to the ncbi website and returns the information of a variance in a dictionary.
    """
    page = requests.get(f'https://www.ncbi.nlm.nih.gov/clinvar/?term={rs_id}')

    soup = BeautifulSoup(page.content)
    ncbi_info = {}
    for selector_k, selector_v in CLINVAR_SELECTORS.items():
        k, v = get_content_from_soup(soup.select(selector_v))
        ncbi_info[k] = v
    return ncbi_info


def get_functional_consequence_from_SNP(rs_id: str) -> str:
    page = requests.get(f'https://www.ncbi.nlm.nih.gov/snp/?term={rs_id}')
    soup = BeautifulSoup(page.content)
    selection = soup.select(SNP_SELECTORS['functional_consequence'])
    functional_consequence = None
    for elem in selection:
        functional_consequence = elem.text

    return functional_consequence


if __name__ == '__main__':
    rs_id = 'rs369122589'
    get_functional_consequence_from_SNP(rs_id)
