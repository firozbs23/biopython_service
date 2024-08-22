import re
from time import sleep
import utility.utils as utils
import pycountry

from Bio import Entrez
from Bio import Medline

MAX_RETRY = 3
SLEEP_TIME = 3


def fetch_pubmed_record(pubmed_id):
    Entrez.email = 'firoz@brainstation-23.com'  # Set your email here
    handle = Entrez.efetch(db='pubmed', id=pubmed_id, rettype='medline', retmode='text')
    record = handle.read()
    handle.close()
    return record


def extract_mesh_terms(record_text):
    mesh_terms = []
    for line in record_text.split('\n'):
        if line.startswith('MH  - '):
            mesh_term = line.replace('MH  - ', '')
            mesh_terms.append(mesh_term)
    return mesh_terms


def get_mesh_terms(pubmed_id):
    record_text = fetch_pubmed_record(pubmed_id)
    mesh_terms = extract_mesh_terms(record_text)
    return mesh_terms


def fetch_pubmed_title(pmid):
    Entrez.email = "firoz@brainstation-23.com"
    retry = 0
    while retry < MAX_RETRY:
        try:
            """
                Fetches the title of a PubMed article given its PubMed ID.

                Args:
                pubmed_id (str): The PubMed ID (PMID) of the article.

                Returns:
                str: The title of the article.
                """

            # Fetch the record using Entrez
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
            record = Medline.read(handle)
            handle.close()

            data = record.get("TI", "")
            return data

        except Exception as e:
            print(f'Exception occurred while extracting pubmed title. Message: {e}')
            retry += 1
            sleep(SLEEP_TIME)

    return ''


def fetch_pubmed_abstract(pmid):
    Entrez.email = "firoz@brainstation-23.com"

    retry = 0
    while retry < MAX_RETRY:
        try:
            """
                Fetches the abstract of a PubMed article given its PubMed ID.

                Args:
                pubmed_id (str): The PubMed ID (PMID) of the article.

                Returns:
                str: The text of the abstract.
                """

            # Provide your email address to NCBI

            # Fetch the abstract using Entrez
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
            record = Medline.read(handle)
            handle.close()

            data = record.get("AB", "")
            return data

        except Exception as e:
            print(f'Exception occurred while extracting pubmed abstract. Message: {e}')
            retry += 1
            sleep(SLEEP_TIME)

    return ''


def fetch_pubmed_pmcid(pmid):
    Entrez.email = "firoz@brainstation-23.com"

    retry = 0
    while retry < MAX_RETRY:
        try:
            """
                Fetches the PMCID of a PubMed article given its PubMed ID.

                Args:
                pubmed_id (str): The PubMed ID (PMID) of the article.

                Returns:
                str: The PMCID of the article if available, else an empty string.
                """

            # Provide your email address to NCBI

            # Fetch the record using Entrez
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
            record = Medline.read(handle)
            handle.close()

            data = record.get('PMC', '')
            return data

        except Exception as e:
            print(f'Exception occurred while collecting PMCID. Mesage: {e}')
            retry += 1
            sleep(SLEEP_TIME)

    return ''


def _fetch_pubmed_doi(pmid):
    Entrez.email = "firoz@brainstation-23.com"

    """
    Fetches the DOI of a PubMed article given its PubMed ID.

    Args:
    pubmed_id (str): The PubMed ID (PMID) of the article.

    Returns:
    str: The DOI of the article if available, else an empty string.
    """

    retry = 0
    while retry < MAX_RETRY:
        try:
            # Fetch the record using Entrez
            handle = Entrez.efetch(db='pubmed', id=pmid, rettype="medline", retmode="text")
            record = Medline.read(handle)
            handle.close()

            data = record.get('AID', '')
            return data

        except Exception as e:
            print(f'Exception occurred while collecting DOI. Mesage: {e}')
            retry += 1
            sleep(SLEEP_TIME)

    return ''


def _extract_doi(text):
    """
    Extracts the DOI number from the provided text.

    Args:
    text (str): The text containing the DOI.

    Returns:
    str: The extracted DOI number.
    """

    # Regular expression pattern to match DOI
    pattern = r'\b10\.\d+\/[^\s]+\b'

    # Find all matches of the pattern in the text
    dois = re.findall(pattern, text)

    # Return the first match, if found
    if dois:
        return dois[0]
    else:
        return ''


def fetch_doi(pmid):
    try:
        pubmed_doi = _fetch_pubmed_doi(pmid)
        doi = ''
        for d in pubmed_doi:
            if 'doi' in d or 'DOI' in d:
                doi = d
                break
        actual_doi = _extract_doi(doi)
        return actual_doi
    except Exception as e:
        print(f'Exception occurred while extracting DOI. Message: {e}')
        return ''


def fetch_publication_date(pmid):
    Entrez.email = "firoz@brainstation-23.com"

    retry = 0
    while retry < MAX_RETRY:
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
            record = Entrez.read(handle)
            handle.close()

            try:
                pub_date = record['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleDate']
                # Check if publication date is available in YYYY-MM-DD format
                if pub_date[0]['Year'] and pub_date[0]['Month'] and pub_date[0]['Day']:
                    publication_date = f"{pub_date[0]['Year']}-{pub_date[0]['Month']}-{pub_date[0]['Day']}"
                else:
                    # If day is not available, return YYYY-MM
                    publication_date = f"{pub_date[0]['Year']}-{pub_date[0]['Month']}"
            except:
                # If publication date is not available in structured format, try to get it from history
                try:
                    pub_date_str = record['PubmedArticle'][0]['MedlineCitation']['DateCompleted']
                    publication_date = f"{pub_date_str['Year']}-{pub_date_str['Month']}-{pub_date_str['Day']}"
                except:
                    publication_date = ''  # Publication date not available

            return publication_date
        except Exception as e:
            print(e)
            retry += 1
            sleep(SLEEP_TIME)

    return ''


def fetch_publication_types(pmid):
    Entrez.email = "firoz@brainstation-23.com"

    retry = 0
    while retry < MAX_RETRY:
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
            record = Entrez.read(handle)
            handle.close()

            publication_types = []
            try:
                pub_types = record['PubmedArticle'][0]['MedlineCitation']['Article']['PublicationTypeList']
                for pub_type in pub_types:
                    publication_types.append(pub_type)
            except Exception as e:
                print(e)
                publication_types = ['']

            p_types = []

            for p_type in publication_types:
                p_types.append(p_type.title())

            types = ' | '.join(publication_types)
            return types
        except Exception as e:
            print(e)
            retry += 1
            sleep(SLEEP_TIME)

    return ''


def fetch_journal(pmid):
    Entrez.email = "firoz@brainstation-23.com"

    retry = 0
    while retry < MAX_RETRY:
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
            record = Entrez.read(handle)
            handle.close()

            try:
                journal_info = record['PubmedArticle'][0]['MedlineCitation']['Article']['Journal']
                journal_name = journal_info['Title']
                journal_abbr = journal_info['ISOAbbreviation']
                journal = f"{journal_name} ({journal_abbr})"
            except:
                journal = ''  # Journal information not available

            return journal
        except Exception as e:
            print(e)
            retry += 1
            sleep(SLEEP_TIME)

    return ''


def fetch_author_affiliations(pmid):
    Entrez.email = "firoz@brainstation-23.com"

    retry = 0
    while retry < MAX_RETRY:
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
            record = Entrez.read(handle)
            handle.close()

            affiliations = {}

            articles = record.get('PubmedArticle', [])
            for article in articles:
                medline_citation = article.get('MedlineCitation', {})
                article_data = medline_citation.get('Article', {})
                author_list = article_data.get('AuthorList', [])
                for author in author_list:
                    first_name = author.get('ForeName', '')
                    last_name = author.get('LastName', '')
                    full_name = ' '.join(filter(None, [first_name, last_name]))
                    try:
                        affiliation_info = author.get('AffiliationInfo', [])
                        affiliations[full_name] = ' | '.join([aff.get('Affiliation', '') for aff in affiliation_info])
                    except Exception as e:
                        print(e)
                        affiliations[full_name] = ''  # Affiliations not available

            return affiliations
        except Exception as e:
            print("Error while fetching or processing PubMed record:", e)
            retry += 1
            sleep(SLEEP_TIME)

    print("Failed to fetch author affiliations after maximum retries.")
    return None


def fetch_affiliations_by_name(pmid, pubmed_name: str):
    Entrez.email = "firoz@brainstation-23.com"

    retry = 0
    while retry < MAX_RETRY:
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
            record = Entrez.read(handle)
            handle.close()

            affiliations = {}

            articles = record.get('PubmedArticle', [])
            for article in articles:
                medline_citation = article.get('MedlineCitation', {})
                article_data = medline_citation.get('Article', {})
                author_list = article_data.get('AuthorList', [])
                for author in author_list:
                    first_name = author.get('ForeName', '')
                    last_name = author.get('LastName', '')
                    full_name = ' '.join(filter(None, [first_name, last_name]))
                    try:
                        affiliation_info = author.get('AffiliationInfo', [])
                        affiliations[full_name] = ' | '.join([aff.get('Affiliation', '') for aff in affiliation_info])
                    except Exception as e:
                        print(e)
                        affiliations[full_name] = ''  # Affiliations not available

            try:
                return affiliations[pubmed_name]
            except Exception as e:
                print(e)
                return ''
        except Exception as e:
            print("Error while fetching or processing PubMed record:", e)
            retry += 1
            sleep(SLEEP_TIME)

    print("Failed to fetch author affiliations after maximum retries.")
    return None


def fetch_pubmed_authors(pmid):
    Entrez.email = "firoz@brainstation-23.com"
    authors_list = []

    retry = 0
    while retry < MAX_RETRY:
        try:
            handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
            records = Entrez.read(handle)
            handle.close()

            authors_info = records['PubmedArticle'][0]['MedlineCitation']['Article']['AuthorList']

            if isinstance(authors_info, list):
                for author in authors_info:
                    first_name = author.get('ForeName', '')
                    middle_name = author.get('Initials', '')
                    last_name = author.get('LastName', '')
                    authors_list.append((first_name, middle_name, last_name))
            else:
                first_name = authors_info.get('ForeName', '')
                middle_name = authors_info.get('Initials', '')
                last_name = authors_info.get('LastName', '')
                authors_list.append((first_name, middle_name, last_name))

                """
                
                authors = fetch_pubmed_authors(32442160)
                for author in authors:
                print("First Name:", author[0])
                print("Middle Name:", author[1])
                print("Last Name:", author[2])
                print()
                
                """
            return authors_list
        except Exception as e:
            print(e)
            retry += 1
            sleep(SLEEP_TIME)
            authors_list.clear()

    return authors_list


def retrieve_pubmed_xml(pmid):
    Entrez.email = "firoz@brainstation-23.com"
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
    xml_record = handle.read()
    handle.close()
    return xml_record


# print(retrieve_pubmed_xml(35711451))

def search_pubmed_by_author(a_name, _retmax=1000):
    Entrez.email = "firoz@brainstation-23.com"
    search_term = a_name + "[Author]"
    handle = Entrez.esearch(db="pubmed", term=search_term, retmax=_retmax)  # Change retmax as needed
    record = Entrez.read(handle)
    handle.close()
    return record['IdList']


def search_pubmed_by_name_and_pmid(author_last_name: str, _pmid: int):
    pmids = search_pubmed_by_author(author_last_name)
    _found_pmid = None
    _found_author = None

    for pmid in pmids:
        if int(pmid) == _pmid:
            _found_pmid = pmid
            break
    if _found_pmid:
        p_authors = fetch_pubmed_authors(_found_pmid)
        ratio = 0.0
        for p_author in p_authors:
            p = f'{p_author[0]} {p_author[2]}'.strip()
            new_ratio = utils.fuzzy_matching_ratio(p_author[2], author_last_name, True)
            if new_ratio > ratio:
                ratio = new_ratio
                _found_author = p_author

    return _found_author, _found_pmid


def search_other_pubmed_by_name(author_last_name: str, _pmid=0):
    pmids = search_pubmed_by_author(author_last_name, 1000)

    other_pmids = []

    index = 1

    for pmid in pmids:
        print(f'Index: {index}')
        index = index + 1

        try:
            if int(pmid) != int(_pmid):
                authors = fetch_pubmed_authors(pmid)

                for author in authors:
                    full_name = f'{author[0]} {author[2]}'.strip()
                    ratio = utils.fuzzy_matching_ratio(author_last_name, full_name, True)
                    if ratio > 90:
                        other_pmids.append({
                            'first_name_in_other_pubmed': author[0],
                            'initials_in_other_pubmed': author[1],
                            'last_name_in_other_pubmed': author[2],
                            'full_name_in_other_pubmed': f'{author[0]} {author[2]}'.strip(),
                            'other_pubmed_id': pmid,
                            'other_pubmed_url': f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/'
                        })
        except:
            pass

    return other_pmids


def fetch_author_type(author_name: str, pmid: int):
    authors = fetch_pubmed_authors(pmid)

    first_author = 1
    last_author = len(authors)

    ratio = 0
    count = 0

    for index, author in enumerate(authors):
        try:
            first_name = author[0]
        except Exception as e:
            print(e)
            first_name = ''

        try:
            last_name = author[2]
        except Exception as e:
            print(e)
            last_name = ''

        full_name = f'{first_name} {last_name}'.strip()

        new_ratio = utils.fuzzy_matching_ratio(author_name, full_name, True)
        if new_ratio > ratio:
            ratio = new_ratio
            count = index + 1

    if count == first_author:
        return 'First Author'
    elif count == last_author:
        return 'Last Author'
    else:
        return 'Co Author'


def get_pubmed_link_by_pmid(pmid):
    _pmid = int(pmid)
    return f'https://pubmed.ncbi.nlm.nih.gov/{_pmid}/'


'''
# Example usage

author_name = "Junjie Yin"  # Replace with the author name you want to search for
pmids = search_pubmed_by_author(author_name)
for pmid in pmids:
    print(pmid)
'''
# print(fetch_publication_date(36497573))
