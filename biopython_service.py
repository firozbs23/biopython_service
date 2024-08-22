import re
from time import sleep
import utility.utils as utils
import pycountry
from fuzzywuzzy import fuzz
from datetime import datetime
import pytz
import concurrent.futures



from Bio import Entrez
from Bio import Medline

MAX_RETRY = 5
SLEEP_TIME = 3
TIMEOUT = 120  # Timeout in seconds



def fuzzy_matching_ratio(str1, str2, ignore_case=False):
    if ignore_case:
        str1 = str1.lower()
        str2 = str2.lower()
    ratio = fuzz.ratio(str1, str2)
    return ratio


# Variable to store the currently cached PubMed record
current_cached_record = None
current_cached_pmid = None



def fetch_pubmed_record(pubmed_id):
    global current_cached_record
    global current_cached_pmid

    # Check if the requested PubMed ID is the same as the one currently cached
    if pubmed_id == current_cached_pmid:
        return current_cached_record

    def fetch():
        Entrez.email = 'firoz@brainstation-23.com'  # Set your email here
        handle = Entrez.efetch(db='pubmed', id=pubmed_id, rettype='medline', retmode='text')
        records = Medline.parse(handle)
        record = next(records)
        handle.close()
        return record

    # If the requested PubMed ID is different, try fetching the record with retries
    retry = 0
    while retry < MAX_RETRY:
        try:
            sleep(SLEEP_TIME)
            with concurrent.futures.ThreadPoolExecutor() as executor:
                future = executor.submit(fetch)
                record = future.result(timeout=TIMEOUT)

            current_cached_record = record
            current_cached_pmid = pubmed_id
            return record

        except concurrent.futures.TimeoutError:
            print(f'Timeout occurred while fetching PubMed record. Retrying... (Attempt {retry + 1}/{MAX_RETRY})')
            retry += 1
            sleep(SLEEP_TIME)

        except Exception as e:
            print(f'Exception occurred while fetching PubMed record. Retrying... (Attempt {retry + 1}/{MAX_RETRY})')
            print(f'Error: {e}')
            retry += 1
            sleep(SLEEP_TIME)

    # If maximum retries reached without success, return None
    print('Failed to fetch PubMed record after maximum retries.')
    return None

def fetch_pubmed_record_old(pubmed_id):
    global current_cached_record
    global current_cached_pmid

    # Check if the requested PubMed ID is the same as the one currently cached
    if pubmed_id == current_cached_pmid:
        return current_cached_record

    # If the requested PubMed ID is different, try fetching the record with retries
    retry = 0
    while retry < MAX_RETRY:
        try:
            sleep(SLEEP_TIME)
            
            Entrez.email = 'firoz@brainstation-23.com'  # Set your email here
            handle = Entrez.efetch(db='pubmed', id=pubmed_id, rettype='medline', retmode='text')
            # record = handle.read()
            records = Medline.parse(handle)
            record = next(records)
            handle.close()

            current_cached_record = record
            current_cached_pmid = pubmed_id
            return record
        
        except TimeoutError:
            print(f'Timeout occurred while fetching PubMed record. Retrying... (Attempt {retry + 1}/{MAX_RETRY})')
            retry += 1
            sleep(SLEEP_TIME)
        
        except Exception as e:
            print(f'Exception occurred while fetching PubMed record. Retrying... (Attempt {retry + 1}/3)')
            print(f'Error: {e}')
            retry += 1
            sleep(SLEEP_TIME)

    # If maximum retries reached without success, return None
    print('Failed to fetch PubMed record after maximum retries.')
    return None


# Variable to store the currently cached PubMed record
current_cached_record_xml = None
current_cached_pmid_xml = None


def fetch_pubmed_record_xml_old(pubmed_id):
    global current_cached_record_xml
    global current_cached_pmid_xml

    # Check if the requested PubMed ID is the same as the one currently cached
    if pubmed_id == current_cached_pmid_xml:
        return current_cached_record_xml

    # If the requested PubMed ID is different, try fetching the record with retries
    retry = 0
    while retry < MAX_RETRY:
        try:
            Entrez.email = "firoz@brainstation-23.com"
          
            handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="xml")
            record = Entrez.read(handle)
            handle.close()

            current_cached_record_xml = record
            current_cached_pmid_xml = pubmed_id
            return record
        
        except TimeoutError:
            print(f'Timeout occurred while fetching PubMed record. Retrying... (Attempt {retry + 1}/{MAX_RETRY})')
            retry += 1
            sleep(SLEEP_TIME)

        except Exception as e:
            print(f'Exception occurred while fetching PubMed record. Retrying... (Attempt {retry + 1}/3)')
            print(f'Error: {e}')
            retry += 1
            sleep(SLEEP_TIME)

    # If maximum retries reached without success, return None
    print('Failed to fetch PubMed record after maximum retries.')
    return None


def fetch_xml(pubmed_id):
    Entrez.email = "firoz@brainstation-23.com"
    handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="xml")
    record = Entrez.read(handle)
    handle.close()
    return record

def fetch_pubmed_record_xml(pubmed_id):
    global current_cached_record_xml
    global current_cached_pmid_xml

    # Check if the requested PubMed ID is the same as the one currently cached
    if pubmed_id == current_cached_pmid_xml:
        return current_cached_record_xml

    # If the requested PubMed ID is different, try fetching the record with retries
    retry = 0
    while retry < MAX_RETRY:
        try:
            sleep(SLEEP_TIME)
            with concurrent.futures.ThreadPoolExecutor() as executor:
                future = executor.submit(fetch_xml, pubmed_id)
                record = future.result(timeout=TIMEOUT)

            current_cached_record_xml = record
            current_cached_pmid_xml = pubmed_id
            return record

        except concurrent.futures.TimeoutError:
            print(f'Timeout occurred while fetching PubMed record. Retrying... (Attempt {retry + 1}/{MAX_RETRY})')
            retry += 1
            sleep(SLEEP_TIME)

        except Exception as e:
            print(f'Exception occurred while fetching PubMed record. Retrying... (Attempt {retry + 1}/{MAX_RETRY})')
            print(f'Error: {e}')
            retry += 1
            sleep(SLEEP_TIME)

    # If maximum retries reached without success, return None
    print('Failed to fetch PubMed record after maximum retries.')
    return None


def get_pubmed_record(pubmed_id):
    return fetch_pubmed_record(pubmed_id)


def get_pubmed_record_xml(pubmed_id):
    return fetch_pubmed_record_xml(pubmed_id)


def get_empty_string():
    return ''


def get_pubmed_title(pubmed_id):
    record = get_pubmed_record(pubmed_id)
    if record is None:
        return get_empty_string()
    return record.get('TI', '')


def get_pubmed_abstract(pubmed_id):
    record = get_pubmed_record(pubmed_id)
    if record is None:
        return get_empty_string()
    return record.get('AB', '')


def get_pubmed_pmcid(pubmed_id):
    record = get_pubmed_record(pubmed_id)
    if record is None:
        return get_empty_string()
    return record.get('PMC', '')


def get_pubmed_doi(pubmed_id):
    record = get_pubmed_record(pubmed_id)
    if record is None:
        return get_empty_string()

    try:
        pubmed_doi = record.get('AID', '')
        doi = ''
        for d in pubmed_doi:
            if 'doi' in d or 'DOI' in d:
                doi = d
                break

        # Regular expression pattern to match DOI
        pattern = r'\b10\.\d+\/[^\s]+\b'
        dois = re.findall(pattern, doi)
        if dois:
            return dois[0]
        else:
            return get_empty_string()
    except Exception as e:
        print(f'Exception: {e}')
        return get_empty_string()


def get_pubmed_publication_date(pubmed_id):
    record = get_pubmed_record_xml(pubmed_id)
    if record is None:
        return get_empty_string()

    publication_date = ''
    try:
        publication_date = ''
        pub_date = record['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleDate']
        if pub_date[0]['Year'] and pub_date[0]['Month'] and pub_date[0]['Day']:
            # Check if publication date is available in YYYY-MM-DD format
            publication_date = f"{pub_date[0]['Year']}-{pub_date[0]['Month']}-{pub_date[0]['Day']}"
        else:
            # If day is not available, return YYYY-MM
            publication_date = f"{pub_date[0]['Year']}-{pub_date[0]['Month']}"
            if len(publication_date) == 7:
                publication_date = publication_date + '-01'
            elif len(publication_date) == 4:
                publication_date = publication_date + '-01-01'

    except Exception as e:
        print(f'Exception: {e}')
        try:
            pub_date_str = record['PubmedArticle'][0]['MedlineCitation']['DateCompleted']
            publication_date = f"{pub_date_str['Year']}-{pub_date_str['Month']}-{pub_date_str['Day']}"
        except Exception as e:
            print(f'Exception: {e}')
            pass

    return publication_date


def get_pubmed_publication_types(pubmed_id):
    record = get_pubmed_record_xml(pubmed_id)
    if record is None:
        return get_empty_string()

    try:
        publication_types = []
        pub_types = record['PubmedArticle'][0]['MedlineCitation']['Article']['PublicationTypeList']
        for pub_type in pub_types:
            publication_types.append(pub_type)
        p_types = []
        for p_type in publication_types:
            p_types.append(p_type.title())
        types = ' | '.join(publication_types)
        return types
    except Exception as e:
        print(f'Exception: {e}')

    return get_empty_string()


def get_pubmed_journal(pubmed_id):
    record = get_pubmed_record_xml(pubmed_id)
    if record is None:
        return get_empty_string()
    journal = ''
    
    try:
        journal_info = record['PubmedArticle'][0]['MedlineCitation']['Article']['Journal']
        journal_name = journal_info['Title']
        journal_abbr = journal_info['ISOAbbreviation']
        journal = f'{journal_name} ({journal_abbr})'
    except Exception as e:
        print(f'Exception: {e}')

    return journal


def get_pubmed_affiliations(pubmed_id):
    record = get_pubmed_record_xml(pubmed_id)
    if record is None:
        return get_empty_string()

    try:
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
        print(f'Exception: {e}')
        return get_empty_string()
    
    

def get_pubmed_issn(pubmed_id):
    record = get_pubmed_record(pubmed_id)
    if record is None:
        return get_empty_string()

    try:
        # Fetch the article details from PubMed
        # handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="medline", retmode="text")
        # records = Medline.parse(handle)
        # record = next(records)
        # handle.close()
        
        # Extract the ISSN number
        issn = record.get('IS', [])
        return issn
        
    except Exception as e:
        print(f'Exception: {e}')
        return get_empty_string()


def get_pubmed_affiliation_by_name(pubmed_id, author_name):
    record = get_pubmed_record_xml(pubmed_id)
    if record is None:
        return get_empty_string()
    try:
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
                    affiliations[full_name] = ''
        return affiliations[author_name]
    except Exception as e:
        print(f'Exception: {e}')
        return get_empty_string()


def get_pubmed_all_authors(pubmed_id):
    record = get_pubmed_record_xml(pubmed_id)
    if record is None:
        return get_empty_string()

    authors_list = []
    try:
        authors_info = record['PubmedArticle'][0]['MedlineCitation']['Article']['AuthorList']
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
    except Exception as e:
        print(f'Exception: {e}')

    return authors_list


def get_pubmed_idlist_by_author_name(a_name, _retmax=5):
    Entrez.email = "firoz@brainstation-23.com"
    search_term = a_name + "[Author]"
    handle = Entrez.esearch(db="pubmed", term=search_term, retmax=_retmax)  # Change retmax as needed
    record = Entrez.read(handle)
    handle.close()
    return record['IdList']


def get_authors_by_name_and_pubmed_id(author_last_name: str, _pmid: int):
    pmids = get_pubmed_idlist_by_author_name(author_last_name)
    _found_pmid = None
    _found_author = None

    for pmid in pmids:
        if int(pmid) == _pmid:
            _found_pmid = pmid
            break
    if _found_pmid:
        p_authors = get_pubmed_all_authors(_found_pmid)
        ratio = 0.0
        for p_author in p_authors:
            p = f'{p_author[0]} {p_author[2]}'.strip()
            new_ratio = utils.fuzzy_matching_ratio(p_author[2], author_last_name, True)
            if new_ratio > ratio:
                ratio = new_ratio
                _found_author = p_author

    return _found_author, _found_pmid


def search_other_pubmed_by_name(author_last_name: str, _pmid=0):
    pmids = get_authors_by_name_and_pubmed_id(author_last_name, 1000)

    other_pmids = []

    index = 1

    for pmid in pmids:
        print(f'Index: {index}')
        index = index + 1

        try:
            if int(pmid) != int(_pmid):
                authors = get_pubmed_all_authors(pmid)

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


def get_pubmed_author_type(author_name: str, pmid: int):
    authors = get_pubmed_all_authors(pmid)

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


def get_author_name(pubmed_id, name):
    all_authors = get_pubmed_all_authors(pubmed_id)

    first_name = ''
    last_name = ''
    initials = ''

    ratio = 0
    for author in all_authors:
        f_name = author[0]
        l_name = author[2]
        full_name = f'{f_name} {l_name}'.strip()

        new_ratio = fuzzy_matching_ratio(full_name, name, True)
        if new_ratio > ratio:
            ratio = new_ratio
            first_name = f_name
            last_name = l_name
            initials = author[1]

    return first_name, last_name, initials


def get_full_name(pubmed_id, name):
    first_name, last_name, initials = get_author_name(pubmed_id, name)
    return f'{first_name} {last_name}'.strip()


def get_mesh_terms(pubmed_id):
    record = get_pubmed_record(pubmed_id)
    if record is None:
        return get_empty_string()

    try:
        # Extract MeSH terms
        mesh_terms = record.get('MH', [])  # Get MeSH terms from the record, default to empty list if not found
        mesh_terms_str = " | ".join(mesh_terms)  # Join MeSH terms into a single string with separator

        return mesh_terms_str
    except Exception as e:
        print(f'Exception: {e}')
        return get_empty_string()


def get_pubmed_all(pubmed_id, name, search_name):
    
    pubmed_title = get_pubmed_title(pubmed_id)
    print(f'pubmed_title: {pubmed_title}')
    
    first_name, last_name, initials = get_author_name(pubmed_id, name)
    print(f'first_name: {first_name}')
    print(f'last_name: {last_name}')
    print(f'initials: {initials}')
    
    full_name = f'{first_name} {last_name}'.strip()
    print(f'full_name: {full_name}')
    
    journal = get_pubmed_journal(pubmed_id)
    print(f'journal: {journal}')
    
    publication_date = get_pubmed_publication_date(pubmed_id)
    print(f'publication_date: {publication_date}')
    
    abstract = get_pubmed_abstract(pubmed_id)
    print(f'abstract: {abstract}')
    
    author_type = get_pubmed_author_type(full_name, pubmed_id)
    print(f'author_type: {author_type}')
    
    publication_types = get_pubmed_publication_types(pubmed_id)
    print(f'publication_types: {publication_types}')
    
    url = get_pubmed_link_by_pmid(pubmed_id)
    print(f'url: {url}')
    
    affiliations = get_pubmed_affiliation_by_name(pubmed_id, full_name)
    print(f'affiliations: {affiliations}')
    
    pmcid = get_pubmed_pmcid(pubmed_id)
    print(f'pmcid: {pmcid}')
    
    doi = get_pubmed_doi(pubmed_id)
    print(f'doi: {doi}')
    
    mesh_terms = get_mesh_terms(pubmed_id)
    print(f'mesh_terms: {mesh_terms}')
    
    issn = get_pubmed_issn(pubmed_id)
    print(f'issn: {issn}')
    
    aff = ''
    try:
        aff = affiliations[0]
    except Exception as e:
        print(f'Exception: {e}')
        aff = ''
        
    print(f'aff: {aff}')

    data = {
        'publication_id': pubmed_id,
        'title': pubmed_title,
        'journal': journal,
        'publication_date': publication_date,
        'abstract': abstract,
        'hcp_role': author_type,
        'publication_type': publication_types,
        'url': url,
        'affiliations': aff,
        'pmcid': pmcid,
        'doi': doi,
        'mesh_terms': mesh_terms,
        
        'first_name': first_name,
        'last_name': last_name,
        'full_name' : full_name,
        'initials' : initials,
        'issn' : issn,
        'publication_platform' : 'PubMed',
        'search_name' : search_name
    }
    
   # 'timestamp': datetime.now(pytz.utc),

    return data


def get_pubmed_data(name):
    
    ids = get_pubmed_idlist_by_author_name(name)

    print(f'Fetching pubmed data for the name : {name}')

    data_list = []
    for pmid in ids:
        print(f'Processing for pubmed id : {pmid}')

        full_name = get_full_name(pmid, name)
        if full_name != name:
            continue

        data = get_pubmed_all(pmid, full_name, name)
        data_list.append(data)

    return data_list


def get_pubmed_data_by_pmid(pubmed_id):

    print(f'Fetching pubmed data for the pubmed_id : {pubmed_id}')
    
    pubmed_title = get_pubmed_title(pubmed_id)
    print(f'pubmed_title: {pubmed_title}')
    
    journal = get_pubmed_journal(pubmed_id)
    print(f'journal: {journal}')
    
    publication_date = get_pubmed_publication_date(pubmed_id)
    print(f'publication_date: {publication_date}')
    
    abstract = get_pubmed_abstract(pubmed_id)
    print(f'abstract: {abstract}')
    
    publication_types = get_pubmed_publication_types(pubmed_id)
    print(f'publication_types: {publication_types}')
    
    url = get_pubmed_link_by_pmid(pubmed_id)
    print(f'url: {url}')
    
    pmcid = get_pubmed_pmcid(pubmed_id)
    print(f'pmcid: {pmcid}')
    
    doi = get_pubmed_doi(pubmed_id)
    print(f'doi: {doi}')
    
    mesh_terms = get_mesh_terms(pubmed_id)
    print(f'mesh_terms: {mesh_terms}')
    
    issn = get_pubmed_issn(pubmed_id)
    print(f'issn: {issn}')

    data = {
        'publication_id': pubmed_id,
        'title': pubmed_title,
        'journal': journal,
        'publication_date': publication_date,
        'abstract': abstract,
        'publication_type': publication_types,
        'url': url,
        'pmcid': pmcid,
        'doi': doi,
        'mesh_terms': mesh_terms,
        'issn' : issn,
        'publication_platform' : 'PubMed',
    }
    
    return data

# print(get_pubmed_affiliation_by_name(27612954, 'Yulia Burakova'))

# print('================================')
# print(get_pubmed_issn(19896784))