import utility.chat_gpt as gpt
import re


def clean_text(text):
    # Define a regular expression pattern to match unwanted characters at the beginning and end
    pattern = r"^[^\w]+|[^\w]+$"
    
    # Use re.sub() to substitute the matched characters with an empty string
    cleaned_text = re.sub(pattern, "", text)
    
    return cleaned_text


def word_count(text):
    words = text.split()
    return len(words)




def find_tag_type(question: str):
    question_pre_text = 'Process below text and classify what type of article is it Disease or Treatment or Diagnosis. I need response in single word which is Disease or Treatment or Diagnosis without fullstop:'
    question_to_ask = f'{question_pre_text}\n{question}'

    response_data, error = gpt.ask_chat_gpt(question_to_ask)
    
    retries = 0
    while retries < 10:
        response_data, error = gpt.ask_chat_gpt(question_to_ask)
        if error:
            print(error)
        
        if response_data and word_count(response_data) == 1 and ( 'disease' == response_data.strip().lower() or 'treatment' == response_data.strip().lower() or 'diagnosis' == response_data.strip().lower() ):
            response_data = response_data.strip().replace('.', '')
            return clean_text(response_data)
        
        retries += 1

    return None



def find_tag_category(question: str):
    question_pre_text = 'Can you describe with one word what is this below text about: (Note: Result never have the words i. Disease, ii. Treatment, iii. Diagnosis and If result not possible in one word, then please use maximum three words)'
    question_to_ask = f'{question_pre_text}\n{question}'

    
    retries = 0
    while retries < 10:
        response_data, error = gpt.ask_chat_gpt(question_to_ask)
        if error:
            print(error)
        
        if response_data and word_count(response_data) <= 3:
            return clean_text(response_data)
        
        retries += 1

    return None


def find_tag_value(question: str):
    question_pre_text = 'please tell us the disease in one word mentioned in the below text: (Note: Result never have the words i. Disease, ii. Treatment, iii. Diagnosis and If result not possible in one word, then please use maximum three words)'
    question_to_ask = f'{question_pre_text}\n{question}'

    retries = 0
    while retries < 10:
        response_data, error = gpt.ask_chat_gpt(question_to_ask)
        if error:
            print(error)
        
        if response_data and word_count(response_data) <= 3:
            return clean_text(response_data)
        
        retries += 1

    return None


def find_tag_type_rational(question: str):
    question_pre_text = 'Process the text below and classify what type of article is it Diagnosis or Disease or Treatment: (Note: Result never have the words i. Disease, ii. Treatment, iii. Diagnosis)'
    question_to_ask = f'{question_pre_text}\n{question}'

    response_data, error = gpt.ask_chat_gpt(question_to_ask)
    
    retries = 0
    while retries < 10:
        response_data, error = gpt.ask_chat_gpt(question_to_ask)
        if error:
            print(error)
        
        if response_data and 'disease' != response_data.strip().lower() and 'treatment' != response_data.strip().lower() and 'diagnosis' != response_data.strip().lower():
            response_data = response_data.strip().replace('.', '')
            return clean_text(response_data)
        
        retries += 1

    return None



def find_chat_analysis(question: str):
    question_pre_text = 'Please write a summary and highlight all positive and negative aspects about this below text:'
    question_to_ask = f'{question_pre_text}\n{question}'

    response_data, error = gpt.ask_chat_gpt(question_to_ask)
    if response_data:
        return response_data
    
    print(error)
    return None




def find_country_and_iso2(affiliation: str):
    if affiliation == '':
        return None
    question_to_ask = f'What are the country name and iso2 code [ISO 3166-1 alpha-2 code] for the text below: \n{affiliation}  Please give me result in json format like : \n {{"country" : "United States", "iso2" : "US"}}'
    response_data, error = gpt.ask_chat_gpt(question_to_ask)
    if response_data:
        return response_data
    
    print(error)
    return None



def find_ai_values(prompt):
    tag_type = find_tag_type(prompt)
    if tag_type is None:
        tag_type = ''
    
    tag__category = find_tag_category(prompt)
    if tag__category is None:
        tag__category = ''
    
    tag_value = find_tag_value(prompt)
    if tag_value is None:
        tag_value = ''
    
    tag_type_rational = find_tag_type_rational(prompt)
    if tag_type_rational is None:
        tag_type_rational = ''
    
    chat_analysis = find_chat_analysis(prompt)
    if chat_analysis is None:
        chat_analysis = ''
    
    return {
        'tag_type': tag_type,
        'tag_category': tag__category,
        'tag_value': tag_value,
        'tag_type_rational': tag_type_rational,
        'chat_analysis': chat_analysis
    }