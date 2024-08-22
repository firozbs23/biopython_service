import utility.chat_gpt as gpt


def find_chat_analysis(question: str):
    question_pre_text = 'Please write a summary and highlight all positive and negative aspects about this below text:'
    question_to_ask = f'{question_pre_text}\n{question}'

    response_data, error = gpt.ask_chat_gpt(question_to_ask)
    if response_data:
        return response_data
    return error


def find_tag_type_rational(question: str):
    question_pre_text = 'Process the text below and classify what type of article is it Diagnosis or Disease or Treatment:'
    question_to_ask = f'{question_pre_text}\n{question}'

    response_data, error = gpt.ask_chat_gpt(question_to_ask)
    if response_data:
        return response_data
    return error


def find_tag_type(question: str):
    question_pre_text = 'Process below text and classify what type of article is it Disease or Treatment or Diagnosis. I need response in single word like Disease or Treatment or Diagnosis without fullstop:'
    question_to_ask = f'{question_pre_text}\n{question}'

    response_data, error = gpt.ask_chat_gpt(question_to_ask)
    if response_data:
        return response_data
    return error


def find_tag_category(question: str):
    question_pre_text = 'Can you describe with one word what is this below text about:'
    question_to_ask = f'{question_pre_text}\n{question}'

    response_data, error = gpt.ask_chat_gpt(question_to_ask)
    if response_data:
        return response_data
    return error


def find_tag_value(question: str):
    question_pre_text = 'please tell us the disease in one word mentioned in the below text:'
    question_to_ask = f'{question_pre_text}\n{question}'

    response_data, error = gpt.ask_chat_gpt(question_to_ask)
    if response_data:
        return response_data
    return error


def find_country_and_iso2(affiliation: str):
    if affiliation == '':
        return None
    question_to_ask = f'What are the country name and iso2 code [ISO 3166-1 alpha-2 code] for the text below: \n{affiliation}  Please give me result in json format like : \n {{"country" : "United States", "iso2" : "US"}}'
    response_data, error = gpt.ask_chat_gpt(question_to_ask)
    if response_data:
        return response_data
    return error
