from time import sleep
from openai import OpenAI

# api_key = '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'
client = OpenAI(api_key=api_key)


def _ask_gpt(input_text):
    chat_completion = client.chat.completions.create(
        messages=[
            {
                "role": "user",
                "content": f"{input_text}",
            }
        ],
        model="gpt-3.5-turbo",
    )
    return chat_completion.choices[0].message.content


def ask_chat_gpt(input_data):
    retry = 0
    max_retry = 5
    sleep_time = 5
    error_msg = ''
    result = None

    while retry < max_retry:
        try:
            result = _ask_gpt(input_data)
            if len(str(result)) > 0:
                retry = max_retry
        except Exception as e:
            error_msg = str(e)
            retry += 1
            sleep(sleep_time)

    if result:
        return result, None
    return None, error_msg


data = """
What are the country name and country iso2 code for this location: ` TEMPLE FACULTY PRACTICE PLAN INC
7600 Central Ave, Philadelphia, PA 19111` ? I expect answer in json formet like `{"country" : "United States", "iso2" : "US"}` and I need only json data not other things.
"""

# print(ask_chat_gpt(data))
