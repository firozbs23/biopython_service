from fuzzywuzzy import fuzz
import sys

from datetime import datetime


def fuzzy_matching_ratio(str1, str2, ignore_case=False):
    if ignore_case:
        str1 = str1.lower()
        str2 = str2.lower()
    ratio = fuzz.ratio(str1, str2)
    return ratio


def read_int_value_from_text_file(file_path=None, default_value=None):
    if file_path is None or len(file_path) == 0:
        file_path = 'sample_file.txt'
    try:
        with open(file_path, 'r') as file:
            value = int(file.read().strip())
            return value
    except FileNotFoundError as e:
        print(f'Unable to read {file_path} file. Exception: {e}', file=sys.stderr)
        return default_value
    except ValueError as e:
        print(f'Error occurred when reading {file_path} file. Exception: {e}', file=sys.stderr)
        return default_value


def write_value_to_file(file_path, value):
    if file_path is None or len(file_path) == 0:
        file_path = 'sample_file.txt'
    try:
        with open(file_path, 'w') as file:
            file.write(str(value))
        return True, 'Saved'
    except Exception as e:
        return False, e


def get_empty_if_null(original_value):
    original_value = str(original_value).strip()
    string_value = original_value if original_value.lower() != 'nan' else ''
    return string_value


def is_empty(input_text):
    new_text = str(input_text).strip()
    new_text = new_text if new_text.lower() != 'nan' else ''
    if len(new_text) > 0:
        return False
    return True


def get_current_time_omnizia_format():
    current_datetime = datetime.utcnow()
    formatted_datetime = current_datetime.strftime("%Y-%m-%d %H:%M:%S.%f+00")
    return formatted_datetime
