import requests


def check_internet_connection():
    try:
        # Try to send a GET request to google.com
        response = requests.get("http://www.google.com", timeout=5)
        # If the request was successful, return True
        return response.status_code == 200
    except requests.ConnectionError:
        # If there was a connection error (no internet), return False
        return False
