import os

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait


# For Selenium wire driver
# from selenium.webdriver.chrome.service import Service as ChromeService
# from seleniumwire import webdriver


# ======================= Start Webdriver Utilities =========================================


def get_driver(proxy_address=None, headless=True, lang_code='en', incognito=False):
    chrome_options = Options()

    if headless:
        chrome_options.add_argument('--headless')
    if incognito:
        chrome_options.add_argument('--incognito')
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument("--window-size=1920x1080")
    chrome_options.add_argument("--ignore-certificate-errors")
    chrome_options.add_argument('--disable-dev-shm-usage')

    if proxy_address is not None:
        chrome_options.add_argument('--proxy-server=' + proxy_address)

    chrome_options.add_argument(f"--lang={lang_code}")  # Fix language setting

    # chrome_options.add_argument(
    #    "--user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/104.0.5112.102 Safari/537.3")

    chrome_options.add_argument(
        "--user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/119.0.6045.105 Safari/537.3")

    # Set a realistic user agent
    # chrome_options.add_argument(
    #   "--user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/94.0.4606.81 Safari/537.36"
    # )

    chrome_options.add_argument("--disable-gpu")
    chrome_options.page_load_strategy = 'normal'

    # driver = webdriver.Chrome(service=Service('/home/apps/chromedriver'), options=chrome_options)
    driver = webdriver.Chrome(service=Service('utility/chromedriver_linux64/chromedriver'), options=chrome_options)
    driver.maximize_window()

    return driver


def get_seleniumwire_driver(proxy_address=None, headless=True, lang_code='en'):
    chrome_options = Options()
    """
    uncomment the following line if you want to use Chrome as headless mode
    there is no gui in headless mode. Selenium will use chrome internally
    """
    if headless:
        chrome_options.add_argument('--headless')
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument("--window-size=1920x1080")  # Set the desired window size
    chrome_options.add_argument("--ignore-certificate-errors")
    chrome_options.add_argument('--disable-dev-shm-usage')
    chrome_options.add_argument(f"--lang={lang_code}")  # Fix language setting

    chrome_options.add_argument(
        "--user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.6099.129 Safari/537.3")
    chrome_options.add_argument("--disable-gpu")
    chrome_options.page_load_strategy = 'normal'  # Wait until page is fully loaded

    seleniumwire_options = {
        'proxy': {
            'http': 'http://' + str(proxy_address),
            'https': 'https://' + str(proxy_address),
            'no_proxy': 'localhost:127.0.0.1'
        },
    }

    driver = webdriver.Chrome(
        service=Service('utility/chromedriver_linux64/chromedriver'),
        options=chrome_options,
        # seleniumwire options
        seleniumwire_options=seleniumwire_options,
    )

    driver.maximize_window()
    return driver


def wait_until_fully_loaded(driver):
    try:
        driver.implicitly_wait(10)
        wait = WebDriverWait(driver, 10)  # 10 seconds timeout
        wait.until(EC.presence_of_element_located((By.TAG_NAME, "body")))
        return driver.execute_script("return document.readyState") == "complete"
    except Exception as e:
        return e


def scroll_to_element(driver, element):
    # Scroll to the element using JavaScript
    try:
        driver.execute_script("arguments[0].scrollIntoView(true);", element)
        return True
    except Exception as e:
        return e


def wait_until_element_loaded(driver, element_id):
    try:
        # Wait for the element with ID "doctors-list" to be present and visible, with a maximum timeout of 10 seconds
        element = WebDriverWait(driver, 10).until(
            EC.visibility_of_element_located((By.ID, element_id))
        )

        # Perform actions with the fully loaded "doctors-list" element
        print("Element 'doctors-list' is loaded:", element.text)

    except Exception as e:
        print(f'Error occurred while waiting for element: {element_id}. Message: {e}')


def scroll_to_element_by_id(driver, element_id):
    try:
        # Scroll to the element with the specified ID using JavaScript
        script = f"document.getElementById('{element_id}').scrollIntoView(true);"
        driver.execute_script(script)
    except Exception as e:
        print(f'Error occurred while scrolling element by id : {element_id}. Exception: {e}')


def wait_until_xpath_loaded(driver, element_xpath):
    try:
        wait = WebDriverWait(driver, 10)  # 10 seconds timeout, you can adjust as needed
        element = wait.until(EC.presence_of_element_located((By.XPATH, element_xpath)))
    except:
        pass


def take_screenshot(driver, filename='sample.png', location=None):
    """
    Takes a screenshot using Selenium with Chrome driver and saves it to a specific directory.

    :param driver: WebDriver instance (e.g., created using webdriver.Chrome())
    :param filename: Name of the screenshot file (default is 'sample.png')
    :param location: Relative path of the directory where the screenshot will be saved
                     (default is the current working directory)
    """
    try:
        # Set the default location to the current working directory if not provided
        if location is None:
            location = os.getcwd()

        # Navigate to the desired directory
        driver.get("file://" + location)
        wait_until_fully_loaded(driver)

        # Take a screenshot
        driver.save_screenshot(os.path.join(location, filename))
        print(f"Screenshot saved successfully at {os.path.join(location, filename)}")

    except Exception as e:
        print(f"Error: {e}")

# ================= End Webdriver Utilities ===================================
