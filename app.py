from flask import Flask, request
from biopython_service import get_pubmed_data, get_pubmed_data_by_pmid

# Create a Flask application
app = Flask(__name__)


# Define a route for the "Hello, World!" endpoint
@app.route('/')
def hello_world():
    return 'Hello, World!'


@app.route('/pubmed')
def get_data():
    name = request.args.get('name')
    if not name:
        return "Parameter 'name' is missing or empty."

    pubmed_data = get_pubmed_data(name)
    return pubmed_data

@app.route('/pubmed/pmid')
def get_data_by_pmid():
    pmid = request.args.get('pmid')
    if not pmid:
        return "Parameter 'pmid' is missing or empty."

    pubmed_data = get_pubmed_data_by_pmid(pmid)
    return pubmed_data


# Run the Flask app
if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
