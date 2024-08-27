from flask import Flask, request
from flask_cors import CORS
from biopython_service import get_pubmed_data, get_pubmed_data_by_pmid
from utility.find_from_chat_gpt import find_ai_values

# Create a Flask application
app = Flask(__name__)
CORS(app)


# Define a route for the "Hello, World!" endpoint
@app.route('/')
def hello_world():
    return 'Hello, World!'


@app.route('/pubmed')
def get_data():
    name = request.args.get('name')
    if not name:
        return "Parameter 'name' is missing or empty."

    data = get_pubmed_data(name)
    response = []
    
    for pubmed_data in data:
        try:
            title = pubmed_data['title']
        except:
            title = ''
    
    
        try:
            abstract = pubmed_data['abstract']
        except:
            abstract = ''
    
        
        prompt = f'{title} {abstract}'
        
        ai_data = find_ai_values(prompt=prompt)
        
        response.append ({
            **pubmed_data, **ai_data
        })
    
    return response

@app.route('/pubmed/pmid')
def get_data_by_pmid():
    pmid = request.args.get('pmid')
    if not pmid:
        return "Parameter 'pmid' is missing or empty."

    data = get_pubmed_data_by_pmid(pmid)
    
    response = []
    
    for pubmed_data in data:
        try:
            title = pubmed_data['title']
        except:
            title = ''
    
    
        try:
            abstract = pubmed_data['abstract']
        except:
            abstract = ''
    
        
        prompt = f'{title} {abstract}'
        
        ai_data = find_ai_values(prompt=prompt)
        
        response.append ({
            **pubmed_data, **ai_data
        })
    
    return response


@app.route('/chat-gpt')
def chat_gpt():
    prompt = request.args.get('prompt')
    if not prompt:
        return "Parameter 'prompt' is missing or empty."
    response = find_ai_values(prompt)
    return response


# Run the Flask app
if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
