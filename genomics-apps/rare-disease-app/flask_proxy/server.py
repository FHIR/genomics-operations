# Import flask and cors handling
from flask import Flask
from flask_cors import CORS, cross_origin
import requests

# Initializing flask app

app = Flask(__name__)
CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'
# Route for getting article information from google search


@app.route('/<url>')
@cross_origin()
def getFHIRResponse(url):
    url = url.replace("@!abab@!", "$")
    url = url.replace("@", "/")
    url = url.replace("!", "?")

    print(url)

    headers = {'Accept': 'application/json'}

    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        return response.json()
    else:
        print('Error: status code='+str(response.status_code))


# Running app
if __name__ == '__main__':
    app.run(port=8000, debug=True)
