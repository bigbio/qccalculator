from flask import Flask, request, render_template, jsonify

app = Flask(__name__, template_folder='serving_static')
#~ app = Flask(__name__, static_url_path='/home/walzer')

@app.route('/')
def root():
    message = "Hello there!"
    return render_template('p0.html', message=message)
    #~ return 'Hello, World!
    
@app.route('/echo', methods=['POST'])
def hello():
    if request.is_json:
        print("whip")
        content = request.json
        print content['mzQC']['creationDate']
    return jsonify(request.json)

if __name__ == "__main__":
    app.run(debug=True)

#~ python flask_p0.py