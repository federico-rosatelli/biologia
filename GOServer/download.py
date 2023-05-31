import requests
from clint.textui import progress

url = 'https://drive.google.com/file/d/1NFQqeo4af6It3ieXOhe1EV5hxh7T3Alp/view?usp=sharing'

r = requests.get(url, stream=True)

with open("MongoDumpAll.zip", "wb") as file:
    total_length = int(r.headers.get('content-length'))
    for ch in progress.bar(r.iter_content(chunk_size = 2391975), expected_size=(total_length/1024) + 1):
        if ch:
            file.write(ch)