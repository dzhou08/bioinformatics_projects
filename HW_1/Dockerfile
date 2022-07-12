FROM python:3.8-slim-buster

WORKDIR /

COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt
COPY . .

CMD [ "python", "HW_1_biopython.py"]