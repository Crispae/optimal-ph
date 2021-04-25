FROM continuumio/miniconda3


ADD environment.yml .
RUN conda env create -f environment.yml

COPY submission/model.pkl .
COPY submission/EnzymePh .

 ENV PATH=/root/.local:$PATH
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "delta", "python", "predictor.py", "--input_csv", "input.csv"]