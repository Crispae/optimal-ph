FROM continuumio/miniconda3


ADD environment.yml .
RUN conda env create -f environment.yml

COPY src/model.pkl .
COPY src/EnzymePh .
COPY src/predict.py .

 ENV PATH=/root/.local:$PATH
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "delta", "python", "predict.py", "--input_csv", "input.csv"]