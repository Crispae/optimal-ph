FROM continuumio/miniconda3
ADD environment.yml .
RUN conda env create -f environment.yml


#COPY src/model.pkl .
#COPY src/predict.py .
#COPY src/EnzymePh .

COPY src .




ENV PATH=/root/.local:$PATH
 
ENTRYPOINT ["conda","run","-n","delta","python", "predict.py", "--input_csv", "input.csv"]