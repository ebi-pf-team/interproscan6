FROM python:3.12
LABEL authors="Matthias Blum (mblum@ebi.ac.uk)"
RUN wget -O idrpred.zip https://github.com/matthiasblum/idrpred/archive/refs/tags/1.0.1.zip
RUN unzip idrpred.zip
RUN pip install idrpred-1.0.1/
RUN rm -r idrpred*
