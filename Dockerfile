FROM public.ecr.aws/lambda/python:3.12

RUN dnf install -y tar gzip libgomp && \
    dnf clean all && \
    rm -rf /var/cache/yum


ARG DEPLOY_TYPE=prod
ARG LAMBDA_TASK_ROOT=/var/task
ENV MPLCONFIGDIR=/tmp

RUN pip install "poetry>=1.1"
RUN pip install pysam
RUN pip install boto3
RUN pip install Bio
RUN pip install matplotlib

# RUN apt-get install minimap2
# RUN apt-get install samtools

WORKDIR ${LAMBDA_TASK_ROOT}
COPY poetry.lock pyproject.toml ${LAMBDA_TASK_ROOT}

COPY speedy/artifacts/speedy-0.1.0-cp312-none-any.whl /tmp

RUN poetry config virtualenvs.create false \
  && poetry install $(test "$DEPLOY_TYPE" == prod && echo "--only main") --no-interaction --no-root --no-ansi \
  && poetry run pip install /tmp/speedy-0.1.0-cp312-none-any.whl

COPY main.py ${LAMBDA_TASK_ROOT}
# COPY cleaner_main.py ${LAMBDA_TASK_ROOT}
COPY sequence_manipulation.py ${LAMBDA_TASK_ROOT}
COPY data_retrieval.py ${LAMBDA_TASK_ROOT}
COPY dataset_creation ${LAMBDA_TASK_ROOT}/dataset_creation
COPY data ${LAMBDA_TASK_ROOT}/data
COPY json_to_csv.py ${LAMBDA_TASK_ROOT}
COPY cred.py ${LAMBDA_TASK_ROOT}
COPY dp_path_finder.py ${LAMBDA_TASK_ROOT}

# Contains the reference gene database for the blastn search.
# It was created from NCBI's Pathogen Detection Reference Gene Catalog (https://www.ncbi.nlm.nih.gov/pathogens/refgene/#)
# db version 2023-09-26.1
# Filtered for Amp, Kan and Chlor resistance genes.

CMD [ "app.acceptance" ]


