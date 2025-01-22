FROM continuumio/anaconda3:2022.05

VOLUME /input
VOLUME /output
VOLUME /data

# set the working directory
WORKDIR /app

# copy the environment file first, to avoid re-running conda install on every code change
COPY environment.yml ./
RUN conda env create -f environment.yml
# mitigate samtools issue
RUN ln -s /opt/conda/envs/sniffles-genomics/lib/libcrypto.so.1.1 /opt/conda/envs/sniffles-genomics/lib/libcrypto.so.1.0.0

# make RUN commands use the new environment
# SHELL ["conda", "run", "-n", "sniffles-genomics", "/bin/bash", "-c"]
# ENV CONDA_ENV="conda run -n sniffles-genomics"
ENV PATH="/opt/conda/envs/sniffles-genomics/bin:$PATH"
# RUN which python && python --version
ENV PYTHONUNBUFFERED=1

# copy the rest of the files
COPY README.md ./
COPY app.py ./

# set environment variables
ENV SERVER_PORT=8050
ENV DEBUG_MODE=False
ENV PROCESS_INTERACTIVE_MODE=False
ENV PROCESS_ONLY_NEW_FILES=True

# expose the server port
EXPOSE ${SERVER_PORT}

### everything above this runs only when building the image
### everything below this runs when the container is started

# run the app
CMD [ "python", "/app/app.py" ]
# CMD [ "python", "/app/watch.py" ]