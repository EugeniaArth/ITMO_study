FROM continuumio/miniconda3:24.1.2-0

COPY envs/all.yml /tmp/all.yml
RUN conda env update -n base -f /tmp/all.yml && conda clean -afy

COPY bin/plot_coverage.py /usr/local/bin/plot_coverage.py
RUN chmod +x /usr/local/bin/plot_coverage.py

WORKDIR /data
