FROM continuumio/miniconda3:24.1.2-0

COPY envs/python_plot.yml /tmp/python_plot.yml
RUN conda env update -n base -f /tmp/python_plot.yml && conda clean -afy

COPY bin/plot_coverage.py /usr/local/bin/plot_coverage.py
RUN chmod +x /usr/local/bin/plot_coverage.py

WORKDIR /data
