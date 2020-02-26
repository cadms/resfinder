FROM python:3.8.1-slim

RUN useradd -u 9999 exec

WORKDIR /app

COPY db db/
COPY src src/
COPY tests tests/
COPY vendors/ncbi-blast-2.10.0+/bin /bin
COPY vendors/kma kma/
COPY requirements.txt .
RUN mkdir /app/results

RUN apt update
RUN apt install -y \
  gcc \
  libz-dev \
  make \
  vim

RUN cd kma && make
RUN chmod +x kma
RUN mv kma/kma /bin
RUN rm -rf /app/kma

RUN pip3 install -r requirements.txt

RUN ln -s /app/src/resfinder.py /bin/resfinder

RUN chown -R exec:exec /app

VOLUME /app/results

USER exec

CMD ["/bin/sh"]
