# Dockerfile
FROM debian:buster

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get clean && apt-get update

RUN apt-get install -y \
    aptly \
    gnupg1

COPY yadedev_pub.gpg /root/

COPY yadedev_sec.gpg /root/

COPY .aptly.conf /root/

RUN gpg1 --import /root/yadedev_pub.gpg && gpg1 --import /root/yadedev_sec.gpg

RUN rm /root/yadedev_*

RUN apt-get autoclean && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/cache/apt/archives/*deb
