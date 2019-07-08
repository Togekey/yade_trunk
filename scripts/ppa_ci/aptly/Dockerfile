# Dockerfile
FROM debian:buster

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get clean && apt-get update

RUN apt-get install -y \
    aptly \
    gnupg1 \
    unzip \
    wget

COPY yadedev_pub.gpg /root/

COPY yadedev_sec.gpg /root/

COPY .aptly.conf /root/

RUN gpg1 --import /root/yadedev_pub.gpg && gpg1 --import /root/yadedev_sec.gpg

RUN rm /root/yadedev_*

RUN aptly repo create -distribution=stretch -component=main yadedaily-stretch &&\
    aptly repo create -distribution=buster -component=main yadedaily-buster &&\
    aptly repo create -distribution=bionic -component=main yadedaily-bionic

RUN mkdir -p /root/deb

RUN cd /root/deb &&\
    wget https://gitlab.com/api/v4/projects/10133144/jobs/artifacts/feature%2Fdailypackages/download?job=deb_bionic -O yade.zip &&\
    unzip yade.zip &&\
    aptly repo add yadedaily-bionic deb/*.deb &&\
    aptly repo add yadedaily-bionic deb/*.dsc &&\
    rm -rf /root/deb/* &&\
    AWS_ACCESS_KEY_ID="asdf"  AWS_SECRET_ACCESS_KEY="ASDF" aptly publish repo yadedaily-bionic s3:yadedaily-repo:/debian/

RUN apt-get autoclean && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/cache/apt/archives/*deb
