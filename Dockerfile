FROM rocker/r-base:4.0.2
MAINTAINER "Man Chen" manchen9005@gmail.com

# Install R packages
RUN install2.r --error \
    tidyr \
    dplyr \
    purrr \
    skellam \
    metafor \
    clubSandwich \
    simhelpers \
    plyr


EXPOSE 1410
ENV HTTR_LOCALHOST 0.0.0.0