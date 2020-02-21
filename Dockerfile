FROM ubuntu:14.04

MAINTAINER Jindrich Vimr <jvimr@softeu.com>

RUN echo "1.565.1" > .lts-version-number

RUN apt-get update && apt-get install -y wget git curl zip vim libgomp1

RUN apt-get update && apt-get install -y perl 

#-------------Marc Code -----------------------
RUN if [ ! -d /opt ]; then mkdir /opt; fi
COPY src /opt/marc
RUN chmod +x /opt/marc/*

#----------------------fastANI dependence ---------
RUN wget -O /opt/fastANI.zip https://github.com/ParBLiSS/FastANI/releases/download/v1.3/fastANI-Linux64-v1.3.zip
RUN unzip /opt/fastANI.zip -d /opt/fastANI && chmod +x /opt/fastANI && rm /opt/fastANI.zip

#----------------------Environment set ---------------
WORKDIR /home
ENV PATH $PATH:/opt/marc:/opt/fasttree:/opt/fastANI
#--------------------------------------------------------

#RUN mkdir /opt/fastANI && tar -C /opt/nw -xzvf /fastANI.zip && cd /opt/nw/newick-utils-1.6 && cp src/nw_* /usr/local/bin
####___________________________________________________________________________________________________________________________________
## Installing perl modules
#RUN apt-get install make && curl -L http://cpanmin.us | perl - App::cpanminus
#RUN cpanm SVG
#RUN cpanm Statistics::Basic
#RUN cpanm IO::Tee
#RUN cpanm Bio::SeqIO

###____________________________________________

# Installing NewickTools
#RUN wget -O /opt/newick-utils-1.6.tar.gz http://cegg.unige.ch/pub/newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz 
#RUN mkdir /opt/nw && tar -C /opt/nw -xzvf /opt/newick-utils-1.6.tar.gz && cd /opt/nw/newick-utils-1.6 && cp src/nw_* /usr/local/bin
#_________________________________________________________________________________________________

## Instaling FastTree
#RUN mkdir /opt/fasttree && wget -O /opt/fasttree/FastTree http://www.microbesonline.org/fasttree/FastTree && chmod +x /opt/fasttree/FastTree
#--------------------------------------------------------------------
#WORKDIR /var/www/
#RUN mkdir /var/www/html
#RUN chmod -R 777 /var/www/html

## EvoMining
#RUN git clone https://github.com/chevrm/evomining3
#--------------------------------------------------------------------
#RUN chmod -R 777 /var/www/html/EvoMining

######### PATHS ENVIRONMENT
#RUN chmod +x /var/www/html/EvoMining/cgi-bin/*pl

#CMD ["perl", "startEvoMining.pl"]
