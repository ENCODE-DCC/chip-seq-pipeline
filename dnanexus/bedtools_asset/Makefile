SHELL=/bin/bash -e

all:
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus
	wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
	tar -zxvf bedtools-2.26.0.tar.gz
	cd bedtools2 && make
	sudo cp -r ~/bedtools2/bin/* /usr/local/bin/
	