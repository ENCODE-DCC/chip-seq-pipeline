#!/bin/bash

main() {

    echo "Hours to live: '$hours_to_live'"

    sudo mv /etc/apt/apt.conf.d/99dnanexus /tmp
    sudo add-apt-repository -y ppa:deadsnakes/ppa
    sudo apt-get update
    sudo apt-get -y install python3.4-dev libfreetype6-dev

    source env3/bin/activate
    mkdir idr_test
    cd idr_test
    idr --version
    idr --samples ~/env3/idr/tests/data/peak1 ~/env3/idr/tests/data/peak2 --plot --output-file IDR_test.narrowPeak --log-output-file IDR_test-out.txt

    sleep $(($hours_to_live * 3600))

}
