 #/bin/bash
port=$1
name=$2
datapath=$3
pass=$4
cpus=$5

docker run --tty --interactive  -d -p $port:8787 --cpus=$cpus --name $name -v $datapath:/home/rstudio/tidycell/data/ -e PASSWORD=$pass tidycell
