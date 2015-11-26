# Communitude

Atsushi Miyauchi and Yasushi Kawase,
What Is a Network Community?: A Novel Quality Function and Detection Algorithms, CIKM2015
http://dl.acm.org/citation.cfm?id=2806555&CFID=733215621&CFTOKEN=42619121



## Compiling
`g++ -std=c++11 generate.cpp -o generate`

`g++ -std=c++11 communitude.cpp -o communitude`


## Usage
### Single community model
    ./generate -n 1000 -c 100 -in 0.4 -out 0.01 -single
    ./communitude data/network.dat -t data/group.dat -o result.txt -r 10    

### l-partition model
    ./generate -n 1000 -c 100 -in 0.4 -out 0.01 -multi
    ./communitude data/network.dat -t data/group.dat -f data/fix.dat -o result.txt -r 10    


