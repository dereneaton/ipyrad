
## simulate test_rad1
simrrls -o sim_rad1 -f rad -dm 10 -ds 2

## simulate test_rad2
simrrls -o sim_rad2 -f rad -dm 5 -ds 5

## simulate test_pairddrad
simrrls -o sim_pairddrad -f pairddrad -dm 10 -ds 2

## simulate test_pairddrad w/ merged reads
simrrls -o sim_mergepairddrad -f pairddrad -dm 10 -ds 2 -i1 -50 -i2 50 -L 50

