# create the required directories
mkdir -p CU-RESULTS
mkdir -p CU-REF-RESULTS
mkdir -p CURH-RESULTS
mkdir -p CURH-REF-RESULTS
mkdir -p PLOTS

# remove all the existing data files to create new ones
# find . -name "*.txt" -type f -delete
# find . -name "*.gif" -type f -delete
# find . -name "*.png" -type f -delete
find . -name "*.exe" -type f -delete
find . -name "*.x" -type f -delete

# compile the source file and run
g++ main.cpp -o prg.x
./prg.x

# plot the results according to the mode user used
python plot.py




