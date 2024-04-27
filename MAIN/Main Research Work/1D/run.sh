# make some required directories
mkdir -p ENV
mkdir -p RESULTS
mkdir -p PLOTS

# remove the old and useless files
find . -name "*.exe" -type f -delete
find . -name "*.x" -type f -delete
find . -name "*.txt" -type f -delete

# run the program 
g++ main.cpp -o prg.x
./prg.x

# plot the results using a python script
python plot.py