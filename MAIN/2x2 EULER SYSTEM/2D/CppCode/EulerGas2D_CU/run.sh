#* Bash file to run the program in a single command
#* @author - Lakshya Bamne
#* @supervisor - Prof. Naveen Kumar Garg

#! run this file using the following command in a bash terminal (linux)
#! -> ./run.sh

# create the required directories
mkdir -p result
mkdir -p plots
mkdir -p env

# remove all the existing data files to create new ones
find . -name "*.txt" -type f -delete
find . -name "*.gif" -type f -delete
find . -name "*.png" -type f -delete
find . -name "*.exe" -type f -delete
find . -name "*.x" -type f -delete

# compile the source file and run
g++ main.cpp -o prg.x
./prg.x

# plot the results according to the mode user used
python plot.py
