pC-SAC Installation Guide
This guide provides instructions for installing the pC-SAC software, which is used for chromatin structure analysis and computation.

Prerequisites
Before installing pC-SAC, make sure the following dependencies are installed on your system:

CMake
Boost
Eigen

These can typically be loaded on a HPC environment using module commands, like so:

module load cmake
module load boost
module load eigen

Installation
Follow these steps to install pC-SAC:

1. Setup
First, set the main directory where pC-SAC will be installed. Replace /path/to/pC_SAC_G2 with the actual path on your system.

export MAIN_DIR="/path/to/pC_SAC_G2"

2. Download and Install Eigen
Eigen is a dependency for pC-SAC. Download and install it using the following commands:

wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar xvfz eigen-3.4.0.tar.gz
rm eigen-3.4.0.tar.gz
cd eigen-3.4.0
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../ -DINCLUDE_INSTALL_DIR=../include/
make install
cd $MAIN_DIR

3. Compile pC-SAC
Edit the CMakeLists.txt in the pC-SAC directory to include the path to the Eigen include directory:
nano ${MAIN_DIR}/CMakeLists.txt
# Modify the include path to /path/to/pC_SAC_G2/eigen-3.4.0/include/

Now, compile pC-SAC:

mkdir build
cd build
cmake ..
make


Testing
To verify that pC-SAC has been installed correctly, you can run a test using provided example files:

Files Description
Within the test_data directory, you will find several example input files required for running pC-SAC. Please explore these files in order to understand each required input:

int_mat_seg_len.txt: Length of this file represents the length of any chain for a given reconstruction.
interaction_matrix.txt: Initial long-range probabilities.
test.conf: Configuration file with pC-SAC parameters.

cd ${MAIN_DIR}/script/test_data

# Set the LD_LIBRARY_PATH to include the Boost library path
export LD_LIBRARY_PATH="/nfs/sw/boost/boost-1.72.0/lib"

# Load the Boost module
module load boost/1.72.0

# Run the test
${MAIN_DIR}/build/bin/chromatin.sis.coarse -conf ${MAIN_DIR}/script/test.conf -prefix test_pCSAC

