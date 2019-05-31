This readme and the code were contributed by Timothee Masquelier
timothee.masquelier@alum.mit.edu
March 2016

This code was used in:

This code was used in: Masquelier T, Portelli G and Kornprobst P (2016). Microsaccades enable efficient synchrony-based coding in the retina: a simulation study. Scientific Reports. 

Feel free to use/modify but please cite us if appropriate.

It is a Matlab code (with mex files)

The retina simulator, developed by Adrien Wohrer, is not included in this archive.
It can be downloaded here:
http://www-sop.inria.fr/neuromathcomp/public/software/virtualretina/
At the bottom of this README you will find help to install it on a Linux machine (provided without warranty).

./src contains the source files
./img contains images
./data contains data files

Here are the following steps to reproduce the paper's main simulation.
The order of magnitude of processing times are given (for 15000 s of biological time).

1) engbert_brownian.m : generates the gaze trajectory
Processing time: minutes

2) generateSaccadicFrames.m : generates the frames from one input image and the gaze trajectory
Processing time: hours

3) generateNameFiles.m : generates n=8 lists of frames for Virtual Retina (this can be useful to launch Virtual Retina on multiple threads)
Processing time: seconds

4) batch_vr.py : a Python script to launch several threads of Virtual Retina. This is useful if multiple cores are available.
For example the command: python batch_vr.py -i 0 -f 8 -s 1 will launch 8 threads.
The output spikes are saved in ./data/fr#/spikes.spk If you don't have Python or multiple cores, you can always process the frame files successively (the command
should be something like: ~/retina_package-2.2.2/VirtualRetina-2.2.2/bin/Retina -i ../data/frame/file_name_#.txt  -ret ./human.parvo.xml -r 5 -outD ../data/fr_# )
Processing time: hours

5) formatSL.m : formats the spikes.spk in mat files adequate for the STDP code.
Processing time: minutes

6) STDP.m : performs the STDP-based learning on the RGC spike trains.
This code involves mex files, which should be compiled first (from Matlab, enter: mex pspKernel.c; mex STDPContinuous.c)
Use plots.m to reconstruct preferred stimuli after learning.
Note: one may need to launch the script several time to reach convergence. The weights are read and saved in ./data/weight*.txt files.
Remove the files to start learning from scratch (i.e. from random synaptic weights).
A convergence index for each neuron is stored in ../data/conv.#.mat (0 means perfect convergence)
Processsing time: minutes (for one execution)


—-
Here are the steps we took to install Virtual Retina on a recent Ubuntu machine:
 * download and extract the retina_package-2.2.2 archive
 * install `libxml2-dev`, `libx11-dev` and `cmake` through the normal package manager
 * after running the download script in the retina_package once, the packages were downloaded, but did not compile successfully
 * in `mvaspike-1.0.16/src/`, in file `stdpf.cc` `NULL` had to be replaced with `0`. (then recompile mvaspike by running `make` in `External_Libraries/mvaspike-1.0.16/`)
 * in `xmlparameters++`, in the file `src/initfunctor.h` `return create(...` had to be replaced with `return this->create(...` (then run `make clean` and `make` in `External_Libraries/xmlparameters++` to recompile xmlparameters)
 * The CImg library had to be downloaded and put in a folder called `External_Libraries/CImg-1.3.4` (since this is the version linked in the VirtualRetina-2.2.2/ext-lib-links folder) 

Running the download_all.bash script again then successfully compiled the virtual retina programs.
