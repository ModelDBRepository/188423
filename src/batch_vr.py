# Jan 2015
#
# timothee.masquelier@alum.mit.edu
#
# This code was used in: Masquelier T, Portelli G and Kornprobst P (2016). Microsaccades enable efficient synchrony-based coding in the retina: a simulation study. Scientific Reports. 
#
# This is a Python script to launch several threads of Virtual Retina. This is useful if multiple cores are available.
# For example batch_vr.py -i 0 -f 8 -s 1 will launch 8 threads, each thread processing the frames in ../data/frame/file_name_# and saving the spikes in ../data/fr#

import subprocess
from numpy import *
import sys
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:f:s:', ['initial=', 'final=','step='])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
#    usage()
    sys.exit(2)
#opts, extraparams = getopt.getopt(sys.argv[1:],'r:w:',['randomSeed='])
# starts at the second element of argv since the first one is the script name
# extraparms are extra arguments passed after all option/keywords are assigned
# opts is a list containing the pair "option"/"value"
for o,p in opts:
  if o in ['-i','--initial']:
     initial = double(p)
  elif o in ['-f','--final']:
     final = double(p)
  elif o in ['-s','--step']:
     step = double(p)

        
import multiprocessing
import subprocess

def calculate(alpha):
    

    subprocess.call('mkdir ../data/fr_' + '%01d' % alpha , shell=True)

    subprocess.call('~/retina_package-2.2.2/VirtualRetina-2.2.2/bin/Retina -i ../img/frame/file_name_' + '%01d' % alpha + '.txt  -ret ./human.parvo.xml -r 5 -outD ../data/fr_' + '%01d' % alpha + '/ -nodisp',shell=True)

    #speed = 16
    #str = 'VirtualRetina-2.2.2/bin/Retina VirtualRetina-2.2.2/test/stim/sliding_edge/*.pgm  -bf ' + '%d'%(256-256/speed) + ' -ef ' +'%d'%(256+256/speed)+ ' -ret VirtualRetina-2.2.2/test/retina_files/human.parvo.1cell.xml -r ' +'%d'%speed+  ' -outD VirtualRetina-2.2.2/tmp/fr_' + '%01d' % alpha + '/ -nodisp';
    #str = 'VirtualRetina-2.2.2/bin/Retina VirtualRetina-2.2.2/test/stim/sliding_edge_infinite_speed/*.pgm -ret VirtualRetina-2.2.2/test/retina_files/human.parvo.1cell.xml -r 256 -outD VirtualRetina-2.2.2/tmp/fr_' + '%01d' % alpha + '/ -nodisp';
    #str = 'VirtualRetina-2.2.2/bin/Retina VirtualRetina-2.2.2/test/stim/flashes/*.bmp -ret VirtualRetina-2.2.2/test/retina_files/human.parvo.xml -r 256 -outD VirtualRetina-2.2.2/tmp/fr_' + '%01d' % alpha + '/ -nodisp';
    #print str
    #if alpha==0:
    #    str = str + ' -saveCP'
    subprocess.call(str,shell=True)
    

if __name__ == '__main__':
    pool = multiprocessing.Pool(None)
    tasks = arange(initial,final,step)
    results = []
    r = pool.map_async(calculate, tasks, callback=results.append)
    r.wait() # Wait on the results
    print results
