import os
import time
import subprocess

class Runner(object):
    ''' General class for carrying out a run

    :arg templatefilename: Name of template file
    :arg outputdir: Directory which will contain output
    :arg executable: name of executable
    :arg paramdict: Dictionary with parameter (to be substituted when
                    parsing the template file)
    '''
    def __init__(self,templatefilename,outputdir,executable,paramdict):
        self.templatefilename = templatefilename
        self.outputdir = outputdir
        self.executable = executable
        self.paramdict = paramdict
        directory = os.path.dirname(outputdir)
        if not os.path.exists(directory):
            os.makedirs(directory)

        self.t_elapsed = None

    def execute(self,outputfilename,force=False):
        ''' Run executable and return output
        
        :arg outputfilename: Name of output file (will be stored in
                             outputdir)
        '''
        if (not os.path.exists(self.outputdir+'/'+outputfilename) or force):
            tmpfilename = self.outputdir+'/'+'parametersTMP.in'
            with open(tmpfilename,'w') as f:
                with open(self.templatefilename) as ftpl:
                    for line in ftpl.readlines():
                        print (line.rstrip() % self.paramdict,file = f)



            with open(self.outputdir+'/'+outputfilename,'w') as fout:
                t_start = time.time()
                run_call = subprocess.run([self.executable,tmpfilename],
                                          stdout=subprocess.PIPE,
                                          universal_newlines=True)
                t_finish = time.time()
                self.t_elapsed = t_finish - t_start
                print(run_call.stdout,file=fout)
        else:
            self.t_elapsed = 0.0

    def elapsedtime(self):
        ''' Return elapsed time of last call'''
        return self.t_elapsed
