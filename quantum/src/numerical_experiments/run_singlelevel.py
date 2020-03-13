from runner import *
import time
import re
import pandas
import numpy as np
from matplotlib import pyplot as plt

executable = '../c/driver.x'


class SingleLevelExperiment(object):
    def __init__(self,Tfinal,m0,epsilon,Mlatlist):
        self.Tfinal = Tfinal
        self.m0 = 0.25
        self.epsilon = epsilon
        self.Mlatlist = Mlatlist
        self.outputdir = './runs_singlelevel/'
        self.templatefilename = 'parameters_template_singlelevel.tpl'

    def execute(self,force=False):
        t_start = time.time()
        print ('Starting run at ',self.show_clock())
        print()
        for Mlat in self.Mlatlist:
            print ('Running Mlat = '+('%6d' % Mlat))
            paramdict = {'TFINAL':self.Tfinal,
                         'MLAT':Mlat,
                         'EPSILON':self.epsilon,
                         'M0':self.m0}
            runner = Runner(self.templatefilename,
                            self.outputdir,
                            executable,
                            paramdict)
            outputfilename = 'output_Mlat_'+str(Mlat)+'.txt'
            runner.execute(outputfilename,force=force)

            t_current = time.time()
            t_elapsed = runner.elapsedtime()
            print('... completed at '+self.show_clock())
            print ('  time spent in this run = '+('%8.3f' % t_elapsed)+' s')
            print ('  total elapsed time = '+('%8.2f' % ((t_current-t_start)/60.))+' m')
            print ('')

    def analyse(self):
        raw_data = []
        for Mlat in self.Mlatlist:
            outputfilename = 'output_Mlat_'+str(Mlat)+'.txt'
            chit, dchit = self.extract_output(self.outputdir+'/'+outputfilename)
            raw_data.append([Mlat,chit,dchit])
        df = pandas.DataFrame(raw_data,columns = ['Mlat' , 'chit', 'dchit'])
        df['bias'] = df['chit']-self.chit_exact()
        df['alat'] = 1./df['Mlat']*self.Tfinal
        print (df)
        plt.clf()
        ax = plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('lattice spacing $a$')
        ax.set_ylabel('bias')
        plt.plot(df['alat'],df['bias'],
                 linewidth=2,
                 color='blue',
                 marker='o',
                 markeredgewidth=0,
                 markersize=8,
                 label='measured')
        plt.errorbar(df['alat'],df['bias'],yerr=df['dchit'],
                      linewidth=2,
                      color='blue',label='')
        c1,c0 = np.polyfit(np.log(df['alat']),np.log(df['bias']),deg=1)
        plt.plot(df['alat'],np.exp(c0)*df['alat']**c1,
                 linewidth=2,
                 color='red',
                 label='fit $'+('%5.2f' % np.exp(c0))+'a^{'+('%5.2f' % c1)+'}$')
        plt.legend(loc='lower right')
        plt.savefig('bias.pdf',bbox_inches='tight')

    def extract_output(self,filename):
        with open(filename) as f:
            for line in f.readlines():
                m = re.match(' *Q: Avg \+\/\- Err = (.*) \+\/\- (.*)',line)
                if m:
                    chit = float(m.group(1))
                    dchit = float(m.group(2))
        return chit, dchit

    def chit_exact(self):
        z = 4.*np.pi**2*self.m0/self.Tfinal
        S_num = 0.
        S_denom = 1.
        for Q in range(1,1000):
            S_num += 2.*Q**2*np.exp(-0.5*z*Q**2)
            S_denom += 2.*np.exp(-0.5*z*Q**2)
        return 1./(4.*np.pi**2*m0)*z*S_num/S_denom

    
    def show_clock(self):
        return time.strftime('%a, %d %b %Y %H:%M:%S +0000',time.gmtime())

    
if (__name__ == '__main__'):
    Tfinal = 4.0
    m0 = 0.25
    epsilon = 1.E-4
    Mlatlist = (32,64,128,256,512,1024)
    experiment = SingleLevelExperiment(Tfinal,m0,epsilon,Mlatlist)
    experiment.execute()
    experiment.analyse()
