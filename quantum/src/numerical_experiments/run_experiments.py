from runner import *
import argparse
import time
import re
import pandas
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

executable = '../c/driver.x'

class Experiment():
    def __init__(self,Tfinal,m0,Mlatlist):
        self.Tfinal = Tfinal
        self.m0 = 0.25
        self.Mlatlist = Mlatlist
        self.outputdir = './runs_singlelevel/'
        self.templatefilename = 'parameters_template_singlelevel.tpl'
        self.paramdict = {'TFINAL':self.Tfinal,
                          'M0':self.m0}


    def execute(self,force=False):
        t_start = time.time()
        print ('Starting run at ',self.show_clock())
        print()
        for Mlat in self.Mlatlist:
            print ('Running Mlat = '+('%6d' % Mlat))
            self.paramdict['MLAT'] = Mlat
            runner = Runner(self.templatefilename,
                            self.outputdir,
                            executable,
                            self.paramdict)
            outputfilename = 'output_Mlat_'+str(Mlat)+'.txt'
            runner.execute(outputfilename,force=force)

            t_current = time.time()
            t_elapsed = runner.elapsedtime()
            print('... completed at '+self.show_clock())
            print ('  time spent in this run = '+('%8.3f' % t_elapsed)+' s')
            print ('  total elapsed time = '+('%8.2f' % ((t_current-t_start)/60.))+' m')
            print ('')
            
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


class SingleLevelExperiment(Experiment):
    def __init__(self,Tfinal,m0,Mlatlist,epsilon):
        super().__init__(Tfinal,m0,Mlatlist)
        self.epsilon=epsilon
        self.outputdir = './runs_singlelevel/'
        self.templatefilename = 'parameters_template_singlelevel.tpl'
        self.paramdict['epsilon'] = epsilon

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
        ax.set_xlim(np.sqrt(0.5)*np.array(df['alat'])[-1],
                    np.sqrt(2.0)*np.array(df['alat'])[0])
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
                 color='blue',
                 linestyle='--',
                 label='fit $'+('%6.3f' % np.exp(c0))+'a^{'+('%6.3f' % c1)+'}$')
        for j in range(len(df['alat'])):
            plt.annotate('   '+str(np.array(df['Mlat'])[j]),
                         (np.array(df['alat'])[j],
                          np.array(df['bias'])[j]),
                         va='top',ha='left')

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

class TwoLevelExperiment(Experiment):
    def __init__(self,Tfinal,m0,Mlatlist,nsamples):
        super().__init__(Tfinal,m0,Mlatlist)
        self.nsamples=nsamples
        self.outputdir = './runs_twolevel/'
        self.templatefilename = 'parameters_template_twolevel.tpl'
        self.paramdict['NSAMPLES'] = nsamples

    def analyse(self):
        raw_data = []
        for Mlat in self.Mlatlist:
            outputfilename = 'output_Mlat_'+str(Mlat)+'.txt'
            var_chit, var_dchit, tau_int = self.extract_output(self.outputdir+'/'+outputfilename)
            raw_data.append([Mlat,var_chit, var_dchit, tau_int])
        df = pandas.DataFrame(raw_data,columns = ['Mlat','var_chit' , 'var_dchit', 'tau_int'])
        df['alat'] = 1./df['Mlat']*self.Tfinal
        print (df)
        plt.clf()
        ax = plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(np.sqrt(0.5)*np.array(df['alat'])[-1],
                    np.sqrt(2.0)*np.array(df['alat'])[0])
        plt.plot(df['alat'],df['var_chit'],
                 linewidth=2,
                 color='red',
                 marker='s',
                 markerfacecolor='red',
                 markeredgewidth=2,
                 markeredgecolor='red',
                 label=r'Var[$\chi_t(a)$]')
        plt.plot(df['alat'],df['var_dchit'],
                 linewidth=2,
                 color='blue',
                 markersize=8,
                 marker='o',
                 markerfacecolor='white',
                 markeredgewidth=2,
                 markeredgecolor='blue',
                 label=r'Var[$\chi_t(a)-\chi_t(2a)$]')
        c1,c0 = np.polyfit(np.log(df['alat']),np.log(df['var_dchit']),deg=1)
        plt.plot(df['alat'],np.exp(c0)*df['alat']**c1,
                 linewidth=2,
                 color='blue',
                 linestyle='--',
                 label='fit $'+('%6.3f' % np.exp(c0))+'a^{'+('%6.3f' % c1)+'}$')
        c0 = np.polyfit(np.array(np.log(df['alat']))[-3:-1],
                        np.array(np.log(df['var_chit']))[-3:-1],deg=0)
        for j in range(len(df['alat'])):
            plt.annotate('   '+str(np.array(df['Mlat'])[j]),
                         (np.array(df['alat'])[j],
                          np.array(df['var_dchit'])[j]),
                         va='top',ha='left')

        plt.plot(df['alat'],0*df['alat']+np.exp(c0),
                 linewidth=2,
                 color='red',
                 linestyle='--',
                 label='fit $'+('%6.3f' % np.exp(c0))+'$')

        ax.set_xlabel('lattice spacing $a$')
        ax.set_ylabel('variance')

        plt.legend(loc='lower right')
        plt.savefig('variance.pdf',bbox_inches='tight')

    def extract_output(self,filename):
        with open(filename) as f:
            for line in f.readlines():
                # variance of QoI itself (on the fine level)
                m = re.match(' *QoI\[fine\]: Var *= (.*) *',line)
                if m:
                    var_chit = float(m.group(1))
                # variance of difference in QoI
                m = re.match(' *delta QoI: Var *= *(.*) *',line)
                if m:
                    var_dchit = float(m.group(1))
                # integrated autocorrelation time
                m = re.match(' *delta QoI: tau_{int} *= *(.*) *',line)
                if m:
                    tau_int = float(m.group(1))

        return var_chit, var_dchit, tau_int

    
if (__name__ == '__main__'):

    parser = argparse.ArgumentParser(description='Numerical experiments.')
    parser.add_argument('--bias', action='store_true',
                        default=False,
                        help='generate bias data by running single level experiment')

    parser.add_argument('--variance', action='store_true',
                        default=False,
                        help='generate variance data by running two-level experiment')

    Tfinal = 4.0
    m0 = 0.25
    Mlatlist = (32,64,128,256,512,1024)
    print ('Tfinal = ',Tfinal)
    print ('m0 = ',m0)
    args = parser.parse_args()
    # Run single-level experiment to calculate bias
    if args.bias:
        epsilon = 1.E-4
        print ('epsilon = ',epsilon)
        experiment = SingleLevelExperiment(Tfinal,m0,Mlatlist,epsilon)
        experiment.execute()
        experiment.analyse()

    # Run two-level experiment to calculate variance
    if args.variance:
        nsamples = 100000
        print ('nsamples = ',nsamples)
        experiment = TwoLevelExperiment(Tfinal,m0,Mlatlist,nsamples)
        experiment.execute()
        experiment.analyse()
