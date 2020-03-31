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
        self.outputdir = './runs_bias/'
        self.templatefilename = None
        self.paramdict = {'TFINAL':self.Tfinal,
                          'M0':self.m0}
        self.varying_paramdict = {'MLAT':{}}
        for Mlat in Mlatlist:
            self.varying_paramdict['MLAT'][Mlat] = Mlat


    def execute(self,force=False):
        t_start = time.time()
        print ('Starting run at ',self.show_clock())
        print()
        for Mlat in self.Mlatlist:
            print ('Running Mlat = '+('%6d' % Mlat))
            for key in self.varying_paramdict.keys():
                self.paramdict[key] = self.varying_paramdict[key][Mlat]
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

    def Sigma_hat(self,xi,p):
        sum_1 = 0.0
        sum_2 = 1.0
        for m in range(1,10):
            sum_1 += 2*m**p*np.exp(-0.5*xi*m**2)
            sum_2 += 2*np.exp(-0.5*xi*m**2)
        return sum_1/sum_2
            
    def show_clock(self):
        return time.strftime('%a, %d %b %Y %H:%M:%S +0000',time.gmtime())


class BiasExperiment(Experiment):
    '''
    Measure bias as a function of the lattice spacing
    '''
    def __init__(self,Tfinal,m0,Mlatlist,epsilon):
        super().__init__(Tfinal,m0,Mlatlist)
        self.epsilon=epsilon
        self.outputdir = './runs_bias/'
        self.templatefilename = 'parameters_template_bias.tpl'
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
        plt.plot(df['alat'],self.bias_slope()*df['alat'],
                 linewidth=2,
                 color='red',
                 linestyle='-',
                 label='theory $'+('%6.3f' % self.bias_slope())+'a $')
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

    def bias_slope(self):
        xi = self.Tfinal/self.m0
        S_hat2 = self.Sigma_hat(xi,2)
        S_hat4 = self.Sigma_hat(xi,4)
        return 1./(4.*np.pi**2*self.m0**2)*(0.5-xi*S_hat2+0.25*xi**2*(S_hat4-S_hat2**2))
    
class VarianceExperiment(Experiment):
    ''' 
    Measure variance decay, i.e. variance of difference estimator
    as a function of the lattice spacing.
    '''
    def __init__(self,Tfinal,m0,Mlatlist,nsamples):
        super().__init__(Tfinal,m0,Mlatlist)
        self.nsamples=nsamples
        self.outputdir = './runs_variance/'
        self.templatefilename = 'parameters_template_variance.tpl'
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

class SingleLevelCostExperiment(Experiment):
    ''' 
    Measure runtime of single level method for different lattice spacings
    '''
    def __init__(self,Tfinal,m0,Mlatlist,B0):
        super().__init__(Tfinal,m0,Mlatlist)
        self.B0=B0
        self.outputdir = './runs_singlelevel/'
        self.templatefilename = 'parameters_template_singlelevel.tpl'
        self.varying_paramdict['EPSILON'] = {}
        for Mlat in Mlatlist:
            self.varying_paramdict['EPSILON'][Mlat] = self.epsilon(Mlat)

    def epsilon(self,Mlat):
        alat = self.Tfinal/Mlat
        return self.B0*alat

    def analyse(self):
        raw_data = []
        for Mlat in self.Mlatlist:
            outputfilename = 'output_Mlat_'+str(Mlat)+'.txt'
            epsilon = self.epsilon(Mlat)
            t_elapsed = self.extract_output(self.outputdir+'/'+outputfilename)
            raw_data.append([Mlat,epsilon,t_elapsed])
        df = pandas.DataFrame(raw_data,columns = ['Mlat','epsilon','t_elapsed'])
        df['alat'] = 1./df['Mlat']*self.Tfinal
        print (df)
        plt.clf()
        ax = plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(np.sqrt(0.5)*np.array(df['epsilon'])[-1],
                    np.sqrt(2.0)*np.array(df['epsilon'])[0])
        plt.plot(df['epsilon'],df['t_elapsed'],
                 linewidth=2,
                 color='blue',
                 marker='s',
                 markerfacecolor='blue',
                 markeredgewidth=2,
                 markeredgecolor='blue',
                 label=r'single level')
        epsilon_ref_1 = 1.E-3
        epsilon_ref_2 = 3.E-3
        C_ref = 10.
        plt.plot([epsilon_ref_1,epsilon_ref_2],[C_ref,C_ref*(epsilon_ref_2/epsilon_ref_1)**(-3)],
                 linewidth=2,
                 color='blue',
                 linestyle='--',
                 label=r'$\propto \epsilon^{-3}$')
        for j in range(len(df['alat'])):
            plt.annotate('   '+str(np.array(df['Mlat'])[j]),
                         (np.array(df['epsilon'])[j],
                          np.array(df['t_elapsed'])[j]),
                         va='bottom',ha='left')

        ax.set_xlabel(r'tolerance $\epsilon$')
        ax.set_ylabel('elapsed time [s]')

        plt.legend(loc='upper right')
        plt.savefig('time_singlelevel.pdf',bbox_inches='tight')

    def extract_output(self,filename):
        with open(filename) as f:
            for line in f.readlines():
                m = re.match(' *\[timer SinglevelMC\] *: *(.*) *s',line)
                if m:
                    t_elapsed = float(m.group(1))
        return t_elapsed

class MultiLevelCostExperiment(Experiment):
    ''' 
    Measure runtime of multilevel method for different lattice spacings
    '''
    def __init__(self,Tfinal,m0,Mlatlist,B0):
        super().__init__(Tfinal,m0,Mlatlist)
        self.B0=B0
        self.outputdir = './runs_multilevel/'
        self.templatefilename = 'parameters_template_multilevel.tpl'
        self.varying_paramdict['EPSILON'] = {}
        self.varying_paramdict['NLEVEL'] = {}
        for Mlat in Mlatlist:
            self.varying_paramdict['EPSILON'][Mlat] = self.epsilon(Mlat)
            self.varying_paramdict['NLEVEL'][Mlat] = self.nlevel(Mlat)

    def epsilon(self,Mlat):
        alat = self.Tfinal/Mlat
        return self.B0*alat

    def nlevel(self,Mlat):        
        return (Mlat.bit_length()-1)-6

    def analyse(self):
        raw_data = []
        for Mlat in self.Mlatlist:
            outputfilename = 'output_Mlat_'+str(Mlat)+'.txt'
            epsilon = self.epsilon(Mlat)
            t_elapsed = self.extract_output(self.outputdir+'/'+outputfilename)
            raw_data.append([Mlat,epsilon,t_elapsed])
        df = pandas.DataFrame(raw_data,columns = ['Mlat','epsilon','t_elapsed'])
        df['alat'] = 1./df['Mlat']*self.Tfinal
        print (df)
        plt.clf()
        ax = plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(np.sqrt(0.5)*np.array(df['epsilon'])[-1],
                    np.sqrt(2.0)*np.array(df['epsilon'])[0])
        plt.plot(df['epsilon'],df['t_elapsed'],
                 linewidth=2,
                 color='blue',
                 marker='s',
                 markerfacecolor='blue',
                 markeredgewidth=2,
                 markeredgecolor='blue',
                 label=r'multi level')
        epsilon_ref_1 = 2.E-3
        epsilon_ref_2 = 6.E-3
        C_ref = 10.
        plt.plot([epsilon_ref_1,epsilon_ref_2],[C_ref,C_ref*(epsilon_ref_2/epsilon_ref_1)**(-3)],
                 linewidth=2,
                 color='blue',
                 linestyle='--',
                 label=r'$\propto \epsilon^{-3}$')
        for j in range(len(df['alat'])):
            plt.annotate('   '+str(np.array(df['Mlat'])[j]),
                         (np.array(df['epsilon'])[j],
                          np.array(df['t_elapsed'])[j]),
                         va='bottom',ha='left')

        ax.set_xlabel(r'tolerance $\epsilon$')
        ax.set_ylabel('elapsed time [s]')

        plt.legend(loc='upper right')
        plt.savefig('time_multilevel.pdf',bbox_inches='tight')

    def extract_output(self,filename):
        with open(filename) as f:
            for line in f.readlines():
                m = re.match(' *\[timer MultilevelMC\] *: *(.*) *s',line)
                if m:
                    t_elapsed = float(m.group(1))
        return t_elapsed

    
if (__name__ == '__main__'):

    parser = argparse.ArgumentParser(description='Numerical experiments.')
    parser.add_argument('--bias', action='store_true',
                        default=False,
                        help='generate bias data by running single level experiment')

    parser.add_argument('--variance', action='store_true',
                        default=False,
                        help='generate variance data by running two-level experiment')
    parser.add_argument('--singlelevel', action='store_true',
                        default=False,
                        help='generate single level runtime data')
    parser.add_argument('--multilevel', action='store_true',
                        default=False,
                        help='generate multi level runtime data')

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
        experiment = BiasExperiment(Tfinal,m0,Mlatlist,epsilon)
        experiment.execute()
        experiment.analyse()

    # Run two-level experiment to calculate variance
    if args.variance:
        nsamples = 100000
        print ('nsamples = ',nsamples)
        experiment = VarianceExperiment(Tfinal,m0,Mlatlist,nsamples)
        experiment.execute()
        experiment.analyse()

    # Run single level cost experiment
    if args.singlelevel:
        B0 = 0.5
        Mlatlist = (32,64,128,256,512,1024,2048,4096)
        print ('B0 = ',B0)
        experiment = SingleLevelCostExperiment(Tfinal,m0,Mlatlist,B0)
        experiment.execute()
        experiment.analyse()

    # Run single level cost experiment
    if args.multilevel:
        B0 = 0.5
        Mlatlist = (128,256,512,1024,2048,4096)
        print ('B0 = ',B0)
        experiment = MultiLevelCostExperiment(Tfinal,m0,Mlatlist,B0)
        experiment.execute()
        experiment.analyse()
