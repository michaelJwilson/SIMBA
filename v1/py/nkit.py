##  salloc -N 1 -C haswell -q interactive -t 01:00:00
##
##  source /global/common/software/m3035/conda-activate.sh 3.7
##

import matplotlib;  matplotlib.use('PDF')

import matplotlib.pyplot            as      plt 
import pylab                        as      pl
import numpy                        as      np
import time

from   utils                        import  latexify


latexify(columns=1, equal=True, fontsize=10, ggplot=False, usetex=True)

compute          = False
hubble           = 0.68

t0               = time.time()

cols             = ['x', 'y', 'z']

for fpath, label in zip(['/global/homes/m/mjwilson/SIMBA/simba_pos_3.00307.fits', '/global/homes/m/mjwilson/SIMBA/simba_dmpos_3.00307.fits'], ['LBG', 'DM']):
  if compute:
    from   nbodykit.lab                 import  *
    from   nbodykit.io.fits             import  FITSFile
    from   nbodykit.io.csv              import  CSVFile
    from   nbodykit.source.catalog.file import  CSVCatalog, FITSCatalog
    from   nbodykit.transform           import  StackColumns
    from   nbodykit                     import  setup_logging


    setup_logging()
    
    ##  Comoving Mpc/h.         
    ##  cat          = CSVCatalog(fpath, names=cols, delim_whitespace=True)
    cat              = FITSCatalog(fpath)
    
    cat['Position']  = StackColumns(cat['x'], cat['y'], cat['z'])

    ##  Box size in Mpc/h.
    mesh             = cat.to_mesh(resampler='CIC', Nmesh=16, compensated=False, BoxSize=1.e2 * hubble)

    r                = FFTPower(mesh, mode='2d', dk=0.05, kmin=0.1, Nmu=120, los=[0,0,1], poles=[0,2,4])

    poles            = r.poles

    t1               = time.time()

    print(t1 - t0)

    print('Number of galaxies: {}'.format(len(cat['x'])))
    print('Shotnoise: {}'.format(poles.attrs['shotnoise']))

    for ell in [0]:    
      if ell == 0:
        P = P - poles.attrs['shotnoise']

      k     = poles['k']
      P     = poles['power_%d' % ell].real
  
      ##   pl.axhline(poles.attrs['shotnoise'], c='k', linestyle='--')

      np.savetxt('pl_{}_3.00307.txt'.format(label), np.c_[k, P], fmt='%.6le')

  else:
    k, P = np.loadtxt('pl_{}_3.00307.txt'.format(label), unpack=True)

    plt.semilogy(k, P, label=label, marker='^')

    
# format the axes
plt.legend(loc=0, frameon=False)

plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$P_0$ [$h^{-3} \ \mathrm{Mpc}^3$]")

plt.xlim(0.1,   0.6)
plt.ylim(1.e1, 1.e4)

plt.tight_layout()

pl.savefig('pk.pdf')
