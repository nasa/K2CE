#!/usr/bin/env python

# Kepler K2 Cadence Events (k2ce)

# file://k2_cadence_envents.py

__version__ = '0.86'

"""

## Notices

Copyright (C) 2019 United States Government as represented by the 
Administrator of the National Aeronautics and Space Administration. 
All Rights Reserved.

NASA acknowledges the SETI Institute's primary role in authoring and
producing the Kepler-K2 Cadence Events application under Cooperative
Agreement Number NNX13AD01A.


## Disclaimers

No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF 
ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED 
TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY
IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR 
FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE 
ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO 
THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN 
ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, 
RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS 
RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY 
DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF 
PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."

Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST 
THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS 
ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN 
ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, 
INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S 
USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE 
UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY 
PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR 
ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS 
AGREEMENT.
"""

import sys

pyver = (sys.version_info.major*10) + (sys.version_info.minor)
if (pyver < 27):
    print('*** ERROR *** This application needs Python 2.7 or higher.')
    sys.exit(1)
pass#if

import numpy as np
import matplotlib.pyplot as plt
import os.path
import sys
import argparse
import ast
from astropy.io import fits

command_line = False

try:
    import lightkurve as lk    
except:
    print('\n***** ERROR *****\n') 
    print('The Python package lightkurve needs to be installed.\n')
    print('This is the installation command for lightkurve using pip:\n')
    print('pip install lightkurve --upgrade\n')
    print('For further installation details see the lightkurve homepage:\n')
    print('http://lightkurve.keplerscience.org/install.html\n')
    sys.exit(1)

def check_file_exists(filename,overwrite):
    """
    Utility function for k2_cadence_events.
    """
    assert(isinstance(filename,str))
    assert(isinstance(overwrite,bool))
    msg = 'Requested output file already exists (overwrite=False):\n'
    if (not overwrite):
        if (os.path.isfile(filename)):
            print('\n***** ERROR *****\n\n%s' % (msg))
            print("new_filename='%s'\n" % filename)
            sys.exit(1)
pass  #} def

def k2_cadence_events(\
  filename=None,
  bitmask=None,
  from_archive=True,
  target=None, 
  cadence=None,
  campaign=None,
  tag=None,
  plotfile=None,
  scatter=True,
  bars_yy=None,
  xlim=None,
  ylim=None,
  SAP_FLUX=True,
  report=False,
  report_filename=None,
  n_before=0,
  n_after=0,
  bitmask_decode=False,
  bitmask_flags=False,
  show_plot=True,
  new_filename=None,
  overwrite=False,
  useTPF=False,
  xcut=None,
  ycut=None):
    """
    Parameters
    ----------
    filename : str  [default: None]
        Filename of the KeplerLightCurveFile to be analyzed
    bitmask : int  [default: None]
        Bitmask value (integer) specifying quality flag bitmask
        of cadences to *show* events.  See Table 2-3 of the 
        Kepler Archive Manual (KDMC-10008-006) for more information.
    from_archive : str  [default: True]
        If True, get the data from the Barbara A. Mikulski Archive for Space
        Telescopes (MAST) at the Space Telescope Science Institute (STScI).
    target : str or int  [default: None]
        Target name or EPIC number.
    cadence : str  [default: 'long']
        Type of Kepler/K2 cadence: 'short' or 'long'
    campaign : int  [default: None]
        The K2 Campaign number.
    tag : str  [default: None]
        String written at the start of the title of the plot.
    plotfile : str  [default: None]
        Filename of the output plotfile (if any).
    scatter: bool  [default: True]
        If True: the data is plotted as a scatter plot.
        If False: the data is plotted as a line plot.
    bars_yy : float  [default: None]
        Used to set the Y axis location of the gray vertical 
        lines showing events.
    xlim : 2-item tuple  [default: None]
        User-defined right and left limits for the X axis. 
        Example: xlim=(3360,3390)
    ylim : 2-item tuple  [default: None]
        User-defined bottom and top limits for the Y axis. 
        Example: ylim=(0.0,1.1)
    SAP_FLUX : bool  [default: True]
        If True: flux is SAP_FLUX
        If False: flux is PDCSAP_FLUX
    report : bool  [default: False]
        If True, print(out the time, flux, cadence number and
        the QUALITY value for each event.
    report_filename : str  [default: None]
        Filename of the report (if any).
    n_before : int [default: 0]
        Number of observations (cadences) before an event to mark as bad.
    n_after : int [default: 0]
        Number of observations (cadences) after an event to mark as bad.
    bitmask_decode : [default: False]
        Decodes (translate) the bitmask value to K2 Quality Flag Events 
    bitmask_flags : [default: False]
        If True, show the QUALITY bit flags. See Table 2-3 of the Kepler
        Archive Manual (KDMC-10008-006) for more information.
    show_plot : bool [default: True]
        If True, show the plot
    new_filename : str  [default: None]
        Filename of the new long (llc) or short (slc) light curve file (or
        target pixel file) with the event and bad cadences removed.
    overwrite : bool [default: False]
        If True and new_filename is not None, 
        overwrite ("clobber") the new_filename if it exists. 
    useTPF: bool [default: False]
        if False, return a KeplerLightCurveFile.  
        If True, return a KeplerTargetPixelFile.
    xcut : 2-item tuple [default: None]
        Cadences with time (X axis) values within the xcut limits will be 
        flagged for removal.
    ycut : 2-item tuple [default: None]
        Cadences with normalized flux (Y axis) values within the ycut limits 
        will be flagged for removal.

    Returns
    -------
    ax : matplotlib.axes._subplots.AxesSubplot
        The matplotlib axes object created by the function
    n_events : int
        The number of events (cadences) with a QUALITY value 
        featuring at least one bit in the bitmask)
    objf: KeplerLightCurveFile object or KeplerTargetPixelFile object
        The KeplerLightCurveFile object or KeplerTargetPixelFile object 
        created by the function.
    idx : numpy boolean array
        Normally: array of events created by the function.  If the keywords 
        n_before or n_after are greater than zero, then the idx array is a
        combination of events and bad cadences.  
    """
    #
    ftypes = ['Light Curve File','Target Pixel File']
    color = ['dodgerblue','red','slategrey','navy']
    #
    print('**********************************************')
    print('%s %s' % ('Kepler K2 Cadence Events (k2ce): Version',__version__))
    print('**********************************************')
    # if no target information given, use this default target:
    if ((filename is None)and(target is None)and(campaign is None)):
        if (command_line):
            print('\n*******************************************************')
            print('***** Use --help to see the command line options. *****')
            print('*******************************************************\n')
        from_archive = True
        target = '212803289'  # exoplanet K2-99b
        campaign = 17
        cadence = 'short'
        print('\nUsing default target (exoplanet K2-99b):\n')
        print('  from_archive=%s' % (from_archive))
        print('  target=%s' % (target))
        print('  campaign=%d' % (campaign))
        print('  cadence=%s' % (cadence))
        bitmask_decode = True
    #
    if (filename is not None):
        from_archive = False
    isLCF = False
    isTPF = False
    keplerData = False
    tessData = False
    objf = None
    ok = False
    if (from_archive):
        msg0 = '***** ERROR *****\n'
        msg1 = 'If the keyword from_archive=True,\n'
        msg2 = 'The keyword target must be an integer or string (KIC/EPIC ID '\
          +'or object name)'
        msg3 = "The keyword cadence must be a string: 'short' or 'long'"
        msg4 = 'The keyword campaign must be an integer (a valid K2 Campaign '\
          +'number)'
        assert (target is not None),msg0+msg1+msg2
        assert (cadence is not None),msg0+msg1+msg3
        assert (campaign is not None),msg0+msg1+msg4
        if (not useTPF):
            try:
                objf = lk.search_lightcurvefile(target=target,cadence=cadence,\
                  campaign=campaign).download(quality_bitmask=0)
            except Exception as e:
                print("[1] Exception raised: {}".format(e))
                sys.exit(1)
            isLCF = True
        else:
            try:
                objf = lk.search_targetpixelfile(target=target,\
                  cadence=cadence,
                  campaign=campaign).download(quality_bitmask=0)
            except Exception as e:
                print("[2] Exception raised: {}".format(e))
                sys.exit(1)
            isTPF = True
        keplerData = True
        ok = True
    else:
        assert (filename is not None),\
          '***** ERROR ***** A filename must be given'
        assert (os.path.isfile(filename)),\
          '***** ERROR ***** This file does not exist:\n %s' % (filename)
        ok = False
        if (not ok):
            try:
                objf = lk.open(filename,quality_bitmask=0)
            except Exception as e:
                print("[2] Exception raised: {}".format(e))
                sys.exit(1)
            if isinstance(objf,lk.lightcurvefile.KeplerLightCurveFile):
                keplerData = True
                isLCF = True
                ok = True
            elif isinstance(objf,lk.targetpixelfile.KeplerTargetPixelFile):
                keplerData = True
                isTPF = True
                ok = True
            elif isinstance(objf,lk.lightcurvefile.TessLightCurveFile):
                tessData = True
                isLCF = True
                ok = True
            elif isinstance(objf,lk.targetpixelfile.TessTargetPixelFile):
                tessData = True
                isTPF = True
                ok = True
            else:
                str_ = type(objf)
                print('***** ERROR *****:  '\
                  'lk.open returned an unknown type of object!'+str_)
                sys.exit(1)
    pass #}  if (from_archive):
    assert(ok)
    pass #} if (from_archive):
    assert (ok)
    if (isLCF):
        lcf = objf
        ftype = ftypes[0]
    if (isTPF): 
        tpf = objf
        ftype = ftypes[1]
    #
    filename = objf.path
    path, fn = os.path.split(filename)
    if (len(path)==0):
        path = os.getcwd()
    print('\nfilename=%s/%s' % (path,fn))
    #
    try:
        telescop = objf.hdu[0].header['TELESCOP']
    except:
        print('\n***** ERROR *****\n\nMissing keyword: TELESCOP')
        sys.exit(1)
    telescop = objf.hdu[0].header['TELESCOP']
    if (telescop=='Kepler'):
        assert(keplerData)
        assert(not tessData)
    elif (telescop=='TESS'):
        assert(not keplerData)
        assert(tessData)
    #
    extname = objf.hdu[1].header['EXTNAME']  # FITS keyword
    if (isLCF):
        assert (extname == 'LIGHTCURVE')
    else:
        if (keplerData):
            assert (extname == 'TARGETTABLES')
        if (tessData):
            assert (extname) == 'PIXELS'
    #
    obsmode = None
    if (keplerData):
        try:
            obsmode = objf.hdu[0].header['OBSMODE']
        except:
            obsmode = ''
    try:
        mission = objf.hdu[0].header['MISSION']
    except:
        mission = None
    if ((mission is not None)and(obsmode is not None)):
        if (mission=='K2'):
            print('\n%s/%s %s %s\n' % (telescop,mission,obsmode,ftype))
        if (mission=='Kepler'):
            print('\n%s %s %s\n' % (telescop,obsmode,ftype))  
    elif (obsmode is not None):
        print('\n%s %s %s\n' % (telescop,obsmode,ftype))
    else:
        print('\n%s %s\n' % (telescop,ftype))
    # 
    time = objf.hdu[1].data['time'].copy()
    cadenceno = objf.hdu[1].data['CADENCENO'].copy()
    if (SAP_FLUX):
        flux_type = 'SAP_FLUX'
    else:
        flux_type = 'PDCSAP_FLUX'
    if (isLCF):
        flux = lcf.hdu[1].data[flux_type].copy()
        if (keplerData):
            quality = lcf.hdu[1].data['SAP_QUALITY'].copy()
        else:
            quality = lcf.hdu[1].data['QUALITY'].copy()
    else:
        flux = tpf.to_lightcurve().flux.copy()
        quality = tpf.hdu[1].data['QUALITY'].copy()
        flux_type = 'SAP_FLUX'
    assert (time.size == cadenceno.size)
    assert (time.size == quality.size)
    assert (time.size == flux.size)
    flux_median = np.nanmedian(flux)
    flux /= flux_median  # Normalized flux
    #
    try: 
        campaign = objf.hdu[0].header['CAMPAIGN']  # FITS keyword
        print('K2 Campaign %d\n' % campaign)
    except:
        campaign = None
    try: 
        quarter = objf.hdu[0].header['QUARTER']  # FITS keyword
        print('Kepler Quarter %d\n' % quarter)
    except:
        quarter = None
    try: 
        sector = objf.hdu[0].header['SECTOR']  # FITS keyword
        print('TESS Sector %d\n' % sector)
    except:
        sector = None
    bjdrefi = objf.hdu[1].header['BJDREFI']  # FITS keyword
    object = objf.hdu[0].header['OBJECT']  # FITS keyword
    print('target: %s\n' % (object))
    print('%d cadences\n' % (len(time)))
    if (bitmask is not None):
        mask = bitmask
    else:
        if (keplerData):
            bit21 = (1<<20)  # thruster firing
            #assert (bit21 == 1048576)
            bit20 = (1<<19)  # possible thruster firing
            #assert (bit20 == 524288)
            bit16 = (1<<15)  # spacecraft not in fine point
            #assert (bit16 == 32768)
            bit03 = (1<<2)  # spacecraft is in coarse point
            #assert (bit03 == 4)
            mask = bit21 + bit20 + bit16 + bit03
            #assert (mask == 1605636)
        else: 
            # TESS data
            bit01 = (1<<0)  # AttitudeTweak
            #assert (bit01 == 1)
            bit02 = (1<<1)  # SafeMode
            #assert (bit02 == 2)
            bit03 = (1<<2)  # CoarsePoint
            #assert (bit03 == 4)
            bit04 = (1<<3)  # EarthPoint
            #assert (bit04 == 8)
            bit05 = (1<<4)  # Argabrightening
            #assert (bit05 == 16)
            bit06 = (1<<5)  # Desat
            #assert (bit06 == 32)
            bit08 = (1<<7)  # ManualExclude
            assert (bit08 == 128)
            bit10 = (1<<9)  # ImpulsiveOutlier
            #assert (bit10 == 512)
            mask = bit01+bit02+bit03+bit04+bit05+bit06+bit08+bit10
            assert (mask == 703)
        bitmask = mask
        print('Using default bitmask value of %d.' % (bitmask))
    #
    if (bitmask_decode):
        if (keplerData):
            bitmask_str = '{0:021b}'.format(bitmask)        
            print('\nThe bitmask value of %d = %s\n' % (bitmask,bitmask_str))
            print('translates as\n')
            print(lk.KeplerQualityFlags.decode(bitmask))
        if (tessData):
            bitmask_str = '{0:012b}'.format(bitmask)        
            print('\nThe bitmask value of %d = %s\n' % (bitmask,bitmask_str))
            print('translates as\n')
            print(lk.TessQualityFlags.decode(bitmask))
        print('')
    #
    if (bitmask_flags):
        if (keplerData):
            d = lk.KeplerQualityFlags.STRINGS
            # supply missing dictionary item:
            d[512] = 'This bit unused by Kepler'
            print('\nName     Value   Explanation (Kepler/K2)')
        if (tessData):
            d = lk.TessQualityFlags.STRINGS
            print('\nName     Value   Explanation (TESS)')
        if (pyver == 27):
            list_sorted = sorted( ((k,v) for k,v in d.iteritems()))
        if (pyver >= 30):
            list_sorted = sorted( ((k,v) for k,v in d.items()))
        for j,(v,k) in enumerate(list_sorted):
            print('Bit%02d  %7d : %s' % ((j+1),v,k))
        print('')
    #
    xx = time
    yy = flux
    cc = cadenceno
    qq = quality
    #campaign_str = 'C'+str(campaign)
    mtag = ''
    if (mission is not None):
        if (mission == 'K2'):
            mtag = 'C'+str(campaign)
        if (mission == 'Kepler'):
            mtag = 'Q'+str(quarter)
    if (keplerData):
        mask_str = '{0:021b}'.format(mask)
    else:
        mask_str = '{0:012b}'.format(mask)
    title = object+'   '
    if (keplerData):
        title += '['+mtag+']'
    title += '   [bitmask: '+mask_str+']'
    if (tag is not None):
        title = tag + title
    kwargs1  = dict(color=color[0],zorder=0)
    fig, ax = plt.subplots(figsize=(14,5))
    if (not scatter):
        ax.plot(xx,yy,**kwargs1)
    else:
        ax.scatter(xx,yy,s=6,**kwargs1)
    xlabel = 'Time  [BJD - '+str(bjdrefi)+']  [days]'
    ax.set_xlabel(xlabel,size='x-large')
    ax.set_ylabel('Normalized '+flux_type,size='x-large')
    ax.set_title(title,size='x-large');
    ax.grid(alpha=0.5)
    #    
    idx = (qq & mask) > 0
    n_events = np.count_nonzero(idx)
    #
    xxx = xx[idx].copy()
    yyy = yy[idx].copy()
    ccc = cc[idx].copy()
    qqq = qq[idx].copy()
    #
    if (report_filename is not None):
        assert(isinstance(report_filename,str))
        report = True
    if (report):
        if (report_filename is not None):
            check_file_exists(report_filename,overwrite)
            f = open(report_filename,'w')
            f.write('# REPORT\n')
            f.write('#\n')
            f.write('# %s\n' % (filename))
            f.write('#\n')
            f.write('# %d events with bitmask value of %d (= %s)\n' %\
              (n_events,mask,mask_str))
            if (n_events > 0):
                f.write('#\n')
                f.write('#  Event         Time  Normalized_Flux  CADENCENO  '\
                  +'===========QUALITY============\n')
                for j, (xxx_, yyy_, ccc_, qqq_) in \
                  enumerate(zip(xxx,yyy,ccc,qqq)):
                    if (keplerData):
                        qqq_bitmask_str = '{0:021b}'.format(qqq_)
                    else:
                        qqq_bitmask_str = '{0:012b}'.format(qqq_)
                    jj = j + 1
                    f.write('%8d %12.6f %16.11f %9d %8d = %s\n' % \
                      (jj, xxx_, yyy_, ccc_, qqq_, qqq_bitmask_str))
            f.close()
            print('')
            print('%s <--- report written  :-)\n' % (report_filename))
            print('')
        else:
            print('# REPORT')
            print('#')
            print('# %s' % (filename))
            print('#')
            print('# %d events with bitmask value of %d (= %s)' %\
              (n_events,mask,mask_str))
            if (n_events > 0):
                print('#')
                print('#  Event         Time  Normalized_Flux  CADENCENO  '\
                  +'=============QUALITY============')
                for j, (xxx_, yyy_, ccc_, qqq_) in \
                  enumerate(zip(xxx,yyy,ccc,qqq)):
                    if (keplerData):
                        qqq_bitmask_str = '{0:021b}'.format(qqq_)
                    else:
                        qqq_bitmask_str = '{0:012b}'.format(qqq_)
                    jj = j + 1
                    print('%8d %12.6f %16.11f %9d %8d = %s' % \
                      (jj, xxx_, yyy_, ccc_, qqq_, qqq_bitmask_str))        
    #
    assert(n_before >= 0)
    assert(n_after >= 0)
    clipit = (n_before > 0) or (n_after > 0)
    if (clipit):
        # mark bad cadences in the index array
        jdx = idx.copy()
        jje = jdx.size
        jj_min = 0
        jj_max = jje - 1
        for jj in range(jje):
            if (idx[jj]):
                ccm = cc[jj]
                ccb = ccm - n_before
                cce = ccm + n_after
                jb = np.clip(jj - n_before,jj_min,jj_max)
                assert( jb >= jj_min )
                je = np.clip(jj + n_after,jj_min,jj_max)
                assert( je <= jj_max )
                jn = je - jb + 1
                for j in range(jn):
                    k = jb + j
                    cck = cc[k]
                    if ((cck >= ccb)and(cck <= cce)):
                        jdx[k] = True
        idx = jdx.copy() 
    #
    clipx = False
    if (xcut is not None):
        kdx = idx.copy()
        xcut_ = sorted(xcut)
        xmin = xcut_[0]
        xmax = xcut_[1]
        je = kdx.size
        for j in range(je):
            x_ = xx[j]
            if (np.isfinite(x_)):
                if ((x_ >= xmin)and(x_ <= xmax)):
                    kdx[j] = True
        idx = kdx.copy()
        clipx = True
    #
    clipy = False
    if (ycut is not None):
        kdx = idx.copy()
        ycut_ = sorted(ycut)
        ymin = ycut_[0]
        ymax = ycut_[1]
        je = kdx.size
        for j in range(je):
            y_ = yy[j]
            if (np.isfinite(y_)):
                if ((y_ >= ymin)and(y_ <= ymax)):
                    kdx[j] = True
        idx = kdx.copy() 
        clipy = True
    #
    if (clipit or clipx or clipy):                       
        xxx = xx[idx].copy()
        yyy = yy[idx].copy()
        ccc = cc[idx].copy()
        qqq = qq[idx].copy()
    #
    kwargs2 = dict(s=10,color=color[1],zorder=10)
    ax.scatter(xxx,yyy,**kwargs2)
    #
    if (xlim is not None):
        ax.set_xlim(xlim)
    if (ylim is not None):
        ax.set_ylim(ylim)
    #
    # mark events near the top of the plot with grey vertical bars
    yyy = yy[idx].copy()
    if (bars_yy is None):
        bottom, top = ax.get_ylim()
        bars_yy = bottom + (0.94*(top - bottom))
    yyy[:] = bars_yy
    kwargs3  = dict(marker='|',alpha=0.5,color=color[2],zorder=0,s=(20**2))
    ax.scatter(xxx,yyy,**kwargs3)
    #
    # show path and filename on the right side of plot
    path, fn = os.path.split(filename)
    if (len(path)==0):
        path = os.getcwd()
    plt.figtext( 0.95, 0.05, path+'/', ha='right', va='bottom', \
      color=color[3], size='small', rotation=90)
    plt.figtext( 0.96, 0.05, fn, ha='right', va='bottom', \
      color=color[3], size='small', rotation=90)
    #
    if (plotfile is not None):
        check_file_exists(plotfile,overwrite)
        plt.savefig(plotfile,dpi=300)
        print('%s <--- plotfile written  :-)\n' % (plotfile))
        #plt.show()
        plt.close()
    if (show_plot):
        plt.ioff()
        plt.show()  
    if (new_filename is not None):
        check_file_exists(new_filename,overwrite)
        hdul = fits.open(filename)
        hdul[1].data = hdul[1].data[~idx]
        hdul.writeto(new_filename,overwrite=overwrite)
        hdul.close()       
        print('')
        print('%s  <--- new FITS file written  :-)\n' % (new_filename))
        print('')
    sys.stdout.flush()
    return (ax, n_events, objf, idx)
pass  #} def

def str2bool(v):
    """
    Utility function for argparse.
    """
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
    pass  #} def

def accept_str_or_int(v):
    """
    Utility function for argparse.
    """
    if isinstance(v, int):
        return str(v)
    elif isinstance( v, str):
        return v
    else:
        raise argparse.ArgumentTypeError('str or int value expected.')
    pass  #} def

if __name__ == "__main__":
    command_line = True
    verbose = False
    #verbose = True
    #
    parser = argparse.ArgumentParser()
    #
    parser.add_argument('--filename', action="store", type=str, default=None,
      help='Filename of the K2 light curve file (ktwo*llc.fits) to be '\
        +'analyzed  [default: None]')
    parser.add_argument('--bitmask', action="store", type=int, default=None,
      help='Bitmask value (integer) specifying quality flag bitmask of '\
        +'cadences to *show* events.  See Table 2-3 of the Kepler Archive '\
        +'Manual  (KDMC-10009-006) for more information [default: None]')
    parser.add_argument('--from_archive', 
      type=str2bool, default=True,
      help='If True, get the data from the Mikulski Archive for Space '\
        +'Telescopes (MAST) at the Space Telescope Science Institute (STScI)'\
        +'  [default=True]')
    parser.add_argument('--target', action="store", 
      type=accept_str_or_int, default=None,
      help='Target name or EPIC number  [default: None]')
    parser.add_argument('--cadence', action="store", type=str, default='long',
      help="Type of K2 cadence: 'short' or 'long' [default: 'long']")
    parser.add_argument('--campaign', action="store", type=int, default=None, 
      help='K2 campaign number [default: None]')
    parser.add_argument('--tag', action="store", type=str, default=None,
      help='String written at the start of the title of the plot '\
        +'[default: None]')
    parser.add_argument('--plotfile', action="store", type=str, default=None,
      help='Filename of the output plotfile (if any)  [default: None]')
    parser.add_argument('--scatter', 
      type=str2bool, default=True,
      help='If True: the data is plotted as a scatter plot.\n'\
        +'If False: the data is plotted as a line plot\n[default=True]')
    parser.add_argument('--bars_yy', action="store", type=float, default=None,
      help='Used to set the Y axis location of the gray vertical lines '\
        +'showing events [default: None]')
    parser.add_argument('--xlim', action="store", type=ast.literal_eval,\
      default=None,help='User-defined right and left limits for the X axis.\n'\
        +"Example: xlim='(3360,3390)' [default: None]")
    parser.add_argument('--ylim', action="store", type=ast.literal_eval,\
      default=None,help='User-defined bottom and top limits for the Y axis.\n'\
        +"Example: ylim='(0.0,1.1)' [default: None]")
    parser.add_argument('--SAP_FLUX',type=str2bool,default=True, 
      help='If True, flux is SAP_FLUX. If False, flux is PDCSAP_FLUX '\
      +'[default: True]')
    parser.add_argument('--report',type=str2bool,default=False,
      help='If True, print(out the time, flux, cadence number and the '\
        +'QUALITY value fore each event [default: False]')
    parser.add_argument('--report_filename', action="store", type=str, 
      default=None,help='Filename of the report (if any) [default: None]')
    parser.add_argument('--n_before', action="store", type=int, default=0, 
      help='Number of observations (cadences) before an event to mark as '\
        +'bad [default: 0]')
    parser.add_argument('--n_after', action="store", type=int, default=0, 
      help='Number of observations (cadences) after an event to mark as '\
        +'bad [default: 0]')
    parser.add_argument('--bitmask_decode',type=str2bool,default=False,
      help='If True, Decodes (translate) the bitmask value to K2 Quality '\
        +'Flag Events [default: False]')
    parser.add_argument('--bitmask_flags',type=str2bool, default=False,
      help='If True, show the QUALITY bit flags. See Table 2-3 of the Kepler '\
        +'Archive Manual (KDMC-10008-006) for more information. '\
        +'[default=False]')
    parser.add_argument('--show_plot', type=str2bool, default=True,
      help='If True, show the plot [default=True]')
    parser.add_argument('--new_filename',action="store",type=str,default=None,
      help='Filename of the new long (llc) or short (slc) file with the '\
        +'event cadences and the bad cadences removed [default: None].')
    parser.add_argument('--overwrite', 
      type=str2bool, default=False,
      help='If True, overwrite ("clobber") an existing output file '\
        +'[default: False].')
    parser.add_argument('--useTPF', type=str2bool, default=False,
      help='If False, return a KeplerLightCurveFile object.  If True, '\
        +'return a KeplerTargetPixelFile object [default: False].')
    parser.add_argument('--xcut', action="store", type=ast.literal_eval,\
      default=None,help='Cadences with time (X axis) values within the '\
      +'xcut limits will be flagged for removal.\n'\
      +"Example: xcut='(3510,3590)' [default: None]")
    parser.add_argument('--ycut', action="store", type=ast.literal_eval,\
      default=None,help='Cadences with normalized flux (Y axis) values '\
      +'within the ycut limits will be flagged for removal.\n'\
      +"Example: ycut='(0.0,0.9)' [default: None]")
    #
    args = parser.parse_args()
    #
    #verbose = True  # DEBUG
    if (verbose):
        print('%s =args.filename' % (args.filename))
        print('%s =args.bitmask' % (args.bitmask))
        print('%s =args.from_archive' % args.from_archive)
        print('%s =args_target' % (args.target))
        print('%s =args_cadence' % (args.target))
        print('%s =args.campaign' % (args.campaign))
        print('%s =args.tag' % (args.tag))
        print('%s =args.plotfile' % (args.plotfile))
        print('%s =args.scatter' % (args.scatter))
        print('%s =args.bars_yy' % (args.bars_yy))
        print('%s =args.xlim' % (args.xlim))
        print('%s =args.ylim' % (args.ylim))
        print('%s =args.SAP_FLUX' % (args.SAP_FLUX))
        print('%s =args.report' % (args.report))
        print('%s =args.report_filename' % (args.report_filename))
        
        #bingo start here 2019SEP15
                
        print('%s =args.n_before' % (args.n_before))
        print('%s =args.n_after' % (args.n_after))
        print('%s =args.bitmask_decode' % (args.bitmask_decode))
        print('%s =args.bitmask_flags' % (args.bitmask_flags))
        print('%s=args.show_plot' % (args.show_plot))
        print('%s =args.new_filename' % (args.new_filename))
        print('%s =args.overwrite' % (args.overwrite))
        print('%s =args.useTPF' % (args.useTPF))
        print('%s =args.xcut' % (args.xcut))
        print('%s =args.ycut' % (args.ycut))
        print('%s =__version__' % (__version__))
    #
    if ((args.plotfile is not None) and (args.show_plot)):
        args.show_plot = False
        print('\n*** NOTE BENE *** : The plot will not be shown '\
          +'because the plot was written to the plotfile.\n')
    #
    _, n_events, _, idx_ = k2_cadence_events(filename=args.filename,
      bitmask=args.bitmask,
      from_archive=args.from_archive,
      target=args.target,
      cadence=args.cadence,
      campaign=args.campaign,
      tag=args.tag,
      plotfile=args.plotfile,
      scatter=args.scatter,
      bars_yy=args.bars_yy,
      xlim=args.xlim,
      ylim=args.ylim,
      SAP_FLUX=args.SAP_FLUX,
      report=args.report,
      report_filename=args.report_filename,
      n_before=args.n_before,
      n_after=args.n_after,
      bitmask_decode=args.bitmask_decode,
      bitmask_flags=args.bitmask_flags,
      show_plot=args.show_plot,
      new_filename=args.new_filename,
      overwrite=args.overwrite,
      useTPF=args.useTPF,
      xcut=args.xcut,
      ycut=args.ycut)
    print('\n%d events\n' % n_events)
    flagged = np.count_nonzero(idx_)
    print('%d cadences flagged' % (flagged))
    #print('\n',lk.__version__,'=lk.__version__'
pass  #} if __name__ == "__main__":
      
#EOF
