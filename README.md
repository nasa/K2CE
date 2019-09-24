# K2CE

## The Kepler-K2 Cadence Events Application

Since early 2018, the Kepler/K2 project has been performing a uniform global reprocessing of data from K2 Campaigns 0 through 14. Subsequent K2 campaigns (C15-C19) are being processed using the same processing pipeline. One of the major benefits of the reprocessing effort is that, for the first time, short-cadence (1-min) light curves are produced in addition to the standard long-cadence (30-min) light curves. Users have been cautioned that the Kepler pipeline detrending module (PDC), developed for use on original Kepler data, has not been tailored for use on short-cadence K2 observations. Systematics due to events on fast timescales,  such as thruster firings, are sometimes poorly corrected for many short-cadence targets. A Python data visualization and manipulation tool, called Kepler-K2 Cadence Events, has been developed that identifies and removes cadences associated with problematic thruster events, thus producing better light curves. Kepler-K2 Cadence Events can be used to visualize and manipulate light curve files and target pixel files from the Kepler, K2, and TESS missions.  We anticipate this software will be available from <http://code.nasa.gov> 

## Prequisites

The Kepler-K2 Cadence Events application should run on any computer with Python 3 .  The application was developed using the open-source Anaconda Distribution of Python 2.7 which is available at <https://www.anaconda.com/distribution>.  K2CE has now been ported to Python 3.7.3.

The application uses Lightkurve (<http://docs.lightkurve.org/>) which is a community-developed, open-source Python package which offers a "beautiful and user-freindly way to analyze astronomical flux time seris data".  It was specifically designed to analyze the pixels and lightcurves obtained by NASA's Kepler, K2, and TESS exoplanet missions. Installation instructions are given at the following webpage: <http://docs.lightkurve.org/about/install.html>.

## Running the Code

### As a Stand-Alone Application from the Command Line

(1) Copy the Python application (k2\_cadence\_events.py) into the local directory.

(2) Make sure that the code is executable by typing (using MacOS X or Linux)

	chmod u+x k2_cadence_events.py

(3) The application can be run without any arguments.  

If you type the following command from the command line, you will see a demonstration of the application:

	./k2_cadence_events.py
	
You should see something like this:

<pre>
./k2_cadence_events.py 
**********************************************
Kepler K2 Cadence Events (k2ce): Version 0.79
**********************************************

*******************************************************
***** Use --help to see the command line options. *****
*******************************************************


Using default target (exoplanet K2-99b):

  from_archive=True
  target=212803289
  campaign=17
  cadence=short

filename=/Users/kmighell/.lightkurve-cache/mastDownload/K2/ktwo212803289-c17_sc/ktwo212803289-c17_slc.fits

Kepler/K2 short cadence Light Curve File

K2 Campaign 17

target: EPIC 212803289

98490 cadences

Using default bitmask value of 1605636.

The bitmask value of 1605636 = 110001000000000000100

translates as

['Coarse point', 'Possible thruster firing', 'Thruster firing', 'No fine point']

</pre>

A plot should appear on your screen like this:

![](README_fig1.png)

Once you are ready to leave the application, close the plot window (MacOS X: click on the red button on the upper-right corner of the plot window).

The application ends by writing a little more information:

<pre>
1635 events

1635 cadences flagged
</pre>

The demo plot shows the short cadence observations of EPIC 212803289
(exoplanet: K2-99b) that were obtained during K2 Campaign 17.

The short-cadence K2 light curve file was downloaded automatically from
the Barbara A. Mikulski Archive for Space Telescopes (MAST) at the Space
Telescope Science Institute (STScI).

The **red points** are the observations that had a SAP_QUALITY value with either Bit21
(1048576 = (2\*\*20)) and/or Bit20 (524288 = (2\*\*19)) and/or Bit16 (32768 =
(2\*\*15)) and/or (4 = (2\*\*2)) set -- which indicates (1) an actual or
(2) probable thruster firing event, or the ***Kepler*** spacecraft was in (3) coarse point or (4) not in fine point **during an individual cadence observation**.

The **blue points** in the scatter plot are the "**good observations**" that ***did not*** have at least one of these bits set.

Of the total of 98490 short-cadence observations, 1635 (1.66%) of the cadences were flagged as being a "cadence event".

The grey vertical bars near the top of the plot show events.
The darker the grey, the greater the number of events in time.

The name of the Kepler/K2 short cadence light curve file analysed is shown
on the right side of the plot. Note that this is where the file is stored
locally.

#### Command-line arguments 

To see the command-line arguments, type

	k2_cadence_events.py --help
	
or 
	
	k2_cadence_events.py -h

You should see something like this:
	
<pre>
./k2_cadence_events.py 
  [-h] [--help]  
  [--filename FILENAME]    
  [--bitmask BITMASK]  
  [--from_archive FROM_ARCHIVE]                            
  [--target TARGET] 
  [--cadence CADENCE]
  [--campaign CAMPAIGN] 
  [--tag TAG]
  [--plotfile PLOTFILE] 
  [--scatter SCATTER]
  [--bars_yy BARS_YY] 
  [--xlim XLIM]
  [--ylim YLIM] 
  [--SAP_FLUX SAP_FLUX]
  [--report REPORT]
  [--report_filename REPORT_FILENAME]
  [--n_before N_BEFORE] 
  [--n_after N_AFTER]
  [--bitmask_decode BITMASK_DECODE]
  [--bitmask_flags BITMASK_FLAGS]
  [--show_plot SHOW_PLOT]
  [--new_filename NEW_FILENAME]
  [--overwrite OVERWRITE] [--useTPF USETPF]
  [--xcut XCUT] 
  [--ycut YCUT]

optional arguments:

  -h, --help            show this help message and exit

  --filename FILENAME   Filename of the K2 light curve file (ktwo*llc.fits) to
                        be analyzed [default: None]

  --bitmask BITMASK     Bitmask value (integer) specifying quality flag
                        bitmask of cadences to *show* events. See Table 2-3 of
                        the Kepler Archive Manual (KDMC-10009-006) for more
                        information [default: None]

  --from_archive FROM_ARCHIVE
                        If True, get the data from the Mikulski Archive for
                        Space Telescopes (MAST) at the Space Telescope Science
                        Institute (STScI) [default=True]

  --target TARGET       Target name or EPIC number [default: None]

  --cadence CADENCE     Type of K2 cadence: 'short' or 'long' [default:
                        'long']

  --campaign CAMPAIGN   K2 campaign number [default: None]

  --tag TAG             String written at the start of the title of the plot
                        [default: None]

  --plotfile PLOTFILE   Filename of the output plotfile (if any) [default:
                        None]

  --scatter SCATTER     If True: the data is plotted as a scatter plot. If
                        False: the data is plotted as a line plot
                        [default=True]

  --bars_yy BARS_YY     Used to set the Y axis location of the gray vertical
                        lines showing events [default: None]

  --xlim XLIM           User-defined right and left limits for the X axis.
                        Example: xlim='(3360,3390)' [default: None]

  --ylim YLIM           User-defined bottom and top limits for the Y axis.
                        Example: ylim='(0.0,1.1)' [default: None]
  --SAP_FLUX SAP_FLUX   If True, flux is SAP_FLUX. If False, flux is
                        PDCSAP_FLUX [default: True]

  --report REPORT       If True, print out the time, flux, cadence number and
                        the QUALITY value fore each event [default: False]

  --report_filename REPORT_FILENAME
                        Filename of the report (if any) [default: None]

  --n_before N_BEFORE   Number of observations (cadences) before an event to
                        mark as bad [default: 0]

  --n_after N_AFTER     Number of observations (cadences) after an event to
                        mark as bad [default: 0]

  --bitmask_decode BITMASK_DECODE
                        If True, Decodes (translate) the bitmask value to K2
                        Quality Flag Events [default: False]

  --bitmask_flags BITMASK_FLAGS
                        If True, show the QUALITY bit flags. See Table 2-3 of
                        the Kepler Archive Manual (KDMC-10008-006) for more
                        information. [default=False]

  --show_plot SHOW_PLOT
                        If True, show the plot [default=True]

  --new_filename NEW_FILENAME
                        Filename of the new long (llc) or short (slc) file
                        with the event cadences and the bad cadences removed
                        [default: None].

  --overwrite OVERWRITE
                        If True, overwrite ("clobber") an existing output file
                        [default: False].

  --useTPF USETPF       If False, return a KeplerLightCurveFile object. If
                        True, return a KeplerTargetPixelFile object [default:
                        False].

  --xcut XCUT           Cadences with time (X axis) values within the xcut
                        limits will be flagged for removal. Example:
                        xcut='(3510,3590)' [default: None]

  --ycut YCUT           Cadences with normalized flux (Y axis) values within
                        the ycut limits will be flagged for removal. Example:
                        ycut='(0.0,0.9)' [default: None]
</pre>

The above was reformatted for easier reading.

### As a part of a Python script

(1) Using your favorite editor, cut-and-paste the following text to create a
Python script called **spud.py**:

<pre>
#!/usr/bin/env python
import sys
import os
import k2_cadence_events as k2ce
nargs = len(sys.argv) - 1
if (nargs == 0):
    # interactive mode: no arguments
    k2ce.k2_cadence_events()
else: 
    # non-interactive mode: first argument is the name of the plotfile (e.g., foo.png)
    plotfile=sys.argv[1]  
    if (os.path.isfile(plotfile)):  # abort if file already exists
        print '***** ERROR ***** Requested file already exists: ',plotfile
    else:
        k2ce.k2_cadence_events(plotfile=plotfile)
        print
        print 'Show the docstring for k2_cadence_events() :'
        print k2ce.k2_cadence_events.__doc__
#EOF
</pre>

(2) Make sure that the script is executable by typing 

	chmod u+x spud.py
	
(3) You can execute the script in its interactive mode (using no arguments):

	./spud.py
	
(4) You can execute the script in its non-interactive mode using the first argument as the name of the output plotfile (e.g., spud_plot.png):

	./spud.py spud_plot.png
	
After the plot is written the Python docstring for k2\_cadence\_envents() is shown.

## Explanation of QUALITY and SAP_QUALITY bit values

<pre>
Name     Value : Explanation (Kepler/K2)
Bit01        1 : Attitude tweak
Bit02        2 : Safe mode
Bit03        4 : Coarse point
Bit04        8 : Earth point
Bit05       16 : Zero crossing
Bit06       32 : Desaturation event
Bit07       64 : Argabrightening
Bit08      128 : Cosmic ray in optimal aperture
Bit09      256 : Manual exlude
Bit10      512 : This bit not used by ***Kepler***
Bit11     1024 : Sudden sensitivity dropout
Bit12     2048 : Impulsive outlier
Bit13     4096 : Argabrightening on CCD
Bit14     8192 : Cosmic ray in collateral data
Bit15    16384 : Detectpr amp,a;u
Bit16    32768 : No fine point
Bit17    65536 : No data
Bit18   131072 : Rolling band in optimal aperture
Bit19   262144 : Rolling band in full mask
Bit20   524288 : Possible thruster Firing
Bit21  1048576 : Thruster firing
</pre>

## Contact

Kepler-K2 Cadence Events was created by Kenneth J. Mighell and supported by the Kepler/K2 Science Office.  You can contact the author at kenneth dot j dot mighell at nasa dot gov or kmighell at seti dot org.


## More information about using the application

You can learn more about the many options of the k2\_cadence\_events application by using running its demo Jupyter notebook called **k2\_cadence\_events.ipynb**.

	jupyter notebook k2_cadence_events.ipynb
	
## Notices

Copyright © 2019 United States Government as represented by the Administrator of the National Aeronautics and Space Administration. All Rights Reserved.

NASA acknowledges the SETI Institute’s primary role in authoring and producing the Kepler-K2 Cadence Events application under Cooperative Agreement Number NNX13AD01A

## Disclaimers

No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."

Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.

##### 2019SEP10
