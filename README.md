# MODEL8478881246: mw5DF4BED6_54E0_401A_9449_98B93AD602C6

## Installation

Download this repository, and install with distutils

`python setup.py install`

Or, install using pip

`pip install git+https://github.com/biomodels/MODEL8478881246.git`

To install a specific version (in this example, from the 2014-09-16 BioModels release)

`pip install git+https://github.com/biomodels/MODEL8478881246.git@20140916`


# Model Notes


This SBML file is a translation of a MatLab model utilized in the following
paper. It describes the interplay of IkB and P100 negative feedback to control
the cellular localization of the NFkB transcription factor in response to
activation of IkB Kinase 1 (IKK1) and IkB Kinase 2 (IKK2).  

Soumen Basak, Ha Na Kim, Jeffrey D. Kearns, Vinay Tergaonkar, Ellen O'Dea,
Shannon L. Werner, Chris A. Benedict, Carl F. Ware, Gourishankar Ghosh, Inder
M. Verma, Alexander Hoffmann. "A Fourth IkB Protein in the NFkB Signaling
Module." Cell (2007)  

Questions concerning the paper should be addressed to the corresponding
author. Alexander Hoffmann (ahoffmann@ucsd.edu)

Technical questions concerning the implementation of the model within SBML
should be addressed to Jeff Kearns (jkearns@ucsd.edu)  

* * *

- Introduction:

The original model was written and simulated within MathWorks MatLab 2006a
using the ode15s (stiff/NDF) solver. It is highly recommended that those
wanting to model this system use the MatLab version which we will freely
provide upon request. As always, simulation results vary according to the
numerical solver used. In addition, there are also shortcomings in the
Stimulation Phase of the model (listed below) that maybe resolved with future
versions of SBML.

* * *

- SBML Generation :  

Translation to SBML Level 2.1 was performed via reconstruction of the model
within MathWorks SimBiology Desktop (version 2.1)   followed by an Export to
SBML.

* * *

- Running the model :  

Simulations require two sequential runs of the model. It is recommended that
you generate a script to pass the final component concentrations from the
first run into the initial concentrations of the second.

  

1\. Equilibrium Phase.

\- In both phases, total NFkB (0.125 micromole), IKK1 (0.1 micromole), and
IKK2 (0.1 micromole) protein remains constant. These exist in free and bound
forms.

\- Synthesis and degradation reactions are included for IkBa, IkBb, IkBe, and
P100 mRNA and protein.

\- Equilibrium state nuclear NFkB (NFkBn) levels for combinations of IkB and
P100 protein knockouts are included in the Cell paper. Use these to validate
whether your simulation time is long enough and whether your solver is
behaving properly.

  

2\. Stimulation Phase.

\- The simulations of this phase have known shortcomings due to insufficient
support by the SBML framework. These include not being able to utilize an
interpolation function and lack of delay functions (resolved in SBML 2.2, but
not yet in the MathWorks SimBiology tool). Resolution of these should enable
recapitulation of the Stimulation Phase results shown in the Cell paper.

\- The model input signal is defined as activity profiles for IKK1 and IKK2
(measured via biochemical assays). Profiles are included below for Lympotoxin-
Beta (LTB) activation of IKK1 and for Tumor Necrosis Factor (TNF) and
Lipopolysaccharide (LPS) activation of IKK2.

\- Combinations of IKK1 and IKK2 inputs are possible. In the paper, a pulse
stimulation with TNF (IKK2) followed later by LTB (IKK1) stimulation was used.

\- The activity profile encodes what fraction of the IKK is able to mediate
IkB and P100 protein degradation by acting as a rate constant multiplier.  

Specifically....  

ikk1_multiplier = IKK1(t) =   IKK1 active at time 't', where 1= 100%

         pd_c_2pq     =   pd_c_2pq   * ikk1_multiplier;

         pd_c_3pqn   =   pd_c_3pqn   * ikk1_multiplier;

  

ikk2_multiplier = IKK2(t) =   IKK2 active at time 't', where 1 = 100%

         pd_c_2ai     =   pd_c_2ai   * ikk2_multiplier;

         pd_c_2bi     =   pd_c_2bi   * ikk2_multiplier;

         pd_c_2ei     =   pd_c_2ei   * ikk2_multiplier;

         pd_c_3ain   =   pd_c_3ain   * ikk2_multiplier;

         pd_c_3bin   =   pd_c_3bin   * ikk2_multiplier;

         pd_c_3ein   =   pd_c_3ein   * ikk2_multiplier;

  

In the equilibrium phase, the ikk1_multiplier is 0.01 and ikk2_multiplier is
0.001.

  

\- IKK1 activation also regulates the translation rate of P100 and the nuclear
IkB/NFkB association rates by unknown mechanisms. We include this in the
Stimulation Phase via a fold IKK1 activation multiplier.

ikk1_multipler_tsl_asn = IKK1(t) / IKK1(t=0)

ps_c_p = ps_c_p * ikk1_multiplier_tsl_asn

a_n_an = a_n_an * ikk1_multiplier_tsl_asn

a_n_bn = a_n_bn * ikk1_multiplier_tsl_asn

a_n_en = a_n_en * ikk1_multiplier_tsl_asn  
  

\- The IKK profile interpolation does NOT work in the SBML version.   Any
ideas on how to make it work are appreciated (contact Jeff Kearns).

\- In brief, the MatLab version takes in two arrays-- IKK profile versus time
-- to create a time-dependent IKK input function via piecewise cubic
interpolation. As this creates N-1 polynomials for arrays with N datapoints,
it is difficult to incorporate this into SBML via Rules or Events.   It might
be possible to construct the interpolation function using the tools present in
MathML.   The MatLab code is included below to help anyone attempt this.

  

%==========================================================

% GENERATE_IKK_INPUTS

%==========================================================  

function ikk1_curve = interpolate_ikk1 (stimulus)  

global SIM_TIME;  

     if (strcmp(stimulus, 'LTB'))

         values   = [1 1.2 1.4 1.9 2.4 2.6 2.6 2.6] / 100;

         time     = [0 180 300 480 600 900 1200 2400]; %minutes

     elseif (strcmp(stimulus, 'LTB_20hr'))

          values   = [1 1 1.2 1.4 1.9 2.4 2.6 2.6] / 100;

          time     = [0 1200 1380 1500 1680 1800 2100 2400];

     elseif (strcmp(stimulus, 'basal_1'))

          values   = [1 1 1 1 1] / 100;

         time     = [0 360 600 1200 2400];

     else

          error('Incorrect IKK1 input');

      end  

     ikk1_curve   = interp1(time, values, 0:SIM_TIME);

\-------------

function ikk2_curve = interpolate_ikk2(stimulus)  

global SIM_TIME;  

     if (strcmp(stimulus, 'TNF_p15'))

          values   = [0.1 60 100 65 50 36 21 16 10 .1 .1 .1 .1] / 100;

          time     = [0 5 10 15 20 25 30 45 59 60 360 1200 2400];

     elseif (strcmp(stimulus, 'LPS_p45'))

          values   = [0.1 3 8 24 25 15 8 5 0.1 0.1] / 100;

          time     = [0 15 30 45 60 120 240 360 600 2400];

     elseif (strcmp(stimulus, 'basal_1'))

          values   = [.1 .1 .1 .1 .1] / 100;

          time     = [0 360 600 1200 2400];  

     else

          error('Incorrect IKK2 input name');

     end

     ikk2_curve   = interp1(time, values, 0:SIM_TIME);

  

\- mRNA synthesis (transcription) is comprised of a constitutive process and
an inducible process. This inducible process is dependent upon the
concentration of nuclear NFkB (NFkBn). We know from biochemical assays that
the inducible transcription of three of the components (IkBbt,   IkBet, and
P100t) are delayed and we included this delay strictly within the Stimulation
Phase of the MatLab mocel (not in the equilibrium phase).   This requires the
value of NFkBn at a previous time point.   As SimBiology only supports SBML
2.1 (as of January 2007), the delay functions are NOT included in the SBML
model.   As tools include support SBML 2.2, you can add these reactions.

IkBet = 45 minute delay

IkBbt = 45 minute delay

P100t = 90 minute delay

  

This model originates from BioModels Database: A Database of Annotated
Published Models (http://www.ebi.ac.uk/biomodels/). It is copyright (c)
2005-2011 The BioModels.net Team.  
To the extent possible under law, all copyright and related or neighbouring
rights to this encoded model have been dedicated to the public domain
worldwide. Please refer to [CC0 Public Domain
Dedication](http://creativecommons.org/publicdomain/zero/1.0/) for more
information.

In summary, you are entitled to use this encoded model in absolutely any
manner you deem suitable, verbatim, or with modification, alone or embedded it
in a larger context, redistribute it, commercially or not, in a restricted way
or not..  
  
To cite BioModels Database, please use: [Li C, Donizelli M, Rodriguez N,
Dharuri H, Endler L, Chelliah V, Li L, He E, Henry A, Stefan MI, Snoep JL,
Hucka M, Le Novère N, Laibe C (2010) BioModels Database: An enhanced, curated
and annotated resource for published quantitative kinetic models. BMC Syst
Biol., 4:92.](http://www.ncbi.nlm.nih.gov/pubmed/20587024)


