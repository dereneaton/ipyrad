#!/usr/bin/env python

"""Globals used for ipa.mrbayes tool.

"""


# template for standard tree inference
NEX_TEMPLATE_1 = """\
#NEXUS

log start filename={outname}.log replace;
execute {nexus};

begin mrbayes;
set autoclose=yes nowarn=yes;

lset nst=6 rates=gamma;

{constraints}

mcmcp ngen={ngen} nrun={nruns} nchains={nchains};
mcmcp relburnin=yes burninfrac=0.25; 
mcmcp samplefreq={samplefreq} printfreq=10000;
mcmcp filename={outname};
mcmc;

sump filename={outname};
sumt filename={outname};
end;
"""


# template for clock model tree inference
# https://www.groundai.com/project/molecular-clock-dating-using-mrbayes/
TEMPLATE_2_DICT = {
    "brlenspr": "clock:birthdeath",
    "clockratepr": "lognorm(-7,0.6)",
    "clockvarpr": "tk02",
    "tk02varpr": "exp(1.0)",
    "samplestrat": "diversity",
    "sampleprob": "0.1",
    "speciationpr": "beta(2, 200)",
    "treeagepr": "offsetexp(1,5)",
}
NEX_TEMPLATE_2 = """\
#NEXUS

log start filename={outname}.log replace;
execute {nexus};

begin mrbayes;
set autoclose=yes nowarn=yes;

lset nst=6 rates=gamma;

prset brlenspr={brlenspr};
prset clockratepr={clockratepr};
prset clockvarpr={tk02varpr};
prset tk02varpr={tk02varpr};

prset samplestrat={samplestrat};
prset sampleprob={sampleprob};
prset speciationpr={speciationpr};
prset extinctionpr={extinctionpr};
prset treeagepr={treeagepr};

{constraints}

mcmcp ngen={ngen} nrun={nruns} nchains={nchains};
mcmcp relburnin=yes burninfrac=0.25;
mcmcp samplefreq={samplefreq};
mcmcp printfreq=10000 diagnfr=5000;
mcmcp filename={outname};
mcmc;

sump filename={outname};
sumt filename={outname};
end;
"""



# template for clock model tree inference
# https://www.groundai.com/project/molecular-clock-dating-using-mrbayes/
TEMPLATE_3_DICT = {
    "brlenspr": "clock:uniform",
    "clockvarpr": "igr",
    "igrvarpr": "exp(10.0)",
    "clockratepr": "normal(0.01,0.005)",
    "topologypr": "uniform"
}
NEX_TEMPLATE_3 = """\
#NEXUS

[log block]
log start filename={outname}.log replace;

[data block]
execute {nexus};

[tree block]
{treeblock}

[mb block]
begin mrbayes;
  set autoclose=yes nowarn=yes;

  lset nst=6 rates=gamma;

  prset brlenspr={brlenspr};
  prset clockvarpr={clockvarpr};
  prset igrvarpr={igrvarpr};
  prset clockratepr={clockratepr};
  prset topologypr={topologypr};

  {constraints}

  mcmcp ngen={ngen} nrun={nruns} nchains={nchains};
  mcmcp relburnin=yes burninfrac=0.25;
  mcmcp samplefreq={samplefreq};
  mcmcp printfreq=10000 diagnfr=5000;
  mcmcp filename={outname};
  mcmc;

  sump filename={outname};
  sumt filename={outname} contype=allcompat;
end;

[log block]
log stop filename={outname}.log append;
"""
