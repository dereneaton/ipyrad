.. include:: global.rst

.. _assembly_methods:

Assembly Methods
================
ipyrad_ has four methods for Assembling data sets. The first and simplest 
is denovo_, which requires no prior information or genomic resources, while the 
remaining three methods all require some kind of reference genome. It is 
important to note, however, that these methods do not require a complete nuclear
genome to be useful, but they can also be informed by data from plastomes,
transcriptomes, and you can even use the genomes of potential contaminants to 
be filtered out of denovo data sets. 

.. _denovo:  
denovo
------
Sequences are assembled without any reference resources. Homology is inferred during
alignment clustering by sequence similarity using the program vsearch_. 

.. _reference:  
reference
---------
Sequences are mapped to a reference genome using the program smalt_ based on 
sequence similarity. 

reference_add
-------------
Sequences are mapped to a reference genome using the program smalt_ based on 
sequence similarity, and reads that do not match to the reference are assembled
using the denovo method. 

reference_sub
--------------
Sequences which map to a reference genome are excluded, and all remaining reads
are assembled using the denovo method. This method can be used to filter out 
data which match to a chloroplast genome in plants, or to a host genome in a 
study of a parasite. 

