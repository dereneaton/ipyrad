.. include:: global.rst

.. _ipyrad_scripts:

**************
ipyrad scripts
**************

The basic workflow of an ipyrad script is to create an Assembly object, 
assign parameter values to this object, link raw data to it, and then 
execute assembly functions. 

Assembly objects
----------------
You can think of an assembly object as a sort of map to your data. It 
stores the parameters used in assembly steps, the location of data files
on disk, and the statistical output from each step of the assembly. 


To create a new Assembly object named *data*:  

.. code:: python

    data1 = ip.Assembly("data1")


.. code:: parsed-literal  

    new object created...


.. _Assembly

Setting parameters
------------------
Use the get_params() call to show the current parameter settings:
    data.get_params()

To change a parameter use the set_params() call:
    data1.set_params(1, "tests")
    data1.set_params(2, "tests/data/*.gz")
    data1.set_params(3, "tests/data/testrad_barcodes.txt")  

To get more info on how a specific parameter influences the assembly you 
can use ip.get_params_info(N), if you are working interactively, otherwise
you can look in the Parameters_ section. 

    ip.get_params_info(1)   

.. _Setting_parameters

Sample objects
--------------


.. _Sample 


Saving Assembly objects
=======================
Assembly objects are auto-saved every time an assembly function (step or run) 
is called, or it can be saved by the user by calling save_object(). This stores
the object to disk so that it can be re-loaded at a later time. Because the 
Assembly object does not contain any of the sequence data, but simply links to 
its location, the saved object is quite small and so can easily be moved between
computers. This way you can execute computationally heavy jobs on an HPC machine, 
and later load the Assembly object on your laptop to examine the results 
in greater detail, to run downstream analyses_, or to explore through plotting_. 

.. _Saving

