



.. _HPCscript:

HPC script
==========

The API requires starting the parallel Engines before executing step 
functions with ipyrad. This provides a lot of flexibility to run ipyrad 
on with complex connection setups. Below we show a fairly common setup 
to connect to multiple nodes on an HPC cluster. 


Starting ipcluster
^^^^^^^^^^^^^^^^^^^

First update ipyparallel to the most recent version. This should upgrade it 
to version >5.0. 

.. code-block:: bash

    pip install -U ipyparallel


Now you can start an ipcluster running. I recommend doing this in a separate
screen created using the ``screen`` unix command. We use the engines=MPI and 
--ip=* commands to tell ipcluster to connect to cores across multiple nodes.


.. code-block:: bash

    screen

    ## once in new screen type the following
    ipcluster start -n 32 --engines=MPI --ip=* 

    ## now disconnect from this screen by typing (ctrl-a, then d)


Running ipyrad interactively
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now open up a second screen which we will use to run the ipyrad API interactively.
The code below could alternatively be saved as a python script and run as 
`python myscript.py`. The ipyrad API will automatically use all available 
Engines from ipcluster. In this case, 32. 


.. code-block:: python

    ##
    ## First open an IPython session by typing `ipython` into a terminal.
    ##

	## imports
	import ipyrad as ip

    ## make Assembly object
    data = ip.Assembly("test")

    ## Fill params. Set the location of your files here.
    data.set_params("project_dir", "test")
    data.set_params("raw_fastq_path", "iptest/sim_rad1_R1_.fastq.gz")
    data.set_params("barcodes_path", "iptest/sim_rad1_barcodes.txt")

    ## set subsampling for step 2
    data._hackersonly["preview_step2"] = 2000

    ## print params and hackers params
    data.get_params()
    print "hackers dict"
    print data._hackersonly

    ## demultiplex without preview mode. This will take a while.
    data.step1()

    ## run step2 with preview mode. This should run very fast.
    data.step2(preview=True)

    ## save the commands from this session to a file
    %save my_ipyrad_script.py 



Now while the code is running you can disconnect from this session 
(again ctrl-a, then d) and watch the cpus working away using the unix 
``top`` command. And you can peek into the output directory. 
Finally, when the job is done you can go back in and look at the 
resulting stats for your assembly by reconnecting to the interactive 
IPython session using ``screen -r``. 

 .. code-block:: python

    ## print stats
    data.stats

    ## run the next step
    data.step3()

    ## etc.




