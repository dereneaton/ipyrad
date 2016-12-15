

Go to https://aws.amazon.com/hpc and click on the link that says ``Create a Free Account``, 
then either sign in with your existing Amazon account, or create a new one. You will then 
be asked to enter your credit card information, and to pass a security verification. Go ahead
and do that now. Finally, you will be asked to choose a Support Plan. Choose the Basic plan, 
which is *free*. 


You will next be greeted with a message welcoming you to Amazon Web Services. Click 
on *Launch Management Console* to get to work. 

Go to the console page, click on your instance, and a bunch of information about it will
appear at the bottom of the page. Find the part that says ``Public DNS`` and copy the
address. It will look something like this: *ec2-54-244-25-91.us-west-2.compute.amazonaws.com*.

Now we will follow instructions very similar to those provided for launching a jupyter-notebook
on an Institutional HPC cluster, but with slight variations for the EC2 cluster. One important 
difference is that because the EC2 instance is *temporary*, you cannot store your data for long-term
on its cloud platform. Therefore we need to transfer the data to the instance when we start, 
assemble it, and then transfer all of the results back to our local machine when finished, otherwise
they will be deleted (or you can pay to keep them stored online, but I'll assume you are cheap like me.)

### Choose an instance type. 
The free tier plan only let's you connect to one CPU with up to 30GB of storage, and
XXGB of RAM. There are many different paid options, however, that offer much greater 
storage and compute power. In general, you will want at least 1GB of RAM per CPU. 
So let's say we want 32 CPUs, >60GB RAM, and >1TB of storage. We can look on the website
http://www.ec2instances.info to find which instance offers the best price for these 
requirements. In this case it looks like it is probably the `Cluster Compute Eight 
Extra Large` Instance, which costs $1.090 per hour. If you had a data set with 96 
individuals of RAD data sequenced on two lanes of Illumina HiSeq I expect you could
finish your assembly in <24 hours using this setup. Therefore, your RAD-seq assembly
would cost you less that $24. Not too bad.  


### Get setup.
Download the aws command-line tools using pip

	pip install awscli  

Then follow the setup instructions here to create credentials for your login:
http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-set-up.html#cli-signup


### CLI remote commands...

### local commands...


### Step 1: login to the EC2 instance. 

