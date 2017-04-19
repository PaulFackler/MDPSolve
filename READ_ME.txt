SETUP Instructions for MDPSolve

To install the MDPSolve package:

1) Create a directory called MDPSOLVE on your disk (this could be a hard disk or
     a removable disk). Unzip all of the files to this directory. Make sure that you 
     maintain the directory structure (for example in WinZip make sure the 
     "Use Folder Names" box is checked).


2) Make sure the directory structure is correct - this file should be in
     a directory called MDPSOLVE and it should contain 6 subdirectories:
     mdpdemos, mdptools, mdputils, probability influence and factoredbeliefs.
     For previous users of MDPSolve the last two directories are new.
     The influence directory and its subdirectories contain new procedures for model
     specification that use a graphical factored approach. The factoredbelief
     directory contains experimental code that was used in a paper describing the
     new xpomdp approach (Paul L. Fackler and Krishna Pacifici. Addressing Structural
     and Observational Uncertainty in Resource Management. Journal of Environmental 
     Management, 133 (2014): 27-36).

         
3) Start MATLAB, make MDPSOLVE the default directory and set the MATLAB path to include
      MDPSOLVE and all of its subdirectories. To do this type
          cd('c:\MDPSOLVE')
          path(path,genpath(cd))
      You will probably need to modify "C:\" this so the complete directory where you put
      MDPSOLVE is given. Also on UNIX systems you may need to use "/" rather than "\".
      These instructions can be put in your startup.m file so they are executed each time
      you load MATLAB.

      An alternative is to set the MATLAB path interactively. Click on 
          Files/Set Path/Add with Subfolders
      and select the MDPSOLVE directory. You may be given the option to save this path
      and if you can then you will not need to repeat these steps. If your machine
      does not give you rights to save, however, you may not be able to. In this
      case simply repeat the procedure each time you start a MATLAB session.

     
4) A number of routines in MDPSOLVE run far faster if MEX files are available.
     MEX files are files written in C that can be called from MATLAB. There are a 
     number of these files in the MDPSOLVE\mdputils subdirectory. The original
     C source code provided and I've provided compiled MEX files for Windows users using
     either a 32 or 64 bit computer. Type "which lookup" - if this returns lookup.mex32 
     or lookup.mex64 you can probably skip this step. 
  
     The MEX files may not work for you however because you have a non-Windows operating
     system or because your version of MATLAB is not compatible with the ones I used
     to complile the MEX files. If this is the case you should compile the MEX files by running
     "mdpmexall" (which is located in the main MDPSLOLVE directory). 

     If this is the first time you have compiled a MEX file you should be given a list of
     C compilers to choose from. The 32 bit version of MATLAB is distributed with the LCC 
     compiler and it will work. It is a good choice if you are new to C compilers. 
     For 64 bit MATLAB users, you will need to install your own compiler. Go to
         http://www.mathworks.com/support/compilers/R2013a
     for information on supported compilers (change the address to get the proper release 
     that you are using). An easy choice is Microsoft Windows SDK 7.1 which can be
     downloaded and installed it at no charge.

     All of the functions of MDPSOLVE run without the MEX files so you can skip this step
     at first anyway. There are functions that run very slowly or are more likely to 
     run into memory problems if the m-file rather than the MEX-file versions are used.
     For this reason I suggest you make an effort to compile the MEX files.
 
       
5) Run runmdpdemos to have a set of the demos in MDPSOLVE\mdpdemos subdirectory run 
     one after the other. If this file does not run without errors or warnings, go  
     back through the last 3 steps and make sure you have followed the instructions 
     carefully. If you are still having problems post a message to the MDPSolve forum
        https://groups.google.com/forum/#!forum/mdpsolve
     If all else fails contact me with details about the nature of the errors that 
     you are getting (paul_fackler@ncsu.edu).


6) If you want to understand what this program does and how it works and to have 
     the demos explained, read the User's Guide, which is called MDPSOLVEDoc.pdf
     and is located in the MDPSOLVE directory. Be forewarned however, the User's Guide 
     is "under construction" and there are known and unknown gaps in its coverage.
     Suggestions for improvement are welcome.

7) This is the seventh step.

Trying GitHub. Things are still not working. 

Okay - let's introduce another conflict.

Eric is working on his version and I'm on mine.
