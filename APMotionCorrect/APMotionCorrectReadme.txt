APMotionCorrect GUI v. 2.3 beta
Andrew Peters
peters.andrew.j@gmail.com
Email me to help me make it better!

----Installation instructions:

Install ImageJ (online)

Place APMotionCorrect folder and subfolders in matlab path

Move all .java and .class files from "MatlabTurboreg\copy to java" into your matlab java folder (i.e. C:\Program Files\MATLAB\R2010b\java)

Move all .jar files from "MatlabTurboreg\copy to jar" into your matlab java\jar folder (i.e. C:\Program Files\MATLAB\R2010b\java\jar)

----To run:

Type "APMotionCorrect" into matlab

----Operating instructions:
Select Save Path, TIFF files (files to be corrected), and Ref file (if you have a reference), click Add File to add to list

[Multiple TIFF and reference files may be selected at once, they will not show up in the edit box but should appear in the file list once "Add File" is pressed]

Three motion correction algorithms: Greenberg (line based, not very efficient), Turboreg (frame based, great), HMM (line based, good)

Options: 
Automatic Reference: uses average images for alignment
Playback: display corrected movie in matlab when finished
Channels: Number of channels (motion correct currently only fixes green channels)
Multiple Sessions: check this if you want to correct multiple imaging sessions sequentially. If checked, puts group number in front of files to be corrected, and then corrects within groups together (if parsed file is selected) and between groups seperately

HMM - maxdx/dy: max shift of each row/column in pixels
lambda: pixel shift/scan line time (movement velocity)
parsed file: selected files are all of one movie, broken up
(this feature causes HMM or turboreg->HMM to use references that are averaged over all files when autoreference is selected)

Settings that have worked the best so far: 
Autoreference
Turboreg -> HMM
(on 2x and above, double for 1x) dx = 6, dy = 2, lambda = 0.01

WARNING On Windows 7: close windows explorer before running correction. If it is not closed, there will likely be an error about write permission while matlab is creating TIFF files

 

