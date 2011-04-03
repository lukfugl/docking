run "python example.py" to see it in action. right now, it's parsing the 2PWS.xbgf file and scoring the chain against itself... it has a divide by zero check, so when an atom is scored against itself, it's zero, not +/- infinity. what we really want is a file with two chains in it (one for the protein, one for the ligand), then we can score the chains against each other, and rotate just the ligand to find an optimal configuration.

performance is not terrible, but not terrific, either... I'll be looking into that a little more.

tested with Python 2.7.1, BioPython 1.56, and NumPy 1.5.1 (all on 32-bit Windows 7)

BioPython can be obtained from http://biopython.org/wiki/Download
NumPy can be obtained from http://new.scipy.org/download.html

find latest versions of the code at https://github.com/lukfugl/docking

you can check this code out using git. You can get git for windows from http://code.google.com/p/msysgit/downloads/detail?name=Git-1.7.4-preview20110204.exe&can=3&q=. I know you can get git for OSX as well, thought I can't remember the details; git for linux/unix is pretty straightforward.