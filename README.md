# pericardium
Laura's Pericardium code

This code uses the 2D Immersed Boundary Method to simulate fluid flow driven by a flexible 'heart' tube through a racetrack circulatory system using a neuromechanical model of pumping. Simulations are run with and without a pericardium, or a stiff boundary that encases the flexible heart to compare differences in pumping. 

The code is set up to run on UNC-CH's Kure cluster. 

In order to remove the pericardium, certain lines need to be commented out. These are:
Lines 826-829 
Lines 871-872
Lines 887-888

These lines calculate pericardium forces, spread them, and move them. 

To put an effective hole in the pericardium, *uncomment* Lines 849-854. This causes the points to act as markers only. 