===================================
Current State Linearization Toolkit
===================================


Abstract
========
This code repository stores the code needed to execute the current state linearization method discussed in Evans and Phillips (2015) '"Linearization about the Current State: A Computational Method for Approximating Nonlinear Policy Functions during Simulation."<https://drive.google.com/file/d/0B6KGaihAO5TJZGJLemE1V1d5bFE/view>'_

This toolkit can be used to solve and simulate using either the CSL method from the paper or to solve and simulate using standard linearization about the steady state (SSL).  

This code makes extensive use of Harald Uhlig's MATLAB code version 2.0 from H. Uhlig (1995).

Currently this software is available only in MATLAB.  A Python version is forthcoming.

Code examples for the model from Brock and Mirman (1978), a simple RBC model, and Hansen (1985) are included.


Sections
========
* MATLAB files
* Python files


Researchers
===========
- Rick Evans
- Kerk Phillips


License
=======

The following copyright license restrictions apply:

Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.  However, you are not allowed to sell this software or otherwise impinge on its free distribution.

Copyright: K. Phillips.  Feel free to copy, modify and use at your own risk.  However, you are not allowed to sell this software or otherwise impinge on its free distribution.


References
==========

Brock, William A. and Leonard Mirman, (1972) "Optimal Economic Growth and Uncertainty: the Discounted Case," Journal of Economic Theory, 4(3), 479-513.

Hansen, Gary, (1985), "Indivisible labor and the business cycle, Journal of Monetary Economics, 16(3), 309-327.

Uhlig, Harald, (1995) "A Toolkit for Solving Nonlinear Dynamic Stochastic Models Easily," Discussion Paper, Institute for Empirical Macroeconomis, Federal Reserve Bank of Minneapolis #101 or Tilburg University, CentER DP 9597.