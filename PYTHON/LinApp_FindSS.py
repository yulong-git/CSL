'''
Version 1.0, written by Kerk Phillips, April 2014

Adapted by Yulong Li, November 2015
'''

def LinApp_FindSS(funcname,param,guessXY,Zbar,nx,ny):
	'''
	Finds the steady state for a DSGE model numerically

	Parameters
    -----------
    funcname: function
        the name of the function which generates a column vector 
        from ny+nx dynamic equations.
		The ny equations to be linearized into the form below in the first 
		ny rows.
		A X(t) + B X(t-1) + C Y(t) + D Z(t) = 0 
		The function must be written so as to evaluate to zero for all rows
		in the steady state.

	param: array, dtype=float
		a vector of parameter values to be passed to funcname.
	
	guessXY: array, dtype=float
		guess for the steady state values of X and Y
	
	Zbar: array, dtype=float
		1-by-nz vector of Z steady state values
	
	nx: number, int
		number of X variables
	
	ny: number, int
		number of Y variables

	Returns
    --------
	XYbar: array, dtype=float
		1-by-(nx+ny) vector of X and Y steady state values, with the X
		values in positions 1 - nx and the Y values in nx+1 - nx+ny.
	'''
	f = @(XYbar) steady(XYbar, Zbar, funcname, param, nx, ny);


	%  use fsolve
	XYbar = fsolve(f,guessXY,options);

	end

	function out = steady(XYbar, Zbar, funcname, param, nx, ny)
	Xbar = XYbar(1:nx);
	Ybar = XYbar(1+nx:nx+ny);
	in = [Xbar; Xbar; Xbar; Ybar; Ybar; Zbar; Zbar];
	out = funcname(in,param);
