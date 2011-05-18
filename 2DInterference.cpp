#include <iostream>
#include <cmath>
#include <string>

#include "../matrix/matrix.h"
#include "../matrix/matrixmath.h"
#include "../bcBMPLib/bcBMPLib.h"

# define PI (double)3.14159

using namespace std;

// DESCRIPTION: This program creates 2D Interference patterns and saves them as
//				8bit Bitmap files.
//
// AUTHOR: Bryan A. Clifford
// LAST MODIFIED: 01.19.2011
int main()
{
	//__________________________________________________________________________
	// Experimental variables (to be changed by the user)
	double	height = 1e-6,			// Real height of config in meters
			width = .5e-6,			// Real width of config in meters

			waveLength = 50e-9,		// Wavelength of incident plane wave in a
									// vacuum in meters

			rfcIdx = 1,				// Refractive index of the medium
			img_width = 2e-6,		// Real width of image in meters
			img_height = 2e-6,		// Real height of image in meters
			vres = .75e-9,			// Real height of a pixel in meters
			hres = .75e-9;			// Real width of a pixel in meters

	// Name of the image file
	string	img_name = "interference.bmp";


	//__________________________________________________________________________
	// Wave # of the incident plane wave
	double k = 0;

	// Important diagonal lengths
	double img_length = 0,
		   pxl_length = 0;
	
	// Array to contain coordinates of the configuration point sources
	int**	points = NULL;

	// Number of point sources
	int		npoints = 0;

	// Important image variables
	matrix<int>		img_int;		// The actual image
	matrix<double>	img;			// The image used for calculations
	int 			nrows = 0,		// The number of rows in the image
					ncols = 0;		// The number of columns in the image

	// The boundry points of the configuration
	int	ctr_row = 0,				// The indices of the center point of the
		ctr_col = 0,				// configuration

		top = 0,					// The rows corresponding to the top and
		bot = 0,					// bottom edges of the configuration

		left = 0,					// The columns corresponding to the left and
		right = 0;					// right edges of the configuration

	// Descriptive statistics about image
	double	max = 0,			// Maximum pixel value
			min = 0,			// Minimum pixel value
			range = 0,			// max - min
			mean = 0,			// Mean pixel value
			stdev = 0,			// The population standard deviation of
								// pixel values
			std_err = 0,		// The standard error of pixel values
			stdev_divRange = 0,	// The stdev divided by range
			stdev_divMean = 0;	// The stdev divided bye mean

	// Threshold values that can be used to amplify the interference pattern
	double	thresh_hi = 0,	// Upper threshold
			thresh_lo = 0;	// Lower threshold
	
	// Temporary coordinate pairs (x,y)
	int	x0 = 0,
		x1 = 0,
		y0 = 0,
		y1 = 0;

	// Temporary distances
	double	dx = 0,
			dy = 0,
			r = 0;

	// Temporary value
	double val = 0;
		


	//__________________________________________________________________________
	// Calculates wave #
	k = 2*PI*rfcIdx/waveLength;

	// Calculates diagonal lengths
	img_length = sqrt(img_height*img_height + img_width*img_width);
	pxl_length = sqrt(vres*vres + hres*hres);

	// Calculates size of image
	nrows = img_height/vres;
	ncols = img_width/hres;

	// Calculates boundry points of the configuration
	ctr_row = nrows/2;
	ctr_col = ncols/2;
	top = ctr_row - .5*height/vres;
	bot = ctr_row + .5*height/vres;
	left = ctr_col - .5*width/hres;
	right = ctr_col + .5*width/hres;

	// Calculates/sets the number of point sources in the image
	npoints = 5;

	// Resizes image  matrices and fills them with 0s
	img_int(0,0) = 0;
	img(0,0) = 0;
	img_int.grow((nrows - 1), 0, (ncols - 1), 0, 0);
	img.grow((nrows - 1), 0, (ncols - 1), 0, 0);

	// Allocates memory for point source coordinates
	points = new int* [npoints];
	for (int i = 0; i < npoints; i++)
	{
		points[i] = new int [2];
	}

	// Calculates/sets the coordinates of the point sources in the image
	points[0][0] = ctr_row;
	points[0][1] = ctr_col;

	points[1][0] = top;
	points[1][1] = ctr_col;

	points[2][0] = bot;
	points[2][1] = ctr_col;

	points[3][0] = ctr_row;
	points[3][1] = left;

	points[4][0] = ctr_row;
	points[4][1] = right;



	//__________________________________________________________________________
	// Creates amplitude interference patterns
	for (int n = 0; n < npoints; n++)
	{
		// Set current center point
		cout << "Creating interference pattern from source " << n << endl;
		x0 = points[n][1];
		y0 = points[n][0];

		for (y1 = 0; y1 < nrows; y1++)
		{
			for (x1 = 0; x1 < ncols; x1++)
			{
					// Calculate the real distance (r) from (x0,y0) to (x1,y1)
					// r = sqrt((x1-x0)^2 + (y1 - y2)^2)
					dx = ((double)(x1 - x0))*hres;
					dy = ((double)(y1 - y0))*vres;
					r = sqrt(dx*dx + dy*dy);

					// Calculate the undamped amplitude at the point (x1,y1)
					// from the current current source (n)
					val = cos(k*r);

					// Apply damping to the pattern
					// (Disappear linearly towards edges)
					val *= pow(
								(
									((double)1 - 
										(double)2*(double)abs(x1 - ctr_col)/
										((double)ncols)
									)*
									((double)1 - 
										(double)2*(double)abs(y1 - ctr_row)/
										((double)nrows)
									)
								),.5);

					// Add value to image
					img(y1,x1) += val;
			}
		}
		delete[] points[n];
	}
	delete[] points;

	// Calculate the intensities
	cout << "Calculating intensity values" << endl;
	for (int row = 0; row < nrows; row++)
	{
		for (int col = 0; col < ncols; col++)
		{
			img(row,col) *= img(row,col);
		}
	}


	//__________________________________________________________________________
	// Calculate interference pattern statistics
	cout << "Calculating interference pattern statistics" << endl;
	max = max2d(img);
	min = min2d(img);
	range = max - min;
	mean = mean2d(img);
	stdev = stdev2d(img, 0);
	std_err = stdev/sqrt((double)(nrows*ncols));
	stdev_divRange = stdev/range;
	stdev_divMean = stdev/mean;

	// Display interference pattern statistics
	cout	<< "Max:                " << max << endl 
			<< "Min:                " << min << endl
			<< "Range:              " << range << endl
			<< "Mean:               " << mean << endl
			<< "Standard Deviation: " << stdev << endl
			<< "Standard Error:     " << std_err << endl
			<< "StDev/Range:        " << stdev_divRange << endl
			<< "StDev/Mean:         " << stdev_divMean << endl;


	//__________________________________________________________________________
	// Apply thresholding to image to make pattern more visible (pretty)
	cout << "Thresholding pattern" << endl;
	thresh_hi = 4*stdev*stdev + mean;
	thresh_lo = mean - 4*stdev*stdev;
	for (int row = 0; row < nrows; row++)
	{
		for (int col = 0; col < ncols; col++)
		{
			if (img(row,col) > thresh_hi)
			{
				img(row,col) = thresh_hi;
			}

			if (img(row,col) < thresh_lo)
			{
				img(row,col) = thresh_lo;
			}
		}
	}
	max = max2d(img);


	//__________________________________________________________________________
	// Equate image matrices for the creation of a Bitmap
	cout << "Normalizing pattern for BMP creation" << endl;
	for (int row = 0; row < nrows; row++)
	{
		for (int col = 0; col < ncols; col++)
		{
			// Must normalize the image so that the max is 255
			img_int(row,col) = 255*img(row,col)/max;
		}
	}

	// Saves image as a 8bit Bitmap image
	cout << "Saving pattern as " << img_name << endl;
	makeBMP(&img_int,1,img_name);
	cout << "Done" << endl;

	return 0;
}
