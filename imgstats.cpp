#include <iostream>
#include <string>

#include "../matrix/matrixmath.h"
#include "../bcBMPLib/bcBMPLib.h"

using namespace std;

int main()
{
	// Image
	matrix<int>* img = NULL;

	// Number of image layers
	int nlayers = 0;

	// Image Name
	string img_name;

	// Pointer to the layer name array
	string* layerNamePtr = NULL;

	// Layer name arrays
	string	layer1 = "Gray",
			layers3[3] = {"Red","Green","Blue"},
			layers4[4] = {"Red","Green","Blue", "Alpha"};

	// Descriptive statistics about image
	double	*max = NULL,				// Maximum pixel value
			*min = NULL,				// Minimum pixel value
			*range = NULL,				// max - min
			*mean = NULL,				// Mean pixel value
			*stdev = NULL,				// The population standard deviation of
										// pixel values
			*std_err = NULL,			// The standard error of pixel values
			*stdev_divRange = NULL,		// The stdev divided by range
			*stdev_divMean = NULL;		// The stdev divided bye mean

	// Prompt user for image name
	cout	<< "================="
			<< endl
			<< "Enter image name: ";
	cin		>> img_name;
	cout	<< "================="
			<< endl;

	// Read in image and display header information
	nlayers = readBMP(img, img_name, 1);

	// Allocate memory for image statistics
	max = new double [nlayers];
	min = new double [nlayers];
	range = new double [nlayers];
	mean = new double [nlayers];
	stdev = new double [nlayers];
	std_err = new double [nlayers];
	stdev_divRange = new double [nlayers];
	stdev_divMean = new double [nlayers];

	// Set layerNamePtr to point to the correct layer name array
	switch (nlayers)
	{
		case 1:
			layerNamePtr = &layer1;
			break;

		case 3:
			layerNamePtr = layers3;
			break;

		case 4:
			layerNamePtr = layers4;
			break;
	}

	// Calculate and display layer statistics
	for (int n = 0; n < nlayers; n++)
	{
		max[n] = max2d(img[n]);
		min[n] = min2d(img[n]);
		range[n] = max[n] - mean[n];
		mean[n] = mean2d(img[n]);
		stdev[n] = vari2d(img[n], 0);
		std_err[n] = stdev[n]/
			sqrt((double)(img[n].height()*img[n].width()));
		stdev_divRange[n] = stdev[n]/range[n];
		stdev_divMean[n] = stdev[n]/mean[n];

		cout	<< "==================="
				<< endl
				<< layerNamePtr[n]
				<< " Layer"
				<< endl
				<< "==================="
				<< endl
				<< "Max:                " << max[n] << endl 
				<< "Min:                " << min[n] << endl
				<< "Range:              " << range[n] << endl
				<< "Mean:               " << mean[n] << endl
				<< "Standard Deviation: " << stdev[n] << endl
				<< "Standard Error:     " << std_err[n] << endl
				<< "StDev/Range:        " << stdev_divRange[n] << endl
				<< "StDev/Mean:         " << stdev_divMean[n] << endl
				<< endl;
	}

	// Free memory of image statistic variables
	delete[] max;
	delete[] min;
	delete[] range;
	delete[] mean;
	delete[] stdev;
	delete[] std_err;
	delete[] stdev_divRange;
	delete[] stdev_divMean;


	return 0;
}
