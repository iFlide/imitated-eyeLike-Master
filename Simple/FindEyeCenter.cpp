#include "FindEyeCenter.h"
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <queue>

int result_left_x{ 0 };
int result_left_y{ 0 };
int result_right_x{ 0 };
int result_right_y{ 0 };

int leftmaxPx{ 0 };            //center in left eye
int leftmaxPy{ 0 };

int rightmaxPx{ 0 };           //center in right eye
int rightmaxPy{ 0 };

int facewidth{ 0 };            //input by user
int faceheight{ 0 };		   //input by user

int eye_region_width{ 0 };     //estimated by user
int eye_region_height{ 0 };	   //estimated by user
int eye_region_top{ 0 };       //estimated by user
int eye_region_leftside{ 0 };  //estimated by user
int eye_region_rightside{ 0 }; //estimated by user

int eyeROIwidth{ 0 };          //width after being scaled
int eyeROIheight{ 0 };         //height after being scaled

double leftgradientThresh{ 0 };     //gradient thresh
double rightgradientThresh{ 0 };

double leftfloodThresh{ 0 };        //image thresh
double rightfloodThresh{ 0 };

int faceROI[400][400] = { { 0 } };         //input face

int leftEyeRegion[300][300] = { { 0 } };  //estimated left eye
int rightEyeRegion[300][300] = { { 0 } };

double lefteyeROI[100][100] = { { 0 } };  //left eye after scaled
double righteyeROI[100][100] = { { 0 } };

double leftgradientX[100][100] = { { 0 } };//left eye's gradient
double leftgradientY[100][100] = { { 0 } };

double rightgradientX[100][100] = { { 0 } };
double rightgradientY[100][100] = { { 0 } };

double leftmags[100][100] = { { 0 } };     //left eye's magnititude
double rightmags[100][100] = { { 0 } };

double leftweight[100][100] = { { 0 } };   //weight of each pixel
double rightweight[100][100] = { { 0 } };

double leftoutSum[100][100] = { { 0 } };   //outcome after weight*mag
double rightoutSum[100][100] = { { 0 } };

double leftout[100][100] = { { 0 } };      //change outsum's scale
double rightout[100][100] = { { 0 } };

double leftfloodClone[100][100] = { { 0 } };  //floodclone is out after threshing
double rightfloodClone[100][100] = { { 0 } };

double leftmask[100][100] = { { 0 } };        //area should be computed in leftout
double rightmask[100][100] = { { 0 } };


using namespace std;
FindEyeCenter::FindEyeCenter()
{
	leftmaxPx = 0;
	leftmaxPy = 0;
	rightmaxPx = 0;
	rightmaxPy = 0;
	result_left_x = 0;
	result_left_y = 0;
	result_right_x = 0;
	result_right_y = 0;
}
FindEyeCenter::~FindEyeCenter()
{

}

void FindEyeCenter::run(int &lx,int &ly,int &rx,int &ry)
{
	//TODO: to do a complete find-eye-center
	ComputeEye();
	resize(eyeROIwidth,eyeROIheight);
	ComputeXGradient();
	ComputeYGradient();
	MatrixMagnitude();
	computeDynamicThreshold();
	Thresh();
	ComputeWeight(eyeROIheight, eyeROIwidth);
	GetleftCenter();
	GetrightCenter();
	Convert();
	leftFloodKillEdges();
	rightFloodKillEdges();
	leftmaxCenter();
	rightmaxCenter();
	unscaleleftcenter();
	unscalerightcenter();
	result_left_x = eye_region_leftside + leftmaxPx;
	result_left_y = eye_region_top + leftmaxPy;
	result_right_x = eye_region_rightside + rightmaxPx;
	result_right_y = eye_region_top + rightmaxPy;
	lx = result_left_x;
	ly = result_left_y;
	rx = result_right_x;
	ry = result_right_y;
}

int FindEyeCenter::NearDoubleNumber(double D)
{
	double NewNum = floor(D);//get near number
	return ((D - NewNum) - 0.5 <= 0) ? NewNum : NewNum + 1;
}

bool FindEyeCenter::shouldpush(int i, int j)
{
	if (i >= 0 && i < eyeROIheight&&j >= 0 && j < eyeROIwidth)
		return true;
	return false;
}

void FindEyeCenter::ComputeEye()
{
	eye_region_width = facewidth * (kEyePercentWidth / 100.0);
	eye_region_height = facewidth * (kEyePercentHeight / 100.0);
	eye_region_top = faceheight * (kEyePercentTop / 100.0);
	eye_region_leftside = facewidth*(kEyePercentSide / 100.0);
	eye_region_rightside = facewidth - eye_region_width - eye_region_leftside;

	//get leftEyeRegion
	for (int i = 0; i <  eye_region_height; i++)
	{
		for (int j = 0; j <  eye_region_width; j++)
		{
			leftEyeRegion[i][j] = faceROI[i + eye_region_top][j + eye_region_leftside];
		}
	}
	//get rightEyeRegion
	for (int i = 0; i < eye_region_height; i++)
	{
		for (int j = 0; j < eye_region_width; j++)
		{
			rightEyeRegion[i][j] = faceROI[i + eye_region_top][j + eye_region_rightside];
		}
	}
	//get the width and height of eyeROI that we want
	eyeROIwidth = kFastEyeWidth;
	eyeROIheight = (kFastEyeWidth / (double)eye_region_width*(double)eye_region_height);

}

void FindEyeCenter::resize(double width,double height)
{
	double ratex = eye_region_width / width;
	double ratey = eye_region_height / height;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int ni = (i*ratey);
			int nj = (j*ratex);
			lefteyeROI[i][j] = leftEyeRegion[ni][nj];
			righteyeROI[i][j] = rightEyeRegion[ni][nj];
		}
	}
}

void FindEyeCenter::ComputeXGradient()
{
	for (int i = 0; i < eyeROIheight; i++)
	{
		leftgradientX[i][0] = lefteyeROI[i][1]-lefteyeROI[i][0];
		rightgradientX[i][0] = righteyeROI[i][1] - righteyeROI[i][0];
		for (int j = 1; j < eyeROIwidth-1; j++)
		{
			leftgradientX[i][j] = (lefteyeROI[i][j + 1] - lefteyeROI[i][j - 1]) / 2.0;
			rightgradientX[i][j] = (righteyeROI[i][j + 1] - righteyeROI[i][j - 1]) / 2.0;
		}
		leftgradientX[i][eyeROIwidth - 1] = lefteyeROI[i][eyeROIwidth - 1] - lefteyeROI[i][eyeROIwidth - 2];
		rightgradientX[i][eyeROIwidth - 1] = righteyeROI[i][eyeROIwidth - 1] - righteyeROI[i][eyeROIwidth - 2];
	}
}

void FindEyeCenter::ComputeYGradient()
{
	for (int j = 0; j < eyeROIwidth; j++)
	{
		leftgradientY[0][j] = lefteyeROI[1][j] - lefteyeROI[0][j];
		rightgradientY[0][j] = righteyeROI[1][j] - righteyeROI[0][j];
		for (int i = 1; i < eyeROIheight - 1; i++)
		{
			leftgradientY[i][j] = (lefteyeROI[i + 1][j] - lefteyeROI[i - 1][j]) / 2.0;
			rightgradientY[i][j] = (righteyeROI[i + 1][j] - righteyeROI[i - 1][j]) / 2.0;
		}
		leftgradientY[eyeROIheight - 1][j] = lefteyeROI[eyeROIheight - 1][j] - lefteyeROI[eyeROIheight - 2][j];
		rightgradientY[eyeROIheight - 1][j] = righteyeROI[eyeROIheight - 1][j] - righteyeROI[eyeROIheight - 2][j];
	}
}

void FindEyeCenter::MatrixMagnitude()
{
	for (int i = 0; i < eyeROIheight; i++)
	{
		for (int j = 0; j < eyeROIwidth; j++)
		{
			leftmags[i][j] = sqrt(leftgradientX[i][j] * leftgradientX[i][j] + leftgradientY[i][j] * leftgradientY[i][j]);
			rightmags[i][j] = sqrt(rightgradientX[i][j] * rightgradientX[i][j] + rightgradientY[i][j] * rightgradientY[i][j]);
		}
	}
}

void FindEyeCenter::computeDynamicThreshold()
{
	double leftmean = 0;
	double rightmean = 0;
	double leftStdDev = 0;
	double rightStdDev = 0;

	for (int i = 0; i < eyeROIheight; i++)
	{
		for (int j = 0; j < eyeROIwidth; j++)
		{
			leftmean += leftmags[i][j];
			rightmean += rightmags[i][j];
		}
	}
	leftmean = leftmean / (eyeROIheight*eyeROIwidth);
	rightmean = rightmean / (eyeROIheight*eyeROIwidth);
	for (int i = 0; i < eyeROIheight; i++)
	{
		for (int j = 0; j < eyeROIwidth; j++)
		{
			leftStdDev += (leftmags[i][j] - leftmean)*(leftmags[i][j] - leftmean);
			rightStdDev += (rightmags[i][j] - rightmean)*(rightmags[i][j] - rightmean);
		}
	}
	leftStdDev = sqrt(leftStdDev / (eyeROIheight*eyeROIwidth));
	rightStdDev = sqrt(rightStdDev/ (eyeROIheight*eyeROIwidth));
	//refree to source 
	leftStdDev = leftStdDev / (sqrt(double(eyeROIheight*eyeROIwidth)));
	rightStdDev = rightStdDev / (sqrt(double(eyeROIheight*eyeROIwidth)));
	leftgradientThresh = kGradientThreshold*leftStdDev + leftmean;
	rightgradientThresh = kGradientThreshold*rightStdDev + rightmean;
}

void FindEyeCenter::Thresh()
{
	for (int i = 0; i < eyeROIheight; i++)
	{
		for (int j = 0; j < eyeROIwidth; j++)
		{
			if (leftmags[i][j]>leftgradientThresh)
			{
				leftgradientX[i][j] = leftgradientX[i][j] / leftmags[i][j];
				leftgradientY[i][j] = leftgradientY[i][j] / leftmags[i][j];
			}
			else
			{
				leftgradientX[i][j] = 0;
				leftgradientY[i][j] = 0;
			}
			if (rightmags[i][j]>rightgradientThresh)
			{
				rightgradientX[i][j] = rightgradientX[i][j] / rightmags[i][j];
				rightgradientY[i][j] = rightgradientY[i][j] / rightmags[i][j];
			}
			else
			{
				rightgradientX[i][j] = 0;
				rightgradientY[i][j] = 0;
			}
		}
	}
}

void FindEyeCenter::ComputeWeight(int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
//TODO: dont forget to get the verges' value
			if (i<2 || i>rows - 3 || j<2 || j>cols - 3)
			{
				leftweight[i][j] = lefteyeROI[i][j];
				rightweight[i][j] = righteyeROI[i][j];
			}
			else
			{
				leftweight[i][j] =
					  lefteyeROI[i - 2][j - 2] * GuassianFilter[0][0] + lefteyeROI[i - 2][j - 1] * GuassianFilter[0][1] + lefteyeROI[i - 2][j] * GuassianFilter[0][2] + lefteyeROI[i - 2][j + 1] * GuassianFilter[0][3] + lefteyeROI[i - 2][j + 2] * GuassianFilter[0][4]
					+ lefteyeROI[i - 1][j - 2] * GuassianFilter[1][0] + lefteyeROI[i - 1][j - 1] * GuassianFilter[1][1] + lefteyeROI[i - 1][j] * GuassianFilter[1][2] + lefteyeROI[i - 1][j + 1] * GuassianFilter[1][3] + lefteyeROI[i - 1][j + 2] * GuassianFilter[1][4]
					+ lefteyeROI[i][j - 2] * GuassianFilter[2][0] + lefteyeROI[i][j - 1] * GuassianFilter[2][1] + lefteyeROI[i][j] * GuassianFilter[2][2] + lefteyeROI[i][j + 1] * GuassianFilter[2][3] + lefteyeROI[i][j + 2] * GuassianFilter[2][4]
					+ lefteyeROI[i + 1][j - 2] * GuassianFilter[3][0] + lefteyeROI[i + 1][j - 1] * GuassianFilter[3][1] + lefteyeROI[i + 1][j] * GuassianFilter[3][2] + lefteyeROI[i + 1][j + 1] * GuassianFilter[3][3] + lefteyeROI[i + 1][j + 2] * GuassianFilter[3][4]
					+ lefteyeROI[i + 2][j - 2] * GuassianFilter[4][0] + lefteyeROI[i + 2][j - 1] * GuassianFilter[4][1] + lefteyeROI[i + 2][j] * GuassianFilter[4][2] + lefteyeROI[i + 2][j + 1] * GuassianFilter[4][3] + lefteyeROI[i + 2][j + 2] * GuassianFilter[4][4];
				leftweight[i][j] = 255 - leftweight[i][j];
				rightweight[i][j] =
					  righteyeROI[i - 2][j - 2] * GuassianFilter[0][0] + righteyeROI[i - 2][j - 1] * GuassianFilter[0][1] + righteyeROI[i - 2][j] * GuassianFilter[0][2] + righteyeROI[i - 2][j + 1] * GuassianFilter[0][3] + righteyeROI[i - 2][j + 2] * GuassianFilter[0][4]
					+ righteyeROI[i - 1][j - 2] * GuassianFilter[1][0] + righteyeROI[i - 1][j - 1] * GuassianFilter[1][1] + righteyeROI[i - 1][j] * GuassianFilter[1][2] + righteyeROI[i - 1][j + 1] * GuassianFilter[1][3] + righteyeROI[i - 1][j + 2] * GuassianFilter[1][4]
					+ righteyeROI[i][j - 2] * GuassianFilter[2][0] + righteyeROI[i][j - 1] * GuassianFilter[2][1] + righteyeROI[i][j] * GuassianFilter[2][2] + righteyeROI[i][j + 1] * GuassianFilter[2][3] + righteyeROI[i][j + 2] * GuassianFilter[2][4]
					+ righteyeROI[i + 1][j - 2] * GuassianFilter[3][0] + righteyeROI[i + 1][j - 1] * GuassianFilter[3][1] + righteyeROI[i + 1][j] * GuassianFilter[3][2] + righteyeROI[i + 1][j + 1] * GuassianFilter[3][3] + righteyeROI[i + 1][j + 2] * GuassianFilter[3][4]
					+ righteyeROI[i + 2][j - 2] * GuassianFilter[4][0] + righteyeROI[i + 2][j - 1] * GuassianFilter[4][1] + righteyeROI[i + 2][j] * GuassianFilter[4][2] + righteyeROI[i + 2][j + 1] * GuassianFilter[4][3] + righteyeROI[i + 2][j + 2] * GuassianFilter[4][4];
				rightweight[i][j] = 255 - rightweight[i][j];
			}
		}
	}
}

void FindEyeCenter::GetleftCenter()
{
	for (int y = 0; y < eyeROIheight; ++y)
	{
		for (int x = 0; x < eyeROIwidth; ++x) 
		{
			double gX = leftgradientX[y][x], gY = leftgradientY[y][x];
			if (gX == 0.0 && gY == 0.0) {
				continue;
			}
			lefttestPossibleCentersFormula(x, y, gX, gY);
		}
	}
}

void FindEyeCenter::lefttestPossibleCentersFormula(int x, int y, double gx, double gy)
{
	for (int cy = 0; cy < eyeROIheight; ++cy) 
	{
		for (int cx = 0; cx < eyeROIwidth; ++cx) 
		{
			if (x == cx && y == cy) 
			{
				continue;
			}
			// create a vector from the possible center to the gradient origin
			double dx = x - cx;
			double dy = y - cy;
			// normalize d
			double magnitude = sqrt((dx * dx) + (dy * dy));
			dx = dx / magnitude;
			dy = dy / magnitude;
			double dotProduct = dx*gx + dy*gy;
			dotProduct = std::max(0.0, dotProduct);
			// square and multiply by the weight
			if (kEnableWeight) 
			{
				leftoutSum[cy][cx] += dotProduct * dotProduct * (leftweight[cy][cx] / kWeightDivisor);
			}
			else {
				leftoutSum[cy][cx] += dotProduct * dotProduct;
			}
		}
	}
}

void FindEyeCenter::GetrightCenter()
{
	for (int y = 0; y < eyeROIheight; ++y)
	{
		for (int x = 0; x < eyeROIwidth; ++x)
		{
			double gX = rightgradientX[y][x], gY = rightgradientY[y][x];
			if (gX == 0.0 && gY == 0.0) {
				continue;
			}
			righttestPossibleCentersFormula(x, y, gX, gY);
		}
	}
}

void FindEyeCenter::righttestPossibleCentersFormula(int x, int y, double gx, double gy)
{
	for (int cy = 0; cy < eyeROIheight; ++cy)
	{
		for (int cx = 0; cx < eyeROIwidth; ++cx)
		{
			if (x == cx && y == cy)
			{
				continue;
			}
			// create a vector from the possible center to the gradient origin
			double dx = x - cx;
			double dy = y - cy;
			// normalize d
			double magnitude = sqrt((dx * dx) + (dy * dy));
			dx = dx / magnitude;
			dy = dy / magnitude;
			double dotProduct = dx*gx + dy*gy;
			dotProduct = std::max(0.0, dotProduct);
			// square and multiply by the weight
			if (kEnableWeight)
			{
				rightoutSum[cy][cx] += dotProduct * dotProduct * (rightweight[cy][cx] / kWeightDivisor);
			}
			else {
				rightoutSum[cy][cx] += dotProduct * dotProduct;
			}
		}
	}
}

void FindEyeCenter::Convert()
{
	int leftmaxi(0), leftmaxj(0), rightmaxi(0), rightmaxj(0);

	double leftmaxval(0), rightmaxval(0);
	double numGradients = eyeROIheight*eyeROIwidth;
	for (int i = 0; i < eyeROIheight; i++)
	{
		for (int j = 0; j < eyeROIwidth; j++)
		{
			leftout[i][j] = leftoutSum[i][j] / numGradients;
			rightout[i][j] = rightoutSum[i][j] / numGradients;
			if (leftout[i][j]>=leftmaxval)
			{
				leftmaxval = leftout[i][j];
				leftmaxi = i;
				leftmaxj = j;
			}
			if (rightout[i][j] >= rightmaxval)
			{
				rightmaxval = rightout[i][j];
				rightmaxi = i;
				rightmaxj = j;
			}
		}
	}
	leftfloodThresh = leftmaxval * kPostProcessThreshold;
	rightfloodThresh = rightmaxval * kPostProcessThreshold;

	for (int i = 0; i < eyeROIheight; i++)
	{
		for (int j = 0; j < eyeROIwidth; j++)
		{
			leftfloodClone[i][j] = leftout[i][j];
			rightfloodClone[i][j] = rightout[i][j];
			if (leftfloodClone[i][j] <leftfloodThresh)
			{
				leftfloodClone[i][j] = 0;
			}
			if (rightfloodClone[i][j] <rightfloodThresh)
			{
				rightfloodClone[i][j] = 0;
			}
		}
	}
}

void FindEyeCenter::leftFloodKillEdges()
{
	for (int i = 0; i < eyeROIheight; i++)
	{
		for (int j = 0; j < eyeROIwidth; j++)
		{
			leftmask[i][j] = 255;
			rightmask[i][j] = 255;
		}
	}
	std::queue<int> toDox;
	std::queue<int> toDoy;
	toDox.push(0);
	toDoy.push(0);
	while (!toDox.empty() && !toDoy.empty())
	{
		int ci = toDox.front();
		int cj = toDoy.front();
		toDox.pop();
		toDoy.pop();
		if (leftfloodClone[ci][cj] == 0.0)
		{
			continue;
		}
		int ni = ci + 1;
		int nj = cj;
		if (shouldpush(ni, nj))
		{
			toDox.push(ni);
			toDoy.push(nj);
		}
		ni = ci - 1; nj = cj;
		if (shouldpush(ni, nj))
		{
			toDox.push(ni);
			toDoy.push(nj);
		}
		ni = ci; nj = cj + 1;
		if (shouldpush(ni, nj))
		{
			toDox.push(ni);
			toDoy.push(nj);
		}
		ni = ci; nj = cj - 1;
		if (shouldpush(ni, nj))
		{
			toDox.push(ni);
			toDoy.push(nj);
		}
		leftfloodClone[ci][cj] = 0.0;
		leftmask[ci][cj] = 0;
	}
}

void FindEyeCenter::rightFloodKillEdges()
{
	std::queue<int> toDox;
	std::queue<int> toDoy;
	toDox.push(0);
	toDoy.push(0);
	while (!toDox.empty() && !toDoy.empty())
	{
		int ci = toDox.front();
		int cj = toDoy.front();
		toDox.pop();
		toDoy.pop();
		if (rightfloodClone[ci][cj] == 0.0)
		{
			continue;
		}
		int ni = ci + 1;
		int nj = cj;
		if (shouldpush(ni, nj))
		{
			toDox.push(ni);
			toDoy.push(nj);
		}
		ni = ci - 1; nj = cj;
		if (shouldpush(ni, nj))
		{
			toDox.push(ni);
			toDoy.push(nj);
		}
		ni = ci; nj = cj + 1;
		if (shouldpush(ni, nj))
		{
			toDox.push(ni);
			toDoy.push(nj);
		}
		ni = ci; nj = cj - 1;
		if (shouldpush(ni, nj))
		{
			toDox.push(ni);
			toDoy.push(nj);
		}
		rightfloodClone[ci][cj] = 0.0;
		rightmask[ci][cj] = 0;
	}
}

void FindEyeCenter::leftmaxCenter()
{
	int max = 0;
	int maxi = 0;
	int maxj = 0;
	for (int i = 0; i < eyeROIheight; i++)
	{
		for (int j = 0; j < eyeROIwidth;j++)
		if (leftout[i][j]>max&&leftmask[i][j]>0)
		{
			leftmaxPx = j;
			leftmaxPy = i;
			max = leftout[i][j];
		}
	}
}

void FindEyeCenter::rightmaxCenter()
{
	int max = 0;
	int maxi = 0;
	int maxj = 0;
	for (int i = 0; i < eyeROIheight; i++)
	{
		for (int j = 0; j < eyeROIwidth; j++)
		if (rightout[i][j]>max&&rightmask[i][j]>0)
		{
			rightmaxPx = j;
			rightmaxPy = i;
			max = rightout[i][j];
		}
	}
}

void FindEyeCenter::unscaleleftcenter()
{
	float ratio = (((float)kFastEyeWidth) / eye_region_width);
	leftmaxPx = round(leftmaxPx / ratio);
	leftmaxPy = round(leftmaxPy / ratio);
}

void FindEyeCenter::unscalerightcenter()
{
	float ratio = (((float)kFastEyeWidth) / eye_region_width);
	rightmaxPx = round(rightmaxPx / ratio);
	rightmaxPy = round(rightmaxPy / ratio);
}