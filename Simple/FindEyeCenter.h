#pragma once

const bool kPlotVectorField = false;

// Size constants
const int kEyePercentTop = 25;
const int kEyePercentSide = 13;
const int kEyePercentHeight = 30;
const int kEyePercentWidth = 35;

// Preprocessing
const bool kSmoothFaceImage = false;
const float kSmoothFaceFactor = 0.005;

// Algorithm Parameters
const int kFastEyeWidth = 50;
const int kWeightBlurSize = 5;
const bool kEnableWeight = true;
const float kWeightDivisor = 1.0;
const double kGradientThreshold = 50.0;

// Postprocessing
const bool kEnablePostProcess = true;
const float kPostProcessThreshold = 0.97;

// Eye Corner
const bool kEnableEyeCorner = false;

const double GuassianFilter[5][5] = {
	0.0050, 0.0173, 0.0262, 0.0173, 0.0050,
	0.0173, 0.0598, 0.0903, 0.0598, 0.0173,
	0.0262, 0.0903, 0.1366, 0.0903, 0.0262,
	0.0173, 0.0598, 0.0903, 0.0598, 0.0173,
	0.0050, 0.0173, 0.0262, 0.0173, 0.0050
};

extern int result_left_x;
extern int result_left_y;
extern int result_right_x;
extern int result_right_y;

extern int leftmaxPx;            //center in left eye
extern int leftmaxPy;

extern int rightmaxPx;           //center in right eye
extern int rightmaxPy;

extern int facewidth;            //input by user
extern int faceheight;		   //input by user

extern int eye_region_width;     //estimated by user
extern int eye_region_height;	   //estimated by user
extern int eye_region_top;       //estimated by user
extern int eye_region_leftside;  //estimated by user
extern int eye_region_rightside; //estimated by user

extern int eyeROIwidth;          //width after being scaled
extern int eyeROIheight;         //height after being scaled

extern double leftgradientThresh;     //gradient thresh
extern double rightgradientThresh;

extern double leftfloodThresh;        //image thresh
extern double rightfloodThresh;

extern int faceROI[400][400];         //input face

extern int leftEyeRegion[300][300] ;  //estimated left eye
extern int rightEyeRegion[300][300] ;

extern double lefteyeROI[100][100] ;  //left eye after scaled
extern double righteyeROI[100][100] ;

extern double leftgradientX[100][100] ;//left eye's gradient
extern double leftgradientY[100][100] ;

extern double rightgradientX[100][100] ;
extern double rightgradientY[100][100] ;

extern double leftmags[100][100] ;     //left eye's magnititude
extern double rightmags[100][100] ;

extern double leftweight[100][100] ;   //weight of each pixel
extern double rightweight[100][100] ;

extern double leftoutSum[100][100] ;   //outcome after weight*mag
extern double rightoutSum[100][100] ;

extern double leftout[100][100] ;      //change outsum's scale
extern double rightout[100][100];

extern double leftfloodClone[100][100] ;  //floodclone is out after threshing
extern double rightfloodClone[100][100] ;

extern double leftmask[100][100] ;        //area should be computed in leftout
extern double rightmask[100][100] ;

class FindEyeCenter
{
public:

	FindEyeCenter();
	~FindEyeCenter(); 

	void run(int &lx, int &ly, int &rx, int &ry);                                      //run it, conjunction of the other functions
	int NearDoubleNumber(double D);                  //get the nearby number
	bool shouldpush(int i, int j);					 //judge whether a point is  illegal or not
	void ComputeEye();                               //estimate eye based on face
	void resize(double width, double height);        //get eyeROI based on eyeregion
	void ComputeXGradient();                         //get gradient
	void ComputeYGradient();
	void MatrixMagnitude();                          //get mag
	void computeDynamicThreshold();                  //get mag's thresh
	void Thresh();                                   //do thresh
	//sigma = 0.3*((ksize-1)*0.5 - 1) + 0.8=1.1
	void ComputeWeight(int rows, int cols);          //get weight
	void GetleftCenter();                            //compute di*gi
	void GetrightCenter();
	void lefttestPossibleCentersFormula(int x, int y, double gx, double gy);  //test
	void righttestPossibleCentersFormula(int x, int y, double gx, double gy);
	void Convert();                                   //convert
	void leftFloodKillEdges(); 
	void rightFloodKillEdges();
	void leftmaxCenter();                             //get max point based on function max
	void rightmaxCenter();
	void unscaleleftcenter();                         //unscale
	void unscalerightcenter();
};

