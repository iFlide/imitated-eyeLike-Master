#include <iostream>
#include <math.h>
#include <opencv2\opencv.hpp>
#include "FindEyeCenter.h"
#include "facechar.h"

using namespace std;
using namespace cv;

int main()
{
	int lx, ly, rx, ry;
	facewidth = 216;
	faceheight = 216;
	int countc = 0;
	for (int i = 0; i < faceheight; i++)
	{
		for (int j = 0; j < facewidth; j++)
			faceROI[i][j] = facechar[countc++];
	}
	FindEyeCenter find;
	find.run(lx,ly,rx,ry);
	cout << "(" << lx << "," << ly << ")" << endl;
	cout << "(" << rx << "," << ry << ")" << endl;
	system("pause");
	return 0;
}