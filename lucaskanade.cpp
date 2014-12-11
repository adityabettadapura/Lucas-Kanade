// Homework.cpp : Defines the entry point for the console application.
//

#include <afxwin.h>  // necessary for MFC to work properly
#include "../../src/blepo.h"
#include "math.h"
#include <stack>
#include <string>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#define SQR(a) ((a)*(a))
using namespace blepo;

void GetGradient(ImgGray img, int width, int halfwidth, ImgFloat *outx, ImgFloat *outy, double dervkern[], double kern[])
{
	double sum = 0;
	int x, y, i, j;
	double *temp;
	temp = (double *)calloc(img.Height()*img.Width(), sizeof(double));
	for (y = 0; y < img.Height(); y++)
	{
		for (x = halfwidth; x < img.Width() - halfwidth - 1; x++)
		{
			sum = 0;
			for (i = 0; i < width; i++)
			{
				sum += ((img(x + halfwidth - i, y))*(dervkern[width -i -1]));
			}
			temp[y*img.Width() + x] = sum;
		}
	}
	for (y = halfwidth; y < img.Height() - halfwidth - 1; y++)
	{
		for (x = 0; x < img.Width(); x++)
		{
			sum = 0;
			for (i = 0; i < width; i++)
			{
				sum += ((temp[(y + halfwidth - i)* img.Width() + x])*(kern[width -i -1]));
			}
			(*outx)(x, y) = sum;
		}
	}

	for (y = 0; y < img.Height(); y++)
	{
		for (x = halfwidth; x < img.Width() - halfwidth - 1; x++)
		{
			sum = 0;
			for (i = 0; i < width; i++)
			{
				sum += ((img(x + halfwidth - i, y))*(kern[width -i -1]));
			}
			temp[y*img.Width() + x] = sum;
		}
	}
	for (y = halfwidth; y < img.Height() - halfwidth - 1; y++)
	{
		for (x = 0; x < img.Width(); x++)
		{
			sum = 0;
			for (i = 0; i < width; i++)
			{
				sum += ((temp[(y + halfwidth - i)*img.Width() + x])*(dervkern[width -i -1]));
			}
			(*outy)(x, y) = sum;
		}
	}

}

void CreateGaussKern(double width, double halfwidth, double gauss[], double gaussDeriv[], float sigma)
{
	double sumGauss = 0;
	double sumGaussDeriv = 0;

	for (int i = 0; i < width; i++)
	{
		gauss[i] = exp((-(i - halfwidth)*(i - halfwidth)) / (2 * sigma*sigma));
		sumGauss += gauss[i];

		gaussDeriv[i] = (i - halfwidth)*gauss[i];
		sumGaussDeriv += gaussDeriv[i] * i;
	}

	for (int j = 0; j < width; j++)
	{
		gauss[j] = gauss[j] / sumGauss;
		gaussDeriv[j] = gaussDeriv[j] / sumGaussDeriv;
	}
}

double InterpolateBilinear(ImgGray frame, double x, double y)
{
	double x0, y0, ax, ay;
	x0 = floor(x);
	y0 = floor(y);
	ax = x - x0;
	ay = y - y0;

	if (x0 < 0) x0 = 0;
	if (y0 < 0) y0 = 0;
	if (x0 >= frame.Width() - 1) x0 = frame.Width() - 2;
	if (y0 >= frame.Height() - 1) y0 = frame.Height() - 2;

	return ((1 - ax)*(1 - ay)*frame(x0, y0) + ax*(1 - ay)*frame(x0 + 1, y0) + (1 - ax)*ay*frame(x0, y0 + 1) + ax*ay*frame(x0 + 1, y0 + 1));
}

double InterpolateBilinear(ImgFloat frame, double x, double y)
{
	double x0, y0, ax, ay;
	x0 = floor(x);
	y0 = floor(y);
	ax = x - x0;
	ay = y - y0;

	if (x0 < 0) x0 = 0;
	if (y0 < 0) y0 = 0;
	if (x0 >= frame.Width() -1) x0 = frame.Width()-2;
	if (y0 >= frame.Height() - 1) y0 = frame.Height() - 2;

	return ((1 - ax)*(1 - ay)*frame(x0, y0) + ax*(1 - ay)*frame(x0 + 1, y0) + (1 - ax)*ay*frame(x0, y0 + 1) + ax*ay*frame(x0 + 1, y0 + 1));
}

void Compute2x2GradMat(int x, int y, double z[], ImgFloat grad_x, ImgFloat grad_y, int size)
{
	for (int j = -(size-1)/2; j < (size+1)/2; j++)
	{
		for (int i = -(size-1)/2; i < (size + 1) / 2; i++)
		{
				z[0] = z[0] + InterpolateBilinear(grad_x, x + i, y + j) * InterpolateBilinear(grad_x, x + i, y + j);	//gxx
				z[1] = z[1] + InterpolateBilinear(grad_x, x + i, y + j) * InterpolateBilinear(grad_y, x + i, y + j);	//gxy
				z[2] = z[2] + InterpolateBilinear(grad_y, x + i, y + j) * InterpolateBilinear(grad_y, x + i, y + j);	//gyy
		}
	}
}

void ComputeZGradMat(int x, int y, double z[], ImgFloat grad_x, ImgFloat grad_y)
{
	for (int j = -1; j < 2; j++)
	{
		for (int i = -1; i < 2; i++)
		{
			z[0] = z[0] + grad_x(x + i, y + j) * grad_x(x + i, y + j);	//gxx
			z[1] = z[1] + grad_x(x + i, y + j) * grad_y(x + i, y + j);	//gxy
			z[2] = z[2] + grad_y(x + i, y + j) * grad_y(x + i, y + j);	//gyy
		}
	}
}

void Compute2x1ErrorVector(double x, double y, double e[], double u[], ImgFloat grad_x, ImgFloat grad_y, ImgGray currentframe, ImgGray nextframe, int size)
{
	for (int j = -(size - 1) / 2; j < (size + 1) / 2; j++)
	{
		for (int i = -(size - 1) / 2; i < (size + 1) / 2; i++)
		{
			e[0] = e[0] + InterpolateBilinear(grad_x,x + i, y + j)*(InterpolateBilinear(currentframe, x + i, y + j) - InterpolateBilinear(nextframe, x + i + u[0], y + j + u[1]));
			e[1] = e[1] + InterpolateBilinear(grad_y,x + i, y + j)*(InterpolateBilinear(currentframe, x + i, y + j) - InterpolateBilinear(nextframe, x + i + u[0], y + j + u[1]));
		}
	}
}

void DrawSquare(int x, int y, ImgBgr *input)
{
	if((x >0) && (y>0) && (x < (*input).Width()-1) && (y < (*input).Height()-1))
	{
		for (int j = -1; j < 2; j++)
		{
			for (int i = -1; i < 2; i++)
			{
				if ((x+i == x) && (y+j == y))
					continue;
				else
					(*input)(x + i, y + j) = Bgr(0, 0, 255);
			}
		}
	}
}

int main(int argc, const char* argv[], const char* envp[])
{
	// Initialize MFC and return if failure
	HMODULE hModule = ::GetModuleHandle(NULL);
	if (hModule == NULL || !AfxWinInit(hModule, NULL, ::GetCommandLine(), 0))
	{
		printf("Fatal Error: MFC initialization failed (hModule = %x)\n", hModule);
		return 1;
	}

	try
	{
		//Input Section
		CString path = "../../images/";
		CString filename;
		int firstframe, lastframe, winsize;
		float sigma;
		
		const char* format = argv[1];
		firstframe = atoi(argv[2]);
		lastframe = atoi(argv[3]);
		sigma = atof(argv[4]);
		winsize = atoi(argv[5]);
		
		filename.Format(format, firstframe);
		filename = path + filename;

		printf("filename = %s\n", filename);

		ImgGray firstimg;
		Load(filename, &firstimg);

		//Calculate gradient
		float tempSigma = 0.8;
		printf("%f\n", tempSigma);
		int mu = round(2.5*tempSigma - 0.5);
		int w = 2 * mu + 1;
		double *gauss;
		double *gaussDeriv;
		gauss = (double*)calloc(w, sizeof(double));
		gaussDeriv = (double*)calloc(w, sizeof(double));

		//create gaussian kernel
		CreateGaussKern(w, mu, gauss, gaussDeriv, tempSigma);
		
		ImgFloat grad_x, grad_y;
		grad_x.Reset(firstimg.Width(), firstimg.Height());
		grad_y.Reset(firstimg.Width(), firstimg.Height());
		Set(&grad_x, 0);
		Set(&grad_y, 0);

		// Calculate Gradient
		GetGradient(firstimg, w, mu, &grad_x, &grad_y, gaussDeriv, gauss);

		ImgFloat corner;
		corner.Reset(firstimg.Width(), firstimg.Height());
		Set(&corner, 0);

		int threshold = 2000;
		if (firstframe >= 30) threshold = 4000; //for flowergarden
		if (firstframe >= 588) threshold = 2000; // for state_seq
		printf("%d\n", threshold);

		//find cornerness
		for (int y = 1; y < firstimg.Height()-1; y++)
		{
			for (int x = 1; x < firstimg.Width()-1; x++)
			{
				double z[3] = { 0, 0, 0 };
				double lambda1 = 0;
				double lambda2 = 0;
				double cornerness = 0;

				ComputeZGradMat(x, y, z, grad_x, grad_y);
				lambda1 = 0.5*((z[0] + z[2]) + sqrt(SQR(z[0] - z[2]) + 4 * z[1] * z[1]));
				lambda2 = 0.5*((z[0] + z[2]) - sqrt(SQR(z[0] - z[2]) + 4 * z[1] * z[1]));
				cornerness = min(lambda1, lambda2);
				corner(x, y) = ((lambda2 > threshold) ? lambda2 : 0);
			}
		}

		int featurecount = 0;
		ImgBgr goodfeature;
		Convert(firstimg, &goodfeature);

		for (int y = winsize/2; y < firstimg.Height() - (winsize/2); y++)
		{
			for (int x = winsize/2; x < firstimg.Width() - (winsize/2); x++)
			{
				if (corner(x, y) < corner(x - 1, y) || corner(x, y) < corner(x + 1, y) ||
					corner(x, y) < corner(x, y - 1) || corner(x, y) < corner(x, y + 1))
				{
					corner(x, y) = 0;
				}
				if (corner(x, y) > threshold)
				{
					DrawSquare(x,y,&goodfeature);
					featurecount++;
				}
			}
		}
		
		//store coordinates of good feature points
		std::vector<double> feature_x(featurecount + 1);
		std::vector<double> feature_y(featurecount + 1);

		int count = 0;
		for (int y = winsize/2; y < firstimg.Height() - (winsize/2); y++)
		{
			for (int x = winsize/2; x < firstimg.Width() - (winsize/2); x++)
			{
				if (corner(x, y) > threshold)
				{
					feature_x[count] = x;
					feature_y[count] = y;
					count++;
				}
			}
		}

		for (int i = 0; i < featurecount; i++)
		{
			printf("feature_x[%d] = %lf\tfeature_[%d] = %lf\n", i, feature_x[i], i, feature_y[i]);
		}

		Figure fig1;
		fig1.SetTitle("Image with Good feature points");
		fig1.Draw(goodfeature);

		//LucasKanade algorithm
		//Calculate gaussian kernel
		mu = round(2.5*sigma - 0.5);
		printf("mu = %d\n",mu);
		w = 2 * mu + 1;
		double *gaussKern;
		double *gaussDerivKern;

		gaussKern = (double*)calloc(w, sizeof(double));
		gaussDerivKern = (double*)calloc(w, sizeof(double));
		CreateGaussKern(w, mu, gaussKern, gaussDerivKern, sigma);

		//Initialize images
		ImgFloat gradient_x, gradient_y;
		gradient_x.Reset(firstimg.Width(), firstimg.Height());
		gradient_y.Reset(firstimg.Width(), firstimg.Height());
		Set(&gradient_x, 0);
		Set(&gradient_y, 0);

		ImgGray currentimg, nextimg;
		currentimg.Reset(firstimg.Width(), firstimg.Height());
		nextimg.Reset(firstimg.Width(), firstimg.Height());
		Set(&currentimg, 0);
		Set(&nextimg, 0);

		CString file;
		ImgBgr finalfeature;
		finalfeature.Reset(firstimg.Width(), firstimg.Height());
		Figure fig2;
		fig2.SetTitle("Final image with Feature Tracking");

		//for each pair of image
		for (int k = firstframe; k < lastframe; k++)
		{
			//Take first frame
			file.Format(format, k);
			file = path + file;
			Load(file, &currentimg);

			//Take second frame
			file.Format(format, k+1);
			file = path + file;
			Load(file, &nextimg);

			//Compute gradient of first frame
			GetGradient(currentimg, w, mu, &gradient_x, &gradient_y, gaussDerivKern, gaussKern);

			for (int i = 0; i < featurecount; i++)
			{
				//LK algorithm
				//initialize Z, U vectors
				double zmat[3] = { 0, 0, 0 };
				double umat[2] = { 0, 0 };
				double emat[2] = { 0, 0 };
				int iteration = 0;

				////Compute gradient of first frame
				//GetGradient(currentimg, w, mu, &gradient_x, &gradient_y, gaussDerivKern, gaussKern);
				if ((feature_x[i] < 1) || (feature_y[i] < 1))
					continue;
				else
				{
					Compute2x2GradMat(feature_x[i], feature_y[i], zmat, gradient_x, gradient_y, winsize);
				}
				
				while (iteration < 1)
				{
					double det = 0;
					double u_delta = 0;
					double v_delta = 0;

					Compute2x1ErrorVector(feature_x[i], feature_y[i], emat, umat, gradient_x, gradient_y, currentimg, nextimg, winsize);

					det = zmat[0] * zmat[2] - SQR(zmat[1]);
					if (det == 0) det = 1;
	
					u_delta = ((zmat[2] * emat[0]) - (zmat[1] * emat[1]))/det;
					v_delta = ((zmat[0] * emat[1]) - (zmat[1] * emat[0]))/det;

					umat[0] += u_delta;
					umat[1] += v_delta;

					iteration++;
				}
				feature_x[i] += umat[0];
				feature_y[i] += umat[1];
			}

			Convert(nextimg, &finalfeature);
			for (int i = 0; i < featurecount; i++)
			{
				DrawSquare(round(feature_x[i]), round(feature_y[i]), &finalfeature);
			}
			fig2.Draw(finalfeature);
		}
		EventLoop();
	}
	catch (const Exception& e)
	{
		e.Display();    // display exception to user in a popup window 
	}
	return 0;
}
