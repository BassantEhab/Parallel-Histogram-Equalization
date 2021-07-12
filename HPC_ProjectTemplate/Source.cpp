#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#include "Source.h"
#pragma once

#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>
using namespace std;
using namespace msclr::interop;
const int histogramLength = 256;

int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
	int* input;
	int OriginalImageWidth, OriginalImageHeight;
	//*********************************************************Read Image and save it to local arrayss*************************	
	//Read Image and save it to local arrayss
	System::Drawing::Bitmap BM(imagePath);
	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	*w = BM.Width;
	*h = BM.Height;
	int* Red = new int[BM.Height * BM.Width];
	int* Green = new int[BM.Height * BM.Width];
	int* Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height * BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			input[i * BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values
		}
	}
	return input;
}

int* inputHistogram(int s, int* image)
{
	int* histogram = new int[histogramLength];
	for (int i = 0; i < histogramLength; i++)
		histogram[i] = 0;
	for (int i = 0; i < s; i++) {
		histogram[image[i]] += 1;
	}
	return histogram;
}

float* probabilityHistogram(int* histogram, int s, int r)
{
	float* prob = new float[s];
	for (int i = 0; i < s; i++)
		prob[i] = (float)histogram[i] / r;
	return prob;
}

float* cumulativeHistogram(float* prob, int s)
{
	float* cumulative = new float[s];
	cumulative[0] = prob[0];
	for (int i = 1; i < s; i++)
		cumulative[i] = prob[i] + cumulative[i - 1];
	return cumulative;
}

int* cdfHistogram(float* cumulative, int s)
{
	int* cdf = new int[s];
	for (int i = 0; i < s; i++)
		cdf[i] = roundf(cumulative[i] * 255);
	return cdf;
}

int* outputImage(int s, int* cdf, int* image)
{
	int* newImg;
	newImg = new int[s];
	for (int i = 0; i < s; i++)
		newImg[i] = cdf[image[i]];
	return newImg;
}

void createImage(int* image, int width, int height, int index)
{
	System::Drawing::Bitmap MyNewImage(width, height);


	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			//i * OriginalImageWidth + j
			if (image[i * width + j] < 0)
			{
				image[i * width + j] = 0;
			}
			if (image[i * width + j] > 255)
			{
				image[i * width + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	MyNewImage.Save("..//Data//Output//Tslem" + ".png");
	cout << "result Image Saved " << index << endl;
}


int main()
{
	int ImageWidth = 4, ImageHeight = 4;
	int start_s, stop_s, TotalTime = 0, s;
	int* imageData, * freq, * cdf, * cdfFinal, * subImg, * newImg, * subNewImg, *finalfreq;
	float* prob, * cumProb;
	MPI_Init(NULL, NULL);
	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	imageData = new int[ImageWidth * ImageHeight];
	if (rank == 0)
	{
		System::String^ imagePath;
		std::string img;
		img = "..//Data//Input//test.png";
		imagePath = marshal_as<System::String^>(img);
		imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);
		start_s = clock();
	}

	MPI_Bcast(&ImageWidth, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ImageHeight, 1, MPI_INT, 0, MPI_COMM_WORLD);
	s = ImageWidth * ImageHeight;
	subImg = new int[s / size];
	MPI_Scatter(imageData, s / size, MPI_INT, subImg, s / size, MPI_INT, 0, MPI_COMM_WORLD);

	freq = inputHistogram(s/size, subImg);
	finalfreq = new int[histogramLength];
	MPI_Reduce(freq, finalfreq, histogramLength, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	cdf = new int[histogramLength];
	if (rank == 0)
	{
		prob = probabilityHistogram(finalfreq, histogramLength, s);
		cumProb = cumulativeHistogram(prob, histogramLength);
		cdf = cdfHistogram(cumProb, histogramLength);
	}

	newImg = new int[s];
	MPI_Bcast(cdf, histogramLength, MPI_INT, 0, MPI_COMM_WORLD);
	subNewImg = new int[s / size];
	subNewImg = outputImage(s / size, cdf, subImg);
	MPI_Gather(subNewImg, s/size, MPI_INT, newImg, s/size, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		//newImg = outputImage(s, cdf, imageData);
		createImage(newImg, ImageWidth, ImageHeight, 0);
		stop_s = clock();
		TotalTime += (int)((stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000);
		cout << "time: " << TotalTime << endl;
		free(imageData);
	}
	MPI_Finalize();

	return 0;

}
