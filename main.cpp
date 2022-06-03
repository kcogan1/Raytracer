#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <stdio.h>
#include <errno.h>
#include <vector>
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <cfloat>
#define MAX DBL_MAX

using namespace std;

//STRUCTS
struct viewPoint {
	Eigen::Vector3d from;
	Eigen::Vector3d at;
	Eigen::Vector3d up;
	float angle;
	float hither;
	int resolution[2];
};

struct color {
	double fillColor[3];
	double Kd;
	double Ks;
	double shine;
	double T;
	double refraction;
};

struct polygon {
	Eigen::Matrix3d points;
	color colorInfo;
	bool shadowed = false;
};

struct ray {
	Eigen::Vector3d origin;
	Eigen::Vector3d direction;
	bool hit = false;
	polygon hitList[20];
	Eigen::Vector3d IntersectionLocation;
	int indexOfFirstHit = 0;
	double minT = 0;
	int bounces = 0;
};


//GLOBAL VARS
float bgColor[3];
viewPoint camera;
vector <polygon> worldShapes;
vector <color> colorVector;
Eigen::Vector3d light1;
Eigen::Vector3d light2;
int hitIterator = 0;
float refraction;
int recursionLimit = 5;
int recursionCounter = 0;
double totalColor = 0;

//FUNCTION DEFINITIONS
void parseFile();
bool intersections(ray*, polygon*, double, double);
double det(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&);
double shade(ray*, polygon*, int);
double reflect(int, ray*, polygon*);

void parseFile()
{
	//indicates what info is being parsed
	char indicator;
	char tempIndicator;
	//Char of arbitrary length to store the line read in
	char input[200];

	//temp vars to store the info being parsed
	char tempStr[20];
	int tempInt1;
	int tempInt2;
	float temp1;
	float temp2;
	float temp3;
	float temp4;
	float temp5;
	float temp6;
	float temp7;
	float temp8;

	int colorCounter = 0;

	//opens file
	FILE * file;
	file = fopen("teapot-3.nff", "r");
	//file = fopen("tetra-3.nff", "r");

	//if file is empty or does not exist
	if (file == NULL)
	{
		cout << "File is NULL" << endl;
		cout << errno << endl;
	}

	else
	{
		//gets first character of the line
		indicator = fgetc(file);

		//loops through file until end of file
		while (tempIndicator != EOF)
		{
			//BACKGROUND COLOR
			if (indicator == 'b')
			{
				//reads and parses line
				fgets(input, 50, file);
				sscanf(input, "%g %g %g", &temp1, &temp2, &temp3);

				//stores values
				bgColor[0] = temp1;
				bgColor[1] = temp2;
				bgColor[2] = temp3;

				//cout << "bgColor:" << temp1 << " " << temp2 << " " << temp3 << endl;
			}

			//VIEWPOINT
			else if (indicator == 'v')
			{
				//FROM
				fgets(input, 50, file);
				sscanf(input, "%s %g %g %g", tempStr, &temp1, &temp2, &temp3);
				camera.from(0) = temp1;
				camera.from(1) = temp2;
				camera.from(2) = temp3;
				cout << "From: " << camera.from[0] << " " << camera.from[1] << " " << camera.from[2] << endl;

				//AT
				fgets(input, 50, file);
				sscanf(input, "%s %g %g %g", tempStr, &temp1, &temp2, &temp3);
				camera.at(0) = temp1;
				camera.at(1) = temp2;
				camera.at(2) = temp3;
				cout << "at: " << camera.at[0] << " " << camera.at[1] << " " << camera.at[2] << endl;

				//UP
				fgets(input, 50, file);
				sscanf(input, "%s %g %g %g", tempStr, &temp1, &temp2, &temp3);
				camera.up(0) = temp1;
				camera.up(1) = temp2;
				camera.up(2) = temp3;
				cout << "UP: " << camera.up[0] << " " << camera.up[1] << " " << camera.up[2] << endl;

				//ANGLE
				fgets(input, 50, file);
				sscanf(input, "%s %g", tempStr, &temp1);
				camera.angle = temp1;
				//cout << "Angle: " << camera.angle << endl;

				//HITHER
				fgets(input, 50, file);
				sscanf(input, "%s %g", tempStr, &temp1);
				camera.hither = temp1;
				//cout << "Hither: " << camera.hither << endl;

				//RESOLUTION
				fgets(input, 50, file);
				sscanf(input, "%s %d %d", tempStr, &tempInt1, &tempInt2);
				camera.resolution[0] = tempInt1;
				camera.resolution[1] = tempInt2;
				//cout << "resolution: " << camera.resolution[0] << " " << camera.resolution[1] << endl;
			}

			//LIGHT
			else if (indicator == 'l')
			{
				fgets(input, 50, file);
				sscanf(input, "%g %g %g", &temp1, &temp2, &temp3);
			
				light1[0] = temp1;
				light1[1] = temp2;
				light1[2] = temp3;

				cout << "light1: " << light1[0] << " " << light1[1] << " " << light1[2] << endl;

				fgets(input, 50, file);
				sscanf(input, "%s %g %g %g", tempStr, &temp1, &temp2, &temp3);

				light2[0] = temp1;
				light2[1] = temp2;
				light2[2] = temp3;
				cout << "light2: " << light2[0] << " " << light2[1] << " " << light2[2] << endl;
			}

			//FILL COLOR
			else if (indicator == 'f')
			{
				color* tempColor = new color;
				//reads and parses the line
				fgets(input, 50, file);
				sscanf(input, "%g %g %g %g %g %g %g %g", &temp1, &temp2, &temp3, &temp4, &temp5, &temp6, &temp7, &temp8);

				//stores values
				tempColor->fillColor[0] = temp1;
				tempColor->fillColor[1] = temp2;
				tempColor->fillColor[2] = temp3;

				tempColor->Kd = temp4;
				tempColor->Ks = temp5;
				tempColor->shine = temp6;
				tempColor->T = temp7;
				tempColor->refraction = temp8;

				cout << "kd: " << tempColor->Kd << endl << "ks: " << tempColor->Ks << endl;
				cout << "Shine: " << tempColor->shine << endl << "T: " << tempColor->T << endl << "refeact: " << tempColor->refraction << endl;

				colorVector.push_back(*tempColor);
				
				colorCounter++;
				//cout << "color counter = " << colorCounter;
				//cout << "Fill Color:" << fillColor[0] << " " << fillColor[1] << " " << fillColor[2] << endl;
			}

			//POLYGON
			else if (indicator == 'p')
			{
				fgets(input, 200, file);
				sscanf(input, "%d", &tempInt1);

				//if the nff give the normals as well
				if (tempIndicator == 'p')
				{
					//if the polygon has 3 verticies
					if (tempInt1 == 3)
					{
						//creates new polygon object
						polygon* tempPoly = new polygon;

						//loops through the 3 verticies and stores the info
						for (int i = 0; i < 3; i++)
						{
							fgets(input, 200, file);
							sscanf(input, "%g %g %g %g %g %g", &temp1, &temp2, &temp3, &temp4, &temp5, &temp6);

							tempPoly->points(i, 0) = temp1;
							tempPoly->points(i, 1) = temp2;
							tempPoly->points(i, 2) = temp3;
						}

						//adds the new polygon to a vector called "shapes"
						tempPoly->colorInfo = colorVector[(colorCounter - 1)];
						//cout << "color counter = " << colorCounter << endl;
						worldShapes.push_back(*tempPoly);
					}

					//if the polygon has 4 verticies
					else if (tempInt1 == 4)
					{
						//creates new polygon object
						polygon* tempPoly1 = new polygon;
						polygon* tempPoly2 = new polygon;

						int counter = 0;
						for (int i = 0; i < 3; i++)
						{
							fgets(input, 200, file);
							sscanf(input, "%g %g %g %g %g %g", &temp1, &temp2, &temp3, &temp4, &temp5, &temp6);
			
							tempPoly1->points(i, 0) = temp1;
							tempPoly1->points(i, 1) = temp2;
							tempPoly1->points(i, 2) = temp3;

							if (i == 0 || i == 2 || i == 3)
							{
								tempPoly2->points(counter, 0) = temp1;
								tempPoly2->points(counter, 1) = temp2;
								tempPoly2->points(counter, 2) = temp3;
								counter++;
							}
						}

						//after the loop, it has to get the fourth vertex becuase the for loop will never hit it.
						fgets(input, 200, file);
						sscanf(input, "%g %g %g %g %g %g", &temp1, &temp2, &temp3, &temp4, &temp5, &temp6);

						tempPoly2->points(2, 0) = temp1;
						tempPoly2->points(2, 1) = temp2;
						tempPoly2->points(2, 2) = temp3;

						//adds the new polygon to a vector called "shapes"
						tempPoly1->colorInfo = colorVector[(colorCounter - 1)];
						tempPoly2->colorInfo = colorVector[(colorCounter - 1)];
						//cout << "color counter = " << colorCounter << endl;
						worldShapes.push_back(*tempPoly1);
						worldShapes.push_back(*tempPoly2);
					}
				}

				//if the polygon has 3 verticies
				else if (tempInt1 == 3)
				{
						//creates new polygon object
						polygon* tempPoly = new polygon;

						//loops through the 3 verticies and stores the info
						for (int i = 0; i < 3; i++)
						{
							fgets(input, 200, file);
							sscanf(input, "%g %g %g", &temp1, &temp2, &temp3);

							tempPoly->points(i, 0) = temp1;
							tempPoly->points(i, 1) = temp2;
							tempPoly->points(i, 2) = temp3;
						}
						//adds the new polygon to a vector called "shapes"
						tempPoly->colorInfo = colorVector[(colorCounter - 1)];
						//cout << "color counter = " << colorCounter << endl;
						worldShapes.push_back(*tempPoly);
				}

				//if the polygon has 4 verticies
				else if (tempInt1 == 4)
				{
						//creates new polygon object
						polygon* tempPoly1 = new polygon;
						polygon* tempPoly2 = new polygon;
					
						//loops through the 3 verticies and stores the info
						int counter = 0;
						for (int i = 0; i < 3; i++)
						{
							fgets(input, 200, file);
							sscanf(input, "%g %g %g", &temp1, &temp2, &temp3);

							tempPoly1->points(i, 0) = temp1;
							tempPoly1->points(i, 1) = temp2;
							tempPoly1->points(i, 2) = temp3;

							if (i == 0 || i == 2 || i == 3)
							{
								tempPoly2->points(counter, 0) = temp1;
								tempPoly2->points(counter, 1) = temp2;
								tempPoly2->points(counter, 2) = temp3;
								counter++;
							}
						}

						//after the loop, it has to get the fourth vertex becuase the for loop will never hit it.
						fgets(input, 200, file);
						sscanf(input, "%g %g %g", &temp1, &temp2, &temp3);

						tempPoly2->points(2, 0) = temp1;
						tempPoly2->points(2, 1) = temp2;
						tempPoly2->points(2, 2) = temp3;

						//adds the new polygon to a vector called "shapes"
						tempPoly1->colorInfo = colorVector[(colorCounter - 1)];
						tempPoly2->colorInfo = colorVector[(colorCounter - 1)];
						//cout << "color counter = " << colorCounter << endl;

						worldShapes.push_back(*tempPoly1);
						worldShapes.push_back(*tempPoly2);
				}			
			}

			//handles lines in the NFF file that do not have an indicator letter
			else
			{
				fgets(input, 50, file);
			}

			//gets new indicator for the next iteration of the loop
			indicator = fgetc(file);
			tempIndicator = fgetc(file);
		}
	}

	if(file !=NULL)
		fclose(file);
}

//calculate determinant
double det(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c) {
	return a[0] * (b[1] * c[2] - c[1] * b[2]) +
		b[0] * (c[1] * a[2] - a[1] * c[2]) +
		c[0] * (a[1] * b[2] - b[1] * a[2]);
}

bool intersections(ray *rayInfo, polygon *polyInfo, double lower, double upper)
{
	Eigen::Vector3d ba; //v0 - v1 
	Eigen::Vector3d ca; //v0 - v2 
	Eigen::Vector3d ea; //v0 - ray.origin

	//ba
	ba(0) = polyInfo->points(0, 0) - polyInfo->points(1, 0);
	ba(1) = polyInfo->points(0, 1) - polyInfo->points(1, 1);
	ba(2) = polyInfo->points(0, 2) - polyInfo->points(1, 2);

	//ca
	ca(0) = polyInfo->points(0, 0) - polyInfo->points(2, 0);
	ca(1) = polyInfo->points(0, 1) - polyInfo->points(2, 1);
	ca(2) = polyInfo->points(0, 2) - polyInfo->points(2, 2);

	//ea
	ea(0) = polyInfo->points(0, 0) - rayInfo->origin(0);
	ea(1) = polyInfo->points(0, 1) - rayInfo->origin(1);
	ea(2) = polyInfo->points(0, 2) - rayInfo->origin(2);

	double detA = det(ba, ca, rayInfo->direction);
	double beta = det(ea, ca, rayInfo->direction) / detA;

	if (beta < 0 || beta > 1) return false;
	double gamma = det(ba, ea, rayInfo->direction) / detA;
	if (gamma < 0.0 || gamma > 1.0 - beta) return false;
	double t = det(ba, ca, ea) / detA;
	if (t < lower || t > upper) return false;


	//if it is a hit, it stores it in the hitList, increments the list iterator, and returns true

	//if this is the first hit, this hit is minT
	rayInfo->hitList[hitIterator] = *polyInfo;

	if (rayInfo->hit == false)
	{
		rayInfo->minT = t;
		rayInfo->indexOfFirstHit = 0;
	}

	//else, if this new t is less than the current mint, make it the new minT
	else 
	{
		if (t < rayInfo->minT)
		{
			rayInfo->minT = t;
			rayInfo->indexOfFirstHit = hitIterator;
		}
	}

	hitIterator++;
	rayInfo->hit = true;
	return true;
}

double shade(ray* rayInfo, polygon* polyInfo, int index)
{
	//calculate surface normal of current polygon
	Eigen::Vector3d U; //v1 - v0
	Eigen::Vector3d V; //v2 - v0
	Eigen::Vector3d normal;

	U(0) = polyInfo->points(1, 0) - polyInfo->points(0, 0);
	U(1) = polyInfo->points(1, 1) - polyInfo->points(0, 1);
	U(2) = polyInfo->points(1, 2) - polyInfo->points(0, 2);

	V(0) = polyInfo->points(2, 0) - polyInfo->points(0, 0);
	V(1) = polyInfo->points(2, 1) - polyInfo->points(0, 1);
	V(2) = polyInfo->points(2, 2) - polyInfo->points(0, 2);

	//normal is a vector
	normal = U.cross(V);
	normal.normalize();
	
	double diffuse;
	double specular;
	double lightIntensity = (1 / sqrt(2));
	Eigen::Vector3d h;
	Eigen::Vector3d l = light1 - rayInfo->IntersectionLocation;
	l.normalize();
	double localColor;

	//do for first light----------------------------------------------------------------------
	h = (rayInfo->direction + l) / (rayInfo->direction + l).norm();
	diffuse = max(0.0, (normal.dot(l)));
	specular = pow(max(0.0, normal.dot(h)), polyInfo->colorInfo.shine);
	localColor = (polyInfo->colorInfo.Kd * polyInfo->colorInfo.fillColor[index] * diffuse + polyInfo->colorInfo.Ks * specular) *lightIntensity;


	//do for second light----------------------------------------------------------------------
	l = light2 - rayInfo->IntersectionLocation;
	l.normalize();
	h = (rayInfo->direction + l) / (rayInfo->direction + l).norm();
	diffuse = max(0.0, (normal.dot(l)));
	specular = pow(max(0.0, normal.dot(h)), polyInfo->colorInfo.shine);

	localColor += (polyInfo->colorInfo.Kd * polyInfo->colorInfo.fillColor[index] * diffuse + polyInfo->colorInfo.Ks * specular) * lightIntensity;
	
	//checks that the color is between 1 and 0
	if (localColor > 1.0)
		localColor = 1;

	if (localColor < 0.0)
		localColor = 0;

	return localColor;
}

double reflect(int localColor, ray* rayInfo, polygon* polyInfo)
{
	//make shadows for EVERY shape
	ray shadowRay;
	shadowRay.origin = rayInfo->IntersectionLocation;
	shadowRay.direction = shadowRay.origin - light1;
	shadowRay.direction.normalize();

	bool intersect;
	bool shadow = false;
	double oldColor;
	hitIterator = 0;
	rayInfo->hit = false;

		//runs this loop to find new minT value, and shadowed areas
		for (int k = 0; k < worldShapes.size(); k++)
		{
			shadow = intersections(&shadowRay, &worldShapes[k], camera.hither, MAX);
			if (shadow == true) {
				worldShapes[k].shadowed = true;
				//cout << "Shadowed" << endl;
			}

			intersect = intersections(rayInfo, &worldShapes[k], camera.hither, MAX);
			if (intersect == true)
			{
				rayInfo->hit = true;
			}
		}

		//if the hit is true, recurse the reflection ray
		if (rayInfo->hit == true)
		{
			rayInfo->IntersectionLocation = rayInfo->origin + (rayInfo->minT * rayInfo->direction);
			
			//calculate surface normal of current polygon
			Eigen::Vector3d U; //v1 - v0
			Eigen::Vector3d V; //v2 - v0
			Eigen::Vector3d normal;

			U(0) = rayInfo->hitList[rayInfo->indexOfFirstHit].points(1, 0) - rayInfo->hitList[rayInfo->indexOfFirstHit].points(0, 0);
			U(1) = rayInfo->hitList[rayInfo->indexOfFirstHit].points(1, 1) - rayInfo->hitList[rayInfo->indexOfFirstHit].points(0, 1);
			U(2) = rayInfo->hitList[rayInfo->indexOfFirstHit].points(1, 2) - rayInfo->hitList[rayInfo->indexOfFirstHit].points(0, 2);

			V(0) = rayInfo->hitList[rayInfo->indexOfFirstHit].points(2, 0) - rayInfo->hitList[rayInfo->indexOfFirstHit].points(0, 0);
			V(1) = rayInfo->hitList[rayInfo->indexOfFirstHit].points(2, 1) - rayInfo->hitList[rayInfo->indexOfFirstHit].points(0, 1);
			V(2) = rayInfo->hitList[rayInfo->indexOfFirstHit].points(2, 2) - rayInfo->hitList[rayInfo->indexOfFirstHit].points(0, 2);

			normal = U.cross(V);
			normal.normalize();

			//computers direction or reflection ray
			Eigen::Vector3d reflectionRay = rayInfo->direction - 2 * rayInfo->direction.dot(normal) * normal;
			reflectionRay.normalize();

			//shades the polygon that the ray just landed on
			if (polyInfo->shadowed == false)
				oldColor = shade(rayInfo, &rayInfo->hitList[rayInfo->indexOfFirstHit], localColor);
			else
				oldColor = rayInfo->hitList[rayInfo->indexOfFirstHit].colorInfo.fillColor[localColor];

			//constructs a new ray to pass into recursive function
			ray newRay;
			newRay.origin = rayInfo->IntersectionLocation;
			newRay.direction = reflectionRay;
			polygon* newPolygon = &rayInfo->hitList[rayInfo->indexOfFirstHit];
			newPolygon->colorInfo.fillColor[localColor] = oldColor;
			newRay.bounces = rayInfo->bounces + 1;

			//recurse
			totalColor = oldColor;
			if (rayInfo->hitList[rayInfo->indexOfFirstHit].colorInfo.Ks > 0 && newRay.bounces < 5)
			{
				//totalColor += polyInfo->colorInfo.Ks * reflect(localColor, &newRay, newPolygon);
				totalColor += rayInfo->hitList[rayInfo->indexOfFirstHit].colorInfo.Ks * reflect(localColor, &newRay, newPolygon);
			}
		}

		//if ray has no hits, return background color
		else if(rayInfo->hit == false)
		{
			//totalColor += bgColor[localColor];
			return bgColor[localColor];
			//local color is rgb index
		}

		if (totalColor > 1.0)
			totalColor = 1;

		else if (totalColor < 0.0)
			totalColor = 0;
		
		return totalColor;		
}

	
int main() 
{
	parseFile();

	//array for pixels
	unsigned char pixels[512][512][3]; //HARD CODED RESOLUTION

	//w
	Eigen::Vector3d w;
	w = camera.from - camera.at;
	w /= w.norm();
	//u
	Eigen::Vector3d u;
	u = (camera.up).cross(w);
	u.normalize();
	//v
	Eigen::Vector3d v;
	v = w.cross(u);

	double d = (camera.from - camera.at).norm();
	double h =  tan((M_PI / 180.0) * (camera.angle / 2.0)) * d;
	double increment = (2 * h) / camera.resolution[0];
	double l = -h + 0.5 * increment;
	double t = h * (((double)camera.resolution[1]) / camera.resolution[0]) - 0.5 * increment;

	bool intersect;
	bool shadow;
	int counter = 0;
	
	//loop through pixels
	for (int i = 0; i < camera.resolution[1]; i++)
	{
		for (int j = 0; j < camera.resolution[0]; j++)
		{
			ray tempRay;
			//for each pixel, compute ray direction
			double x = l + i*increment;
			double y = t - j*increment;

			Eigen::Vector3d dir;
			dir = x * u + y * v - d * w;
			Eigen::Vector3d origin = camera.from;
			Eigen::Vector3d imagept = camera.from + dir;
			tempRay.direction = imagept - origin;
			tempRay.direction.normalize();

			tempRay.origin = camera.from;

			//for each pixel, loop through each triangle (polygon) and see if there is an intersection
			hitIterator = 0;
			tempRay.hit = false;
			int indexOfPolygon;

			//loops through all polygons
			for (int k = 0; k < worldShapes.size(); k++)
			{
				//every time another t value is added, an polygon is added as well
				intersect = intersections(&tempRay, &worldShapes[k], camera.hither, MAX);

				if (intersect == true)
				{
					tempRay.hit = true;	
					indexOfPolygon = k;
				}
			}
	

			//if ray has a hit
			if (tempRay.hit == true)
			{
				//reset hit flag
				tempRay.hit = false;
				double color;

				double r;
				double g;
				double b;

				tempRay.IntersectionLocation = tempRay.origin + (tempRay.minT * tempRay.direction);

				/*r = shade(&tempRay, &tempRay.hitList[tempRay.indexOfFirstHit], 0);
				g = shade(&tempRay, &tempRay.hitList[tempRay.indexOfFirstHit], 1);
				b = shade(&tempRay, &tempRay.hitList[tempRay.indexOfFirstHit], 2);

				pixels[j][i][0] = r * 255;
				pixels[j][i][1] = g * 255;
				pixels[j][i][2] = b * 255;*/


				//call recursive function
				double newColorr = reflect(0, &tempRay, &tempRay.hitList[tempRay.indexOfFirstHit]); //r
				//reset variables for next recursion
				recursionCounter = 0;
				tempRay.bounces = 0;
				totalColor = 0;

				double newColorg = reflect(1, &tempRay, &tempRay.hitList[tempRay.indexOfFirstHit]);//g	
				recursionCounter = 0;
				tempRay.bounces = 0;
				totalColor = 0;

				double newColorb = reflect(2, &tempRay, &tempRay.hitList[tempRay.indexOfFirstHit]); //b
				
				//assign color
				pixels[j][i][0] = newColorr *255;
				pixels[j][i][1] = newColorg *255;
				pixels[j][i][2] = newColorb *255;
			}

			//if ray had no hits
			else
			{
				pixels[j][i][0] = bgColor[0] * 255;
				pixels[j][i][1] = bgColor[1] * 255;
				pixels[j][i][2] = bgColor[2] * 255;
			}
		}
	}


	//writes pixels to ppm file
	FILE* k = fopen("C:\\Users\\bluej\\OneDrive\\Documents\\UMBC\\SENIOR+\\435\\proj1435\\hide.ppm", "wb");
	fprintf(k, "P6\n%d %d\n%d\n", 512, 512, 255); //HARD CODED RESOLIUTION
	fwrite(pixels, 1, 512 * 512 * 3, k);
	fclose(k);

	return 0;
}