// gga_diff.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "nemagga.h"

using namespace gnssimu_lib;

void print_kml_heder(FILE* fKML)
{
	// write header for KML 
	if (fKML) {
		fprintf(fKML, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
		fprintf(fKML, "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n");
		fprintf(fKML, "<Document>\n");
		// fprintf(fKML, "<Placemark>\n");    
		// fprintf(fKML, "<name>extruded</name>\n");
		// fprintf(fKML, "<LineString>\n");
		// fprintf(fKML, "<extrude>1</extrude>\n");
		// fprintf(fKML, "<tessellate>1</tessellate>\n");
		// fprintf(fKML, "<altitudeMode>relativeToGround</altitudeMode>\n");
		// fprintf(fKML, "<coordinates>\n"); 
		fprintf(fKML, "<Style id=\"solu1\">\n");
		fprintf(fKML, "<IconStyle>\n");
		fprintf(fKML, "<color>ffff00ff</color>\n");
		fprintf(fKML, "<scale>0.300</scale>\n");
		fprintf(fKML, "<Icon>\n");
		fprintf(fKML, "<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>\n");
		fprintf(fKML, "</Icon>\n");
		fprintf(fKML, "</IconStyle>\n");
		fprintf(fKML, "</Style>\n");
		fprintf(fKML, "<Style id=\"solu11\">\n");
		fprintf(fKML, "<IconStyle>\n");
		fprintf(fKML, "<color>ffff00ff</color>\n");
		fprintf(fKML, "<scale>0.500</scale>\n");
		fprintf(fKML, "<Icon>\n");
		fprintf(fKML, "<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>\n");
		fprintf(fKML, "</Icon>\n");
		fprintf(fKML, "</IconStyle>\n");
		fprintf(fKML, "</Style>\n");
		fprintf(fKML, "<Style id=\"solu4\">\n");
		fprintf(fKML, "<IconStyle>\n");
		fprintf(fKML, "<color>ff008800</color>\n");
		fprintf(fKML, "<scale>0.300</scale>\n");
		fprintf(fKML, "<Icon>\n");
		fprintf(fKML, "<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>\n");
		fprintf(fKML, "</Icon>\n");
		fprintf(fKML, "</IconStyle>\n");
		fprintf(fKML, "</Style>\n");
		fprintf(fKML, "<Style id=\"solu5\">\n");
		fprintf(fKML, "<IconStyle>\n");
		fprintf(fKML, "<color>ff00aaff</color>\n");
		fprintf(fKML, "<scale>0.300</scale>\n");
		fprintf(fKML, "<Icon>\n");
		fprintf(fKML, "<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>\n");
		fprintf(fKML, "</Icon>\n");
		fprintf(fKML, "</IconStyle>\n");
		fprintf(fKML, "</Style>\n");
		fprintf(fKML, "<Style id=\"solu15\">\n");
		fprintf(fKML, "<IconStyle>\n");
		fprintf(fKML, "<color>ff00aaff</color>\n");
		fprintf(fKML, "<scale>0.500</scale>\n");
		fprintf(fKML, "<Icon>\n");
		fprintf(fKML, "<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>\n");
		fprintf(fKML, "</Icon>\n");
		fprintf(fKML, "</IconStyle>\n");
		fprintf(fKML, "</Style>\n");
	}
	return;
}

void print_kml_gga(FILE* fKML, double lat, double lon, double ht, int type)
{
	if (fKML == NULL) return;
	if (fKML) {
		fprintf(fKML, "<Placemark>\n");
		if (type == 1)
		{
			fprintf(fKML, "<styleUrl>#solu1</styleUrl>\n");
		}
		else if (type == 4)
		{
			fprintf(fKML, "<styleUrl>#solu4</styleUrl>\n");
		}
		else if (type == 5)
		{
			fprintf(fKML, "<styleUrl>#solu5</styleUrl>\n");
		}
		else if (type == 11)
		{
			fprintf(fKML, "<styleUrl>#solu11</styleUrl>\n");
		}
		else if (type == 15)
		{
			fprintf(fKML, "<styleUrl>#solu15</styleUrl>\n");
		}
		else
		{
			fprintf(fKML, "<styleUrl>#solu1</styleUrl>\n");
		}
		//fprintf(fKML, "<ExtendedData>\n");
		//fprintf(fKML, "<Data name=\"time\">\n");
		//fprintf(fKML, "<value>%.3lf</value>\n", gga.time);
		//fprintf(fKML, "</Data>\n");
		//fprintf(fKML, "</ExtendedData>\n");
		fprintf(fKML, "<Point>\n");
		fprintf(fKML, "<coordinates>%14.9f,%14.9f,%14.4f</coordinates>\n", lon, lat, ht);
		fprintf(fKML, "</Point>\n");
		fprintf(fKML, "</Placemark>\n");
	}
	return;
}

void print_kml_eof(FILE* fKML)
{
	if (fKML)
	{
		// fprintf(fKML, "</coordinates>\n");    
		// fprintf(fKML, "</LineString>\n");
		// fprintf(fKML, "</Placemark>\n");
		fprintf(fKML, "</Document>\n");
		fprintf(fKML, "</kml>\n");

	}
}


void	att2C_nb(double* att, double C_nb[3][3])
{
	// attitude: roll, pitch and heading
	// attitude => C_nb
	double R = att[0], P = att[1], H = att[2];
	C_nb[0][0] = cos(H) * cos(P);
	C_nb[1][0] = sin(H) * cos(P);
	C_nb[2][0] = -sin(P);
	C_nb[0][1] = -sin(H) * cos(R) + cos(H) * sin(P) * sin(R);
	C_nb[1][1] = cos(H) * cos(R) + sin(H) * sin(P) * sin(R);
	C_nb[2][1] = cos(P) * sin(R);
	C_nb[0][2] = sin(H) * sin(R) + cos(H) * sin(P) * cos(R);
	C_nb[1][2] = -cos(H) * sin(R) + sin(H) * sin(P) * cos(R);
	C_nb[2][2] = cos(P) * cos(R);
	return;
}

void rotate_vector1(double C[3][3], double* vec, double *vec_out, int isTranspose)
{
	if (isTranspose == 0)
	{
		vec_out[0] = C[0][0] * vec[0] + C[0][1] * vec[1] + C[0][2] * vec[2];
		vec_out[1] = C[1][0] * vec[0] + C[1][1] * vec[1] + C[1][2] * vec[2];
		vec_out[2] = C[2][0] * vec[0] + C[2][1] * vec[1] + C[2][2] * vec[2];
	}
	else
	{
		vec_out[0] = C[0][0] * vec[0] + C[1][0] * vec[1] + C[2][0] * vec[2];
		vec_out[1] = C[0][1] * vec[0] + C[1][1] * vec[1] + C[2][1] * vec[2];
		vec_out[2] = C[0][2] * vec[0] + C[1][2] * vec[1] + C[2][2] * vec[2];
	}
	return;
}

void rotate_vector(double C[3][3], double* vec, int isTranspose)
{
	double temp[3] = { 0.0 };
	rotate_vector1(C, vec, temp, isTranspose);
	vec[0] = temp[0];
	vec[1] = temp[1];
	vec[2] = temp[2];
	return;
}


typedef struct
{
	int wn;
	double ws;
	double rpy[3]; /* Roll, Pitch and Yaw/Heading */
	double vxyz[3];
	double vENU[3];
	double wxyz[3];
	double fxyz[3];
	double xyz[3];
	double rmsRPY[3];
	double rmsENU[3];
	double rmsVenu[3];
	double blh[3];
	int ns;
	double hdop;
}solu1_t;

typedef struct
{
	int wn;
	double ws;
	double blh[3];
	double hdop;
	double pdop;
	int ns;
	double rmsENU[3];
	int type;
}solu2_t;

typedef struct
{
	double ws;
	double blh[3];
	double m[3];
	double s[3];
}solu_stat_t;

typedef struct
{
	double ned[3];
}coord_t;

#define MAX_TIME 120

void get_stat_data(coord_t* coord, double* m, double* s)
{
	int i, j, k = 0;
	for (k = 0; k < 3; ++k)
	{
		m[k] = 0.0;
		for (i = 0; i < MAX_TIME; ++i)
		{
			m[k] += coord[i].ned[k];
		}
		m[k] /= MAX_TIME;
	}
	for (k = 0; k < 3; ++k)
	{
		s[k] = 0.0;
		for (i = 0; i < MAX_TIME; ++i)
		{
			s[k] += (coord[i].ned[k] - m[k]) * (coord[i].ned[k] - m[k]);
		}
		s[k] = sqrt(s[k] / (MAX_TIME - 1));
	}
	return;
}

unsigned long solu_diff(const char* fname1, const char* fname2, char *key, int type, double *lao)
{
	/* compare the soltion from NovAtel SPAN solution and NMEA GGA */
#if 0
	/* NMEA GGA */
	Week, GPSTime, Roll, Pitch, Heading, VX - ECEF, VY - ECEF, VZ - ECEF, VEast, VNorth, VUp, AngRateX, AngRateY, AngRateZ, AccBdyX, AccBdyY, AccBdyZ, X - ECEF, Y - ECEF, Z - ECEF, RollSD, PitchSD, HdngSD, SDEast, SDNorth, SDHeight, SD - VE, SD - VN, SD - VH, Latitude, Longitude, H - Ell, NS, HDOP
	(weeks, (sec), (deg), (deg), (deg), (m / s), (m / s), (m / s), (m / s), (m / s), (m / s), (deg / s), (deg / s), (deg / s), (m / s ^ 2), (m / s ^ 2), (m / s ^ 2), (m), (m), (m), (deg), (deg), (deg), (m), (m), (m), (m / s), (m / s), (m / s), (deg), (deg), (m), , (dop)
		2069, 423470, 0.391943239, 2.563304313, 60.98594284, 0, 0, 0, 0, 0, 0, 0.0041, -0.0044, 0.0023, 0.023, -0.038, -0.023, -2687728.289, -4281324.497, 3876272.875, 0.007077488, 0.007051909, 0.010994121, 0.011, 0.011, 0.021, 0.003, 0.003, 0.003, 37.66736958, -122.1197648, -19.819, 7, 1.4
#endif
	FILE * fSOL[5] = { NULL };
	char buffer[1024] = { 0 }, outfname[255] = { 0 }, kmlfname[255] = { 0 }, stafname[255] = { 0 };
	fSOL[0] = fopen(fname1, "r");
	fSOL[1] = fopen(fname2, "r");

	if (fSOL[0] == NULL || fSOL[1] == NULL)
	{
		if (fSOL[0] != NULL) fclose(fSOL[0]);
		if (fSOL[1] != NULL) fclose(fSOL[1]);
		return 0;
	}
	memcpy(buffer, fname1, strlen(fname1));

	char* result1 = strrchr(buffer, '.');
	if (result1 != NULL) result1[0] = '\0';

	result1 = strrchr(key, '\n');
	if (result1 != NULL) result1[0] = '\0';

	sprintf(outfname, "%s_%s.csv", buffer, key);
	sprintf(kmlfname, "%s_%s.kml", buffer, key);
	sprintf(stafname, "%s_%s.sta", buffer, key);
	fSOL[2] = fopen(outfname, "w");
	fSOL[3] = fopen(kmlfname, "w");
	fSOL[4] = fopen(stafname, "w");

	print_kml_heder(fSOL[3]);

	/* skip two lines */
	fgets(buffer, sizeof(buffer), fSOL[1]);
	fgets(buffer, sizeof(buffer), fSOL[1]);
	solu1_t sol[2] = { 0 };
	solu2_t sol2 = { 0 };
	double xyz_offset[] = { 0.1612,       -0.0036,        0.0859 };
	unsigned long numofepoch = 0;
	double preTime = 0.0;

	coord_t* coord = new coord_t[MAX_TIME];
	memset(coord, 0, sizeof(coord_t));
	int index = 0;
	bool isFirst = true;

	while (!feof(fSOL[0]))
	{
		memset(&sol2, 0, sizeof(sol2));
		memset(buffer, 0, sizeof(buffer));
		fgets(buffer, sizeof(buffer), fSOL[0]);
		if (sscanf(buffer, "%i,%lf,%lf,%lf,%lf,%lf,%lf,%i,%lf,%lf,%lf,%i"
			, &sol2.wn, &sol2.ws
			, sol2.blh + 0, sol2.blh + 1, sol2.blh + 2
			, &sol2.hdop, &sol2.pdop
			, &sol2.ns
			, sol2.rmsENU + 0, sol2.rmsENU + 1, sol2.rmsENU + 2
			, &sol2.type) < 12)
		{
			continue;
		}
		sol2.ws /= 1000.0;
		sol2.blh[0] *= PI / 180.0;
		sol2.blh[1] *= PI / 180.0;
		while (!feof(fSOL[1]))
		{
			if (sol[0].ws <= sol2.ws && sol[1].ws >= sol2.ws)
			{
				/* compare the solution */
				int ii = 0;
				break;
			}
			solu1_t sol1 = { 0 };
			memset(buffer, 0, sizeof(buffer));
			fgets(buffer, sizeof(buffer), fSOL[1]);
			if (sscanf(buffer, "%i,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%i,%lf"
				, &sol1.wn, &sol1.ws
				, sol1.rpy + 0, sol1.rpy + 1, sol1.rpy + 2
				, sol1.vxyz + 0, sol1.vxyz + 1, sol1.vxyz + 2
				, sol1.vENU + 0, sol1.vENU + 1, sol1.vENU + 2
				, sol1.wxyz + 0, sol1.wxyz + 1, sol1.wxyz + 2
				, sol1.fxyz + 0, sol1.fxyz + 1, sol1.fxyz + 2
				, sol1.xyz + 0, sol1.xyz + 1, sol1.xyz + 2
				, sol1.rmsRPY + 0, sol1.rmsRPY + 1, sol1.rmsRPY + 2
				, sol1.rmsENU + 0, sol1.rmsENU + 1, sol1.rmsENU + 2
				, sol1.rmsVenu + 0, sol1.rmsVenu + 1, sol1.rmsVenu + 2
				, sol1.blh + 0, sol1.blh + 1, sol1.blh + 2
				, &sol1.ns, &sol1.hdop) < 34)
			{
				continue;
			}
			sol1.blh[0] *= PI / 180.0;
			sol1.blh[1] *= PI / 180.0;
			sol[0] = sol[1];
			sol[1] = sol1;
		}
		if (sol[0].ws == 0.0 && sol[1].ws > sol2.ws)
		{
			printf("no overlap start at %10.3f,%s,%s\n", sol2.ws, fname1, fname2);
			break;
		}
		if (sol[0].ws <= sol2.ws && sol[1].ws >= sol2.ws)
		{
			double xyz[3] = { 0.0 };
			double blh[3] = { sol2.blh[0], sol2.blh[1], sol2.blh[2] };
			blh2xyz(blh, xyz);
			/* interpolate data to reference time using reference velocity */
			double dt1 = sol[0].ws - sol2.ws;
			double dt2 = sol[1].ws - sol2.ws;
			double dt = 0.0;
			double xyz_ref[3] = { 0 };
			double vxyz_ref[3] = { 0 };
			double att_ref[3] = { 0 };
			double ws = 0.0;
			if (fabs(dt1) > fabs(dt2))
			{
				for (int i = 0; i < 3; ++i)
				{
					xyz_ref[i] = sol[1].xyz[i];
					vxyz_ref[i] = sol[1].vxyz[i];
					att_ref[i] = sol[1].rpy[i] * PI / 180.0;
					ws = sol[1].ws;
				}
				dt = dt2;
			}
			else
			{
				for (int i = 0; i < 3; ++i)
				{
					xyz_ref[i] = sol[0].xyz[i];
					vxyz_ref[i] = sol[0].vxyz[i];
					att_ref[i] = sol[0].rpy[i] * PI / 180.0;
					ws = sol[0].ws;
				}
				dt = dt1;
			}
			double blh_ref[3] = { 0 };
			double C_en[3][3] = { 0 };
			xyz2blh(xyz_ref, blh_ref);
			blh2C_en(blh_ref, C_en);
			double xyz_[3] = { 0.0 };
			double blh_[3] = { 0.0 };
			for (int i = 0; i < 3; ++i)
				xyz_[i] = xyz[i] + vxyz_ref[i] * dt - xyz_offset[i];
			double C_nb[3][3] = { 0 };
			att2C_nb(att_ref, C_nb);
			double dned[3] = { 0.0 };
			double lao_b[3] = { 0.0 };
			double lao_n[3] = { 0.0 };
			double lao_e[3] = { 0.0 };
			if (lao != NULL)
			{
				lao_b[0] = lao[0];
				lao_b[1] = lao[1];
				lao_b[2] = lao[2];
			}
			rotate_vector1(C_nb, lao_b, lao_n, 0);
			rotate_vector1(C_en, lao_n, lao_e, 0);
			xyz_[0] += lao_e[0];
			xyz_[1] += lao_e[1];
			xyz_[2] += lao_e[2];
			xyz2blh(xyz_, blh_);

			double dxyz[3] = { xyz_[0] - xyz_ref[0], xyz_[1] - xyz_ref[1], xyz_[2] - xyz_ref[2] };

			rotate_vector1(C_en, dxyz, dned, 1); /* xyz => ned */

			double diffH = sqrt(dned[0] * dned[0] + dned[1] * dned[1]);
			double diff3 = sqrt(dned[0] * dned[0] + dned[1] * dned[1] + dned[2] * dned[2]);

			fprintf(fSOL[2], "%10.3f,%14.9f,%14.9f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%02i,%02i,%7.3f,%7.3f,%7.3f,%10.4f,%10.4f,%10.4f\n"
				, ws, blh_[0] * 180.0 / PI, blh_[1] * 180.0 / PI, blh_[2]
				, dned[0], dned[1], dned[2], diffH, diff3, sol2.ns, sol2.type, sol2.rmsENU[0], sol2.rmsENU[1], sol2.rmsENU[2], att_ref[0], att_ref[1], att_ref[2]
				);

			if (floor(sol2.ws + 0.5) != floor(preTime + 0.5))
			{
				print_kml_gga(fSOL[3], blh_[0] * 180.0 / PI, blh_[1] * 180.0 / PI, blh_[2], type);
				print_kml_gga(fSOL[3], blh_ref[0] * 180.0 / PI, blh_ref[1] * 180.0 / PI, blh_ref[2], 4);
			}
			if (index >= MAX_TIME)
			{
				index = 0;
				isFirst = false;
			}
			coord[index].ned[0] = dned[0];
			coord[index].ned[1] = dned[1];
			coord[index].ned[2] = dned[2];
			++index;
			if (!isFirst)
			{
				double m[3] = { 0.0 };
				double s[3] = { 0.0 };
				get_stat_data(coord, m, s);
				fprintf(fSOL[4], "%10.3f,%14.9f,%14.9f,%10.4f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f\n"
					, ws, blh_[0] * 180.0 / PI, blh_[1] * 180.0 / PI, blh_[2]
					, m[0], m[1], m[2], s[0], s[1], s[2]
				);
			}
			++numofepoch;
		}
		preTime = sol2.ws;
		int jj = 0;
	}

	print_kml_eof(fSOL[3]);
	if (fSOL[0] != NULL) fclose(fSOL[0]);
	if (fSOL[1] != NULL) fclose(fSOL[1]);
	if (fSOL[2] != NULL) fclose(fSOL[2]);
	if (fSOL[3] != NULL) fclose(fSOL[3]);
	if (fSOL[4] != NULL) fclose(fSOL[4]);

	delete[]coord;
	return numofepoch;
}

unsigned long solution_match_time(const char* fname1, const char* fname2, char *key)
{
	/* compare the soltion from NovAtel SPAN solution and NMEA GGA */
#if 0
	/* NMEA GGA */
	Week, GPSTime, Roll, Pitch, Heading, VX - ECEF, VY - ECEF, VZ - ECEF, VEast, VNorth, VUp, AngRateX, AngRateY, AngRateZ, AccBdyX, AccBdyY, AccBdyZ, X - ECEF, Y - ECEF, Z - ECEF, RollSD, PitchSD, HdngSD, SDEast, SDNorth, SDHeight, SD - VE, SD - VN, SD - VH, Latitude, Longitude, H - Ell, NS, HDOP
	(weeks, (sec), (deg), (deg), (deg), (m / s), (m / s), (m / s), (m / s), (m / s), (m / s), (deg / s), (deg / s), (deg / s), (m / s ^ 2), (m / s ^ 2), (m / s ^ 2), (m), (m), (m), (deg), (deg), (deg), (m), (m), (m), (m / s), (m / s), (m / s), (deg), (deg), (m), , (dop)
		2069, 423470, 0.391943239, 2.563304313, 60.98594284, 0, 0, 0, 0, 0, 0, 0.0041, -0.0044, 0.0023, 0.023, -0.038, -0.023, -2687728.289, -4281324.497, 3876272.875, 0.007077488, 0.007051909, 0.010994121, 0.011, 0.011, 0.021, 0.003, 0.003, 0.003, 37.66736958, -122.1197648, -19.819, 7, 1.4
#endif
		FILE * fSOL[3] = { NULL };
	char buffer[1024] = { 0 }, outfname[255] = { 0 };
	fSOL[0] = fopen(fname1, "r");
	fSOL[1] = fopen(fname2, "r");

	if (fSOL[0] == NULL || fSOL[1] == NULL)
	{
		if (fSOL[0] != NULL) fclose(fSOL[0]);
		if (fSOL[1] != NULL) fclose(fSOL[1]);
		return 0;
	}

	memcpy(buffer, fname1, strlen(fname1));

	char* result1 = strrchr(buffer, '.');
	if (result1 != NULL) result1[0] = '\0';

	result1 = strrchr(key, '\n');
	if (result1 != NULL) result1[0] = '\0';

	sprintf(outfname, "%s_%s.csv", buffer, key);
	fSOL[2] = fopen(outfname, "w");
	solu_stat_t sol[2] = { 0 };
	solu_stat_t sol2 = { 0 };
	solu_stat_t sol1 = { 0 };
	unsigned long numofepoch = 0;

	while (!feof(fSOL[0]))
	{
		memset(&sol2, 0, sizeof(sol2));
		memset(buffer, 0, sizeof(buffer));
		fgets(buffer, sizeof(buffer), fSOL[0]);
		if (sscanf(buffer, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf"
			, &sol2.ws
			, sol2.blh + 0, sol2.blh + 1, sol2.blh + 2
			, sol2.m + 0, sol2.m + 1, sol2.m + 2
			, sol2.s + 0, sol2.s + 1, sol2.s + 2) < 10)
		{
			continue;
		}
		sol2.blh[0] *= PI / 180.0;
		sol2.blh[1] *= PI / 180.0;
		if (sol2.ws == 422821.000)
			int kk = 0;
		while (!feof(fSOL[1]))
		{
			if (sol[0].ws <= sol2.ws && sol[1].ws >= sol2.ws)
			{
				/* compare the solution */
				int ii = 0;
				break;
			}
			memset(&sol1, 0, sizeof(sol1));
			memset(buffer, 0, sizeof(buffer));
			fgets(buffer, sizeof(buffer), fSOL[1]);
			if (sscanf(buffer, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf"
				, &sol1.ws
				, sol1.blh + 0, sol1.blh + 1, sol1.blh + 2
				, sol1.m + 0, sol1.m + 1, sol1.m + 2
				, sol1.s + 0, sol1.s + 1, sol1.s + 2) < 10)
			{
				continue;
			}
			sol1.blh[0] *= PI / 180.0;
			sol1.blh[1] *= PI / 180.0;
			sol[0] = sol[1];
			sol[1] = sol1;
		}
		if (sol[0].ws == 0.0 && sol[1].ws > sol2.ws)
		{
			printf("no overlap start at %10.3f,%s,%s\n", sol2.ws, fname1,fname2);
			break;
		}
		if (sol[0].ws <= sol2.ws && sol[1].ws >= sol2.ws)
		{
			double dt1 = sol[0].ws - sol2.ws;
			double dt2 = sol[1].ws - sol2.ws;
			if (fabs(dt1) > fabs(dt2))
			{
				sol1 = sol[1];
			}
			else
			{
				sol1 = sol[0];
			}

			fprintf(fSOL[2], "%10.3f,%14.9f,%14.9f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n"
				, sol2.ws, sol2.blh[0] * 180.0 / PI, sol2.blh[1] * 180.0 / PI, sol2.blh[2]
				, sol2.m[0], sol2.m[1], sol2.m[2], sol2.s[0], sol2.s[1], sol2.s[2]
				, sol1.m[0], sol1.m[1], sol1.m[2], sol1.s[0], sol1.s[1], sol1.s[2]
			);
		}
	}

	if (fSOL[0] != NULL) fclose(fSOL[0]);
	if (fSOL[1] != NULL) fclose(fSOL[1]);
	if (fSOL[2] != NULL) fclose(fSOL[2]);

	return numofepoch;
}

void process(const char* fname)
{
	FILE* fINI = fopen(fname, "r");
	char buffer[255];
	char fname1[255] = { 0 };
	char fname2[255] = { 0 };
	char fname3[255] = { 0 };
	char inp_dir[255] = { 0 };
	char out_dir[255] = { 0 };
	int soluType = 0, type = 0;
	double lao[3] = { 0.0 };

	while (fINI!=NULL && !feof(fINI))
	{
		memset(buffer, 0, sizeof(buffer));
		fgets(buffer, sizeof(buffer), fINI);
		if (strlen(buffer) < 2) continue;
		if (buffer[0] == '#') continue;

		memset(fname1, 0, sizeof(fname1));
		memset(fname2, 0, sizeof(fname2));
		memset(fname3, 0, sizeof(fname3));
		soluType = 0;
		type = 0;
		sscanf(buffer, "%i", &type);
		if (type == 1)
		{
			memset(lao, 0, sizeof(lao));
			strncpy(fname1, inp_dir, strlen(inp_dir));
			strncpy(fname2, inp_dir, strlen(inp_dir));
			sscanf(buffer, "%i,%[^\,],%[^\,],%[^\,],%i,%lf,%lf,%lf", &type, fname1 + strlen(inp_dir), fname2 + strlen(inp_dir), fname3, &soluType, lao + 0, lao + 1, lao + 2);
			solu_diff(fname1, fname2, fname3, soluType, lao);
		}
		if (type == 2)
		{
			strncpy(fname1, inp_dir, strlen(inp_dir));
			strncpy(fname2, inp_dir, strlen(inp_dir));
			sscanf(buffer, "%i,%[^\,],%[^\,],%[^\,]", &type, fname1 + strlen(inp_dir), fname2 + strlen(inp_dir), fname3);
			solution_match_time(fname1, fname2, fname3);
		}
		else if (type == 0)
		{
			/* input directory, effective after this command */
			sscanf(buffer, "%i,%[^\,]", &type, inp_dir);
			char* temp = strchr(inp_dir, '\n');
			if (temp != NULL) temp[0] = '\0';
			if (strlen(inp_dir) > 0)
			{
				if (inp_dir[strlen(inp_dir) - 1] != '\\')
				{
					inp_dir[strlen(inp_dir)] = '\\';
				}
			}
		}
	}
	if (fINI != NULL) fclose(fINI);
}

void process_batch(const char* fname)
{
	FILE* fINI = fopen(fname, "r");
	char buffer[255];
	char fname1[255] = { 0 };
	char fname2[255] = { 0 };
	char fname3[255] = { 0 };
	char inp_dir[255] = { 0 };
	char out_dir[255] = { 0 };
	int soluType = 0, type = 0;
	double lao[3] = { 0.0 };

	FILE* fLOG = fopen("rtk_dif.log", "w");
	FILE* fSUM = fopen("rtk_sum.log", "w");

	while (fINI && !feof(fINI))
	{
		memset(buffer, 0, sizeof(buffer));
		fgets(buffer, sizeof(buffer), fINI);
		if (strlen(buffer) < 2) continue;

		char *temp = strrchr(buffer, '\n');
		if (temp) temp[0] = '\0';

		temp = strrchr(buffer, '\r');
		if (temp) temp[0] = '\0';

		printf("%s\n", buffer);

		char old_filename[255] = { 0 };
		strcpy(old_filename, buffer);

		char filename[255] = { 0 };
		strcpy(filename, buffer);
		temp = strrchr(filename, '.');
		if (temp) temp[0] = '\0';

		char outfname_rts[255] = { 0 }, outfname_rtk[255] = { 0 }, outfname_dif[255] = { 0 };

		sprintf(outfname_rtk, "%s_RTK3_pos.nmea", filename);

		temp = strrchr(buffer, '.');
		if (temp) temp[0] = '-';

		sprintf(outfname_rts, "%s-rts.nmea", buffer);

		sprintf(outfname_dif, "%s_sol_diff.csv", filename);

		if (gga_diff(outfname_rtk, outfname_rts, outfname_dif, fLOG) && fSUM)
		{
			fprintf(fSUM, "%s\n", old_filename);
		}
	}
	if (fINI) fclose(fINI);
	if (fLOG) fclose(fLOG);
	if (fSUM) fclose(fSUM);
}

int main(int argc, char* argv[])
{
	//process_batch("D:\\rtk\\data\\RTK_Data\\data.ini");

	if (argc<2)
	{
	}
	else if (argc==2)
	{
		TNEMAGGAReader nmeaReader;
		nmeaReader.ReadGGA(argv[1]);
	}

	return 0;
}
